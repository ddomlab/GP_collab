"""Build the consolidated CPU and GPU performance-result datasets.

Importing this module does not build datasets. Running it executes the
selection in the ``__main__`` block; call ``build_master_performance_data``
directly when using it from another module.
"""

import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, cast

import numpy as np
import pandas as pd
from sklearn.metrics import r2_score

from visualization_setting import ensure_long_path
from visualize_gp import (
    FEATURE_STABILITY_COLUMNS,
    FP_MIXING_METHODS,
    MGK_MIXING_METHODS,
    MODELS,
    PAPER,
    TREE_MODELS,
    count_kernels,
    fp_bit_kernels,
    fp_sk_kernels,
    kendalls_w,
)


HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent.parent / "results"
TRAINING = HERE.parent / "training"
if str(TRAINING) not in sys.path:
    sys.path.insert(0, str(TRAINING))

from utils_uncertainty_calibration import (  # noqa: E402
    compute_Cv,
    compute_RUSC,
    compute_cdf_ama,
    compute_cvpp_ama,
    compute_ece,
    compute_sharpness,
    gaussian_nll,
)


DEFAULT_SCORE_METRICS = [
    "rmse",
    "r2",
    "mae",
    "cvpp_ama",
    "nll",
    "ece",
    "Cv",
    "RUSC",
]
BASE_MASTER_RESULT_COLUMNS = [
    "paper",
    "target",
    "model",
    "fp kernel",
    "count kernel",
    "mixing method",
]
OOF_SCORE_METRICS = [
    "r2",
    "ece",
    "cvpp_ama",
    "cdf_ama",
    "RUSC",
    "nll",
    "Cv",
    "sharpness",
]
OOF_SCORE_COLUMNS = [
    "OOF_R2" if metric == "r2" else f"OOF_{metric}"
    for metric in OOF_SCORE_METRICS
]
TIME_COLUMNS = ["Running time (GPU)", "Running time (CPU)"]
FeatureSet = Literal[
    "count_only",
    "count_and_fingerprint",
    "random_fp_permutation",
]
COUNT_ONLY_MODELS = ["RF", "XGBR", "NGB", "GPytorchMAP"]
COUNT_ONLY_MIXING_METHODS = ["sum", "product"]
RANDOM_FP_PERMUTATION_MODELS = ["GPytorchMAP"]
RANDOM_FP_PERMUTATION_CONFIGS = [
    ("TanimotoMatern32", "Matern32", "sum"),
    ("TanimotoMatern32", "Matern32", "product"),
    ("TanimotoMatern32", "Matern32", "(count:x)+(fp:x)"),
    ("TanimotoMatern32", "Matern32", "(count:+)x(fp:+)"),
]


def _validate_feature_set(feature_set: str) -> FeatureSet:
    valid_feature_sets = {
        "count_only",
        "count_and_fingerprint",
        "random_fp_permutation",
    }
    if feature_set not in valid_feature_sets:
        choices = ", ".join(sorted(valid_feature_sets))
        raise ValueError(f"feature_set must be one of: {choices}.")
    return cast(FeatureSet, feature_set)


def _is_tree_model(model: str) -> bool:
    return model in TREE_MODELS


def _count_and_fingerprint_kernel_configs() -> List[tuple[str, str, str]]:
    return [
        (fp_kernel, count_kernel, mixing_method)
        for fp_kernel in fp_sk_kernels + fp_bit_kernels
        for count_kernel in count_kernels
        for mixing_method in FP_MIXING_METHODS
    ]


def _score_file_stems(
    model: str,
    feature_set: FeatureSet,
    fp_kernel: Optional[str] = None,
    count_kernel: Optional[str] = None,
    mixing_method: Optional[str] = None,
    use_gpu: bool = False,
) -> List[str]:
    if feature_set == "random_fp_permutation":
        configuration = (fp_kernel, count_kernel, mixing_method)
        if (
            not use_gpu
            or model not in RANDOM_FP_PERMUTATION_MODELS
            or configuration not in RANDOM_FP_PERMUTATION_CONFIGS
        ):
            return []
        return [
            (
                f"(ECFP3_count_512-COUNT)_"
                f"({model}_{fp_kernel}-{count_kernel}_{mixing_method})"
                "_hypOFF_Standard_Standard_random_fp_permutation_GPU_scores"
            )
        ]

    if _is_tree_model(model):
        if use_gpu:
            return []
        input_name = (
            "COUNT"
            if feature_set == "count_only"
            else "ECFP3_count_512-COUNT"
        )
        return [f"({input_name})_{model}_hypOFF_Standard_Standard_scores"]

    if feature_set == "count_only":
        if (
            model != "GPytorchMAP"
            or count_kernel not in count_kernels
            or mixing_method not in COUNT_ONLY_MIXING_METHODS
        ):
            return []
        device_suffix = "_GPU" if use_gpu else ""
        return [
            (
                f"(COUNT)_({model}_{count_kernel}_{mixing_method})"
                f"_hypOFF_Standard_Standard{device_suffix}_scores"
            )
        ]

    if model == "MGK":
        if (
            not use_gpu
            or count_kernel not in count_kernels
            or mixing_method not in MGK_MIXING_METHODS
        ):
            return []
        return [
            (
                f"(MG-COUNT)_(MGK_Graph-{count_kernel}_{mixing_method})"
                f"{suffix}_hypOFF_Standard_Standard_GPU_scores"
            )
            for suffix in ("", "_mean")
        ]

    device_suffix = "_GPU" if use_gpu else ""
    return [
        (
            f"(ECFP3_count_512-COUNT)_"
            f"({model}_{fp_kernel}-{count_kernel}_{mixing_method})"
            f"{suffix}_hypOFF_Standard_Standard{device_suffix}_scores"
        )
        for suffix in ("", "_mean")
    ]


def _find_score_path(
    paper_dir: Path,
    model: str,
    feature_set: FeatureSet,
    fp_kernel: Optional[str] = None,
    count_kernel: Optional[str] = None,
    mixing_method: Optional[str] = None,
    use_gpu: bool = False,
) -> Optional[Path]:
    for stem in _score_file_stems(
        model,
        feature_set,
        fp_kernel,
        count_kernel,
        mixing_method,
        use_gpu,
    ):
        path = ensure_long_path(paper_dir / f"{stem}.json")
        if path.exists():
            return path
    return None


def _load_score_file(
    paper_dir: Path,
    model: str,
    feature_set: FeatureSet,
    fp_kernel: Optional[str] = None,
    count_kernel: Optional[str] = None,
    mixing_method: Optional[str] = None,
    use_gpu: bool = False,
) -> tuple[Optional[Dict[str, Any]], Optional[Path]]:
    score_path = _find_score_path(
        paper_dir,
        model,
        feature_set,
        fp_kernel,
        count_kernel,
        mixing_method,
        use_gpu,
    )
    if score_path is None:
        return None, None

    with score_path.open("r", encoding="utf-8") as file:
        return json.load(file), score_path


def _prediction_path_from_score_path(score_path: Optional[Path]) -> Optional[Path]:
    if score_path is None or not score_path.name.endswith("_scores.json"):
        return None

    prediction_stem = score_path.name.removesuffix("_scores.json") + "_predictions"
    for suffix in (".csv", ".json"):
        prediction_path = ensure_long_path(
            score_path.with_name(f"{prediction_stem}{suffix}")
        )
        if prediction_path.exists():
            return prediction_path
    return None


def _prediction_target_column(
    columns: List[str],
    target: Optional[str] = None,
) -> Optional[str]:
    if target is not None:
        for target_name in (target, target.removeprefix("target_")):
            if target_name in columns:
                return target_name

    non_prediction_columns = [
        column
        for column in columns
        if not re.match(r"^seed_.+_y_(pred|std)$", str(column))
    ]
    return non_prediction_columns[0] if len(non_prediction_columns) == 1 else None


def _prediction_seed_pairs(columns: List[str]) -> List[tuple[str, Optional[str]]]:
    pairs = []
    column_set = set(columns)
    for prediction_column in columns:
        match = re.match(r"^(seed_.+)_y_pred$", str(prediction_column))
        if match is None:
            continue
        standard_deviation_column = f"{match.group(1)}_y_std"
        pairs.append(
            (
                prediction_column,
                standard_deviation_column
                if standard_deviation_column in column_set
                else None,
            )
        )
    return pairs


def _pooled_oof_arrays(
    prediction_path: Optional[Path],
    target: Optional[str] = None,
) -> tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray]]:
    if prediction_path is None or prediction_path.suffix.lower() != ".csv":
        return None, None, None

    prediction_df = pd.read_csv(prediction_path)
    target_column = _prediction_target_column(
        prediction_df.columns.tolist(),
        target,
    )
    seed_pairs = _prediction_seed_pairs(prediction_df.columns.tolist())
    if target_column is None or not seed_pairs:
        return None, None, None

    y_true = pd.to_numeric(prediction_df[target_column], errors="coerce").to_numpy()
    pooled_true = []
    pooled_prediction = []
    pooled_standard_deviation = []
    has_standard_deviation = False

    for prediction_column, standard_deviation_column in seed_pairs:
        y_prediction = pd.to_numeric(
            prediction_df[prediction_column],
            errors="coerce",
        ).to_numpy()
        if standard_deviation_column is None:
            y_standard_deviation = np.full_like(
                y_prediction,
                np.nan,
                dtype=float,
            )
        else:
            y_standard_deviation = pd.to_numeric(
                prediction_df[standard_deviation_column],
                errors="coerce",
            ).to_numpy()
            has_standard_deviation = True

        valid_mask = np.isfinite(y_true) & np.isfinite(y_prediction)
        pooled_true.extend(y_true[valid_mask])
        pooled_prediction.extend(y_prediction[valid_mask])
        pooled_standard_deviation.extend(y_standard_deviation[valid_mask])

    if len(pooled_true) < 2:
        return None, None, None

    return (
        np.asarray(pooled_true, dtype=float),
        np.asarray(pooled_prediction, dtype=float),
        (
            np.asarray(pooled_standard_deviation, dtype=float)
            if has_standard_deviation
            else None
        ),
    )


def _safe_oof_metric(
    metric: str,
    y_true: Optional[np.ndarray],
    y_prediction: Optional[np.ndarray],
    y_standard_deviation: Optional[np.ndarray],
) -> Optional[float]:
    if y_true is None or y_prediction is None:
        return None

    try:
        if metric == "r2":
            value = r2_score(y_true, y_prediction)
        elif metric in {"Cv", "sharpness"}:
            if y_standard_deviation is None:
                return None
            valid_mask = np.isfinite(y_standard_deviation)
            if valid_mask.sum() < 2:
                return None
            value = (
                compute_Cv(y_standard_deviation[valid_mask])
                if metric == "Cv"
                else compute_sharpness(y_standard_deviation[valid_mask])
            )
        else:
            if y_standard_deviation is None:
                return None
            valid_mask = (
                np.isfinite(y_true)
                & np.isfinite(y_prediction)
                & np.isfinite(y_standard_deviation)
            )
            if valid_mask.sum() < 2:
                return None
            metric_functions = {
                "ece": compute_ece,
                "cvpp_ama": compute_cvpp_ama,
                "cdf_ama": compute_cdf_ama,
                "RUSC": compute_RUSC,
                "nll": gaussian_nll,
            }
            value = metric_functions[metric](
                y_true[valid_mask],
                y_prediction[valid_mask],
                y_standard_deviation[valid_mask],
            )
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        return None

    return None if pd.isna(value) else float(value)


def _pooled_oof_scores(
    prediction_path: Optional[Path],
    target: Optional[str] = None,
) -> Dict[str, Optional[float]]:
    y_true, y_prediction, y_standard_deviation = _pooled_oof_arrays(
        prediction_path,
        target,
    )
    return {
        ("OOF_R2" if metric == "r2" else f"OOF_{metric}"): _safe_oof_metric(
            metric,
            y_true,
            y_prediction,
            y_standard_deviation,
        )
        for metric in OOF_SCORE_METRICS
    }


def _seed_sort_key(value: Any) -> tuple:
    return (0, int(value)) if str(value).isdigit() else (1, str(value))


def _seed_items(data: Optional[Dict[str, Any]]) -> List[tuple[str, Dict[str, Any]]]:
    if data is None:
        return []
    return [
        (key, data[key])
        for key in sorted(
            (key for key, value in data.items() if isinstance(value, dict)),
            key=_seed_sort_key,
        )
    ]


def _seed_fold_scores(
    data: Optional[Dict[str, Any]],
    metric: str,
) -> Optional[List[Any]]:
    scores = []
    for _, seed_data in _seed_items(data):
        values = seed_data.get(f"test_{metric}")
        if isinstance(values, list):
            scores.extend(values)
    return scores or None


def _seed_fold_feature_records(
    data: Optional[Dict[str, Any]],
    key: str,
) -> Optional[List[Dict[str, Any]]]:
    if data is None:
        return None

    records = []
    top_level_values = data.get(key)
    if isinstance(top_level_values, list):
        records.extend(value for value in top_level_values if isinstance(value, dict))
    elif isinstance(top_level_values, dict):
        records.append(top_level_values)

    for _, seed_data in _seed_items(data):
        values = seed_data.get(key)
        if isinstance(values, list):
            records.extend(value for value in values if isinstance(value, dict))
        elif isinstance(values, dict):
            records.append(values)
    return records or None


def _feature_kendalls_w(records: Any) -> Optional[float]:
    if records is None or not records:
        return None
    value = kendalls_w(records)["Kendall's W"]
    return None if pd.isna(value) else value


def _feature_row_values(data: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    feature_keys = {
        "lengthscale_kendalls_w": "test_lengthscale",
        "feature_importance_MDI_kendalls_w": "test_feature_importance_MDI",
        "feature_importance_SHAP_kendalls_w": "test_feature_importance_SHAP",
    }
    return {
        output_column: _feature_kendalls_w(
            _seed_fold_feature_records(data, source_key)
        )
        for output_column, source_key in feature_keys.items()
    }


def _score_row_values(
    data: Optional[Dict[str, Any]],
    score_metrics: List[str],
    target: Optional[str],
    score_path: Optional[Path],
    use_gpu: bool,
) -> Dict[str, Any]:
    scores = {}
    for metric in score_metrics:
        scores.update(
            {
                f"{metric}_avg": None if data is None else data.get(f"{metric}_avg"),
                f"{metric}_stdev": (
                    None if data is None else data.get(f"{metric}_stdev")
                ),
                f"{metric}_seed_fold_scores": _seed_fold_scores(data, metric),
            }
        )

    scores.update(
        _pooled_oof_scores(
            _prediction_path_from_score_path(score_path),
            target,
        )
    )
    scores["Running time (GPU)" if use_gpu else "Running time (CPU)"] = (
        None if data is None else data.get("run_time_sec")
    )
    scores["Running time (CPU)" if use_gpu else "Running time (GPU)"] = None
    return scores


def _master_result_columns(score_metrics: List[str]) -> List[str]:
    metric_columns = [
        column
        for metric in score_metrics
        for column in (
            f"{metric}_avg",
            f"{metric}_stdev",
            f"{metric}_seed_fold_scores",
        )
    ]
    return (
        BASE_MASTER_RESULT_COLUMNS
        + FEATURE_STABILITY_COLUMNS
        + metric_columns
        + OOF_SCORE_COLUMNS
        + TIME_COLUMNS
    )


def _save_master_performance_data(data: pd.DataFrame, save_path: Path) -> None:
    save_path = Path(save_path)
    ensure_long_path(save_path.parent).mkdir(parents=True, exist_ok=True)
    data.to_pickle(ensure_long_path(Path(f"{save_path}.pkl")))

    csv_data = data.copy()
    for column in csv_data.filter(like="_seed_fold_scores").columns:
        csv_data[column] = csv_data[column].apply(
            lambda value: json.dumps(value) if value is not None else None
        )
    csv_data.to_csv(ensure_long_path(Path(f"{save_path}.csv")), index=False)


def _master_row(
    paper: str,
    target: str,
    model: str,
    fp_kernel: Optional[str],
    count_kernel: Optional[str],
    mixing_method: Optional[str],
    score_metrics: List[str],
    data: Optional[Dict[str, Any]],
    score_path: Optional[Path],
    use_gpu: bool,
) -> Dict[str, Any]:
    return {
        "paper": paper,
        "target": target,
        "model": model,
        "fp kernel": fp_kernel,
        "count kernel": count_kernel,
        "mixing method": mixing_method,
        **_feature_row_values(data),
        **_score_row_values(data, score_metrics, target, score_path, use_gpu),
    }


def build_master_performance_data(
    save_path: Optional[Path] = (
        RESULTS / "master_performance_data" / "Tree_and_GP"
    ),
    score_metrics: Optional[List[str]] = None,
    feature_set: FeatureSet = "count_and_fingerprint",
) -> None:
    """Build CPU and GPU master datasets for the selected input features.

    ``feature_set="count_and_fingerprint"`` includes tree, GPytorchMAP,
    GpyroHMC, and MGK results trained on COUNT plus fingerprints.
    ``feature_set="count_only"`` includes COUNT-only tree results and
    GPytorchMAP results for the Matern32, Matern52, and RBF count kernels with
    the ``sum`` and ``product`` mixing methods.
    ``feature_set="random_fp_permutation"`` includes the GPU GPytorchMAP
    results trained after randomly permuting fingerprint features. It uses the
    TanimotoMatern32 fingerprint kernel, Matern32 count kernel, and the four
    mixing methods present in those result files.

    Tree-model results are included in both device datasets. GP rows use only
    the requested device, with no cross-device fallback. Only precomputed
    Kendall's W values are retained for feature stability.

    When ``save_path`` is provided, standard modes write
    ``<save_path>_CPU.{pkl,csv}`` and ``<save_path>_GPU.{pkl,csv}``. Random
    fingerprint permutation results are GPU-only, so that mode writes only
    ``<save_path>_GPU.{pkl,csv}``.
    """
    feature_set = _validate_feature_set(feature_set)
    metrics = list(
        dict.fromkeys(
            DEFAULT_SCORE_METRICS if score_metrics is None else score_metrics
        )
    )

    devices = (
        ("GPU",)
        if feature_set == "random_fp_permutation"
        else ("CPU", "GPU")
    )
    if feature_set == "count_only":
        models = COUNT_ONLY_MODELS
    elif feature_set == "random_fp_permutation":
        models = RANDOM_FP_PERMUTATION_MODELS
    else:
        models = MODELS

    for device in devices:
        use_gpu = device == "GPU"
        rows = []

        for paper_name, paper_info in PAPER.items():
            for target in paper_info["target"]:
                paper_dir = RESULTS / paper_name / target

                for model in models:
                    if _is_tree_model(model):
                        data, score_path = _load_score_file(
                            paper_dir,
                            model,
                            feature_set,
                        )
                        rows.append(
                            _master_row(
                                paper_name,
                                target,
                                model,
                                None,
                                None,
                                None,
                                metrics,
                                data,
                                score_path,
                                use_gpu=False,
                            )
                        )
                        continue

                    if feature_set == "random_fp_permutation":
                        configurations = RANDOM_FP_PERMUTATION_CONFIGS
                    elif feature_set == "count_only":
                        configurations = (
                            (None, count_kernel, mixing_method)
                            for count_kernel in count_kernels
                            for mixing_method in COUNT_ONLY_MIXING_METHODS
                        )
                    elif model == "MGK":
                        configurations = (
                            ("Graph", count_kernel, mixing_method)
                            for count_kernel in count_kernels
                            for mixing_method in MGK_MIXING_METHODS
                        )
                    else:
                        configurations = _count_and_fingerprint_kernel_configs()

                    for fp_kernel, count_kernel, mixing_method in configurations:
                        data, score_path = _load_score_file(
                            paper_dir,
                            model,
                            feature_set,
                            fp_kernel,
                            count_kernel,
                            mixing_method,
                            use_gpu,
                        )
                        rows.append(
                            _master_row(
                                paper_name,
                                target,
                                model,
                                fp_kernel,
                                count_kernel,
                                mixing_method,
                                metrics,
                                data,
                                score_path,
                                use_gpu,
                            )
                        )

        dataframe = pd.DataFrame(
            rows,
            columns=_master_result_columns(metrics),
            dtype=object,
        )
        if save_path is not None:
            _save_master_performance_data(
                dataframe,
                Path(f"{save_path}_{device}"),
            )


if __name__ == "__main__":
    # Choose one feature set and a distinct output name:
    feature_set: FeatureSet = "random_fp_permutation"
    output_name = {
        "count_and_fingerprint": "Tree_and_GP_count_and_fingerprint",
        "count_only": "Tree_and_GP_COUNT_only",
        "random_fp_permutation": "GP_random_fp_permutation",
    }[feature_set]
    build_master_performance_data(
        save_path=RESULTS / "master_performance_data" / output_name,
        score_metrics=DEFAULT_SCORE_METRICS,
        feature_set=feature_set,
    )
