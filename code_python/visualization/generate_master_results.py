"""Build consolidated regular and cluster-OOD result datasets.

Importing this module does not build datasets. Running it executes the
selection in the ``__main__`` block; call ``build_master_performance_data``,
``build_master_ood_performance_data``, or
``extract_feature_importance_stability`` directly when using it from another
module.
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
DEFAULT_OOD_SCORE_METRICS = [
    "rmse",
    "r2",
    "mae",
    "cvpp_ama",
    "cdf_ama",
    "nll",
    "ece",
    "Cv",
    "RUSC",
    "sharpness",
    "run_time_sec",
]
OOD_DATA_SPLITS = ("ood", "iid")
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
    "HDF",
    "count_only",
    "count_and_fingerprint",
    "structure_and_count",
    "count_and_graph",
    "random_count_permutation",
    "random_fp_permutation",
    "unique_fp_permutation",
]
OODType = Literal[
    "structure_cluster_OOD",
    "continuous_cluster_OOD",
]
COUNT_ONLY_MODELS = ["RF", "XGBR", "NGB", "GPytorchMAP"]
COUNT_ONLY_MIXING_METHODS = ["sum", "product"]
RANDOM_PERMUTATION_MODELS = ["GPytorchMAP"]
RANDOM_PERMUTATION_CONFIGS = [
    ("TanimotoMatern32", "Matern32", "sum"),
    ("TanimotoMatern32", "Matern32", "product"),
    ("TanimotoMatern32", "Matern32", "(count:x)+(fp:x)"),
    ("TanimotoMatern32", "Matern32", "(count:+)x(fp:+)"),
]
RANDOM_PERMUTATION_FEATURE_SETS = {
    "random_count_permutation",
    "random_fp_permutation",
    "unique_fp_permutation",
}
HDF_MODELS = ("RF", "XGBR", "NGB")
HDF_DEPTHS = (2, 4, 6)
HDF_BITS = (128, 256, 512, 1024)
COUNT_AND_GRAPH_MODELS = {
    "GNN-GP": "(GRAPH)_GNN-GP",
    "GNN-GP-SN": "(GRAPH)_GNN-GP-SN",
    "GNN-MLP": "(GRAPH)_GNN-MLP",
    "GNN-MLP-Ens": "(GRAPH)_GNN-MLP-Ens",
    "MLP-GP": "(TABULAR)_MLP-GP",
    "MLP-GP-SN": "(TABULAR)_MLP-GP-SN",
}
ATTRIBUTION_MODEL_NAMES = {
    "gnn_gp": "GNN-GP",
    "gnn_gpsn": "GNN-GP-SN",
    "gnn_mlp": "GNN-MLP",
    "mlp_gp": "MLP-GP",
    "mlp_gpsn": "MLP-GP-SN",
}
ATTRIBUTION_METHOD_COLUMNS = {
    "SHAP": "feature_importance_SHAP_kendalls_w",
    "GP-SHAP": "feature_importance_GP_SHAP_kendalls_w",
    "IG": "feature_importance_IG_kendalls_w",
}
FEATURE_IMPORTANCE_STABILITY_COLUMNS = [
    *FEATURE_STABILITY_COLUMNS,
    "feature_importance_GP_SHAP_kendalls_w",
    "feature_importance_IG_kendalls_w",
]
FEATURE_IMPORTANCE_STABILITY_RESULT_COLUMNS = [
    "paper",
    "target",
    "model",
    "device",
    "fp kernel",
    "count kernel",
    "mixing method",
    *FEATURE_IMPORTANCE_STABILITY_COLUMNS,
]
OOD_MODELS = [*MODELS, *COUNT_AND_GRAPH_MODELS]
StructureRepresentation = Literal["ECFP", "Graph", "LE", "Mordred"]
STRUCTURE_REPRESENTATIONS = ("ECFP", "Graph", "LE", "Mordred")


def _validate_feature_set(feature_set: str) -> FeatureSet:
    valid_feature_sets = {
        "HDF",
        "count_only",
        "count_and_fingerprint",
        "structure_and_count",
        "count_and_graph",
        "random_count_permutation",
        "random_fp_permutation",
        "unique_fp_permutation",
    }
    if feature_set not in valid_feature_sets:
        choices = ", ".join(sorted(valid_feature_sets))
        raise ValueError(f"feature_set must be one of: {choices}.")
    return cast(FeatureSet, feature_set)


def _validate_ood_type(ood_type: str) -> OODType:
    valid_ood_types = {"structure_cluster_OOD", "continuous_cluster_OOD"}
    if ood_type not in valid_ood_types:
        choices = ", ".join(sorted(valid_ood_types))
        raise ValueError(f"ood_type must be one of: {choices}.")
    return cast(OODType, ood_type)


def _is_tree_model(model: str) -> bool:
    return model in TREE_MODELS


def _count_and_fingerprint_kernel_configs() -> List[tuple[str, str, str]]:
    return [
        (fp_kernel, count_kernel, mixing_method)
        for fp_kernel in fp_sk_kernels + fp_bit_kernels
        for count_kernel in count_kernels
        for mixing_method in FP_MIXING_METHODS
    ]


def _structure_and_count_configurations(
    model: str,
) -> List[
    tuple[
        StructureRepresentation,
        Optional[str],
        Optional[str],
        Optional[str],
    ]
]:
    """Return valid structure/count configurations for one model family."""
    if _is_tree_model(model):
        return [
            (representation, None, None, None)
            for representation in ("ECFP", "LE", "Mordred")
        ]
    if model == "MGK":
        return [
            ("Graph", "Graph", count_kernel, mixing_method)
            for count_kernel in count_kernels
            for mixing_method in MGK_MIXING_METHODS
        ]
    return [
        ("ECFP", fp_kernel, count_kernel, mixing_method)
        for fp_kernel, count_kernel, mixing_method
        in _count_and_fingerprint_kernel_configs()
    ]


def _score_file_stems(
    model: str,
    feature_set: FeatureSet,
    fp_kernel: Optional[str] = None,
    count_kernel: Optional[str] = None,
    mixing_method: Optional[str] = None,
    use_gpu: bool = False,
    structure_representation: Optional[StructureRepresentation] = None,
    hdf_depth: Optional[int] = None,
    hdf_bits: Optional[int] = None,
) -> List[str]:
    if feature_set == "HDF":
        if (
            use_gpu
            or model not in HDF_MODELS
            or hdf_depth not in HDF_DEPTHS
            or hdf_bits not in HDF_BITS
        ):
            return []
        return [
            f"(HDF{hdf_depth}_{hdf_bits}bit-COUNT)_{model}"
            "_hypOFF_Standard_Standard_scores"
        ]

    if feature_set == "structure_and_count":
        if structure_representation not in STRUCTURE_REPRESENTATIONS:
            return []
        if structure_representation in {"LE", "Mordred"}:
            if not _is_tree_model(model) or use_gpu:
                return []
            return [
                f"({structure_representation}-COUNT)_{model}"
                "_hypOFF_Standard_Standard_scores"
            ]
        if structure_representation == "Graph":
            if model != "MGK":
                return []
        elif model == "MGK":
            return []
        feature_set = "count_and_fingerprint"

    if feature_set == "count_and_graph":
        if model not in COUNT_AND_GRAPH_MODELS:
            return []
        device_suffix = "_GPU" if use_gpu else ""
        return [
            f"{COUNT_AND_GRAPH_MODELS[model]}_hypOFF_Standard_Standard"
            f"{device_suffix}_scores"
        ]

    if feature_set in RANDOM_PERMUTATION_FEATURE_SETS:
        configuration = (fp_kernel, count_kernel, mixing_method)
        if (
            not use_gpu
            or model not in RANDOM_PERMUTATION_MODELS
            or configuration not in RANDOM_PERMUTATION_CONFIGS
        ):
            return []
        return [
            (
                f"(ECFP3_count_512-COUNT)_"
                f"({model}_{fp_kernel}-{count_kernel}_{mixing_method})"
                f"_hypOFF_Standard_Standard_{feature_set}_GPU_scores"
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
    structure_representation: Optional[StructureRepresentation] = None,
    hdf_depth: Optional[int] = None,
    hdf_bits: Optional[int] = None,
) -> Optional[Path]:
    for stem in _score_file_stems(
        model,
        feature_set,
        fp_kernel,
        count_kernel,
        mixing_method,
        use_gpu,
        structure_representation,
        hdf_depth,
        hdf_bits,
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
    structure_representation: Optional[StructureRepresentation] = None,
    hdf_depth: Optional[int] = None,
    hdf_bits: Optional[int] = None,
) -> tuple[Optional[Dict[str, Any]], Optional[Path]]:
    score_path = _find_score_path(
        paper_dir,
        model,
        feature_set,
        fp_kernel,
        count_kernel,
        mixing_method,
        use_gpu,
        structure_representation,
        hdf_depth,
        hdf_bits,
    )
    if score_path is None:
        return None, None

    with score_path.open("r", encoding="utf-8") as file:
        return json.load(file), score_path


def _load_feature_stability_score_file(
    paper_dir: Path,
    model: str,
    fp_kernel: Optional[str] = None,
    count_kernel: Optional[str] = None,
    mixing_method: Optional[str] = None,
    use_gpu: bool = False,
) -> tuple[Optional[Dict[str, Any]], Optional[Path]]:
    """Load a COUNT-plus-fingerprint score file unless it is an LFS pointer."""
    score_path = _find_score_path(
        paper_dir,
        model,
        "count_and_fingerprint",
        fp_kernel,
        count_kernel,
        mixing_method,
        use_gpu,
    )
    if score_path is None:
        return None, None

    with score_path.open("r", encoding="utf-8") as file:
        first_line = file.readline()
        if first_line.startswith("version https://git-lfs.github.com/spec/"):
            return None, score_path
        file.seek(0)
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
    if not pairs:
        pairs.extend(
            (str(column), None)
            for column in columns
            if str(column).isdigit()
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


def _has_feature_importance_stability(values: Dict[str, Any]) -> bool:
    return any(
        value is not None and not pd.isna(value)
        for value in values.values()
    )


def _feature_importance_stability_row(
    paper: str,
    target: str,
    model: str,
    device: str,
    fp_kernel: Optional[str],
    count_kernel: Optional[str],
    mixing_method: Optional[str],
    values: Dict[str, Any],
) -> Dict[str, Any]:
    return {
        "paper": paper,
        "target": target,
        "model": model,
        "device": device,
        "fp kernel": fp_kernel,
        "count kernel": count_kernel,
        "mixing method": mixing_method,
        **dict.fromkeys(FEATURE_IMPORTANCE_STABILITY_COLUMNS),
        **values,
    }


def _attribution_feature_importance_stability_rows(
    paper: str,
    target: str,
    paper_dir: Path,
) -> List[Dict[str, Any]]:
    consistency_path = ensure_long_path(
        paper_dir / "attribution" / "summary" / "seed_consistency.csv"
    )
    if not consistency_path.exists():
        return []

    consistency = pd.read_csv(consistency_path)
    required_columns = {"method", "model", "kendall_w"}
    missing_columns = required_columns.difference(consistency.columns)
    if missing_columns:
        raise ValueError(
            f"{consistency_path} is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    consistency = consistency[
        consistency["model"].isin(ATTRIBUTION_MODEL_NAMES)
        & consistency["method"].isin(ATTRIBUTION_METHOD_COLUMNS)
    ].copy()
    duplicate_rows = consistency.duplicated(["model", "method"], keep=False)
    if duplicate_rows.any():
        duplicates = consistency.loc[
            duplicate_rows,
            ["model", "method"],
        ].drop_duplicates()
        duplicate_names = ", ".join(
            f"{row.model}/{row.method}"
            for row in duplicates.itertuples(index=False)
        )
        raise ValueError(
            f"{consistency_path} contains duplicate model/method rows: "
            f"{duplicate_names}"
        )

    rows = []
    for source_model, tidy_model in ATTRIBUTION_MODEL_NAMES.items():
        model_data = consistency[consistency["model"] == source_model]
        values = {}
        for record in model_data.itertuples(index=False):
            value = pd.to_numeric(record.kendall_w, errors="coerce")
            values[ATTRIBUTION_METHOD_COLUMNS[record.method]] = (
                None if pd.isna(value) else float(value)
            )
        if values:
            rows.append(
                _feature_importance_stability_row(
                    paper,
                    target,
                    tidy_model,
                    "CPU",
                    None,
                    None,
                    None,
                    values,
                )
            )
    return rows


def extract_feature_importance_stability(
    save_path: Optional[Path] = (
        RESULTS
        / "master_performance_data"
        / "feature_importance_stability"
    ),
) -> pd.DataFrame:
    """Extract Kendall's W feature-stability results into a separate table.

    Tree-model MDI and SHAP stability is calculated from CPU score files.
    GpyroHMC and GPytorchMAP lengthscale stability is calculated only from GPU
    COUNT-plus-fingerprint score files. Neural and deep-kernel model SHAP,
    GP-SHAP, and integrated-gradients stability is read from each task's
    ``attribution/summary/seed_consistency.csv`` file. The attribution model
    aliases are converted to their tidy display names.

    MGK is intentionally excluded because it has no corresponding feature
    importance. Only configurations with at least one stability value are
    returned. If ``save_path`` is provided, matching ``.pkl`` and ``.csv``
    files are written.
    """
    rows = []

    for paper_name, paper_info in PAPER.items():
        for target in paper_info["target"]:
            paper_dir = RESULTS / paper_name / target

            for model in (model for model in MODELS if _is_tree_model(model)):
                data, _ = _load_feature_stability_score_file(
                    paper_dir,
                    model,
                )
                values = _feature_row_values(data)
                if _has_feature_importance_stability(values):
                    rows.append(
                        _feature_importance_stability_row(
                            paper_name,
                            target,
                            model,
                            "CPU",
                            None,
                            None,
                            None,
                            values,
                        )
                    )

            gp_models = [
                model
                for model in MODELS
                if not _is_tree_model(model) and model != "MGK"
            ]
            for model in gp_models:
                for (
                    fp_kernel,
                    count_kernel,
                    mixing_method,
                ) in _count_and_fingerprint_kernel_configs():
                    data, _ = _load_feature_stability_score_file(
                        paper_dir,
                        model,
                        fp_kernel,
                        count_kernel,
                        mixing_method,
                        use_gpu=True,
                    )
                    values = _feature_row_values(data)
                    if _has_feature_importance_stability(values):
                        rows.append(
                            _feature_importance_stability_row(
                                paper_name,
                                target,
                                model,
                                "GPU",
                                fp_kernel,
                                count_kernel,
                                mixing_method,
                                values,
                            )
                        )

            rows.extend(
                _attribution_feature_importance_stability_rows(
                    paper_name,
                    target,
                    paper_dir,
                )
            )

    dataframe = pd.DataFrame(
        rows,
        columns=FEATURE_IMPORTANCE_STABILITY_RESULT_COLUMNS,
        dtype=object,
    )
    if save_path is not None:
        _save_master_performance_data(dataframe, Path(save_path))
    return dataframe


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
    running_time = None if data is None else data.get("run_time_sec")
    if running_time is None:
        seed_running_times = []
        for _, seed_data in _seed_items(data):
            fit_times = seed_data.get("fit_time")
            score_times = seed_data.get("score_time")
            if isinstance(fit_times, list):
                seed_running_times.append(
                    sum(fit_times)
                    + (sum(score_times) if isinstance(score_times, list) else 0)
                )
        if seed_running_times:
            running_time = float(np.mean(seed_running_times))

    scores["Running time (GPU)" if use_gpu else "Running time (CPU)"] = running_time
    scores["Running time (CPU)" if use_gpu else "Running time (GPU)"] = None
    return scores


def _master_result_columns(
    score_metrics: List[str],
    include_structure_representation: bool = False,
    include_hdf_parameters: bool = False,
) -> List[str]:
    metric_columns = [
        column
        for metric in score_metrics
        for column in (
            f"{metric}_avg",
            f"{metric}_stdev",
            f"{metric}_seed_fold_scores",
        )
    ]
    identifying_columns = BASE_MASTER_RESULT_COLUMNS.copy()
    if include_structure_representation:
        identifying_columns.insert(3, "structure_representation")
    if include_hdf_parameters:
        identifying_columns[3:3] = ["HDF depth", "HDF bits"]
    return (
        identifying_columns
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
    structure_representation: Optional[StructureRepresentation] = None,
    hdf_depth: Optional[int] = None,
    hdf_bits: Optional[int] = None,
) -> Dict[str, Any]:
    return {
        "paper": paper,
        "target": target,
        "model": model,
        "HDF depth": hdf_depth,
        "HDF bits": hdf_bits,
        "structure_representation": structure_representation,
        "fp kernel": fp_kernel,
        "count kernel": count_kernel,
        "mixing method": mixing_method,
        **_score_row_values(data, score_metrics, target, score_path, use_gpu),
    }


def _ood_score_row_values(
    data: Optional[Dict[str, Any]],
    score_metrics: List[str],
) -> Dict[str, Any]:
    """Return the raw cluster dictionaries for each OOD/IID metric.

    Current structure-cluster files contain one top-level training seed. If a
    future file contains multiple seeds, the original cluster dictionaries are
    retained below an additional seed-keyed dictionary rather than discarded.
    """
    values = {
        f"{split}_test_{metric}": None
        for split in OOD_DATA_SPLITS
        for metric in score_metrics
    }

    for split in OOD_DATA_SPLITS:
        for metric in score_metrics:
            seed_values = []
            for seed, seed_data in _seed_items(data):
                split_data = seed_data.get(split)
                if not isinstance(split_data, dict):
                    continue
                metric_key = f"test_{metric}"
                if metric_key in split_data:
                    seed_values.append((seed, split_data[metric_key]))

            output_column = f"{split}_test_{metric}"
            if len(seed_values) == 1:
                values[output_column] = seed_values[0][1]
            elif seed_values:
                values[output_column] = dict(seed_values)

    return values


def _ood_master_result_columns(score_metrics: List[str]) -> List[str]:
    return BASE_MASTER_RESULT_COLUMNS + [
        f"{split}_test_{metric}"
        for split in OOD_DATA_SPLITS
        for metric in score_metrics
    ]


def _ood_model_configurations(
    model: str,
) -> List[tuple[Optional[str], Optional[str], Optional[str]]]:
    if _is_tree_model(model) or model in COUNT_AND_GRAPH_MODELS:
        return [(None, None, None)]
    if model == "MGK":
        return [
            ("Graph", count_kernel, mixing_method)
            for count_kernel in count_kernels
            for mixing_method in MGK_MIXING_METHODS
        ]
    return _count_and_fingerprint_kernel_configs()


def build_master_ood_performance_data(
    save_path: Optional[Path] = (
        RESULTS
        / "master_ood_performance_data"
        / "Tree_GP_MLP_structure_cluster_OOD"
    ),
    score_metrics: Optional[List[str]] = None,
    ood_type: OODType = "structure_cluster_OOD",
) -> Dict[str, pd.DataFrame]:
    """Build CPU and GPU master datasets for the selected OOD split type.

    Each metric cell retains the JSON representation for that validation
    split: OOD cells map cluster IDs to scalar scores, while IID cells map
    cluster IDs to lists of scores. In addition to the tree, GP, and MGK
    configurations, the count-and-graph GNN/MLP models contribute one row per
    task with empty kernel metadata. Tree-model results are included in both
    device datasets from their CPU files. Every non-tree row uses only the
    requested device, with no cross-device fallback. Point-prediction GNN-MLP
    rows retain their available prediction metrics while uncertainty and
    calibration cells are left empty.

    When ``save_path`` is provided, the function writes
    ``<save_path>_CPU.{pkl,csv}`` and ``<save_path>_GPU.{pkl,csv}``, matching
    :func:`build_master_performance_data`. Lengthscales and feature importances
    are intentionally excluded. Missing result configurations are represented
    by rows whose metric cells are ``None``.
    """
    ood_type = _validate_ood_type(ood_type)

    metrics = list(
        dict.fromkeys(
            metric.removeprefix("test_")
            for metric in (
                DEFAULT_OOD_SCORE_METRICS
                if score_metrics is None
                else score_metrics
            )
        )
    )
    dataframes = {}
    for device in ("CPU", "GPU"):
        use_gpu = device == "GPU"
        rows = []

        for paper_name, paper_info in PAPER.items():
            for target in paper_info["target"]:
                ood_dir = RESULTS / paper_name / target / ood_type

                for model in OOD_MODELS:
                    for fp_kernel, count_kernel, mixing_method in (
                        _ood_model_configurations(model)
                    ):
                        is_count_and_graph = model in COUNT_AND_GRAPH_MODELS
                        data, _ = _load_score_file(
                            ood_dir,
                            model,
                            (
                                "count_and_graph"
                                if is_count_and_graph
                                else "count_and_fingerprint"
                            ),
                            fp_kernel,
                            count_kernel,
                            mixing_method,
                            use_gpu=use_gpu and not _is_tree_model(model),
                        )
                        rows.append(
                            {
                                "paper": paper_name,
                                "target": target,
                                "model": model,
                                "fp kernel": fp_kernel,
                                "count kernel": count_kernel,
                                "mixing method": mixing_method,
                                **_ood_score_row_values(data, metrics),
                            }
                        )

        dataframe = pd.DataFrame(
            rows,
            columns=_ood_master_result_columns(metrics),
            dtype=object,
        )
        dataframes[device] = dataframe

        if save_path is not None:
            output_path = Path(f"{save_path}_{device}")
            ensure_long_path(output_path.parent).mkdir(
                parents=True,
                exist_ok=True,
            )
            dataframe.to_pickle(
                ensure_long_path(Path(f"{output_path}.pkl"))
            )

            csv_data = dataframe.copy()
            metric_columns = _ood_master_result_columns(metrics)[
                len(BASE_MASTER_RESULT_COLUMNS):
            ]
            for column in metric_columns:
                csv_data[column] = csv_data[column].apply(
                    lambda value: (
                        json.dumps(value) if value is not None else None
                    )
                )
            csv_data.to_csv(
                ensure_long_path(Path(f"{output_path}.csv")),
                index=False,
            )

    return dataframes


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
    ``feature_set="structure_and_count"`` extends that comparison with LE and
    Mordred tree-model results. Its ``structure_representation`` column marks
    ECFP, Graph (MGK), LE, or Mordred for every row.
    ``feature_set="count_only"`` includes COUNT-only tree results and
    GPytorchMAP results for the Matern32, Matern52, and RBF count kernels with
    the ``sum`` and ``product`` mixing methods.
    ``feature_set="count_and_graph"`` includes the CPU and GPU GNN-GP,
    GNN-GP-SN, GNN-MLP, GNN-MLP-Ens, MLP-GP, and MLP-GP-SN results. The
    ``*-GP-SN`` models are the spectrally normalized SNGP variants, and
    GNN-MLP-Ens is the uncertainty-aware GNN-MLP ensemble. The point-prediction
    GNN-MLP is also retained; its uncertainty and calibration fields are left
    empty because its score and prediction files do not contain predictive
    standard deviations. These files have no kernel configuration metadata,
    so those master-data columns are left empty.
    ``feature_set="HDF"`` includes CPU tree-model results for every HDF depth
    and bit count listed in ``HDF_DEPTHS`` and ``HDF_BITS``. Those parameters
    are retained in the ``HDF depth`` and ``HDF bits`` columns.
    ``feature_set="random_fp_permutation"`` includes the GPU GPytorchMAP
    results trained after randomly permuting fingerprint features. It uses the
    TanimotoMatern32 fingerprint kernel, Matern32 count kernel, and the four
    mixing methods present in those result files.
    ``feature_set="random_count_permutation"`` uses the same GPU GPytorchMAP
    configurations for results trained after randomly permuting COUNT
    features.
    ``feature_set="unique_fp_permutation"`` uses the same GPU GPytorchMAP
    configurations for results trained after permuting unique fingerprint
    feature vectors.

    Tree-model results are included in both device datasets. GP and
    count-and-graph rows use only the requested device, with no cross-device
    fallback. HDF results are CPU-only. Feature-importance stability is
    excluded; use ``extract_feature_importance_stability`` to build its
    separate table.

    When ``save_path`` is provided, standard modes write
    ``<save_path>_CPU.{pkl,csv}`` and ``<save_path>_GPU.{pkl,csv}``. Random
    fingerprint and COUNT permutation results are GPU-only, so those modes
    write only ``<save_path>_GPU.{pkl,csv}``. Count-and-graph results write
    both CPU and GPU outputs. HDF results are written only to the CPU output.
    """
    feature_set = _validate_feature_set(feature_set)
    metrics = list(
        dict.fromkeys(
            DEFAULT_SCORE_METRICS if score_metrics is None else score_metrics
        )
    )

    devices = (
        ("GPU",)
        if feature_set in RANDOM_PERMUTATION_FEATURE_SETS
        else ("CPU",) if feature_set == "HDF"
        else ("CPU", "GPU")
    )
    if feature_set == "count_only":
        models = COUNT_ONLY_MODELS
    elif feature_set == "count_and_graph":
        models = list(COUNT_AND_GRAPH_MODELS)
    elif feature_set == "HDF":
        models = list(HDF_MODELS)
    elif feature_set in RANDOM_PERMUTATION_FEATURE_SETS:
        models = RANDOM_PERMUTATION_MODELS
    else:
        models = MODELS

    for device in devices:
        use_gpu = device == "GPU"
        rows = []

        for paper_name, paper_info in PAPER.items():
            for target in paper_info["target"]:
                paper_dir = RESULTS / paper_name / (
                    f"{target}_HDF" if feature_set == "HDF" else target
                )

                for model in models:
                    if feature_set == "HDF":
                        for hdf_depth in HDF_DEPTHS:
                            for hdf_bits in HDF_BITS:
                                data, score_path = _load_score_file(
                                    paper_dir,
                                    model,
                                    feature_set,
                                    hdf_depth=hdf_depth,
                                    hdf_bits=hdf_bits,
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
                                        hdf_depth=hdf_depth,
                                        hdf_bits=hdf_bits,
                                    )
                                )
                        continue

                    if feature_set == "structure_and_count":
                        for (
                            structure_representation,
                            fp_kernel,
                            count_kernel,
                            mixing_method,
                        ) in _structure_and_count_configurations(model):
                            row_uses_gpu = use_gpu and not _is_tree_model(model)
                            data, score_path = _load_score_file(
                                paper_dir,
                                model,
                                feature_set,
                                fp_kernel,
                                count_kernel,
                                mixing_method,
                                row_uses_gpu,
                                structure_representation,
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
                                    row_uses_gpu,
                                    structure_representation,
                                )
                            )
                        continue

                    if feature_set == "count_and_graph":
                        data, score_path = _load_score_file(
                            paper_dir,
                            model,
                            feature_set,
                            use_gpu=use_gpu,
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
                                use_gpu=use_gpu,
                            )
                        )
                        continue

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

                    if feature_set in RANDOM_PERMUTATION_FEATURE_SETS:
                        configurations = RANDOM_PERMUTATION_CONFIGS
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
            columns=_master_result_columns(
                metrics,
                include_structure_representation=(
                    feature_set == "structure_and_count"
                ),
                include_hdf_parameters=(feature_set == "HDF"),
            ),
            dtype=object,
        )
        if save_path is not None:
            _save_master_performance_data(
                dataframe,
                Path(f"{save_path}_{device}"),
            )


if __name__ == "__main__":
    # Extract feature stability separately from the performance tables:
    # extract_feature_importance_stability()

    # Choose one feature set and a distinct output name:
    # feature_set: FeatureSet = "count_and_graph"

    # output_name = {
    #     "HDF": "Tree_HDF",
    #     "structure_and_count": "Tree_and_GP_structure_and_count",
    #     "count_and_fingerprint": "Tree_and_GP_count_and_fingerprint",
    #     "count_and_graph": "MLP_GNN_GP_SNGP_count_and_graph",
    #     "count_only": "Tree_and_GP_COUNT_only",
    #     "random_count_permutation": "GP_random_count_permutation",
    #     "random_fp_permutation": "GP_random_fp_permutation",
    #     "unique_fp_permutation": "GP_unique_fp_permutation",
    # }[feature_set]

    # build_master_performance_data(
    #     save_path=RESULTS / "master_performance_data" / output_name,
    #     score_metrics=DEFAULT_SCORE_METRICS,
    #     feature_set=feature_set,
    # )

    ood_type: OODType = "structure_cluster_OOD"
    output_name = {
        "structure_cluster_OOD": "Tree_GP_MLP_structure_cluster_OOD",
        "continuous_cluster_OOD": "Tree_GP_MLP_continuous_cluster_OOD",
    }[ood_type]

    build_master_ood_performance_data(
        save_path=RESULTS / "master_ood_performance_data" / output_name,
        score_metrics=DEFAULT_OOD_SCORE_METRICS,
        ood_type=ood_type,
    )
