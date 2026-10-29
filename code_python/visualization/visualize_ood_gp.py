"""Compare model performance on cluster-OOD and size-matched IID tests.

Each device-specific master OOD pickle stores one scalar score per OOD cluster
and five IID random-seed scores per matching cluster. This module expands those
nested values into performance-profile cases and reports each model's profile
AUC separately for OOD and IID evaluation.
"""

from __future__ import annotations

import ast
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional
import warnings

import matplotlib.pyplot as plt
from matplotlib.colors import to_hex, to_rgb
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

try:  # Support both direct execution and package-style imports.
    from .visualization_setting import (
        apply_vertical_bar_gradients,
        display_model_name as _display_model_name,
        ensure_long_path,
        filter_selection as _filter_selection,
        lighten_color as _lighten_color,
        metric_higher_is_better as _higher_is_better,
        model_config_label as _configuration_label,
        model_config_sort_key as _configuration_sort_key,
        model_type_color as _model_type_color,
        selection_values as _selection_values,
        set_plot_style,
    )
except ImportError:
    from visualization_setting import (
        apply_vertical_bar_gradients,
        display_model_name as _display_model_name,
        ensure_long_path,
        filter_selection as _filter_selection,
        lighten_color as _lighten_color,
        metric_higher_is_better as _higher_is_better,
        model_config_label as _configuration_label,
        model_config_sort_key as _configuration_sort_key,
        model_type_color as _model_type_color,
        selection_values as _selection_values,
        set_plot_style,
    )


set_plot_style()

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent.parent / "results"
DEFAULT_RESULTS_FILE = (
    RESULTS
    / "master_ood_performance_data"
)
DEFAULT_SAVE_DIR = HERE / "ood_result_analysis"

ID_COLUMNS = ["paper", "target", "model", "fp kernel", "count kernel", "mixing method"]
CASE_COLUMNS = ["paper", "target", "cluster"]
EVALUATIONS = ("OOD", "IID")
HIGHEST_BAR_HATCH = "xxxxx"
# OOD-specific model palette. Keep this local so the shared visualization
# palette—and figures produced by other modules—remain unchanged.
OOD_MODEL_COLORS = {
    "tree": "#225075",
    "gp": "#f25d71",
    "gnn": "#4D4E4E",
    "GNN-GP": "#7C1F44",
    "GNN-GP-SN": "#7C1F44",
    "GNN-MLP": "#7BAA82",
    "GNN-MLP-Ens": "#7BAA82",
    "MLP-GP": "#FF7551",
    "MLP-GP-SN": "#FF7551",
}
# Match the hybridization palette and notation used by visualize_gp.py.
HYBRIDIZATION_METHOD_COLORS = {
    "sum": "#2A9C9D",
    "product": "#C77B7B",
    "averageProduct": "#EDC525",
    "(count:+)x(fp:x)": "#5B8FF9",
    "(count:+)x(graph:x)": "#5B8FF9",
    "(count:+)x(fp:+)": "#5B8FF9",
    "(count:+)x(graph:+)": "#61DDAA",
    "(count:x)+(fp:x)": "#F6BD16",
    "(count:x)+(graph:x)": "#F6BD16",
}
HYBRIDIZATION_METHOD_LABELS = {
    "sum": "ΣC + ΣF",
    "product": "ΠC × ΠF",
    "averageProduct": "Av(C) × F",
    "(count:+)x(fp:x)": "ΣC × ΠF",
    "(count:+)x(graph:x)": "ΣC × G",
    "(count:+)x(fp:+)": "ΣC × ΣF",
    "(count:+)x(graph:+)": "ΣC × ΣG",
    "(count:x)+(fp:x)": "ΠC + ΠF",
    "(count:x)+(graph:x)": "ΠC + G",
}


def _is_missing(value: Any) -> bool:
    return value is None or (
        not isinstance(value, (list, tuple, dict)) and bool(pd.isna(value))
    )


def _kernel_free_configuration_mask(df: pd.DataFrame) -> pd.Series:
    return df[["fp kernel", "count kernel", "mixing method"]].isna().all(axis=1)


def _normalize_kernel_triples(kernel_triples: Any) -> Optional[list[tuple[Any, Any, Any]]]:
    if kernel_triples is None:
        return None
    if (
        isinstance(kernel_triples, tuple)
        and len(kernel_triples) == 3
        and not isinstance(kernel_triples[0], (tuple, list))
    ):
        kernel_triples = [kernel_triples]

    triples = [tuple(values) for values in kernel_triples]
    if any(len(values) != 3 for values in triples):
        raise ValueError(
            "Each kernel_triples entry must contain "
            "(fp_kernel, count_kernel, mixing_method)."
        )
    return triples


def _filter_configurations(
    df: pd.DataFrame,
    model: Any = None,
    fp_kernels: Any = None,
    count_kernels: Any = None,
    mixing_methods: Any = None,
    kernel_triples: Any = None,
) -> pd.DataFrame:
    """Apply global or per-model kernel filters and retain kernel-free models."""
    filtered = _filter_selection(df, "model", model)
    kernel_free_mask = _kernel_free_configuration_mask(filtered)

    if isinstance(kernel_triples, Mapping):
        selected_by_model = {}
        for model_name, model_triples in kernel_triples.items():
            model_key = str(model_name).strip().lower()
            if model_key in selected_by_model:
                raise ValueError(
                    f"Duplicate model key in kernel_triples: {model_name!r}."
                )
            triples = _normalize_kernel_triples(model_triples)
            if triples is None:
                raise ValueError(
                    "Each model in kernel_triples must map to a triple or a "
                    "list of triples."
                )
            selected_by_model[model_key] = {
                (str(fp).lower(), str(count).lower(), str(mix).lower())
                for fp, count, mix in triples
            }

        row_triples = pd.Series(
            list(
                zip(
                    filtered["fp kernel"].astype(str).str.lower(),
                    filtered["count kernel"].astype(str).str.lower(),
                    filtered["mixing method"].astype(str).str.lower(),
                )
            ),
            index=filtered.index,
        )
        raw_model_keys = filtered["model"].astype(str).str.strip().str.lower()
        display_model_keys = (
            filtered["model"]
            .map(_display_model_name)
            .astype(str)
            .str.strip()
            .str.lower()
        )
        triple_mask = pd.Series(False, index=filtered.index)
        for model_key, selected_triples in selected_by_model.items():
            model_mask = raw_model_keys.eq(model_key) | display_model_keys.eq(
                model_key
            )
            triple_mask |= model_mask & row_triples.isin(selected_triples)
        filtered = filtered[kernel_free_mask | triple_mask].copy()
    else:
        triples = _normalize_kernel_triples(kernel_triples)
        if triples is None:
            triples = []

    if not isinstance(kernel_triples, Mapping) and triples:
        triple_mask = pd.Series(False, index=filtered.index)
        for fp_kernel, count_kernel, mixing_method in triples:
            triple_mask |= (
                filtered["fp kernel"].eq(fp_kernel)
                & filtered["count kernel"].eq(count_kernel)
                & filtered["mixing method"].eq(mixing_method)
            )
        filtered = filtered[kernel_free_mask | triple_mask].copy()

    for column, selection in [
        ("fp kernel", fp_kernels),
        ("count kernel", count_kernels),
        ("mixing method", mixing_methods),
    ]:
        selected = _selection_values(selection)
        if selected is None:
            continue
        kernel_free_mask = _kernel_free_configuration_mask(filtered)
        selected_lower = {str(value).lower() for value in selected}
        matches = filtered[column].astype(str).str.lower().isin(selected_lower)
        filtered = filtered[kernel_free_mask | matches].copy()

    return filtered


def _metric_name(metric: str) -> str:
    metric = str(metric).strip()
    for prefix in ("ood_test_", "iid_test_"):
        if metric.lower().startswith(prefix):
            metric = metric[len(prefix) :]
    return metric


def _metric_columns(df: pd.DataFrame, metric: str) -> tuple[str, str, str]:
    metric = _metric_name(metric)
    lower_lookup = {column.lower(): column for column in df.columns}
    requested = [f"ood_test_{metric}", f"iid_test_{metric}"]
    try:
        ood_column, iid_column = [lower_lookup[column.lower()] for column in requested]
    except KeyError as exc:
        available = sorted(
            column.removeprefix("ood_test_")
            for column in df.columns
            if column.startswith("ood_test_")
        )
        raise ValueError(
            f"Metric {metric!r} is unavailable. Available metrics: {available}"
        ) from exc
    return metric, ood_column, iid_column


def _normalize_evaluations(ood_and_iid: Optional[str]) -> tuple[str, ...]:
    """Normalize a requested evaluation split; ``None`` selects both splits."""
    if ood_and_iid is None:
        return EVALUATIONS
    evaluation = str(ood_and_iid).strip().upper()
    if evaluation not in EVALUATIONS:
        raise ValueError("ood_and_iid must be either 'OOD' or 'IID'.")
    return (evaluation,)


def _as_mapping(value: Any) -> Optional[Mapping[Any, Any]]:
    if isinstance(value, Mapping):
        return value
    if isinstance(value, str):
        try:
            parsed = ast.literal_eval(value)
        except (SyntaxError, ValueError):
            return None
        return parsed if isinstance(parsed, Mapping) else None
    return None


def _numeric_values(value: Any) -> list[float]:
    values: Iterable[Any]
    if isinstance(value, (list, tuple, np.ndarray, pd.Series)):
        values = value
    else:
        values = [value]

    numeric_values = []
    for item in values:
        try:
            numeric = float(item)
        except (TypeError, ValueError):
            continue
        if np.isfinite(numeric):
            numeric_values.append(numeric)
    return numeric_values


def _configuration_key(row: pd.Series) -> tuple[Any, Any, Any, Any]:
    return tuple(
        None if _is_missing(row[column]) else row[column]
        for column in ID_COLUMNS[2:]
    )


def _expand_scores(
    df: pd.DataFrame,
    metric: str,
    include_kernel_config: bool,
    evaluations: tuple[str, ...] = EVALUATIONS,
) -> pd.DataFrame:
    """Expand nested OOD cluster scalars and IID seed lists to long form."""
    metric, ood_column, iid_column = _metric_columns(df, metric)
    records: list[dict[str, Any]] = []

    for _, row in df.iterrows():
        configuration_key = _configuration_key(row)
        metadata = {
            "metric": metric,
            "configuration key": configuration_key,
            "configuration": _configuration_label(row, include_kernel_config),
            "model": row["model"],
            "fp kernel": row["fp kernel"],
            "count kernel": row["count kernel"],
            "mixing method": row["mixing method"],
        }

        score_columns = {"OOD": ood_column, "IID": iid_column}
        for evaluation in evaluations:
            score_column = score_columns[evaluation]
            cluster_scores = _as_mapping(row[score_column])
            if cluster_scores is None:
                continue
            for cluster, value in cluster_scores.items():
                for repeat, score in enumerate(_numeric_values(value)):
                    records.append(
                        {
                            **metadata,
                            "paper": row["paper"],
                            "target": row["target"],
                            "cluster": str(cluster),
                            "evaluation": evaluation,
                            "repeat": repeat,
                            "score": score,
                        }
                    )

    columns = [
        "paper", "target", "cluster", "evaluation", "repeat", "score",
        "metric", "configuration key", "configuration", "model", "fp kernel",
        "count kernel", "mixing method",
    ]
    return pd.DataFrame.from_records(records, columns=columns)


def _keep_common_complete_clusters(
    scores: pd.DataFrame,
    evaluations: tuple[str, ...] = EVALUATIONS,
) -> pd.DataFrame:
    """Keep cluster cases present for every configuration in selected splits."""
    number_of_configurations = scores["configuration key"].nunique()
    availability = (
        scores.groupby(["evaluation", *CASE_COLUMNS], dropna=False)["configuration key"]
        .nunique()
        .eq(number_of_configurations)
        .unstack("evaluation", fill_value=False)
    )
    for evaluation in evaluations:
        if evaluation not in availability:
            availability[evaluation] = False

    common_index = availability.index[availability[list(evaluations)].all(axis=1)]
    score_index = pd.MultiIndex.from_frame(scores[CASE_COLUMNS])
    return scores[score_index.isin(common_index)].copy()


def _select_scores(
    df: pd.DataFrame,
    metric: str,
    model: Any,
    fp_kernels: Any,
    count_kernels: Any,
    mixing_methods: Any,
    kernel_triples: Any,
    include_kernel_config: bool,
    complete_cases: bool,
    ood_and_iid: Optional[str] = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Filter, expand, and validate scores for the selected evaluation splits.

    Returns the long-form scores and one metadata row per configuration. Only
    configurations with scores for every requested evaluation are kept.
    """
    evaluations = _normalize_evaluations(ood_and_iid)
    if not isinstance(df, pd.DataFrame):
        raise TypeError(f"df must be a pandas DataFrame, got {type(df).__name__}.")
    missing = set(ID_COLUMNS).difference(df.columns)
    if missing:
        raise ValueError(f"Input data is missing columns: {sorted(missing)}")
    filtered = _filter_configurations(
        df,
        model=model,
        fp_kernels=fp_kernels,
        count_kernels=count_kernels,
        mixing_methods=mixing_methods,
        kernel_triples=kernel_triples,
    )
    if filtered.empty:
        raise ValueError("No rows match the requested model and kernel filters.")
    scores = _expand_scores(
        filtered,
        metric,
        include_kernel_config,
        evaluations=evaluations,
    )
    if scores.empty:
        raise ValueError("No OOD/IID scores match the requested model and kernel filters.")

    requested_models = _selection_values(model)
    if requested_models is not None:
        found = {str(value).lower() for value in scores["model"]}
        absent = [value for value in requested_models if str(value).lower() not in found]
        if absent:
            warnings.warn(
                f"Requested models have no {_metric_name(metric)!r} scores after "
                f"filtering: {absent}",
                RuntimeWarning,
                stacklevel=3,
            )

    selected_configurations = {
        _configuration_key(row): _configuration_label(row, include_kernel_config)
        for _, row in filtered.iterrows()
    }
    evaluation_coverage = (
        scores.groupby(["configuration key", "evaluation"])
        .size()
        .unstack("evaluation", fill_value=0)
    )
    for evaluation in evaluations:
        if evaluation not in evaluation_coverage:
            evaluation_coverage[evaluation] = 0
    eligible_configurations = set(
        evaluation_coverage.index[
            evaluation_coverage[list(evaluations)].gt(0).all(axis=1)
        ]
    )
    unavailable = [
        label
        for key, label in selected_configurations.items()
        if key not in eligible_configurations
    ]
    if unavailable:
        warnings.warn(
            f"Skipping configurations without {' and '.join(evaluations)} "
            f"{_metric_name(metric)!r} scores: {unavailable}",
            RuntimeWarning,
            stacklevel=3,
        )
    scores = scores[
        scores["configuration key"].map(
            lambda key: key in eligible_configurations
        )
    ].copy()
    if scores.empty:
        raise ValueError(
            f"No selected configuration has {' and '.join(evaluations)} "
            f"{_metric_name(metric)!r} scores."
        )

    metadata_columns = [
        "configuration key", "configuration", "model", "fp kernel",
        "count kernel", "mixing method",
    ]
    metadata = scores[metadata_columns].drop_duplicates("configuration key")
    if not include_kernel_config and metadata["model"].duplicated().any():
        duplicated = metadata.loc[metadata["model"].duplicated(False), "model"].unique()
        raise ValueError(
            "include_kernel_config=False requires one selected configuration per model; "
            f"multiple configurations were found for {duplicated.tolist()}."
        )

    duplicate_columns = ["evaluation", *CASE_COLUMNS, "repeat", "configuration key"]
    if scores.duplicated(duplicate_columns).any():
        raise ValueError(
            "Duplicate rows were found for the same evaluation/case/configuration. "
            "Narrow the model or kernel filters."
        )

    if complete_cases:
        scores = _keep_common_complete_clusters(scores, evaluations=evaluations)
        if scores.empty:
            raise ValueError(
                "No cluster cases are complete for every selected configuration in "
                f"{' and '.join(evaluations)}. Select fewer models/configurations or "
                "use complete_cases=False."
            )
    return scores, metadata.reset_index(drop=True)


def _ordered_metadata(metadata: pd.DataFrame, model: Any) -> pd.DataFrame:
    """Sort configurations for the x-axis and require unique display labels."""
    if metadata["configuration"].duplicated().any():
        duplicated = metadata.loc[
            metadata["configuration"].duplicated(False), "configuration"
        ].unique()
        raise ValueError(
            "The compact x-axis labels are not unique for the selected kernel "
            f"configurations: {duplicated.tolist()}. Select one configuration "
            "per GP family."
        )
    selected_models = _selection_values(model) or metadata["model"].tolist()
    model_order = {
        str(model_name).lower(): position
        for position, model_name in enumerate(selected_models)
    }
    ordered = metadata.assign(
        _sort_key=metadata.apply(
            lambda row: _configuration_sort_key(row, model_order),
            axis=1,
        )
    )
    return (
        ordered.sort_values("_sort_key", kind="stable")
        .drop(columns="_sort_key")
        .reset_index(drop=True)
    )


def _evaluation_colors(metadata: pd.DataFrame) -> dict[str, dict[str, str]]:
    """Return model-aware OOD colors and lighter matching IID colors."""
    ood_colors = {
        row["configuration"]: _model_type_color(
            row["model"], palette=OOD_MODEL_COLORS
        )
        for _, row in metadata.iterrows()
    }
    return {
        "OOD": ood_colors,
        "IID": {label: _lighten_color(color) for label, color in ood_colors.items()},
    }


def _evaluation_legend(ax: plt.Axes, **kwargs: Any) -> None:
    legend_color = "#666666"
    ax.legend(
        handles=[
            Patch(facecolor=legend_color, label="OOD"),
            Patch(facecolor=_lighten_color(legend_color), label="IID"),
        ],
        frameon=False,
        title="Evaluation",
        **kwargs,
    )


def _darken_color(color: Any, amount: float = 0.35) -> str:
    """Blend a color toward black while preserving its hue."""
    rgb = np.asarray(to_rgb(color), dtype=float)
    return to_hex(rgb * (1.0 - amount))


def _hatch_highest_bars(
    bars: Iterable[Any],
    values: Iterable[float],
    colors: Iterable[Any],
) -> None:
    """Hatch the highest finite bar, including all bars tied for highest."""
    bar_list = list(bars)
    numeric_values = np.asarray(list(values), dtype=float)
    bar_colors = list(colors)
    finite_values = numeric_values[np.isfinite(numeric_values)]
    if not bar_list or finite_values.size == 0:
        return

    highest_value = finite_values.max()
    for bar, value, color in zip(bar_list, numeric_values, bar_colors):
        if np.isfinite(value) and np.isclose(
            value, highest_value, rtol=1e-12, atol=1e-12
        ):
            bar.set_hatch(HIGHEST_BAR_HATCH)
            bar.set_edgecolor(_darken_color(color))
            bar.set_linewidth(0)


def _finalize_figure(
    fig: plt.Figure,
    save_dir: Optional[Path],
    file_name: str,
    high_quality: bool,
    show: bool,
) -> None:
    fig.tight_layout()
    if save_dir is not None:
        save_dir = ensure_long_path(Path(save_dir))
        save_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            ensure_long_path(save_dir / file_name),
            bbox_inches="tight",
            dpi=900 if high_quality else 150,
        )
    if show:
        plt.show()
    else:
        plt.close(fig)


def _performance_profile_results(
    df: pd.DataFrame,
    metric: str = "r2",
    model: Any = None,
    fp_kernels: Any = None,
    count_kernels: Any = None,
    mixing_methods: Any = None,
    kernel_triples: Any = None,
    include_kernel_config: bool = True,
    higher_is_better: Optional[bool] = None,
    complete_cases: bool = True,
    ood_and_iid: Optional[str] = None,
) -> dict[str, Any]:
    evaluations = _normalize_evaluations(ood_and_iid)
    scores, metadata = _select_scores(
        df,
        metric=metric,
        model=model,
        fp_kernels=fp_kernels,
        count_kernels=count_kernels,
        mixing_methods=mixing_methods,
        kernel_triples=kernel_triples,
        include_kernel_config=include_kernel_config,
        complete_cases=complete_cases,
        ood_and_iid=ood_and_iid,
    )
    if len(metadata) < 2:
        raise ValueError("A performance profile requires at least two model configurations.")

    metric = _metric_name(metric)
    if higher_is_better is None:
        higher_is_better = _higher_is_better(metric)
    auc_rows: list[dict[str, Any]] = []

    for evaluation in evaluations:
        evaluation_scores = scores[scores["evaluation"].eq(evaluation)]
        matrix = evaluation_scores.pivot(
            index=[*CASE_COLUMNS, "repeat"],
            columns="configuration key",
            values="score",
        )
        if complete_cases:
            matrix = matrix.dropna(axis=0, how="any")
        else:
            matrix = matrix[matrix.notna().sum(axis=1) >= 2]
        if matrix.empty:
            raise ValueError(f"No usable {evaluation} cases remain after filtering.")

        ranks = matrix.rank(axis=1, ascending=not higher_is_better, method="min")
        percentiles = ranks.div(matrix.notna().sum(axis=1), axis=0)
        for configuration_key in metadata["configuration key"]:
            if configuration_key not in percentiles:
                continue
            values = percentiles[configuration_key].dropna().to_numpy(dtype=float)
            if values.size == 0:
                continue
            meta = metadata.loc[
                metadata["configuration key"].map(
                    lambda value: value == configuration_key
                )
            ].iloc[0]
            auc_rows.append(
                {
                    "evaluation": evaluation,
                    "metric": metric,
                    "configuration": meta["configuration"],
                    "model": meta["model"],
                    "fp kernel": meta["fp kernel"],
                    "count kernel": meta["count kernel"],
                    "mixing method": meta["mixing method"],
                    "direction": (
                        "higher is better" if higher_is_better else "lower is better"
                    ),
                    # Integral of an empirical CDF on [0, 1].
                    "auc": float(np.mean(1.0 - values)),
                    "n_cases": int(values.size),
                    "n_cluster_cases": int(
                        evaluation_scores[CASE_COLUMNS].drop_duplicates().shape[0]
                    ),
                }
            )

    auc_df = pd.DataFrame(auc_rows)
    auc_df["rank"] = auc_df.groupby("evaluation")["auc"].rank(
        ascending=False, method="min"
    ).astype(int)
    return {
        "auc_df": auc_df,
        "scores": scores,
        "metadata": metadata,
    }


def calculate_performance_profile_auc(
    df: pd.DataFrame,
    metric: str = "r2",
    model: Any = None,
    fp_kernels: Any = None,
    count_kernels: Any = None,
    mixing_methods: Any = None,
    kernel_triples: Any = None,
    include_kernel_config: bool = True,
    higher_is_better: Optional[bool] = None,
    complete_cases: bool = True,
    ood_and_iid: Optional[str] = None,
) -> pd.DataFrame:
    """Return performance-profile AUCs for each requested evaluation split.

    OOD cases are ``paper x target x cluster``. IID cases additionally include
    the random-seed position in each cluster's score list. By default, only
    cluster cases available for every selected configuration are used. Pass
    ``ood_and_iid="OOD"`` or ``"IID"`` to calculate only one split; the default
    of ``None`` retains both splits for callers that need the numerical table.
    """
    return _performance_profile_results(
        df=df,
        metric=metric,
        model=model,
        fp_kernels=fp_kernels,
        count_kernels=count_kernels,
        mixing_methods=mixing_methods,
        kernel_triples=kernel_triples,
        include_kernel_config=include_kernel_config,
        higher_is_better=higher_is_better,
        complete_cases=complete_cases,
        ood_and_iid=ood_and_iid,
    )["auc_df"].copy()


def plot_model_profile_comparison(
    df: pd.DataFrame,
    model: Any = None,
    fp_kernels: Any = None,
    count_kernels: Any = None,
    mixing_methods: Any = None,
    kernel_triples: Any = None,
    metric: str = "r2",
    ood_and_iid: str = "OOD",
    include_kernel_config: bool = True,
    higher_is_better: Optional[bool] = None,
    complete_cases: bool = True,
    figsize: tuple[float, float] = (7, 5),
    fontsize: int = 12,
    title: Optional[str] = None,
    x_label: str = "Model",
    y_label: str = "Performance profile AUC",
    x_tick_rotation: int = 35,
    y_lim: tuple[float, float] = (0.0, 1.05),
    show_values: bool = True,
    show: bool = True,
    high_quality: bool = True,
    save_dir: Optional[Path] = DEFAULT_SAVE_DIR,
    file_name: Optional[str] = None,
) -> pd.DataFrame:
    """Plot model profile AUC for either the OOD or matched-IID evaluation.

    Exact GP configurations can be selected with ``kernel_triples``. Pass one
    triple, a list of triples, or a dictionary mapping each raw/display model
    name to its own triple selection. Kernel-free tree, GNN, and MLP models are
    retained automatically. ``ood_and_iid`` is case-insensitive and must be
    either ``"OOD"`` or ``"IID"``. The returned DataFrame contains only the
    plotted split, including its AUC, rank, and number of contributing cases.
    """
    evaluation = _normalize_evaluations(ood_and_iid)[0]
    results = _performance_profile_results(
        df=df,
        metric=metric,
        model=model,
        fp_kernels=fp_kernels,
        count_kernels=count_kernels,
        mixing_methods=mixing_methods,
        kernel_triples=kernel_triples,
        include_kernel_config=include_kernel_config,
        higher_is_better=higher_is_better,
        complete_cases=complete_cases,
        ood_and_iid=evaluation,
    )
    auc_df = results["auc_df"]
    metadata = _ordered_metadata(results["metadata"], model)
    order = metadata["configuration"].tolist()
    values = auc_df.set_index("configuration")["auc"].reindex(order)
    bar_colors = {
        row["configuration"]: _model_type_color(
            row["model"], palette=OOD_MODEL_COLORS
        )
        for _, row in metadata.iterrows()
    }

    fig, ax = plt.subplots(figsize=figsize)
    x_positions = np.arange(len(order), dtype=float)
    heights = values.to_numpy(dtype=float)
    bars = ax.bar(
        x_positions,
        heights,
        width=0.67,
        color=[bar_colors[label] for label in order],
    )
    apply_vertical_bar_gradients(ax, bars=bars, edge_linewidth=0)
    _hatch_highest_bars(
        bars,
        heights,
        [bar_colors[label] for label in order],
    )
    if show_values:
        ax.bar_label(
            bars,
            labels=["" if np.isnan(value) else f"{value:.2f}" for value in heights],
            padding=3,
            fontsize=fontsize - 2,
        )

    ax.set_xlabel(x_label, fontsize=fontsize, fontweight="bold")
    ax.set_ylabel(y_label, fontsize=fontsize - 2, fontweight="bold")
    if title:
        ax.set_title(title, fontsize=fontsize + 2)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(
        order,
        rotation=x_tick_rotation,
        ha="right" if x_tick_rotation else "center",
    )
    ax.set_ylim(*y_lim)
    ax.tick_params(axis="both", labelsize=fontsize - 2)

    _finalize_figure(
        fig,
        save_dir=save_dir,
        file_name=(
            file_name
            or f"{_metric_name(metric)}_{evaluation.lower()}_model_profile_auc.png"
        ),
        high_quality=high_quality,
        show=show,
    )
    return auc_df


def _draw_half_violin(
    ax: plt.Axes,
    values: np.ndarray,
    position: float,
    side: str,
    color: str,
    width: float,
    edge_color: str,
    inner_color: str,
) -> None:
    """Draw one half violin (``side='low'`` left, ``'high'`` right) with an inner box.

    The density is evaluated between the observed minimum and maximum (no
    extrapolation past the data). The inner box marks the interquartile range,
    whiskers at 1.5 IQR clipped to the data, and a white median dot.
    """
    if side not in {"low", "high"}:
        raise ValueError(f"side must be 'low' or 'high', got {side!r}.")
    values = np.asarray(values, dtype=float)
    if values.size == 0 or not np.isfinite(values).all():
        raise ValueError("Half violins require a non-empty array of finite values.")

    direction = -1.0 if side == "low" else 1.0
    half_width = width / 2
    if np.unique(values).size >= 2:
        parts = ax.violinplot(
            values,
            positions=[position],
            widths=width,
            side=side,
            showextrema=False,
            showmedians=False,
            points=200,
        )
        for body in parts["bodies"]:
            body.set_facecolor(color)
            body.set_edgecolor(edge_color)
            body.set_linewidth(1)
            body.set_alpha(1)
    else:
        # A KDE is undefined for constant data; mark the single value instead.
        ax.hlines(
            values[0],
            position,
            position + direction * half_width,
            color=color,
            linewidth=3,
        )

    q1, median, q3 = np.percentile(values, [25, 50, 75])
    iqr = q3 - q1
    lower_whisker = values[values >= q1 - 1.5 * iqr].min()
    upper_whisker = values[values <= q3 + 1.5 * iqr].max()
    inner_x = position + direction * half_width * 0.12
    ax.vlines(
        inner_x, lower_whisker, upper_whisker,
        color=inner_color, linewidth=1.5, capstyle="round", zorder=3,
    )
    ax.vlines(
        inner_x, q1, q3,
        color=inner_color, linewidth=5, capstyle="round", zorder=3,
    )
    ax.scatter(
        [inner_x], [median],
        s=12, color="white", edgecolors="none", zorder=4,
    )


def plot_model_comparison(
    df: pd.DataFrame,
    metric: str = "r2",
    model: Any = None,
    fp_kernels: Any = None,
    count_kernels: Any = None,
    mixing_methods: Any = None,
    kernel_triples: Any = None,
    include_kernel_config: bool = True,
    complete_cases: bool = True,
    average_iid_repeats: bool = False,
    figsize: tuple[float, float] = (9, 5),
    fontsize: int = 12,
    title: Optional[str] = None,
    x_label: str = "Model",
    y_label: Optional[str] = None,
    x_tick_rotation: int = 35,
    y_lim: Optional[tuple[float, float]] = None,
    log_y: bool = False,
    violin_width: float = 0.8,
    show: bool = True,
    high_quality: bool = True,
    save_dir: Optional[Path] = DEFAULT_SAVE_DIR,
    file_name: Optional[str] = None,
) -> pd.DataFrame:
    """Split-violin comparison of raw OOD (left) and IID (right) scores per model.

    Each OOD point is one ``paper x target x cluster`` score. Each IID point is
    one random-seed score of the size-matched IID split for that cluster; set
    ``average_iid_repeats=True`` to average the seeds so both halves contain
    one point per cluster. OOD halves use the shared model-aware color and IID
    halves a lighter shade of it. Set ``log_y=True`` for a log-scaled y-axis (all
    scores must be positive). ``kernel_triples`` accepts a global triple/list
    or a dictionary mapping each model to separate triples. Returns the
    long-form scores that were plotted.
    """
    if not 0 < violin_width <= 1:
        raise ValueError("violin_width must be in the interval (0, 1].")
    if y_lim is not None:
        if len(y_lim) != 2 or not y_lim[0] < y_lim[1]:
            raise ValueError(
                "y_lim must be a (lower, upper) pair with lower < upper."
            )

    scores, metadata = _select_scores(
        df,
        metric=metric,
        model=model,
        fp_kernels=fp_kernels,
        count_kernels=count_kernels,
        mixing_methods=mixing_methods,
        kernel_triples=kernel_triples,
        include_kernel_config=include_kernel_config,
        complete_cases=complete_cases,
    )
    metric = _metric_name(metric)
    if log_y and scores["score"].le(0).any():
        raise ValueError(f"log_y=True requires all {metric!r} scores to be positive.")
    # Validates unique labels, which the per-label grouping below relies on.
    metadata = _ordered_metadata(metadata, model)
    if average_iid_repeats:
        group_columns = ["evaluation", *CASE_COLUMNS, "configuration"]
        averaged = scores.groupby(group_columns, as_index=False, dropna=False)[
            "score"
        ].mean()
        averaged["repeat"] = 0
        scores = averaged.merge(
            scores.drop(columns=["score", "repeat"]).drop_duplicates(group_columns),
            on=group_columns,
            how="left",
            validate="one_to_one",
        )[scores.columns]

    order = metadata["configuration"].tolist()
    colors = _evaluation_colors(metadata)
    grouped = {
        key: group["score"].to_numpy(dtype=float)
        for key, group in scores.groupby(["configuration", "evaluation"])
    }
    missing_halves = [
        f"{label} ({evaluation})"
        for label in order
        for evaluation in ("OOD", "IID")
        if grouped.get((label, evaluation), np.empty(0)).size == 0
    ]
    if missing_halves:
        raise ValueError(f"No scores remain to plot for: {missing_halves}")

    fig, ax = plt.subplots(figsize=figsize)
    for position, label in enumerate(order):
        for evaluation, side in (("OOD", "low"), ("IID", "high")):
            _draw_half_violin(
                ax,
                grouped[(label, evaluation)],
                position=float(position),
                side=side,
                color=colors[evaluation][label],
                width=violin_width,
                edge_color="0.3",
                inner_color="0.4",
            )

    ax.set_xlabel(x_label, fontsize=fontsize, fontweight="bold")
    ax.set_ylabel(y_label or metric.upper(), fontsize=fontsize, fontweight="bold")
    if title:
        ax.set_title(title, fontsize=fontsize + 2)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(
        order,
        rotation=x_tick_rotation,
        ha="right" if x_tick_rotation else "center",
    )
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.tick_params(axis="both", labelsize=fontsize - 2)
    if log_y:
        ax.set_yscale("log")
    if y_lim is not None:
        if log_y and y_lim[0] <= 0:
            raise ValueError("log_y=True requires a positive lower y_lim.")
        outside = scores["score"].lt(y_lim[0]) | scores["score"].gt(y_lim[1])
        if outside.any():
            warnings.warn(
                f"{int(outside.sum())} of {len(scores)} {metric!r} scores fall outside "
                f"y_lim={tuple(y_lim)} and are clipped from view.",
                RuntimeWarning,
                stacklevel=2,
            )
        ax.set_ylim(*y_lim)
    _evaluation_legend(ax)

    _finalize_figure(
        fig,
        save_dir=save_dir,
        file_name=file_name or f"{metric}_ood_vs_iid_model_comparison.png",
        high_quality=high_quality,
        show=show,
    )
    return scores.drop(columns="configuration key").reset_index(drop=True)


def plot_model_mean_std_comparison(
    df: pd.DataFrame,
    model: Any = None,
    fp_kernels: Any = None,
    count_kernels: Any = None,
    mixing_methods: Any = None,
    kernel_triples: Any = None,
    metric: str = "r2",
    include_kernel_config: bool = True,
    complete_cases: bool = True,
    figsize: tuple[float, float] = (9, 5),
    fontsize: int = 12,
    title: Optional[str] = None,
    x_label: str = "Model",
    y_label: Optional[str] = None,
    x_tick_rotation: int = 35,
    y_lim: Optional[tuple[float, float]] = None,
    gradient_bars: bool = True,
    show: bool = True,
    high_quality: bool = True,
    save_dir: Optional[Path] = DEFAULT_SAVE_DIR,
    file_name: Optional[str] = None,
) -> pd.DataFrame:
    """Plot paired OOD/IID task-level means and standard deviations.

    Scores are aggregated in two stages for each model configuration and
    evaluation split. Repeats are first averaged within each cluster, then the
    cluster means are averaged within each ``paper x target`` task. Bar heights
    are the mean of those equally weighted task means, and error bars are their
    sample standard deviation. Thus IID repeats and datasets with more clusters
    do not receive extra weight relative to their matched OOD cases.

    ``kernel_triples`` accepts a global triple/list or a dictionary mapping
    each model to separate triples. Kernel filters retain kernel-free tree,
    GNN, and MLP models. The returned DataFrame contains one OOD row and one
    IID row per plotted configuration.
    """
    if y_lim is not None:
        if len(y_lim) != 2 or not y_lim[0] < y_lim[1]:
            raise ValueError("y_lim must be a (lower, upper) pair with lower < upper.")

    scores, metadata = _select_scores(
        df,
        metric=metric,
        model=model,
        fp_kernels=fp_kernels,
        count_kernels=count_kernels,
        mixing_methods=mixing_methods,
        kernel_triples=kernel_triples,
        include_kernel_config=include_kernel_config,
        complete_cases=complete_cases,
    )
    metric = _metric_name(metric)
    metadata = _ordered_metadata(metadata, model)
    configuration_columns = [
        "configuration",
        "model",
        "fp kernel",
        "count kernel",
        "mixing method",
    ]

    case_values = (
        scores.groupby(
            ["evaluation", *CASE_COLUMNS, *configuration_columns],
            dropna=False,
            as_index=False,
        )["score"]
        .mean()
    )
    task_values = (
        case_values.groupby(
            ["evaluation", "paper", "target", *configuration_columns],
            dropna=False,
            as_index=False,
        )["score"]
        .mean()
    )
    summary_df = (
        task_values.groupby(
            ["evaluation", *configuration_columns],
            dropna=False,
            as_index=False,
        )["score"]
        .agg(mean="mean", std="std", n_tasks="count")
    )
    summary_df["std"] = summary_df["std"].fillna(0.0)

    case_counts = (
        case_values.groupby(
            ["evaluation", *configuration_columns],
            dropna=False,
        )
        .size()
        .rename("n_cases")
        .reset_index()
    )
    score_counts = (
        scores.groupby(
            ["evaluation", *configuration_columns],
            dropna=False,
        )
        .size()
        .rename("n_scores")
        .reset_index()
    )
    summary_df = summary_df.merge(
        case_counts,
        on=["evaluation", *configuration_columns],
        how="left",
        validate="one_to_one",
    ).merge(
        score_counts,
        on=["evaluation", *configuration_columns],
        how="left",
        validate="one_to_one",
    )

    order = metadata["configuration"].tolist()
    order_lookup = {label: index for index, label in enumerate(order)}
    evaluation_lookup = {
        evaluation: index
        for index, evaluation in enumerate(EVALUATIONS)
    }
    summary_df = (
        summary_df.assign(
            _configuration_order=summary_df["configuration"].map(order_lookup),
            _evaluation_order=summary_df["evaluation"].map(evaluation_lookup),
        )
        .sort_values(
            ["_configuration_order", "_evaluation_order"],
            kind="stable",
        )
        .drop(columns=["_configuration_order", "_evaluation_order"])
        .reset_index(drop=True)
    )

    colors = _evaluation_colors(metadata)
    fig, ax = plt.subplots(figsize=figsize)
    x_positions = np.arange(len(order), dtype=float)
    width = 0.33
    offsets = {"OOD": -width / 2, "IID": width / 2}
    plotted_bars = []
    plotted_values = []
    plotted_colors = []

    for evaluation in EVALUATIONS:
        evaluation_df = (
            summary_df[summary_df["evaluation"].eq(evaluation)]
            .set_index("configuration")
            .reindex(order)
        )
        if evaluation_df["mean"].isna().any():
            missing = evaluation_df.index[evaluation_df["mean"].isna()].tolist()
            raise ValueError(f"No {evaluation} summary values remain for: {missing}")

        bar_colors = [colors[evaluation][label] for label in order]
        bars = ax.bar(
            x_positions + offsets[evaluation],
            evaluation_df["mean"].to_numpy(dtype=float),
            width=width,
            yerr=evaluation_df["std"].to_numpy(dtype=float),
            capsize=4,
            color=bar_colors,
            edgecolor="none",
            linewidth=0,
            error_kw={"elinewidth": 1.2, "capthick": 1.2},
        )
        if gradient_bars:
            apply_vertical_bar_gradients(
                ax,
                bars=list(bars),
                colors=bar_colors,
                edge_linewidth=0,
            )
        plotted_bars.extend(bars)
        plotted_values.extend(evaluation_df["mean"].to_numpy(dtype=float))
        plotted_colors.extend(bar_colors)

    _hatch_highest_bars(plotted_bars, plotted_values, plotted_colors)

    ax.set_xticks(x_positions)
    ax.set_xticklabels(
        order,
        rotation=x_tick_rotation,
        ha="right" if x_tick_rotation else "center",
    )
    ax.set_xlabel(x_label, fontsize=fontsize, fontweight="bold")
    ax.set_ylabel(
        y_label or f"Mean {metric.upper()}",
        fontsize=fontsize,
        fontweight="bold",
    )
    if title:
        ax.set_title(title, fontsize=fontsize + 2)
    ax.tick_params(axis="both", labelsize=fontsize - 2)
    if y_lim is not None:
        ax.set_ylim(*y_lim)
    _evaluation_legend(ax)

    _finalize_figure(
        fig,
        save_dir=save_dir,
        file_name=(
            file_name
            or f"{metric}_ood_vs_iid_model_mean_std_comparison.png"
        ),
        high_quality=high_quality,
        show=show,
    )
    return summary_df


def _hybridization_metadata(
    metadata: pd.DataFrame,
    model: Any,
    mixing_methods: Any,
) -> pd.DataFrame:
    """Label configurations without collapsing their calculation identities."""
    metadata = metadata.copy()
    include_model = metadata["model"].nunique() > 1

    def label(row: pd.Series) -> str:
        method = str(row["mixing method"])
        value = HYBRIDIZATION_METHOD_LABELS.get(method, method)
        structure = "G" if str(row["fp kernel"]).lower() == "graph" else "FP"
        value = (
            value.replace("Av(C) × F", f"ΣC × {structure}")
            .replace("ΣF", structure)
            .replace("ΠF", structure)
        )
        if include_model:
            value = f"{_configuration_label(row, False)}\n{value}"
        return value

    metadata["hybridization method"] = metadata.apply(label, axis=1)
    # Legacy averageProduct and newer methods can share the same math label.
    # Keep their bars separate even when the compact notation is identical.
    collisions = metadata.groupby("hybridization method")["mixing method"].nunique()
    ambiguous = metadata["hybridization method"].isin(collisions.index[collisions.gt(1)])
    metadata.loc[ambiguous, "hybridization method"] += (
        "\n[" + metadata.loc[ambiguous, "mixing method"].astype(str) + "]"
    )
    selected_models = _selection_values(model) or metadata["model"].tolist()
    selected_methods = _selection_values(mixing_methods) or list(HYBRIDIZATION_METHOD_COLORS)
    model_order = {str(value).lower(): i for i, value in enumerate(selected_models)}
    method_order = {str(value).lower(): i for i, value in enumerate(selected_methods)}
    metadata["_sort_key"] = metadata.apply(
        lambda row: (
            model_order.get(str(row["model"]).lower(), len(model_order)),
            method_order.get(str(row["mixing method"]).lower(), len(method_order)),
            str(row["mixing method"]),
        ),
        axis=1,
    )
    return metadata.sort_values("_sort_key", kind="stable").drop(columns="_sort_key")


def _hybridization_input(df: pd.DataFrame) -> pd.DataFrame:
    """Exclude models with no hybridization method from the comparison."""
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")
    missing = set(ID_COLUMNS).difference(df.columns)
    if missing:
        raise ValueError(f"Input data is missing columns: {sorted(missing)}")
    return df[df["mixing method"].notna()].copy()


def _draw_hybridization_summary(
    summary: pd.DataFrame,
    metadata: pd.DataFrame,
    metric: str,
    profile: bool,
    figsize: tuple[float, float],
    fontsize: int,
    title: Optional[str],
    x_label: str,
    y_label: Optional[str],
    x_tick_rotation: int,
    y_lim: Optional[tuple[float, float]],
    gradient_bars: bool,
    show_values: bool,
    show: bool,
    high_quality: bool,
    save_dir: Optional[Path],
    file_name: str,
) -> pd.DataFrame:
    """Draw one split's profile or paired raw-score summaries."""
    if y_lim is not None and (len(y_lim) != 2 or not y_lim[0] < y_lim[1]):
        raise ValueError("y_lim must be a (lower, upper) pair with lower < upper.")
    order = metadata["hybridization method"].drop_duplicates().tolist()
    order_lookup = {label: i for i, label in enumerate(order)}
    summary = summary.assign(
        _order=summary["hybridization method"].map(order_lookup),
        _evaluation=summary["evaluation"].map({"OOD": 0, "IID": 1}),
    ).sort_values(["_order", "_evaluation"], kind="stable").drop(
        columns=["_order", "_evaluation"]
    ).reset_index(drop=True)
    palette = {
        row["hybridization method"]: HYBRIDIZATION_METHOD_COLORS.get(
            str(row["mixing method"]), "#808080"
        )
        for _, row in metadata.iterrows()
    }
    fig, ax = plt.subplots(figsize=figsize)
    positions = np.arange(len(order), dtype=float)
    evaluations = summary["evaluation"].drop_duplicates().tolist()
    for evaluation in evaluations:
        values = summary[summary["evaluation"].eq(evaluation)].set_index(
            "hybridization method"
        ).reindex(order)
        if values["mean"].isna().any():
            missing = values.index[values["mean"].isna()].tolist()
            raise ValueError(f"No usable {evaluation} summary values for: {missing}")
        colors = [palette[label] for label in order]
        if not profile and evaluation == "IID":
            colors = [_lighten_color(color) for color in colors]
        offset = 0 if profile else (-0.165 if evaluation == "OOD" else 0.165)
        bars = ax.bar(
            positions + offset,
            values["mean"].to_numpy(dtype=float),
            width=0.67 if profile else 0.33,
            yerr=values["std"].to_numpy(dtype=float),
            capsize=4,
            color=colors,
            edgecolor="none",
            error_kw={"elinewidth": 1.2, "capthick": 1.2},
        )
        if gradient_bars:
            apply_vertical_bar_gradients(ax, bars=list(bars), colors=colors, edge_linewidth=0)
        if show_values:
            heights = values["mean"].to_numpy(dtype=float)
            ax.bar_label(
                bars,
                labels=[
                    "" if np.isnan(value) else f"{value:.2f}"
                    for value in heights
                ],
                padding=3,
                fontsize=max(fontsize - 2, 6),
            )
    ax.set_xticks(positions)
    ax.set_xticklabels(order, rotation=x_tick_rotation,
                       ha="right" if x_tick_rotation else "center")
    ax.set_xlabel(x_label, fontsize=fontsize, fontweight="bold")
    ax.set_ylabel(y_label or ("Performance profile AUC" if profile else f"Mean {metric.upper()}"),
                  fontsize=fontsize, fontweight="bold")
    ax.tick_params(axis="both", labelsize=fontsize - 2)
    if title:
        ax.set_title(title, fontsize=fontsize + 2)
    if y_lim is not None:
        ax.set_ylim(*y_lim)
    if not profile:
        _evaluation_legend(ax)
    _finalize_figure(fig, save_dir, file_name, high_quality, show)
    return summary


def plot_hybridization_profile_comparison(
    df: pd.DataFrame,
    model: Any = "GPytorchMAP",
    fp_kernels: Any = None,
    count_kernels: Any = None,
    mixing_methods: Any = None,
    kernel_triples: Any = None,
    metric: str = "r2",
    ood_and_iid: str = "OOD",
    higher_is_better: Optional[bool] = None,
    complete_cases: bool = True,
    figsize: tuple[float, float] = (7, 5),
    fontsize: int = 12,
    title: Optional[str] = None,
    x_label: str = "Hybridization method",
    y_label: str = "Performance profile AUC",
    x_tick_rotation: int = 30,
    y_lim: Optional[tuple[float, float]] = (0, 1.05),
    gradient_bars: bool = False,
    show_values: bool = True,
    show: bool = True,
    high_quality: bool = False,
    save_dir: Optional[Path] = DEFAULT_SAVE_DIR,
    file_name: Optional[str] = None,
) -> pd.DataFrame:
    """Compare hybridization profile AUCs for one evaluation split.

    Pass ``ood_and_iid='OOD'`` or ``'IID'``. Each selected model/kernel/method
    configuration competes on the same cases in that split by default. AUCs
    are then summarized by model and hybridization label. Selecting one FP/count
    pair gives one AUC per method; selecting several pairs gives mean and sample
    std across their configuration AUCs. ``n_configurations`` records that count.
    Kernel-free models are excluded. Set ``show_values=False`` to omit the AUC
    annotation above each bar.
    """
    if ood_and_iid is None:
        raise ValueError("ood_and_iid must be either 'OOD' or 'IID'.")
    evaluation = _normalize_evaluations(ood_and_iid)[0]
    results = _performance_profile_results(
        _hybridization_input(df), metric=metric, model=model,
        fp_kernels=fp_kernels, count_kernels=count_kernels,
        mixing_methods=mixing_methods, kernel_triples=kernel_triples,
        higher_is_better=higher_is_better, complete_cases=complete_cases,
        ood_and_iid=evaluation,
    )
    metadata = _hybridization_metadata(results["metadata"], model, mixing_methods)
    auc_df = results["auc_df"].merge(
        metadata[["model", "fp kernel", "count kernel", "mixing method", "hybridization method"]],
        on=["model", "fp kernel", "count kernel", "mixing method"],
        validate="one_to_one",
    )
    summary = auc_df.groupby(
        ["evaluation", "model", "mixing method", "hybridization method"],
        as_index=False, dropna=False,
    )["auc"].agg(mean="mean", std="std", n_configurations="count")
    summary["std"] = summary["std"].fillna(0.0)
    summary["metric"] = _metric_name(metric)
    return _draw_hybridization_summary(
        summary, metadata, _metric_name(metric), True, figsize, fontsize, title,
        x_label, y_label, x_tick_rotation, y_lim, gradient_bars, show_values,
        show, high_quality, save_dir,
        file_name or f"{_metric_name(metric)}_{evaluation.lower()}_hybridization_profile_auc.png",
    )


def plot_hybridization_mean_std_comparison(
    df: pd.DataFrame,
    model: Any = "GPytorchMAP",
    fp_kernels: Any = None,
    count_kernels: Any = None,
    mixing_methods: Any = None,
    kernel_triples: Any = None,
    metric: str = "r2",
    complete_cases: bool = True,
    figsize: tuple[float, float] = (7, 5),
    fontsize: int = 12,
    title: Optional[str] = None,
    x_label: str = "Hybridization method",
    y_label: Optional[str] = None,
    x_tick_rotation: int = 30,
    y_lim: Optional[tuple[float, float]] = None,
    gradient_bars: bool = False,
    show: bool = True,
    high_quality: bool = False,
    save_dir: Optional[Path] = DEFAULT_SAVE_DIR,
    file_name: Optional[str] = None,
) -> pd.DataFrame:
    """Compare hybridizations with paired OOD/IID mean and sample std bars.

    Average repeats within clusters, then clusters within each paper/target and
    kernel pair. Average selected kernel pairs within each paper/target and
    hybridization before computing mean/std across equally weighted tasks.
    By default only clusters complete in both splits for every selected
    configuration are retained. Kernel-free models and annotations are omitted.
    """
    scores, metadata = _select_scores(
        _hybridization_input(df), metric=metric, model=model,
        fp_kernels=fp_kernels, count_kernels=count_kernels,
        mixing_methods=mixing_methods, kernel_triples=kernel_triples,
        include_kernel_config=True, complete_cases=complete_cases,
    )
    metadata = _hybridization_metadata(metadata, model, mixing_methods)
    scores = scores.merge(
        metadata[["configuration key", "hybridization method"]],
        on="configuration key", validate="many_to_one",
    )
    groups = ["evaluation", "model", "mixing method", "hybridization method"]
    cases = scores.groupby(
        [*groups, "paper", "target", "fp kernel", "count kernel", "cluster"],
        as_index=False, dropna=False,
    )["score"].mean()
    kernel_tasks = cases.groupby(
        [*groups, "paper", "target", "fp kernel", "count kernel"],
        as_index=False, dropna=False,
    )["score"].mean()
    tasks = kernel_tasks.groupby(
        [*groups, "paper", "target"], as_index=False, dropna=False,
    )["score"].mean()
    summary = tasks.groupby(groups, as_index=False, dropna=False)["score"].agg(
        mean="mean", std="std", n_tasks="count",
    )
    summary["std"] = summary["std"].fillna(0.0)
    summary["metric"] = _metric_name(metric)
    return _draw_hybridization_summary(
        summary, metadata, _metric_name(metric), False, figsize, fontsize, title,
        x_label, y_label, x_tick_rotation, y_lim, gradient_bars, False,
        show, high_quality, save_dir,
        file_name or f"{_metric_name(metric)}_ood_vs_iid_hybridization_mean_std.png",
    )


def load_ood_model_results(
    master_device: str = "GPU",
) -> dict[str, pd.DataFrame]:
    """Load structure- and continuous-cluster OOD masters for one device."""
    if master_device not in {"CPU", "GPU"}:
        raise ValueError("master_device must be either 'CPU' or 'GPU'.")

    return {
        ood_type: pd.read_pickle(
            DEFAULT_RESULTS_FILE
            / f"Tree_GP_MLP_{ood_type}_cluster_OOD_{master_device}.pkl"
        )
        for ood_type in ("structure", "continuous")
    }


if __name__ == "__main__":
    master_device = "GPU"
    ood_results = load_ood_model_results(master_device)

    models = [
        "RF",
        "XGBR",
        "NGB",
        "GPytorchMAP",
        "GpyroHMC",
        "MGK",
        "GNN-GP",
        "GNN-GP-SN",
        "GNN-MLP",
        "GNN-MLP-Ens",
        "MLP-GP",
        "MLP-GP-SN",
    ]
    selected_kernel_triples = {
        "MGK": ("Graph", "Matern32", "(count:+)x(graph:x)"),
        "GP-MAP": [
            ("Matern32", "Matern32", "(count:x)+(fp:x)"),
            ("TanimotoMatern32", "Matern32", "(count:x)+(fp:x)"),
        ],
        "GP-HMC": [
            ("Matern32", "Matern32", "(count:x)+(fp:x)"),
            ("TanimotoMatern32", "Matern32", "(count:x)+(fp:x)"),
        ],
    }
    for ood_type, result_df in ood_results.items():
        for ood_and_iid in ["OOD", "IID"]:
            plot_model_profile_comparison(
                df=result_df,
                model=models,
                kernel_triples=selected_kernel_triples,
                metric="cvpp",
                ood_and_iid=ood_and_iid,
                fontsize=15,
                y_label="Profile AUC: R²",
                figsize=(7.5, 5),
                save_dir=(
                    DEFAULT_SAVE_DIR / "performance_profile" / "model_comparison"
                ),
                show=False,
                high_quality=True,
                file_name=f"R2_{ood_type}_{ood_and_iid.lower()}.png",
            )

    # nll_models = [
    #     "RF",
    #     "XGBR",
    #     "NGB",
    #     "GPytorchMAP",
    #     "GNN-GP",
    #     "GNN-GP-SN",
    #     "GNN-MLP",
    #     "GNN-MLP-Ens",
    #     "MLP-GP",
    #     "MLP-GP-SN",
    # ]
    # for ood_type, result_df in ood_results.items():
        # plot_model_comparison(
        #     df=result_df,
        #     model=nll_models,
        #     kernel_triples=kernel_triples,
        #     metric="r2",
        #     y_label="R²",
        #     y_lim=(-2, 1),
        #     fontsize=15,
        #     # log_y=False,
        #     figsize=(7.5, 5),
        #     save_dir=(DEFAULT_SAVE_DIR / "score_distribution" / "model_comparison"),
        #     file_name=f"R2_{ood_type}_ood_vs_iid.png",
        # )

        # plot_model_mean_std_comparison(
        #     df=result_df,
        #     model=models,
        #     kernel_triples=kernel_triples,
        #     metric="RMSE",
        #     y_label="RMSE",
        #     # y_lim=(-1, 1),
        #     high_quality=False,
        #     gradient_bars=False,
        #     fontsize=15,
        #     figsize=(9, 5),
        #     show=False,
        #     save_dir=(DEFAULT_SAVE_DIR / "average_performance" / "model_comparison"),
        #     file_name=f"RMSE_{ood_type}_ood_vs_iid.png",
        # )

    # Hybridization comparisons: select kernel pairs rather than the model
    # comparison's kernel_triples, and compare the four requested methods.
    # for ood_type, result_df in ood_results.items():
    #     for evaluation in EVALUATIONS:
    #         plot_hybridization_profile_comparison(
    #             df=result_df,
    #             model="GPytorchMAP",
    #             fp_kernels="TanimotoMatern32",
    #             count_kernels="Matern32",
    #             mixing_methods=[
    #                 "sum",
    #                 "product",
    #                 "(count:+)x(fp:+)",
    #                 "(count:x)+(fp:x)",
    #             ],
    #             metric="r2",
    #             ood_and_iid=evaluation,
    #             fontsize=15,
    #             show_values=True,
    #             show=False,
    #             save_dir=DEFAULT_SAVE_DIR / "performance_profile" / "hybridization_comparison",
    #             file_name=f"R2_{ood_type}_{evaluation.lower()}.png",
    #         )
    #     plot_hybridization_mean_std_comparison(
    #         df=result_df,
    #         model="GPytorchMAP",
    #         fp_kernels="TanimotoMatern32",
    #         count_kernels="Matern32",
    #         mixing_methods=[
    #             # "sum",
    #             "product",
    #             "(count:+)x(fp:+)",
    #             "(count:x)+(fp:x)",
    #         ],
    #         metric="r2",
    #         y_label="Mean R²",
    #         # y_lim=(-2, 1),
    #         fontsize=15,
    #         show=False,
    #         save_dir=DEFAULT_SAVE_DIR / "average_performance" / "hybridization_comparison",
    #         file_name=f"R2_{ood_type}_ood_vs_iid.png",
    #     )
