"""Compare model performance on cluster-OOD and size-matched IID tests.

The master OOD pickle stores one scalar score per OOD cluster and five IID
random-seed scores per matching cluster. This module expands those nested
values into performance-profile cases and reports each model's profile AUC
separately for OOD and IID evaluation.
"""

from __future__ import annotations

import ast
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional
import warnings

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

try:  # Support both direct execution and package-style imports.
    from .visualization_setting import (
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
    / "Tree_and_GP_structure_cluster_OOD.pkl"
)
DEFAULT_SAVE_DIR = HERE / "ood_result_analysis"

ID_COLUMNS = ["paper", "target", "model", "fp kernel", "count kernel", "mixing method"]
CASE_COLUMNS = ["paper", "target", "cluster"]
OOD_MODEL_TYPE_COLORS = {
    "tree": "#468FC6",
    "gp": "#FF5757",
    "mgk": "#FF5757",
    "gnn": "#72B7B2",
}


def _is_missing(value: Any) -> bool:
    return value is None or (
        not isinstance(value, (list, tuple, dict)) and bool(pd.isna(value))
    )


def _tree_configuration_mask(df: pd.DataFrame) -> pd.Series:
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
    """Apply model/kernel filters while retaining kernel-free tree models."""
    filtered = _filter_selection(df, "model", model)
    tree_mask = _tree_configuration_mask(filtered)

    triples = _normalize_kernel_triples(kernel_triples)
    if triples is not None:
        triple_mask = pd.Series(False, index=filtered.index)
        for fp_kernel, count_kernel, mixing_method in triples:
            triple_mask |= (
                filtered["fp kernel"].eq(fp_kernel)
                & filtered["count kernel"].eq(count_kernel)
                & filtered["mixing method"].eq(mixing_method)
            )
        filtered = filtered[tree_mask | triple_mask].copy()

    for column, selection in [
        ("fp kernel", fp_kernels),
        ("count kernel", count_kernels),
        ("mixing method", mixing_methods),
    ]:
        selected = _selection_values(selection)
        if selected is None:
            continue
        tree_mask = _tree_configuration_mask(filtered)
        selected_lower = {str(value).lower() for value in selected}
        matches = filtered[column].astype(str).str.lower().isin(selected_lower)
        filtered = filtered[tree_mask | matches].copy()

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

        for evaluation, score_column in [("OOD", ood_column), ("IID", iid_column)]:
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


def _keep_common_complete_clusters(scores: pd.DataFrame) -> pd.DataFrame:
    """Keep base cluster cases present for every configuration in both splits."""
    number_of_configurations = scores["configuration key"].nunique()
    availability = (
        scores.groupby(["evaluation", *CASE_COLUMNS], dropna=False)["configuration key"]
        .nunique()
        .eq(number_of_configurations)
        .unstack("evaluation", fill_value=False)
    )
    for evaluation in ("OOD", "IID"):
        if evaluation not in availability:
            availability[evaluation] = False

    common_index = availability.index[availability["OOD"] & availability["IID"]]
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
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Filter, expand, and validate OOD/IID scores for the selected models.

    Returns the long-form scores and one metadata row per configuration. Only
    configurations with both OOD and IID scores are kept.
    """
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
    scores = _expand_scores(filtered, metric, include_kernel_config)
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
    for evaluation in ("OOD", "IID"):
        if evaluation not in evaluation_coverage:
            evaluation_coverage[evaluation] = 0
    eligible_configurations = set(
        evaluation_coverage.index[
            evaluation_coverage["OOD"].gt(0)
            & evaluation_coverage["IID"].gt(0)
        ]
    )
    unavailable = [
        label
        for key, label in selected_configurations.items()
        if key not in eligible_configurations
    ]
    if unavailable:
        warnings.warn(
            f"Skipping configurations without both OOD and IID "
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
            f"No selected configuration has both OOD and IID "
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
        scores = _keep_common_complete_clusters(scores)
        if scores.empty:
            raise ValueError(
                "No cluster cases are complete for every selected configuration in both "
                "OOD and IID. Select fewer models/configurations or use complete_cases=False."
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
    """Return ``{evaluation: {configuration: color}}``; IID is the lighter shade."""
    ood_colors = {
        row["configuration"]: _model_type_color(
            row["model"],
            palette=OOD_MODEL_TYPE_COLORS,
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
) -> dict[str, Any]:
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
    if len(metadata) < 2:
        raise ValueError("A performance profile requires at least two model configurations.")

    metric = _metric_name(metric)
    if higher_is_better is None:
        higher_is_better = _higher_is_better(metric)
    auc_rows: list[dict[str, Any]] = []

    for evaluation in ("OOD", "IID"):
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
) -> pd.DataFrame:
    """Return OOD and IID performance-profile AUCs for every configuration.

    OOD cases are ``paper x target x cluster``. IID cases additionally include
    the random-seed position in each cluster's score list. By default, only
    cluster cases available for every selected configuration in both OOD and
    IID are used.
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
    )["auc_df"].copy()


def plot_model_profile_comparison(
    df: pd.DataFrame,
    model: Any = None,
    fp_kernels: Any = None,
    count_kernels: Any = None,
    mixing_methods: Any = None,
    kernel_triples: Any = None,
    metric: str = "r2",
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
    """Plot OOD and matched-IID AUC beside one another for each model.

    Exact GP configurations can be selected with ``kernel_triples``. Kernel-free
    tree models are retained automatically. The returned long-form DataFrame
    contains the plotted AUC, rank, and number of contributing cases.
    """
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
    )
    auc_df = results["auc_df"]
    metadata = _ordered_metadata(results["metadata"], model)
    order = metadata["configuration"].tolist()
    values = auc_df.pivot(index="configuration", columns="evaluation", values="auc")
    values = values.reindex(order)
    colors = _evaluation_colors(metadata)

    fig, ax = plt.subplots(figsize=figsize)
    x_positions = np.arange(len(order), dtype=float)
    width = 0.33
    offsets = {"OOD": -width / 2, "IID": width / 2}
    for evaluation in ("OOD", "IID"):
        heights = values[evaluation].to_numpy(dtype=float)
        bars = ax.bar(
            x_positions + offsets[evaluation],
            heights,
            width=width,
            color=[colors[evaluation][label] for label in order],
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
    _evaluation_legend(ax)

    _finalize_figure(
        fig,
        save_dir=save_dir,
        file_name=file_name or f"{_metric_name(metric)}_ood_vs_iid_model_profile_auc.png",
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
    one point per cluster. OOD halves use the model-type color and IID halves
    a lighter shade of it. Set ``log_y=True`` for a log-scaled y-axis (all
    scores must be positive). Returns the long-form scores that were plotted.
    """
    if not 0 < violin_width <= 1:
        raise ValueError("violin_width must be in the interval (0, 1].")
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


if __name__ == "__main__":
    result_df = pd.read_pickle(DEFAULT_RESULTS_FILE)
    models = ["RF", "XGBR", "NGB", "GPytorchMAP", "MGK"]
    kernel_triples = [
        ("TanimotoMatern32", "Matern32", "product"),
        ("Graph", "Matern32", "product"),
    ]
    # profile_scores = plot_model_profile_comparison(
    #     df=result_df,
    #     model=models,
    #     kernel_triples=kernel_triples,
    #     metric="rmse",
    #     fontsize=15,
    #     y_label= "Performance profile AUC (RMSE)",
    #     figsize=(7.5, 5),
    #     save_dir=DEFAULT_SAVE_DIR/"performance_profile"/"model_comparison",
    #     file_name="RMSE_ood_vs_iid.png",
    # )
    score_distribution = plot_model_comparison(
        df=result_df,
        model=models,
        kernel_triples=kernel_triples,
        metric="nll",
        y_label= "NLL",
        y_lim=(-2, 50),
        fontsize=15,
        figsize=(7.5, 5),
        save_dir=DEFAULT_SAVE_DIR/"score_distribution"/"model_comparison",
        file_name="NLL_ood_vs_iid.png",
    )
    # print(profile_scores.to_string(index=False))
