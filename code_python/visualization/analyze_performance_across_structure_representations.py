"""Compare molecular structure representations with performance-profile AUC."""

import os
from pathlib import Path
from typing import Any, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import auc

from visualization_setting import apply_vertical_bar_gradients, ensure_long_path
from visualize_gp import (
    _add_datapoint_counts,
    _display_model_name,
    _expand_master_scores_for_profile,
    _filter_kernel_triples,
    _metric_higher_is_better,
    _selection_values,
)


HERE = Path(__file__).resolve().parent
MASTER_DATA = HERE.parent.parent / "results" / "master_performance_data"

DEFAULT_STRUCTURE_REPRESENTATION_ORDER = [
    "No Fingerprint",
    "ECFP",
    "Graph",
    "LE",
    "Mordred",
    "HDF2_256bit",
]
STRUCTURE_REPRESENTATION_COLORS = {
    "LE": "#59a89c",           # teal
    "ECFP": "#f0c571",         # gold
    "Mordred": "#e02b35",      # red
    "HDF2_256bit": "#082a54",  # dark blue
    "No Fingerprint": "#808080",  # gray
}
STRUCTURE_REPRESENTATION_DISPLAY_NAMES = {
    "HDF2_256bit": "HDF",
}


def _display_structure_representation(value: Any) -> str:
    """Return a concise representation label for plot text."""
    return STRUCTURE_REPRESENTATION_DISPLAY_NAMES.get(str(value), str(value))


def add_hdf_configuration_to_structure_results(
    structure_df: pd.DataFrame,
    hdf_df: pd.DataFrame,
    depth: int,
    bits: int,
    representation_label: Optional[str] = None,
) -> pd.DataFrame:
    """Append one HDF configuration as a structure representation."""
    required_hdf_columns = {"HDF depth", "HDF bits"}
    missing_columns = required_hdf_columns.difference(hdf_df.columns)
    if missing_columns:
        raise ValueError(
            "hdf_df is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )
    if "structure_representation" not in structure_df.columns:
        raise ValueError(
            "structure_df is missing required column: structure_representation"
        )

    label = representation_label or f"HDF{int(depth)}_{int(bits)}bit"
    hdf_depth = pd.to_numeric(hdf_df["HDF depth"], errors="coerce")
    hdf_bits = pd.to_numeric(hdf_df["HDF bits"], errors="coerce")
    selected_hdf = hdf_df.loc[
        hdf_depth.eq(int(depth)) & hdf_bits.eq(int(bits))
    ].copy()
    if selected_hdf.empty:
        raise ValueError(
            f"No HDF results were found for depth={depth}, bits={bits}."
        )
    selected_hdf["structure_representation"] = label

    # Replacing an existing row set with the same label keeps this helper
    # idempotent when it is used repeatedly in an interactive session.
    existing_labels = structure_df["structure_representation"].astype(str)
    structure_without_selected_hdf = structure_df.loc[
        ~existing_labels.str.casefold().eq(label.casefold())
    ]
    return pd.concat(
        [structure_without_selected_hdf, selected_hdf],
        ignore_index=True,
        sort=False,
    )


def add_no_fingerprint_to_structure_results(
    structure_df: pd.DataFrame,
    count_only_df: pd.DataFrame,
    representation_label: str = "No Fingerprint",
) -> pd.DataFrame:
    """Append COUNT-only results as a molecular-representation baseline."""
    column = "structure_representation"
    if column not in structure_df.columns:
        raise ValueError(f"structure_df is missing required column: {column}")
    if count_only_df.empty:
        raise ValueError("count_only_df does not contain any results.")

    selected_count_only = count_only_df.copy()
    selected_count_only[column] = representation_label

    # Replace an existing baseline with the same label so repeated notebook or
    # interactive calls cannot duplicate its rows.
    existing_labels = structure_df[column].astype(str)
    structure_without_count_only = structure_df.loc[
        ~existing_labels.str.casefold().eq(representation_label.casefold())
    ]
    return pd.concat(
        [structure_without_count_only, selected_count_only],
        ignore_index=True,
        sort=False,
    )


def _filter_structure_representations(
    df: pd.DataFrame,
    structure_representations: Any,
) -> tuple[pd.DataFrame, list[str]]:
    """Filter representations case-insensitively and return their plot order."""
    column = "structure_representation"
    if column not in df.columns:
        raise ValueError(f"df is missing required column: {column}")

    available = [
        str(value)
        for value in df[column].dropna().drop_duplicates().tolist()
    ]
    available_by_lower = {value.lower(): value for value in available}
    selected = _selection_values(structure_representations)

    if selected is None:
        ordered = [
            available_by_lower[representation.lower()]
            for representation in DEFAULT_STRUCTURE_REPRESENTATION_ORDER
            if representation.lower() in available_by_lower
        ]
        ordered.extend(
            representation
            for representation in available
            if representation.lower()
            not in {value.lower() for value in ordered}
        )
    else:
        ordered = list(
            dict.fromkeys(
                available_by_lower[str(representation).lower()]
                for representation in selected
                if str(representation).lower() in available_by_lower
            )
        )

    if not ordered:
        raise ValueError(
            "None of the requested structure representations are available. "
            f"Available values: {available}."
        )

    selected_lower = {value.lower() for value in ordered}
    mask = df[column].astype(str).str.lower().isin(selected_lower)
    return df.loc[mask].copy(), ordered


def _configuration_task_label(row: pd.Series) -> str:
    """Identify a kernel configuration without using its representation."""
    values = [
        str(row[column])
        for column in ["fp kernel", "count kernel", "mixing method"]
        if pd.notna(row[column])
    ]
    return " | ".join(values) if values else "default"


def _calculate_structure_profile_auc(
    df: pd.DataFrame,
    metric: str,
    model: Any,
    structure_representations: Any,
    fp_kernels: Any,
    count_kernels: Any,
    mixing_methods: Any,
    kernel_triples: Any,
    tree_feature_importance: str,
    dataset_as_experiment_points: bool,
    compare_models_as_methods: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Calculate representation AUCs within models or across model methods."""
    df = _filter_kernel_triples(df, kernel_triples, keep_missing=True)
    filtered_df, representation_order = _filter_structure_representations(
        df,
        structure_representations,
    )

    profile_frames = []
    for representation in representation_order:
        representation_df = filtered_df.loc[
            filtered_df["structure_representation"].astype(str).str.lower()
            == representation.lower()
        ]
        expanded = _expand_master_scores_for_profile(
            representation_df,
            metric=metric,
            model=model,
            fp_kernels=fp_kernels,
            count_kernels=count_kernels,
            mixing_methods=mixing_methods,
            tree_feature_importance=tree_feature_importance,
        )
        if expanded.empty:
            continue
        expanded["structure representation"] = representation
        profile_frames.append(expanded)

    if not profile_frames:
        raise ValueError(
            f"No {metric} values were found for model={model}, "
            f"structure_representations={structure_representations}, "
            f"fp_kernels={fp_kernels}, count_kernels={count_kernels}, "
            f"mixing_methods={mixing_methods}, and "
            f"kernel_triples={kernel_triples}."
        )

    profile_df = pd.concat(profile_frames, ignore_index=True)
    source_model_order = profile_df["model"].drop_duplicates().tolist()
    profile_df["configuration task"] = profile_df.apply(
        _configuration_task_label,
        axis=1,
    )

    task_columns = ["dataset", "target"]
    if not compare_models_as_methods:
        task_columns.extend(["model", "configuration task"])
    if not dataset_as_experiment_points:
        task_columns.append("repeat")

    # In cross-model mode, every model/configuration/representation combination
    # is one competing method. Otherwise, models and configurations define the
    # task and only representations compete within that task.
    method_columns = ["structure representation"]
    if compare_models_as_methods:
        method_columns = [
            "model",
            "configuration task",
            "structure representation",
        ]

    # Reduce repeat/fold scores or duplicate records to one score per method.
    profile_rank_df = (
        profile_df.groupby(
            task_columns + method_columns,
            dropna=False,
            as_index=False,
        )[metric]
        .mean()
    )

    if compare_models_as_methods:
        comparable_counts = profile_rank_df.groupby(
            task_columns,
            dropna=False,
        ).size()
    else:
        comparable_counts = profile_rank_df.groupby(
            task_columns,
            dropna=False,
        )["structure representation"].nunique()
    comparable_tasks = (
        comparable_counts.loc[lambda counts: counts >= 2]
        .reset_index()[task_columns]
    )
    profile_rank_df = profile_rank_df.merge(
        comparable_tasks,
        on=task_columns,
        how="inner",
    )
    if profile_rank_df.empty:
        raise ValueError(
            "No comparable tasks contain at least two selected methods."
        )

    higher_is_better = _metric_higher_is_better(metric)
    profile_rank_df["rank"] = profile_rank_df.groupby(
        task_columns,
        dropna=False,
    )[metric].rank(
        ascending=not higher_is_better,
        method="min",
    )
    maximum_rank = profile_rank_df.groupby(
        task_columns,
        dropna=False,
    )["rank"].transform("max")
    profile_rank_df["percentile rank"] = profile_rank_df["rank"].div(
        maximum_rank.replace(0, np.nan)
    )

    n_tasks = profile_rank_df[task_columns].drop_duplicates().shape[0]
    tau_grid = np.linspace(0, 1, max(n_tasks, 2))

    def profile_curve(values: pd.Series) -> np.ndarray:
        percentiles = values.dropna().to_numpy(dtype=float)
        return np.array(
            [(percentiles <= tau).mean() for tau in tau_grid],
            dtype=float,
        )

    representation_curves = {
        representation: profile_curve(group["percentile rank"])
        for representation, group in profile_rank_df.groupby(
            "structure representation",
            sort=False,
        )
    }
    representation_auc_rows = []
    for representation, curve in representation_curves.items():
        representation_rows = profile_rank_df.loc[
            profile_rank_df["structure representation"] == representation
        ]
        representation_auc_rows.append(
            {
                "metric": metric,
                "structure representation": representation,
                "auc": auc(tau_grid, curve),
                "n_profile_points": len(representation_rows),
                "n_models": representation_rows["model"].nunique(),
            }
        )
    representation_auc_df = pd.DataFrame(representation_auc_rows)
    representation_auc_df = representation_auc_df.sort_values(
        "auc",
        ascending=False,
        kind="stable",
    ).reset_index(drop=True)
    representation_auc_df.insert(
        0,
        "rank",
        np.arange(1, len(representation_auc_df) + 1),
    )

    model_auc_rows = []
    model_curves = {}
    for (model_name, representation), group in profile_rank_df.groupby(
        ["model", "structure representation"],
        sort=False,
    ):
        curve = profile_curve(group["percentile rank"])
        model_curves[(model_name, representation)] = curve
        model_auc_rows.append(
            {
                "metric": metric,
                "model": model_name,
                "structure representation": representation,
                "auc": auc(tau_grid, curve),
                "n_profile_points": len(group),
            }
        )

    return representation_auc_df, pd.DataFrame(model_auc_rows), {
        "tau_grid": tau_grid,
        "curves": representation_curves,
        "model_curves": model_curves,
        "representation_order": representation_order,
        "model_order": source_model_order,
    }


def _plot_grouped_profiles(
    auc_df: pd.DataFrame,
    profile_results: dict[str, Any],
    representation_order: list[str],
    options: dict[str, Any],
) -> None:
    """Optionally plot the generalized profile for each representation."""
    allowed_options = {
        "figsize",
        "fontsize",
        "title",
        "legend_ncols",
        "show",
        "high_quality",
        "save_dir",
        "file_name",
    }
    unknown_options = set(options).difference(allowed_options)
    if unknown_options:
        raise ValueError(
            "Unsupported profile_options: "
            + ", ".join(sorted(unknown_options))
        )

    figsize = options.get("figsize", (7, 5))
    fontsize = options.get("fontsize", 12)
    show = options.get("show", True)
    high_quality = options.get("high_quality", True)
    save_dir = options.get("save_dir", HERE / "result_analysis")
    file_name = options.get(
        "file_name",
        "structure_representation_performance_profile.png",
    )

    fig, ax = plt.subplots(figsize=figsize)
    tau_grid = profile_results["tau_grid"]
    curves = profile_results["curves"]
    for representation in representation_order:
        if representation not in curves:
            continue
        representation_auc = float(
            auc_df.loc[
                auc_df["structure representation"] == representation,
                "auc",
            ].iloc[0]
        )
        ax.plot(
            tau_grid,
            curves[representation],
            color=STRUCTURE_REPRESENTATION_COLORS.get(representation),
            label=(
                f"{_display_structure_representation(representation)}, "
                f"{representation_auc:.2f}"
            ),
        )

    legend_ncols = options.get("legend_ncols")
    if legend_ncols is None:
        legend_ncols = max(1, min(3, len(representation_order)))
    elif legend_ncols < 1:
        raise ValueError("profile_options['legend_ncols'] must be at least 1.")
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=legend_ncols,
        frameon=False,
        fontsize=fontsize,
    )
    ax.set_xlabel(r"$\tau$", fontsize=fontsize, fontweight="bold")
    ax.set_ylabel(r"$p_i(\tau)$", fontsize=fontsize, fontweight="bold")
    if options.get("title") is not None:
        ax.set_title(options["title"], fontsize=fontsize + 2)
    ax.tick_params(axis="both", labelsize=fontsize - 2)
    fig.tight_layout()

    if save_dir is not None:
        save_path = ensure_long_path(Path(save_dir))
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(
            ensure_long_path(save_path / file_name),
            bbox_inches="tight",
            format="png",
            dpi=900 if high_quality else 100,
        )
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_structure_representation_profile_comparison(
    df: pd.DataFrame,
    model: Any = None,
    structure_representations: Any = None,
    fp_kernels: Any = None,
    count_kernels: Any = None,
    mixing_methods: Any = None,
    kernel_triples: Any = None,
    metric: str = "r2",
    auc_column: str = "auc",
    tree_feature_importance: str = "MDI",
    dataset_as_experiment_points: bool = False,
    compare_models_as_methods: bool = False,
    plot_profile: bool = False,
    profile_options: Optional[dict[str, Any]] = None,
    figsize: tuple[float, float] = (7, 5),
    fontsize: int = 12,
    title: Optional[str] = None,
    x_label: str = "Structure representation",
    y_label: str = "Performance profile AUC",
    x_tick_rotation: int = 20,
    y_lim: tuple[float, float] = (0, 1.05),
    show_values: bool = True,
    gradient_bars: bool = True,
    show: bool = True,
    high_quality: bool = True,
    save_dir: Optional[Path] = HERE / "result_analysis",
    file_name: Optional[str] = None,
) -> pd.DataFrame:
    """Compare representations using generalized performance-profile AUC.

    By default, models are part of the task definition and representations are
    ranked within the same model. Set ``compare_models_as_methods=True`` to
    treat every selected model/kernel/representation combination as a distinct
    competing method on the same dataset-target-repeat tasks. This allows a
    selected MGK Graph configuration to compete with tree-model methods.

    Only tasks containing at least two selected representations contribute to
    the profile. Set ``dataset_as_experiment_points=True`` to average
    repeat/fold scores within each dataset-target-model task before ranking.
    Set ``plot_profile=True`` to additionally draw the profile curves.
    """
    calculated_auc_df, _, profile_results = _calculate_structure_profile_auc(
        df=df,
        metric=metric,
        model=model,
        structure_representations=structure_representations,
        fp_kernels=fp_kernels,
        count_kernels=count_kernels,
        mixing_methods=mixing_methods,
        kernel_triples=kernel_triples,
        tree_feature_importance=tree_feature_importance,
        dataset_as_experiment_points=dataset_as_experiment_points,
        compare_models_as_methods=compare_models_as_methods,
    )
    calculated_auc_df[auc_column] = calculated_auc_df["auc"]
    representation_order = [
        representation
        for representation in profile_results["representation_order"]
        if representation
        in set(calculated_auc_df["structure representation"])
    ]
    order_lookup = {
        representation: index
        for index, representation in enumerate(representation_order)
    }
    calculated_auc_df["_plot_order"] = calculated_auc_df[
        "structure representation"
    ].map(order_lookup)
    calculated_auc_df = (
        calculated_auc_df.sort_values("_plot_order", kind="stable")
        .drop(columns="_plot_order")
        .copy()
    )
    bar_palette = {
        representation: STRUCTURE_REPRESENTATION_COLORS.get(
            representation,
            sns.color_palette("Set2", n_colors=len(representation_order))[
                index
            ],
        )
        for index, representation in enumerate(representation_order)
    }

    if plot_profile:
        _plot_grouped_profiles(
            calculated_auc_df,
            profile_results,
            representation_order,
            {
                "show": show,
                "high_quality": high_quality,
                "save_dir": save_dir,
                **dict(profile_options or {}),
            },
        )

    fig, ax = plt.subplots(figsize=figsize)
    sns.barplot(
        data=calculated_auc_df,
        x="structure representation",
        y=auc_column,
        hue="structure representation",
        order=representation_order,
        hue_order=representation_order,
        palette=bar_palette,
        width=0.67,
        dodge=False,
        legend=False,
        errorbar=None,
        ax=ax,
    )
    if gradient_bars:
        apply_vertical_bar_gradients(
            ax,
            list(ax.patches[: len(representation_order)]),
            [bar_palette[value] for value in representation_order],
            gradient_fraction=0.35,
            edge_linewidth=0,
        )

    auc_by_representation = calculated_auc_df.set_index(
        "structure representation"
    )[auc_column]
    max_label_y = None
    for patch, representation in zip(ax.patches, representation_order):
        auc_value = float(auc_by_representation.at[representation])
        x_position = patch.get_x() + patch.get_width() / 2
        label_y = auc_value + 0.015
        max_label_y = label_y if max_label_y is None else max(
            max_label_y,
            label_y,
        )
        if show_values:
            ax.text(
                x_position,
                label_y,
                f"{auc_value:.2f}",
                ha="center",
                va="bottom",
                fontsize=fontsize - 2,
            )

    ax.set_xlabel(x_label, fontsize=fontsize, fontweight="bold")
    ax.set_ylabel(y_label, fontsize=fontsize, fontweight="bold")
    if title is not None:
        ax.set_title(title, fontsize=fontsize + 2)
    ax.tick_params(axis="both", labelsize=fontsize - 2)
    ax.set_xticks(range(len(representation_order)))
    ax.set_xticklabels(
        [
            _display_structure_representation(representation)
            for representation in representation_order
        ],
        rotation=x_tick_rotation,
        ha="right" if x_tick_rotation else "center",
    )
    ax.set_ylim(*y_lim)
    if max_label_y is not None:
        bottom, top = ax.get_ylim()
        ax.set_ylim(bottom, max(top, max_label_y + 0.08))
    fig.tight_layout()

    if save_dir is not None:
        save_path = ensure_long_path(Path(save_dir))
        os.makedirs(save_path, exist_ok=True)
        if file_name is None:
            selected_models = _selection_values(model)
            model_name = "_".join(
                str(value) for value in (selected_models or ["all"])
            )
            file_name = (
                f"{model_name}_{metric}_structure_representation_"
                "profile_auc_barplot.png"
            )
        fig.savefig(
            ensure_long_path(save_path / file_name),
            bbox_inches="tight",
            format="png",
            dpi=900 if high_quality else 100,
        )
    if show:
        plt.show()
    else:
        plt.close(fig)
    return calculated_auc_df


def plot_structure_representation_mean_std_comparison(
    df: pd.DataFrame,
    model: Any = None,
    structure_representations: Any = None,
    fp_kernels: Any = None,
    count_kernels: Any = None,
    mixing_methods: Any = None,
    kernel_triples: Any = None,
    metric: str = "OOF_R2",
    figsize: tuple[float, float] = (7, 5),
    fontsize: int = 12,
    title: Optional[str] = None,
    x_label: str = "Structure representation",
    y_label: Optional[str] = None,
    x_tick_rotation: int = 0,
    y_lim: Optional[tuple[float, float]] = None,
    show_values: bool = True,
    gradient_bars: bool = True,
    show: bool = True,
    high_quality: bool = True,
    save_dir: Optional[Path] = HERE / "result_analysis",
    file_name: Optional[str] = None,
) -> pd.DataFrame:
    """Plot each representation's task-level metric mean and standard deviation.

    A task is one dataset, target, model, kernel configuration, and repeat/fold.
    Thus, a selected MGK kernel triple contributes as one model method, while
    each selected tree model contributes separately for every representation it
    supports. Error bars show the sample standard deviation across those tasks.
    """
    df = _filter_kernel_triples(df, kernel_triples, keep_missing=True)
    filtered_df, representation_order = _filter_structure_representations(
        df,
        structure_representations,
    )

    task_frames = []
    for representation in representation_order:
        representation_df = filtered_df.loc[
            filtered_df["structure_representation"].astype(str).str.lower()
            == representation.lower()
        ]
        expanded = _expand_master_scores_for_profile(
            representation_df,
            metric=metric,
            model=model,
            fp_kernels=fp_kernels,
            count_kernels=count_kernels,
            mixing_methods=mixing_methods,
        )
        if expanded.empty:
            continue
        expanded["structure representation"] = representation
        expanded["configuration task"] = expanded.apply(
            _configuration_task_label,
            axis=1,
        )
        task_frames.append(expanded)

    if not task_frames:
        raise ValueError(
            f"No {metric} values were found for the requested models, "
            "representations, and kernel configurations."
        )

    task_values = pd.concat(task_frames, ignore_index=True)
    task_columns = [
        "dataset",
        "target",
        "model",
        "configuration task",
        "repeat",
        "structure representation",
    ]
    task_values = (
        task_values.groupby(task_columns, dropna=False, as_index=False)[metric]
        .mean()
    )
    summary = (
        task_values.groupby("structure representation", as_index=False)[metric]
        .agg(mean="mean", std="std", n_tasks="count")
    )
    summary["std"] = summary["std"].fillna(0.0)
    order_lookup = {
        representation: position
        for position, representation in enumerate(representation_order)
    }
    summary["_plot_order"] = summary["structure representation"].map(
        order_lookup
    )
    summary = (
        summary.sort_values("_plot_order", kind="stable")
        .drop(columns="_plot_order")
        .reset_index(drop=True)
    )

    plotted_representations = summary["structure representation"].tolist()
    colors = [
        STRUCTURE_REPRESENTATION_COLORS.get(
            representation,
            sns.color_palette("Set2", n_colors=len(summary))[position],
        )
        for position, representation in enumerate(plotted_representations)
    ]
    x_positions = np.arange(len(summary))
    fig, ax = plt.subplots(figsize=figsize)
    bars = ax.bar(
        x_positions,
        summary["mean"].to_numpy(dtype=float),
        yerr=summary["std"].to_numpy(dtype=float),
        capsize=5,
        color=colors,
        edgecolor="none",
        linewidth=0,
        error_kw={"elinewidth": 1.2, "capthick": 1.2},
    )
    if gradient_bars:
        apply_vertical_bar_gradients(
            ax,
            list(bars),
            colors,
            gradient_fraction=0.35,
            edge_linewidth=0,
        )

    ax.set_xticks(x_positions)
    ax.set_xticklabels(
        [
            _display_structure_representation(representation)
            for representation in plotted_representations
        ],
        rotation=x_tick_rotation,
        ha="right" if x_tick_rotation else "center",
    )
    ax.set_xlabel(x_label, fontsize=fontsize, fontweight="bold")
    ax.set_ylabel(
        y_label or f"Mean {metric}",
        fontsize=fontsize,
        fontweight="bold",
    )
    if title is not None:
        ax.set_title(title, fontsize=fontsize + 2)
    ax.tick_params(axis="both", labelsize=fontsize - 2)
    if y_lim is not None:
        ax.set_ylim(*y_lim)

    if show_values:
        lower, upper = ax.get_ylim()
        label_offset = 0.02 * (upper - lower)
        for bar, row in zip(bars, summary.itertuples(index=False)):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                row.mean + row.std + label_offset,
                f"{row.mean:.2f}",
                ha="center",
                va="bottom",
                fontsize=fontsize - 2,
            )
        _, updated_upper = ax.get_ylim()
        required_upper = max(
            summary["mean"] + summary["std"] + 3 * label_offset
        )
        ax.set_ylim(lower, max(updated_upper, required_upper))

    fig.tight_layout()
    if save_dir is not None:
        save_path = ensure_long_path(Path(save_dir))
        os.makedirs(save_path, exist_ok=True)
        if file_name is None:
            file_name = (
                f"{metric}_structure_representation_mean_std_barplot.png"
            )
        fig.savefig(
            ensure_long_path(save_path / file_name),
            bbox_inches="tight",
            format="png",
            dpi=900 if high_quality else 100,
        )
    if show:
        plt.show()
    else:
        plt.close(fig)
    return summary


def plot_structure_representation_performance_by_target_grid(
    df: pd.DataFrame,
    model: Any = None,
    structure_representations: Any = None,
    fp_kernels: Any = None,
    count_kernels: Any = None,
    mixing_methods: Any = None,
    kernel_triples: Any = None,
    metric: str = "OOF_R2",
    figsize: tuple[float, float] = (7, 5),
    ncols: int = 3,
    fontsize: int = 12,
    x_label: str = "Structure representation",
    y_label: Optional[str] = None,
    x_tick_rotation: int = 20,
    y_lim: Optional[tuple[float, float]] = None,
    show_values: bool = True,
    gradient_bars: bool = False,
    show: bool = False,
    high_quality: bool = False,
    save_dir: Optional[Path] = HERE / "result_analysis",
    file_name: Optional[str] = None,
) -> pd.DataFrame:
    """Plot representation performance for every dataset-target in one grid.

    Each subplot contains one bar per selected molecular representation and is
    ordered by the full dataset-target size. Scores are first reduced to one
    value per model, kernel configuration, and repeat/fold. Bars then show the
    mean +/- sample standard deviation across those task values, matching the
    aggregation used by :func:`plot_structure_representation_mean_std_comparison`.

    ``figsize`` specifies the size of one subplot, rather than the size of the
    complete figure. The returned table contains the plotted summary values,
    datapoint counts, and the shared output path.
    """
    if ncols < 1:
        raise ValueError("ncols must be at least 1.")
    if y_lim is not None:
        if (
            len(y_lim) != 2
            or not np.isfinite(y_lim).all()
            or not y_lim[0] < y_lim[1]
        ):
            raise ValueError(
                "y_lim must contain two finite values with lower < upper."
            )

    df = _filter_kernel_triples(df, kernel_triples, keep_missing=True)
    filtered_df, representation_order = _filter_structure_representations(
        df,
        structure_representations,
    )

    task_frames = []
    for representation in representation_order:
        representation_df = filtered_df.loc[
            filtered_df["structure_representation"].astype(str).str.lower()
            == representation.lower()
        ]
        expanded = _expand_master_scores_for_profile(
            representation_df,
            metric=metric,
            model=model,
            fp_kernels=fp_kernels,
            count_kernels=count_kernels,
            mixing_methods=mixing_methods,
        )
        if expanded.empty:
            continue
        expanded["structure representation"] = representation
        expanded["configuration task"] = expanded.apply(
            _configuration_task_label,
            axis=1,
        )
        task_frames.append(expanded)

    if not task_frames:
        raise ValueError(
            f"No {metric} values were found for the requested models, "
            "representations, and kernel configurations."
        )

    task_values = pd.concat(task_frames, ignore_index=True)
    task_columns = [
        "dataset",
        "target",
        "model",
        "configuration task",
        "repeat",
        "structure representation",
    ]
    task_values = task_values.groupby(
        task_columns,
        dropna=False,
        as_index=False,
    )[metric].mean()
    summary = task_values.groupby(
        ["dataset", "target", "structure representation"],
        as_index=False,
        dropna=False,
    )[metric].agg(mean="mean", std="std", n_tasks="count")
    summary["std"] = summary["std"].fillna(0.0)
    summary["metric"] = metric
    summary = _add_datapoint_counts(summary)

    representation_rank = {
        representation: position
        for position, representation in enumerate(representation_order)
    }
    summary["_representation order"] = summary[
        "structure representation"
    ].map(representation_rank)
    summary = summary.sort_values(
        [
            "n datapoints",
            "dataset",
            "target",
            "_representation order",
        ],
        kind="stable",
    ).reset_index(drop=True)
    summary["output path"] = None

    target_groups = list(
        summary.groupby(["dataset", "target"], sort=False, dropna=False)
    )
    nrows = int(np.ceil(len(target_groups) / ncols))
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(figsize[0] * ncols, figsize[1] * nrows),
        squeeze=False,
    )
    axes_flat = axes.ravel()
    fallback_colors = sns.color_palette(
        "Set2",
        n_colors=len(representation_order),
    )
    representation_colors = {
        representation: STRUCTURE_REPRESENTATION_COLORS.get(
            representation,
            fallback_colors[position],
        )
        for position, representation in enumerate(representation_order)
    }

    try:
        for ax, ((dataset, target), target_df) in zip(
            axes_flat,
            target_groups,
        ):
            target_df = target_df.sort_values(
                "_representation order",
                kind="stable",
            )
            plotted_representations = target_df[
                "structure representation"
            ].tolist()
            positions = np.arange(len(plotted_representations), dtype=float)
            colors = [
                representation_colors[representation]
                for representation in plotted_representations
            ]
            bars = ax.bar(
                positions,
                target_df["mean"].to_numpy(dtype=float),
                width=0.67,
                yerr=target_df["std"].to_numpy(dtype=float),
                capsize=3,
                color=colors,
                edgecolor="none",
                error_kw={"elinewidth": 1.0, "capthick": 1.0},
            )
            if gradient_bars:
                apply_vertical_bar_gradients(
                    ax,
                    list(bars),
                    colors,
                    gradient_fraction=0.35,
                    edge_linewidth=0,
                )
            if show_values:
                means = target_df["mean"].to_numpy(dtype=float)
                ax.bar_label(
                    bars,
                    labels=[
                        f"{value:.2f}" if value >= 0 else ""
                        for value in means
                    ],
                    padding=2,
                    fontsize=max(fontsize - 3, 8),
                    clip_on=True,
                )
                for bar, value in zip(bars, means):
                    if value >= 0:
                        continue
                    # Keep negative values readable when the displayed axis
                    # begins at zero by placing their labels above the baseline.
                    ax.annotate(
                        f"{value:.2f}",
                        xy=(bar.get_x() + bar.get_width() / 2, 0.0),
                        xytext=(0, 2),
                        textcoords="offset points",
                        ha="center",
                        va="bottom",
                        fontsize=max(fontsize - 3, 8),
                        clip_on=True,
                    )

            n_datapoints = float(target_df["n datapoints"].iloc[0])
            n_label = (
                f"{int(n_datapoints)}"
                if n_datapoints.is_integer()
                else f"{n_datapoints:g}"
            )
            target_label = str(target).removeprefix("target_")
            ax.set_title(
                f"{target_label} (n={n_label})",
                fontsize=fontsize + 1,
            )
            ax.set_xticks(positions)
            ax.set_xticklabels(
                [
                    _display_structure_representation(representation)
                    for representation in plotted_representations
                ],
                rotation=x_tick_rotation,
                ha="right" if x_tick_rotation else "center",
            )
            ax.set_xlabel(x_label, fontsize=fontsize)
            ax.set_ylabel(y_label or f"Mean {metric}", fontsize=fontsize)
            ax.tick_params(axis="both", labelsize=fontsize - 2)
            if y_lim is not None:
                ax.set_ylim(*y_lim)

        for ax in axes_flat[len(target_groups):]:
            ax.set_visible(False)

        fig.tight_layout()
        if save_dir is not None:
            output_dir = ensure_long_path(Path(save_dir))
            output_dir.mkdir(parents=True, exist_ok=True)
            output_name = file_name or (
                f"{metric}_structure_representation_by_target.png"
            )
            output_path = ensure_long_path(output_dir / output_name)
            fig.savefig(
                output_path,
                bbox_inches="tight",
                format="png",
                dpi=900 if high_quality else 150,
            )
            summary["output path"] = str(output_path)
        if show:
            plt.show()
    finally:
        plt.close(fig)

    return summary.drop(columns="_representation order")


def plot_structure_representation_profile_heatmap(
    df: pd.DataFrame,
    model: Any = None,
    structure_representations: Any = None,
    fp_kernels: Any = None,
    count_kernels: Any = None,
    mixing_methods: Any = None,
    kernel_triples: Any = None,
    metric: str = "OOF_R2",
    tree_feature_importance: str = "MDI",
    dataset_as_experiment_points: bool = False,
    compare_models_as_methods: bool = False,
    figsize: tuple[float, float] = (7, 5),
    fontsize: int = 12,
    title: Optional[str] = None,
    x_label: str = "Structure representation",
    y_label: str = "Model",
    cbar_label: str = "Performance profile AUC",
    cmap: str = "viridis",
    vmin: float = 0,
    vmax: float = 1,
    annotate: bool = True,
    show: bool = True,
    high_quality: bool = True,
    plot_heatmap: bool = True,
    save_dir: Optional[Path] = HERE / "result_analysis",
    file_name: Optional[str] = None,
) -> pd.DataFrame:
    """Plot model-by-representation performance-profile AUCs as a heatmap.

    By default, each cell uses percentile ranks obtained by comparing
    representations within the same model and kernel configuration. Set
    ``compare_models_as_methods=True`` to rank selected model/configuration/
    representation methods against one another on each dataset task. Missing
    combinations are left blank. Set ``plot_heatmap=False`` to return the
    matrix without creating a figure.
    """
    _, model_auc_df, profile_results = _calculate_structure_profile_auc(
        df=df,
        metric=metric,
        model=model,
        structure_representations=structure_representations,
        fp_kernels=fp_kernels,
        count_kernels=count_kernels,
        mixing_methods=mixing_methods,
        kernel_triples=kernel_triples,
        tree_feature_importance=tree_feature_importance,
        dataset_as_experiment_points=dataset_as_experiment_points,
        compare_models_as_methods=compare_models_as_methods,
    )

    heatmap_matrix = model_auc_df.pivot_table(
        index="model",
        columns="structure representation",
        values="auc",
        aggfunc="mean",
    )
    selected_models = _selection_values(model)
    source_model_order = profile_results["model_order"]
    if selected_models is None:
        model_order = source_model_order
    else:
        source_by_lower = {
            str(value).lower(): value for value in source_model_order
        }
        model_order = [
            source_by_lower[str(value).lower()]
            for value in selected_models
            if str(value).lower() in source_by_lower
        ]
    representation_order = profile_results["representation_order"]
    heatmap_matrix = heatmap_matrix.reindex(
        index=model_order,
        columns=representation_order,
    )

    if not plot_heatmap:
        return heatmap_matrix

    annotations = heatmap_matrix.copy().astype(object)
    for model_name in heatmap_matrix.index:
        for representation in heatmap_matrix.columns:
            value = heatmap_matrix.loc[model_name, representation]
            annotations.loc[model_name, representation] = (
                "" if pd.isna(value) else f"{value:.2f}"
            )

    custom_cmap = sns.color_palette(cmap, as_cmap=True)
    if hasattr(custom_cmap, "copy"):
        custom_cmap = custom_cmap.copy()
    custom_cmap.set_bad(color="lightgray")

    fig, ax = plt.subplots(figsize=figsize)
    heatmap = sns.heatmap(
        heatmap_matrix,
        annot=annotations if annotate else False,
        fmt="",
        cmap=custom_cmap,
        cbar=True,
        vmin=vmin,
        vmax=vmax,
        mask=heatmap_matrix.isna(),
        linewidths=0.5,
        linecolor="white",
        annot_kws={"fontsize": fontsize},
        ax=ax,
    )
    colorbar = heatmap.collections[0].colorbar
    if colorbar is not None:
        colorbar.set_label(
            cbar_label,
            fontsize=fontsize,
            fontweight="bold",
        )
        colorbar.ax.tick_params(labelsize=fontsize - 1)

    ax.set_xlabel(x_label, fontsize=fontsize, fontweight="bold")
    ax.set_ylabel(y_label, fontsize=fontsize, fontweight="bold")
    if title is not None:
        ax.set_title(title, fontsize=fontsize + 2)
    ax.set_xticklabels(
        heatmap_matrix.columns,
        rotation=0,
        fontsize=fontsize,
    )
    ax.set_yticklabels(
        [_display_model_name(value) for value in heatmap_matrix.index],
        rotation=0,
        fontsize=fontsize,
    )
    fig.tight_layout()

    if save_dir is not None:
        save_path = ensure_long_path(Path(save_dir))
        os.makedirs(save_path, exist_ok=True)
        if file_name is None:
            file_name = (
                f"{metric}_model_structure_representation_"
                "profile_auc_heatmap.png"
            )
        fig.savefig(
            ensure_long_path(save_path / file_name),
            bbox_inches="tight",
            format="png",
            dpi=900 if high_quality else 100,
        )
    if show:
        plt.show()
    else:
        plt.close(fig)

    return heatmap_matrix


def _expand_hdf_depth_bits_scores(
    df: pd.DataFrame,
    metric: str,
    model: Any,
    hdf_depths: Any,
    hdf_bits: Any,
    tree_feature_importance: str,
) -> tuple[pd.DataFrame, list[int], list[int]]:
    """Expand selected HDF results to task-level metric values."""
    required_columns = {"HDF depth", "HDF bits"}
    missing_columns = required_columns.difference(df.columns)
    if missing_columns:
        raise ValueError(
            "df is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    hdf_df = df.copy()
    hdf_df["HDF depth"] = pd.to_numeric(
        hdf_df["HDF depth"], errors="coerce"
    )
    hdf_df["HDF bits"] = pd.to_numeric(
        hdf_df["HDF bits"], errors="coerce"
    )

    available_depths = sorted(
        int(value) for value in hdf_df["HDF depth"].dropna().unique()
    )
    available_bits = sorted(
        int(value) for value in hdf_df["HDF bits"].dropna().unique()
    )
    selected_depths = _selection_values(hdf_depths)
    selected_bits = _selection_values(hdf_bits)
    depth_order = (
        available_depths
        if selected_depths is None
        else [
            int(value)
            for value in selected_depths
            if int(value) in available_depths
        ]
    )
    bit_order = (
        available_bits
        if selected_bits is None
        else [
            int(value)
            for value in selected_bits
            if int(value) in available_bits
        ]
    )
    depth_order = list(dict.fromkeys(depth_order))
    bit_order = list(dict.fromkeys(bit_order))
    if not depth_order or not bit_order:
        raise ValueError(
            "None of the requested HDF depths/bits are available. "
            f"Available depths: {available_depths}; "
            f"available bits: {available_bits}."
        )

    profile_frames = []
    for depth in depth_order:
        for bits in bit_order:
            combination_df = hdf_df.loc[
                hdf_df["HDF depth"].eq(depth)
                & hdf_df["HDF bits"].eq(bits)
            ]
            expanded = _expand_master_scores_for_profile(
                combination_df,
                metric=metric,
                model=model,
                tree_feature_importance=tree_feature_importance,
            )
            if expanded.empty:
                continue
            expanded["HDF depth"] = depth
            expanded["HDF bits"] = bits
            profile_frames.append(expanded)

    if not profile_frames:
        raise ValueError(
            f"No {metric} values were found for model={model}, "
            f"HDF depths={hdf_depths}, and HDF bits={hdf_bits}."
        )

    return pd.concat(profile_frames, ignore_index=True), depth_order, bit_order


def _hdf_task_metric_values(
    profile_df: pd.DataFrame,
    metric: str,
    dataset_as_experiment_points: bool,
) -> tuple[pd.DataFrame, list[str], list[str]]:
    """Reduce duplicate records to one metric value per HDF task/method."""
    task_columns = ["dataset", "target", "model"]
    if not dataset_as_experiment_points:
        task_columns.append("repeat")
    method_columns = ["HDF depth", "HDF bits"]

    # OOF metrics contribute one value per dataset-target-model task. For
    # seed/fold metrics, repeat remains part of the task unless explicitly
    # averaged with dataset_as_experiment_points=True.
    task_values = (
        profile_df.groupby(
            task_columns + method_columns,
            dropna=False,
            as_index=False,
        )[metric]
        .mean()
    )
    return task_values, task_columns, method_columns


def _calculate_hdf_depth_bits_profile_auc(
    df: pd.DataFrame,
    metric: str,
    model: Any,
    hdf_depths: Any,
    hdf_bits: Any,
    tree_feature_importance: str,
    dataset_as_experiment_points: bool,
) -> tuple[pd.DataFrame, list[int], list[int]]:
    """Calculate profile AUC for each HDF depth/bit combination."""
    profile_df, depth_order, bit_order = _expand_hdf_depth_bits_scores(
        df=df,
        metric=metric,
        model=model,
        hdf_depths=hdf_depths,
        hdf_bits=hdf_bits,
        tree_feature_importance=tree_feature_importance,
    )
    profile_rank_df, task_columns, method_columns = _hdf_task_metric_values(
        profile_df,
        metric,
        dataset_as_experiment_points,
    )
    comparable_counts = profile_rank_df.groupby(
        task_columns,
        dropna=False,
    ).size()
    comparable_tasks = (
        comparable_counts.loc[lambda counts: counts >= 2]
        .reset_index()[task_columns]
    )
    profile_rank_df = profile_rank_df.merge(
        comparable_tasks,
        on=task_columns,
        how="inner",
    )
    if profile_rank_df.empty:
        raise ValueError(
            "No comparable tasks contain at least two HDF configurations."
        )

    higher_is_better = _metric_higher_is_better(metric)
    profile_rank_df["rank"] = profile_rank_df.groupby(
        task_columns,
        dropna=False,
    )[metric].rank(
        ascending=not higher_is_better,
        method="min",
    )
    maximum_rank = profile_rank_df.groupby(
        task_columns,
        dropna=False,
    )["rank"].transform("max")
    profile_rank_df["percentile rank"] = profile_rank_df["rank"].div(
        maximum_rank.replace(0, np.nan)
    )

    n_tasks = profile_rank_df[task_columns].drop_duplicates().shape[0]
    tau_grid = np.linspace(0, 1, max(n_tasks, 2))
    auc_rows = []
    for (depth, bits), group in profile_rank_df.groupby(
        method_columns,
        sort=False,
    ):
        percentiles = group["percentile rank"].dropna().to_numpy(dtype=float)
        curve = np.array(
            [(percentiles <= tau).mean() for tau in tau_grid],
            dtype=float,
        )
        auc_rows.append(
            {
                "metric": metric,
                "HDF depth": int(depth),
                "HDF bits": int(bits),
                "auc": auc(tau_grid, curve),
                "n_profile_points": len(percentiles),
                "n_models": group["model"].nunique(),
            }
        )

    return pd.DataFrame(auc_rows), depth_order, bit_order


def plot_hdf_depth_bits_profile_heatmap(
    df: pd.DataFrame,
    model: Any = None,
    hdf_depths: Any = None,
    hdf_bits: Any = None,
    metric: str = "OOF_R2",
    tree_feature_importance: str = "MDI",
    dataset_as_experiment_points: bool = False,
    figsize: tuple[float, float] = (7, 5),
    fontsize: int = 12,
    title: Optional[str] = None,
    x_label: str = "HDF bits",
    y_label: str = "HDF depth",
    cbar_label: str = "Performance profile AUC",
    cmap: str = "viridis",
    vmin: float = 0,
    vmax: float = 1,
    annotate: bool = True,
    show: bool = True,
    high_quality: bool = True,
    plot_heatmap: bool = True,
    save_dir: Optional[Path] = HERE / "result_analysis",
    file_name: Optional[str] = None,
) -> pd.DataFrame:
    """Plot HDF depth-by-bits performance-profile AUCs as a heatmap.

    HDF depth/bit combinations are the competing methods. By default, each
    dataset-target-model-repeat is a task, so configurations are compared
    within the same tree model. OOF metrics contain one pooled value per
    dataset-target-model and therefore effectively have no repeated folds.
    Non-OOF metrics use each stored seed/fold score as a separate task. Set
    ``dataset_as_experiment_points=True`` to average those repeated scores
    before ranking the HDF configurations.

    Missing depth/bit combinations are left blank. Set ``plot_heatmap=False``
    to calculate and return the matrix without creating a figure.
    """
    auc_df, depth_order, bit_order = _calculate_hdf_depth_bits_profile_auc(
        df=df,
        metric=metric,
        model=model,
        hdf_depths=hdf_depths,
        hdf_bits=hdf_bits,
        tree_feature_importance=tree_feature_importance,
        dataset_as_experiment_points=dataset_as_experiment_points,
    )
    heatmap_matrix = auc_df.pivot_table(
        index="HDF depth",
        columns="HDF bits",
        values="auc",
        aggfunc="mean",
    ).reindex(index=depth_order, columns=bit_order)

    if not plot_heatmap:
        return heatmap_matrix

    annotations = heatmap_matrix.copy().astype(object)
    for depth in heatmap_matrix.index:
        for bits in heatmap_matrix.columns:
            value = heatmap_matrix.loc[depth, bits]
            annotations.loc[depth, bits] = (
                "" if pd.isna(value) else f"{value:.2f}"
            )
    custom_cmap = sns.color_palette(cmap, as_cmap=True)
    if hasattr(custom_cmap, "copy"):
        custom_cmap = custom_cmap.copy()
    custom_cmap.set_bad(color="lightgray")

    fig, ax = plt.subplots(figsize=figsize)
    heatmap = sns.heatmap(
        heatmap_matrix,
        annot=annotations if annotate else False,
        fmt="",
        cmap=custom_cmap,
        cbar=True,
        vmin=vmin,
        vmax=vmax,
        mask=heatmap_matrix.isna(),
        linewidths=0.5,
        linecolor="white",
        annot_kws={"fontsize": fontsize},
        ax=ax,
    )
    colorbar = heatmap.collections[0].colorbar
    if colorbar is not None:
        colorbar.set_label(
            cbar_label,
            fontsize=fontsize,
            fontweight="bold",
        )
        colorbar.ax.tick_params(labelsize=fontsize - 1)

    ax.set_xlabel(x_label, fontsize=fontsize, fontweight="bold")
    ax.set_ylabel(y_label, fontsize=fontsize, fontweight="bold")
    if title is not None:
        ax.set_title(title, fontsize=fontsize + 2)
    ax.set_xticklabels(heatmap_matrix.columns, rotation=0, fontsize=fontsize)
    ax.set_yticklabels(heatmap_matrix.index, rotation=0, fontsize=fontsize)
    fig.tight_layout()

    if save_dir is not None:
        save_path = ensure_long_path(Path(save_dir))
        os.makedirs(save_path, exist_ok=True)
        if file_name is None:
            file_name = f"{metric}_HDF_depth_bits_profile_auc_heatmap.png"
        fig.savefig(
            ensure_long_path(save_path / file_name),
            bbox_inches="tight",
            format="png",
            dpi=900 if high_quality else 100,
        )
    if show:
        plt.show()
    else:
        plt.close(fig)

    return heatmap_matrix


def plot_hdf_depth_bits_mean_std_heatmap(
    df: pd.DataFrame,
    model: Any = None,
    hdf_depths: Any = None,
    hdf_bits: Any = None,
    metric: str = "OOF_R2",
    tree_feature_importance: str = "MDI",
    dataset_as_experiment_points: bool = False,
    figsize: tuple[float, float] = (7, 5),
    fontsize: int = 12,
    title: Optional[str] = None,
    x_label: str = "HDF bits",
    y_label: str = "HDF depth",
    cbar_label: Optional[str] = None,
    cmap: str = "viridis",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    annotate: bool = True,
    show: bool = True,
    high_quality: bool = True,
    plot_heatmap: bool = True,
    save_dir: Optional[Path] = HERE / "result_analysis",
    file_name: Optional[str] = None,
) -> pd.DataFrame:
    """Plot mean ± standard deviation across HDF evaluation tasks.

    For an OOF metric, each dataset-target-model contributes one task value.
    For a non-OOF metric, every stored seed/fold value is a separate task by
    default. Set ``dataset_as_experiment_points=True`` to first average the
    seed/fold values within each dataset-target-model and HDF configuration.

    Cell colors represent the mean and annotations show ``mean ± std``. The
    returned long-form dataframe also includes the number of task values used
    for each depth/bit combination.
    """
    expanded, depth_order, bit_order = _expand_hdf_depth_bits_scores(
        df=df,
        metric=metric,
        model=model,
        hdf_depths=hdf_depths,
        hdf_bits=hdf_bits,
        tree_feature_importance=tree_feature_importance,
    )
    task_values, _, method_columns = _hdf_task_metric_values(
        expanded,
        metric,
        dataset_as_experiment_points,
    )
    summary = (
        task_values.groupby(method_columns, as_index=False)[metric]
        .agg(mean="mean", std="std", n_tasks="count")
    )
    summary["std"] = summary["std"].fillna(0.0)
    mean_matrix = summary.pivot_table(
        index="HDF depth",
        columns="HDF bits",
        values="mean",
        aggfunc="mean",
    ).reindex(index=depth_order, columns=bit_order)
    std_matrix = summary.pivot_table(
        index="HDF depth",
        columns="HDF bits",
        values="std",
        aggfunc="mean",
    ).reindex(index=depth_order, columns=bit_order)

    if not plot_heatmap:
        return summary

    annotations = mean_matrix.copy().astype(object)
    for depth in mean_matrix.index:
        for bits in mean_matrix.columns:
            mean_value = mean_matrix.loc[depth, bits]
            std_value = std_matrix.loc[depth, bits]
            annotations.loc[depth, bits] = (
                ""
                if pd.isna(mean_value)
                else f"{mean_value:.2f} ± {std_value:.2f}"
            )

    custom_cmap = sns.color_palette(cmap, as_cmap=True)
    if hasattr(custom_cmap, "copy"):
        custom_cmap = custom_cmap.copy()
    custom_cmap.set_bad(color="lightgray")

    fig, ax = plt.subplots(figsize=figsize)
    heatmap = sns.heatmap(
        mean_matrix,
        annot=annotations if annotate else False,
        fmt="",
        cmap=custom_cmap,
        cbar=True,
        vmin=vmin,
        vmax=vmax,
        mask=mean_matrix.isna(),
        linewidths=0.5,
        linecolor="white",
        annot_kws={"fontsize": fontsize - 1},
        ax=ax,
    )
    colorbar = heatmap.collections[0].colorbar
    if colorbar is not None:
        colorbar.set_label(
            cbar_label or f"Mean {metric}",
            fontsize=fontsize,
            fontweight="bold",
        )
        colorbar.ax.tick_params(labelsize=fontsize - 1)

    ax.set_xlabel(x_label, fontsize=fontsize, fontweight="bold")
    ax.set_ylabel(y_label, fontsize=fontsize, fontweight="bold")
    if title is not None:
        ax.set_title(title, fontsize=fontsize + 2)
    ax.set_xticklabels(mean_matrix.columns, rotation=0, fontsize=fontsize)
    ax.set_yticklabels(mean_matrix.index, rotation=0, fontsize=fontsize)
    fig.tight_layout()

    if save_dir is not None:
        save_path = ensure_long_path(Path(save_dir))
        os.makedirs(save_path, exist_ok=True)
        if file_name is None:
            file_name = f"{metric}_HDF_depth_bits_mean_std_heatmap.png"
        fig.savefig(
            ensure_long_path(save_path / file_name),
            bbox_inches="tight",
            format="png",
            dpi=900 if high_quality else 100,
        )
    if show:
        plt.show()
    else:
        plt.close(fig)

    return summary


if __name__ == "__main__":
    master_device = "GPU"
    structure_results = pd.read_pickle(
        MASTER_DATA / f"Tree_and_GP_structure_and_count_{master_device}.pkl"
    )

    HDF_results = pd.read_pickle(MASTER_DATA / "Tree_HDF_CPU.pkl")

    Count_results = pd.read_pickle(MASTER_DATA / "Tree_and_GP_COUNT_only_GPU.pkl")
    structure_results = add_hdf_configuration_to_structure_results(
        structure_results,
        HDF_results,
        depth=2,
        bits=256,
    )
    structure_results = add_no_fingerprint_to_structure_results(
        structure_results,
        Count_results,
    )
    RESULTS = (
        HERE
        / "result_analysis"
    )

    STRUCTURE_REPRESENTATION_ORDER = [
        "No Fingerprint",
        "LE",
        "ECFP",
        "Mordred",
        "HDF2_256bit",
    ]
    # plot_structure_representation_profile_comparison(
    #     structure_results,
    #     metric="OOF_R2",
    #     model=["RF", "XGBR", "NGB"],
    #     # kernel_triples=("Graph", "Matern32", "(count:+)x(graph:x)"),
    #     compare_models_as_methods=True,
    #     structure_representations=STRUCTURE_REPRESENTATION_ORDER,
    #     y_label="Profile AUC of R² (OOF)",
    #     fontsize=16,
    #     figsize=(7, 5),
    #     save_dir=RESULTS/"performance_profile"/"strcuture_representation_comparison",
    #     file_name=f"R2OOF.png",
    # )
    # for x in ["RF", "XGBR", "NGB"]:
    #     plot_structure_representation_mean_std_comparison(
    #         structure_results,
    #         metric="OOF_R2",
    #         model=[x],
    #         structure_representations=STRUCTURE_REPRESENTATION_ORDER,
    #         y_label="R² (OOF)",
    #         y_lim=(0, 1.05),
    #         fontsize=17,
    #         figsize=(7, 5),
    #         show=False,
    #         save_dir=(
    #             RESULTS
    #             / "average_performance"
    #             / "strcuture_representation_comparison"
    #         ),
    #         file_name=f"R2OOF_{x}.png",
    #     )
    # for x in ["RF", "XGBR", "NGB"]:
    #     plot_structure_representation_performance_by_target_grid(
    #         structure_results,
    #         metric="OOF_R2",
    #         model=[x],
    #         structure_representations=STRUCTURE_REPRESENTATION_ORDER,
    #         y_label="R² (OOF)",
    #         y_lim=(0, 1),
    #         ncols=3,
    #         fontsize=16,
    #         figsize=(7, 5),
    #         show=False,
    #         save_dir=(
    #             RESULTS
    #             / "average_performance"
    #             / "strcuture_representation_comparison"
    #             / "all_targets"
    #         ),
    #         file_name=f"R2OOF_{x}.png",
    #     )
    plot_structure_representation_performance_by_target_grid(
            structure_results,
            metric="OOF_R2",
            model=["RF", "XGBR", "NGB"],
            structure_representations=STRUCTURE_REPRESENTATION_ORDER,
            y_label="R² (OOF)",
            y_lim=(0, 1.05),
            ncols=3,
            fontsize=16,
            figsize=(7, 5),
            show=False,
            save_dir=(
                RESULTS
                / "average_performance"
                / "strcuture_representation_comparison"
                / "all_targets"
            ),
            file_name=f"R2OOF.png",
        )
    # plot_structure_representation_profile_heatmap(
    #     structure_results,
    #     metric="OOF_R2",
    #     model=["RF", "XGBR", "NGB"],
    #     structure_representations=STRUCTURE_REPRESENTATION_ORDER,
    #     save_dir=(
    #         RESULTS
    #         / "performance_profile"
    #         / "strcuture_representation_comparison"
    #     ),
    #     file_name="R2OOF_heatmap.png",
    # )


    # Compare HDFs with different depths and bit counts.
    # plot_hdf_depth_bits_profile_heatmap(
    #     HDF_results,
    #     metric="OOF_R2",
    #     model=["RF", "XGBR", "NGB"],
    #     hdf_depths=[2, 4, 6],
    #     hdf_bits=[128, 256, 512, 1024],
    #     cbar_label="Profile AUC of R² (OOF)",
    #     fontsize=16,
    #     figsize=(7, 5),
    #     save_dir=RESULTS / "performance_profile" / "HDF_depth_bits_comparison",
    #     file_name="R2OOF_heatmap.png",
    # )

    # plot_hdf_depth_bits_mean_std_heatmap(
    #     HDF_results,
    #     metric="OOF_R2",
    #     model=["RF", "XGBR", "NGB"],
    #     hdf_depths=[2, 4, 6],
    #     hdf_bits=[128, 256, 512, 1024],
    #     cbar_label="Mean R² (OOF)",
    #     fontsize=16,
    #     figsize=(8, 5),
    #     save_dir=RESULTS / "average_performance" / "HDF_depth_bits_comparison",
    #     file_name="R2OOF_mean_std_heatmap.png",
    # )
