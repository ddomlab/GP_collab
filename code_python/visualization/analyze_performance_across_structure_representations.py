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
    _display_model_name,
    _expand_master_scores_for_profile,
    _filter_kernel_triples,
    _metric_higher_is_better,
    _selection_values,
)


HERE = Path(__file__).resolve().parent
MASTER_DATA = HERE.parent.parent / "results" / "master_performance_data"

DEFAULT_STRUCTURE_REPRESENTATION_ORDER = ["ECFP", "Graph", "LE", "Mordred"]
STRUCTURE_REPRESENTATION_COLORS = {
    "ECFP": "#900101",
    "Graph": "#F58518",
    "LE": "#298C8C",
    "Mordred": "#716F6F",
}


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
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Calculate representation AUCs with model included in each task."""
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

    task_columns = [
        "dataset",
        "target",
        "model",
        "configuration task",
    ]
    if not dataset_as_experiment_points:
        task_columns.append("repeat")

    # A task can contain several raw rows when repeat/fold scores or duplicate
    # result records are present. Reduce it to one score per representation.
    profile_rank_df = (
        profile_df.groupby(
            task_columns + ["structure representation"],
            dropna=False,
            as_index=False,
        )[metric]
        .mean()
    )

    # A within-task rank is meaningful only when at least two representations
    # are available for the same model and configuration. This intentionally
    # leaves unsupported model/representation cells blank in the heatmap.
    comparable_tasks = (
        profile_rank_df.groupby(task_columns, dropna=False)[
            "structure representation"
        ]
        .nunique()
        .loc[lambda counts: counts >= 2]
        .reset_index()[task_columns]
    )
    profile_rank_df = profile_rank_df.merge(
        comparable_tasks,
        on=task_columns,
        how="inner",
    )
    if profile_rank_df.empty:
        raise ValueError(
            "No comparable tasks contain at least two structure "
            "representations for the same model and kernel configuration."
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
            label=f"{representation}, {representation_auc:.2f}",
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
    plot_profile: bool = False,
    profile_options: Optional[dict[str, Any]] = None,
    figsize: tuple[float, float] = (7, 5),
    fontsize: int = 12,
    title: Optional[str] = None,
    x_label: str = "Structure representation",
    y_label: str = "Performance profile AUC",
    x_tick_rotation: int = 0,
    y_lim: tuple[float, float] = (0, 1.05),
    show_values: bool = True,
    gradient_bars: bool = True,
    show: bool = True,
    high_quality: bool = True,
    save_dir: Optional[Path] = HERE / "result_analysis",
    file_name: Optional[str] = None,
) -> pd.DataFrame:
    """Compare representations using generalized performance-profile AUC.

    Models are part of the task definition together with dataset, target,
    kernel configuration, and (by default) repeat/fold. Representations are
    therefore ranked within the same model rather than treating different
    models as competing alternatives. The resulting bar for each
    representation summarizes all comparable model tasks.

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
        representation_order,
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

    Each cell uses the percentile ranks obtained by comparing representations
    within the same dataset, target, model, kernel configuration, and
    repeat/fold. Missing or non-comparable model/representation combinations
    are left blank. Set ``plot_heatmap=False`` to return the matrix without
    creating a figure.
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


if __name__ == "__main__":
    master_device = "GPU"
    structure_results = pd.read_pickle(
        MASTER_DATA / f"Tree_and_GP_structure_and_count_{master_device}.pkl"
    )
    RESULTS = (
        HERE
        / "result_analysis"
    )

    STRUCTURE_REPRESENTATION_ORDER = ["LE", "ECFP", "Mordred", "Graph"]
    plot_structure_representation_profile_comparison(
        structure_results,
        metric="OOF_R2",
        model=["RF", "XGBR", "NGB"],
        structure_representations=STRUCTURE_REPRESENTATION_ORDER,
        y_label="Profile AUC of R² (OOF)",
        fontsize=16,
        figsize=(6, 5),
        save_dir=RESULTS/"performance_profile"/"strcuture_representation_comparison",
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
