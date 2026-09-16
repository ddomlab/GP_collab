"""Analyze model performance across datasets of different sizes."""

import os
import re
from hashlib import sha1
from pathlib import Path
from typing import Any, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

from visualization_setting import ensure_long_path
from visualize_gp import (
    HYBRIDIZATION_METHOD_COLORS,
    _add_datapoint_counts,
    _display_model_name,
    _expand_master_scores_for_profile,
    _filter_kernel_triples,
    _metric_higher_is_better,
    _mixing_method_label,
    _model_order_sort_key,
    _selection_values,
    mixing_methods as DEFAULT_MIXING_METHODS,
)


HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent.parent / "results"
MASTER_DATA = RESULTS / "master_performance_data"



def plot_hybridization_performance_vs_data_number(
    df: pd.DataFrame,
    metric: str = "r2",
    model: Any = "GPytorchMAP",
    fp_kernels: Any = None,
    count_kernels: Any = None,
    mixing_methods: Any = None,
    kernel_triples: Any = None,
    dataset_as_experiment_points: bool = True,
    show_all_targets: bool = True,
    figsize: tuple = (10, 5),
    fontsize: int = 12,
    title: Optional[str] = None,
    x_label: str = "Number of datapoints",
    y_label: Optional[str] = None,
    x_tick_rotation: int = 0,
    y_lim: Optional[tuple] = None,
    log_y: bool = False,
    show_values: bool = True,
    show: bool = True,
    high_quality: bool = True,
    save_dir: Optional[Path] = HERE / "result_analysis",
    file_name: Optional[str] = None,
) -> pd.DataFrame:
    """Plot grouped bars of hybridization performance versus dataset size.

    Each color represents a hybridization method. When ``show_all_targets`` is
    True, every dataset-target is a separate x-axis group. Repeat/fold scores
    are first averaged for each selected fingerprint/count kernel
    configuration; bars and error bars then show the mean and standard
    deviation across those configurations.

    Set ``show_all_targets=False`` for the previous behavior: targets with the
    same datapoint count are combined, and error bars show the standard
    deviation across their dataset-target means.
    """
    df = _filter_kernel_triples(df, kernel_triples, keep_missing=False)
    is_scalar_metric = (
        metric in df.columns
        and f"{metric}_seed_fold_scores" not in df.columns
    )
    plot_df = _expand_master_scores_for_profile(
        df,
        metric=metric,
        model=model,
        fp_kernels=fp_kernels,
        count_kernels=count_kernels,
        mixing_methods=mixing_methods,
    )
    plot_df = plot_df.dropna(subset=["mixing method"]).copy()
    if plot_df.empty:
        raise ValueError(
            f"No {metric} values found for model={model} with "
            f"fp_kernels={fp_kernels}, count_kernels={count_kernels}, "
            f"mixing_methods={mixing_methods}, and kernel_triples={kernel_triples}."
        )

    plot_df = _add_datapoint_counts(plot_df)
    plot_df["hybridization method"] = plot_df["mixing method"].map(
        lambda method: _mixing_method_label(method)
        .replace("ΣF", "F")
        .replace("ΠF", "F")
    )
    selected_models = _selection_values(model)
    include_model = selected_models is None or len(selected_models) > 1
    hue_col = (
        "model hybridization method"
        if include_model
        else "hybridization method"
    )
    if include_model:
        plot_df[hue_col] = (
            plot_df["model"].map(_display_model_name)
            + "\n"
            + plot_df["hybridization method"].astype(str)
        )

    if show_all_targets:
        # Give every selected kernel configuration equal weight, regardless of
        # the number of repeat/fold scores stored for it.
        configuration_group_cols = [
            "dataset",
            "target",
            "n datapoints",
            hue_col,
            "hybridization method",
            "model",
            "fp kernel",
            "count kernel",
            "mixing method",
        ]
        configuration_group_cols = list(dict.fromkeys(configuration_group_cols))
        configuration_df = (
            plot_df.groupby(
                configuration_group_cols,
                dropna=False,
                as_index=False,
            )[metric]
            .mean()
            .copy()
        )

        summary_group_cols = list(dict.fromkeys([
            "dataset",
            "target",
            "n datapoints",
            hue_col,
            "hybridization method",
            "model",
            "mixing method",
        ]))
        summary_df = (
            configuration_df.groupby(summary_group_cols, dropna=False)[metric]
            .agg(["mean", "std", "count"])
            .reset_index()
            .rename(columns={
                "mean": f"{metric}_mean",
                "std": f"{metric}_std",
                "count": f"{metric}_count",
            })
        )
        x_col = "dataset target"
        summary_df[x_col] = (
            summary_df["dataset"].astype(str)
            + "\n"
            + summary_df["target"].astype(str)
        )
    elif dataset_as_experiment_points:
        dataset_group_cols = [
            "dataset",
            "target",
            "n datapoints",
            hue_col,
            "hybridization method",
            "model",
            "mixing method",
        ]
        dataset_group_cols = list(dict.fromkeys(dataset_group_cols))
        plot_df = (
            plot_df.groupby(
                dataset_group_cols,
                dropna=False,
                as_index=False,
            )[metric]
            .mean()
            .copy()
        )

    if not show_all_targets:
        summary_group_cols = list(dict.fromkeys([
            "n datapoints",
            hue_col,
            "hybridization method",
            "model",
            "mixing method",
        ]))
        summary_df = (
            plot_df.groupby(summary_group_cols, dropna=False)[metric]
            .agg(["mean", "std", "count"])
            .reset_index()
            .rename(columns={
                "mean": f"{metric}_mean",
                "std": f"{metric}_std",
                "count": f"{metric}_count",
            })
        )
        x_col = "n datapoints"
    summary_df[f"{metric}_std"] = summary_df[f"{metric}_std"].fillna(0.0)

    selected_mixes = _selection_values(mixing_methods) or DEFAULT_MIXING_METHODS
    selected_models = selected_models or plot_df["model"].drop_duplicates().tolist()
    method_order = {
        str(method).lower(): index
        for index, method in enumerate(selected_mixes)
    }
    model_order = {
        str(model_name).lower(): index
        for index, model_name in enumerate(selected_models)
    }
    config_order = summary_df[
        [hue_col, "model", "mixing method"]
    ].drop_duplicates(subset=[hue_col])
    config_order["sort_key"] = config_order.apply(
        lambda row: (
            *_model_order_sort_key(row["model"], model_order),
            method_order.get(
                str(row["mixing method"]).lower(),
                len(method_order),
            ),
            str(row["mixing method"]),
        ),
        axis=1,
    )
    config_order = config_order.sort_values("sort_key", kind="mergesort")
    hue_order = config_order[hue_col].tolist()
    if show_all_targets:
        target_order_df = (
            summary_df[[x_col, "dataset", "target", "n datapoints"]]
            .drop_duplicates(subset=[x_col])
            .sort_values(
                ["n datapoints", "dataset", "target"],
                kind="mergesort",
            )
        )
        x_order = target_order_df[x_col].tolist()
        x_tick_labels = []
        for _, row in target_order_df.iterrows():
            n_datapoints = row["n datapoints"]
            n_label = (
                f"{int(n_datapoints)}"
                if float(n_datapoints).is_integer()
                else f"{n_datapoints:g}"
            )
            target_label = str(row["target"]).removeprefix("target_")
            x_tick_labels.append(f"{target_label}\n(n={n_label})")
    else:
        x_order = sorted(summary_df["n datapoints"].dropna().unique())
        x_tick_labels = [
            f"{int(value)}" if float(value).is_integer() else f"{value:g}"
            for value in x_order
        ]

    order_lookup = {label: index for index, label in enumerate(hue_order)}
    summary_df["_plot_order"] = summary_df[hue_col].map(order_lookup)
    summary_df = (
        summary_df.sort_values(
            ["n datapoints", "dataset", "target", "_plot_order"]
            if show_all_targets
            else ["n datapoints", "_plot_order"],
            kind="stable",
        )
        .drop(columns="_plot_order")
        .copy()
    )

    fallback_colors = sns.color_palette("Set2", n_colors=max(len(hue_order), 1))
    palette = {
        row[hue_col]: HYBRIDIZATION_METHOD_COLORS.get(
            str(row["mixing method"]),
            fallback_colors[index % len(fallback_colors)],
        )
        for index, (_, row) in enumerate(config_order.iterrows())
    }

    plot_figsize = figsize
    if show_all_targets:
        plot_figsize = (
            figsize[0],
            max(figsize[1], 0.8 * len(x_order) + 2),
        )

    fig, ax = plt.subplots(figsize=plot_figsize)
    if show_all_targets:
        sns.barplot(
            data=summary_df,
            x=f"{metric}_mean",
            y=x_col,
            hue=hue_col,
            order=x_order,
            hue_order=hue_order,
            palette=palette,
            width=0.72,
            errorbar=None,
            orient="h",
            ax=ax,
        )
    else:
        sns.barplot(
            data=summary_df,
            x=x_col,
            y=f"{metric}_mean",
            hue=hue_col,
            order=x_order,
            hue_order=hue_order,
            palette=palette,
            width=0.9,
            errorbar=None,
            ax=ax,
        )

    summary_lookup = summary_df.set_index([x_col, hue_col])[
        [f"{metric}_mean", f"{metric}_std"]
    ]
    max_label_value = None
    for hue_index, container in enumerate(ax.containers[:len(hue_order)]):
        hue_value = hue_order[hue_index]
        for data_index, bar in enumerate(container.patches):
            if data_index >= len(x_order):
                continue
            x_value = x_order[data_index]
            if (x_value, hue_value) not in summary_lookup.index:
                continue

            values = summary_lookup.loc[(x_value, hue_value)]
            mean_value = float(values[f"{metric}_mean"])
            std_value = float(values[f"{metric}_std"])
            if np.isnan(mean_value):
                continue

            if show_all_targets:
                score_position = bar.get_x() + bar.get_width()
                target_position = bar.get_y() + bar.get_height() / 2
                if std_value > 0:
                    ax.errorbar(
                        score_position,
                        target_position,
                        xerr=std_value,
                        fmt="none",
                        ecolor="black",
                        elinewidth=1.1,
                        capsize=3,
                        capthick=1.1,
                        zorder=4,
                    )
                label_value = mean_value + (
                    std_value if std_value > 0 else 0
                ) + 0.01
                max_label_value = (
                    label_value
                    if max_label_value is None
                    else max(max_label_value, label_value)
                )
                if show_values:
                    left, right = ax.get_xlim()
                    value_x = left + 0.02 * (right - left)
                    ax.text(
                        value_x,
                        target_position,
                        f"{mean_value:.2f}",
                        ha="left",
                        va="center",
                        fontsize=max(fontsize - 6, 6),
                        fontweight="bold",
                        color="black" if mean_value < 0.1 else "white",
                    )
            else:
                data_position = bar.get_x() + bar.get_width() / 2
                if std_value > 0:
                    ax.errorbar(
                        data_position,
                        mean_value,
                        yerr=std_value,
                        fmt="none",
                        ecolor="black",
                        elinewidth=1.1,
                        capsize=3,
                        capthick=1.1,
                        zorder=4,
                    )
                label_value = mean_value + (
                    std_value if std_value > 0 else 0
                ) + 0.01
                max_label_value = (
                    label_value
                    if max_label_value is None
                    else max(max_label_value, label_value)
                )
                if show_values:
                    bottom, top = ax.get_ylim()
                    value_y = bottom + 0.02 * (top - bottom)
                    ax.text(
                        data_position,
                        value_y,
                        f"{mean_value:.2f}",
                        ha="center",
                        va="bottom",
                        fontsize=max(fontsize - 6, 6),
                        fontweight="bold",
                        color="black" if mean_value < 0.1 else "white",
                        rotation=90,
                    )

    if show_all_targets:
        ax.set_xlabel(
            y_label or f"Mean {metric}",
            fontsize=fontsize,
            fontweight="bold",
        )
        target_axis_label = (
            "Target (number of datapoints)"
            if x_label == "Number of datapoints"
            else x_label
        )
        ax.set_ylabel(target_axis_label, fontsize=fontsize, fontweight="bold")
    else:
        ax.set_xlabel(x_label, fontsize=fontsize, fontweight="bold")
        ax.set_ylabel(
            y_label or f"Mean {metric}",
            fontsize=fontsize,
            fontweight="bold",
        )
    if title is not None:
        ax.set_title(title, fontsize=fontsize + 2)
    ax.tick_params(axis="both", labelsize=fontsize - 2)
    if show_all_targets:
        ax.set_yticks(range(len(x_order)))
        ax.set_yticklabels(x_tick_labels)
    else:
        ax.set_xticks(range(len(x_order)))
        ax.set_xticklabels(
            x_tick_labels,
            rotation=x_tick_rotation,
            ha="right" if x_tick_rotation else "center",
        )

    if log_y:
        positive_values = pd.to_numeric(
            summary_df[f"{metric}_mean"],
            errors="coerce",
        )
        if (positive_values.dropna() <= 0).any():
            raise ValueError("log_y=True requires all plotted values to be positive.")
        if show_all_targets:
            ax.set_xscale("log")
        else:
            ax.set_yscale("log")

    if y_lim is not None:
        if show_all_targets:
            ax.set_xlim(*y_lim)
        else:
            ax.set_ylim(*y_lim)
    elif log_y:
        metric_minimum = pd.to_numeric(
            summary_df[f"{metric}_mean"],
            errors="coerce",
        ).min() * 0.8
        if show_all_targets:
            ax.set_xlim(left=metric_minimum)
        else:
            ax.set_ylim(bottom=metric_minimum)
    elif str(metric).strip().lower() in {"r2", "oof_r2"}:
        if show_all_targets:
            ax.set_xlim(0, 1.05)
        else:
            ax.set_ylim(0, 1.05)
    elif is_scalar_metric or not _metric_higher_is_better(metric):
        if show_all_targets:
            ax.set_xlim(left=0)
        else:
            ax.set_ylim(bottom=0)
    if max_label_value is not None:
        if show_all_targets:
            left, right = ax.get_xlim()
            ax.set_xlim(left, max(right, max_label_value + 0.05))
        else:
            bottom, top = ax.get_ylim()
            ax.set_ylim(bottom, max(top, max_label_value + 0.05))

    legend = ax.get_legend()
    if legend is not None:
        handles, labels = ax.get_legend_handles_labels()
        legend.remove()
        legend = ax.legend(
            handles,
            labels,
            loc="lower center",
            bbox_to_anchor=(0.5, 1.02),
            ncol=max(1, len(labels)),
            frameon=False,
            title=None,
        )
        for text in legend.get_texts():
            text.set_fontsize(fontsize - 3)

    plt.tight_layout(rect=(0, 0, 1, 0.92))

    if save_dir is not None:
        save_dir = ensure_long_path(Path(save_dir))
        os.makedirs(save_dir, exist_ok=True)
        if file_name is None:
            model_name = "_".join(
                str(value)
                for value in (_selection_values(model) or ["all"])
            )
            file_name = (
                f"{model_name}_{metric}_hybridization_performance_"
                "vs_data_number.png"
            )
        fig.savefig(
            ensure_long_path(save_dir / file_name),
            bbox_inches="tight",
            format="png",
            dpi=900 if high_quality else 100,
        )

    if show:
        plt.show()
    else:
        plt.close(fig)

    return summary_df


def _safe_filename_component(value: Any) -> str:
    """Return a filesystem-safe, readable filename component."""
    component = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._")
    return component or "unnamed"


def _compact_filename_component(value: Any, max_length: int = 60) -> str:
    """Shorten long components while retaining a collision-resistant suffix."""
    component = _safe_filename_component(value)
    if len(component) <= max_length:
        return component
    digest = sha1(component.encode("utf-8")).hexdigest()[:8]
    return f"{component[:max_length - len(digest) - 1]}_{digest}"


def plot_regular_vs_permuted_fp_performance_by_target(
    regular_df: pd.DataFrame,
    permuted_df: pd.DataFrame,
    count_only_df: pd.DataFrame,
    metric: str = "OOF_R2",
    model: str = "GPytorchMAP",
    fp_kernels: Any = None,
    count_kernels: Any = None,
    mixing_methods: Any = None,
    figsize: tuple = (7, 5),
    ncols: int = 3,
    fontsize: int = 13,
    y_label: Optional[str] = None,
    y_lim: Optional[tuple] = None,
    log_y: bool = False,
    show_values: bool = True,
    show: bool = False,
    high_quality: bool = True,
    save_dir: Optional[Path] = None,
) -> pd.DataFrame:
    """Compare count-only, regular, and permuted-FP performance by target.

    One subplot is created per dataset-target pair, with all subplots combined
    in a single figure. The number of datapoints is shown in each subplot
    title. ``figsize`` specifies the size of one subplot and ``ncols`` controls
    the number of subplot columns. Mixing methods are shown on the x-axis and
    the training conditions are placed next to each other. Regular and
    permuted-fingerprint configurations are included only when both master
    datasets contain a valid score for them. Count-only bars are added for
    ``sum`` and ``product`` when matching results are available.

    ``fp_kernels`` and ``count_kernels`` accept either one kernel or a list.
    When multiple matching kernel configurations are selected, their values
    are averaged once for each mixing method and training condition. With one
    matching configuration, its value is retained.

    The returned DataFrame contains the plotted mean, standard deviation,
    value count, and output path for every condition and mixing method.
    """
    selected_fp_kernels = _selection_values(fp_kernels) or ["TanimotoMatern32"]
    selected_count_kernels = _selection_values(count_kernels) or ["Matern32"]
    selected_methods = _selection_values(mixing_methods) or list(
        DEFAULT_MIXING_METHODS
    )
    condition_frames = []
    for condition, source_df in (
        ("C + FP", regular_df),
        ("C + Permuted FP", permuted_df),
    ):
        expanded = _expand_master_scores_for_profile(
            source_df,
            metric=metric,
            model=model,
            fp_kernels=selected_fp_kernels,
            count_kernels=selected_count_kernels,
            mixing_methods=selected_methods,
        )
        if expanded.empty:
            continue
        expanded = expanded.dropna(subset=[metric, "mixing method"]).copy()
        expanded["training condition"] = condition
        condition_frames.append(expanded)

    if len(condition_frames) != 2:
        raise ValueError(
            "Both regular_df and permuted_df must contain valid "
            f"{metric} values for model={model}, "
            f"fp_kernels={selected_fp_kernels}, and "
            f"count_kernels={selected_count_kernels}."
        )

    paired_comparison_df = pd.concat(condition_frames, ignore_index=True)
    paired_comparison_df = _add_datapoint_counts(paired_comparison_df)
    pair_columns = [
        "dataset",
        "target",
        "fp kernel",
        "count kernel",
        "mixing method",
    ]
    paired_methods = (
        paired_comparison_df.groupby(
            pair_columns,
            dropna=False,
        )["training condition"]
        .nunique()
        .loc[lambda counts: counts == 2]
        .index
    )
    paired_method_df = pd.DataFrame(
        paired_methods.tolist(),
        columns=pair_columns,
    )
    paired_comparison_df = paired_comparison_df.merge(
        paired_method_df,
        on=pair_columns,
        how="inner",
    )
    if paired_comparison_df.empty:
        raise ValueError(
            "No dataset-target and mixing-method combinations have valid "
            "scores in both the regular and permuted master datasets."
        )

    count_only_methods = [
        method
        for method in selected_methods
        if str(method).lower() in {"sum", "product"}
    ]
    count_only_expanded = _expand_master_scores_for_profile(
        count_only_df,
        metric=metric,
        model=model,
        count_kernels=selected_count_kernels,
        mixing_methods=count_only_methods,
    )
    count_only_expanded = count_only_expanded.dropna(
        subset=[metric, "mixing method"]
    ).copy()
    if count_only_expanded.empty:
        raise ValueError(
            "count_only_df contains no valid sum or product values for "
            f"metric={metric}, model={model}, and "
            f"count_kernels={selected_count_kernels}."
        )
    count_only_expanded["training condition"] = "C only"
    count_only_expanded = _add_datapoint_counts(count_only_expanded)
    count_match_columns = [
        "dataset",
        "target",
        "count kernel",
        "mixing method",
    ]
    paired_count_keys = paired_comparison_df[
        count_match_columns
    ].drop_duplicates()
    count_only_expanded = count_only_expanded.merge(
        paired_count_keys,
        on=count_match_columns,
        how="inner",
    )
    if count_only_expanded.empty:
        raise ValueError(
            "No count-only sum or product results match the paired regular "
            "and permuted dataset-target configurations."
        )
    comparison_df = pd.concat(
        [paired_comparison_df, count_only_expanded],
        ignore_index=True,
    )

    comparison_df["hybridization method"] = comparison_df["mixing method"].map(
        lambda method: _mixing_method_label(method)
        .replace("ΣF", "F")
        .replace("ΠF", "F")
    )
    summary_df = (
        comparison_df.groupby(
            [
                "dataset",
                "target",
                "n datapoints",
                "mixing method",
                "hybridization method",
                "training condition",
            ],
            dropna=False,
        )[metric]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(
            columns={
                "mean": f"{metric}_mean",
                "std": f"{metric}_std",
                "count": f"{metric}_count",
            }
        )
    )
    summary_df[f"{metric}_std"] = summary_df[f"{metric}_std"].fillna(0.0)
    summary_df["output path"] = None

    method_rank = {
        str(method).lower(): index
        for index, method in enumerate(selected_methods)
    }
    condition_order = ["C only", "C + FP", "C + Permuted FP"]
    condition_palette = {
        "C only": "#d7191c",
        "C + FP": "#fdae61",
        "C + Permuted FP": "#2c7bb6",
    }
    if ncols < 1:
        raise ValueError("ncols must be at least 1.")

    if save_dir is not None:
        save_dir = ensure_long_path(Path(save_dir))
        os.makedirs(save_dir, exist_ok=True)

    summary_df = summary_df.sort_values(
        ["n datapoints", "dataset", "target"],
        kind="stable",
    ).reset_index(drop=True)
    group_columns = ["dataset", "target"]
    target_groups = list(summary_df.groupby(group_columns, sort=False))
    nrows = int(np.ceil(len(target_groups) / ncols))
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(figsize[0] * ncols, figsize[1] * nrows),
        sharey=False,
        squeeze=False,
    )
    axes_flat = axes.ravel()
    legend_handles = [
        Patch(facecolor=condition_palette[condition], label=condition)
        for condition in condition_order
    ]

    for ax, ((dataset, target), target_df) in zip(axes_flat, target_groups):
        target_df = target_df.copy()
        target_df["_method_rank"] = target_df["mixing method"].map(
            lambda method: method_rank.get(
                str(method).lower(),
                len(method_rank),
            )
        )
        target_df = target_df.sort_values(
            ["_method_rank", "mixing method", "training condition"],
            kind="stable",
        )
        method_order = (
            target_df[["hybridization method", "_method_rank"]]
            .drop_duplicates()
            .sort_values("_method_rank", kind="stable")["hybridization method"]
            .tolist()
        )

        sns.barplot(
            data=target_df,
            x="hybridization method",
            y=f"{metric}_mean",
            hue="training condition",
            order=method_order,
            hue_order=condition_order,
            palette=condition_palette,
            errorbar=None,
            ax=ax,
        )

        summary_lookup = target_df.set_index(
            ["hybridization method", "training condition"]
        )[[f"{metric}_mean", f"{metric}_std"]]
        max_label_value = None
        for condition_index, container in enumerate(
            ax.containers[: len(condition_order)]
        ):
            condition = condition_order[condition_index]
            condition_methods = [
                method
                for method in method_order
                if (method, condition) in summary_lookup.index
            ]
            for bar, method in zip(container.patches, condition_methods):
                values = summary_lookup.loc[(method, condition)]
                mean_value = float(values[f"{metric}_mean"])
                std_value = float(values[f"{metric}_std"])
                x_position = bar.get_x() + bar.get_width() / 2
                if std_value > 0:
                    ax.errorbar(
                        x_position,
                        mean_value,
                        yerr=std_value,
                        fmt="none",
                        ecolor="black",
                        elinewidth=1.0,
                        capsize=3,
                        capthick=1.0,
                        zorder=4,
                    )
                label_value = mean_value + max(std_value, 0) + 0.01
                max_label_value = (
                    label_value
                    if max_label_value is None
                    else max(max_label_value, label_value)
                )
                if show_values:
                    ax.text(
                        x_position,
                        label_value,
                        f"{mean_value:.2f}",
                        ha="center",
                        va="bottom",
                        fontsize=max(fontsize - 3, 7),
                    )

        target_label = str(target).removeprefix("target_")
        n_datapoints = target_df["n datapoints"].iloc[0]
        n_label = (
            f"{int(n_datapoints)}"
            if float(n_datapoints).is_integer()
            else f"{n_datapoints:g}"
        )
        ax.set_title(
            f"{target_label} (n={n_label})",
            fontsize=fontsize + 1,
        )
        ax.set_xlabel("Hybridization method", fontsize=fontsize)
        ax.set_ylabel(y_label or f"Mean {metric}", fontsize=fontsize)
        ax.tick_params(axis="both", labelsize=fontsize - 2)
        ax.tick_params(axis="x", rotation=20)

        if log_y:
            values = pd.to_numeric(target_df[f"{metric}_mean"], errors="coerce")
            if (values.dropna() <= 0).any():
                plt.close(fig)
                raise ValueError(
                    "log_y=True requires all plotted values to be positive."
                )
            ax.set_yscale("log")
        if y_lim is not None:
            ax.set_ylim(*y_lim)
        elif str(metric).strip().lower() in {"r2", "oof_r2"}:
            ax.set_ylim(0, 1.05)
        elif not _metric_higher_is_better(metric):
            ax.set_ylim(bottom=0)
        if max_label_value is not None and not log_y:
            bottom, top = ax.get_ylim()
            ax.set_ylim(bottom, max(top, max_label_value + 0.04))

        legend = ax.get_legend()
        if legend is not None:
            legend.remove()

    for ax in axes_flat[len(target_groups):]:
        ax.set_visible(False)

    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.995),
        ncol=len(condition_order),
        frameon=False,
        fontsize=fontsize + 5,
        handlelength=2.5,
        handleheight=1.5,
        columnspacing=2.0,
    )

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    if save_dir is not None:
        filename_parts = [
            model,
            "-".join(map(str, selected_fp_kernels)),
            "-".join(map(str, selected_count_kernels)),
            metric,
            "count_only_vs_regular_vs_permuted_fp_by_target",
        ]
        file_name = "_".join(
            _compact_filename_component(part) for part in filename_parts
        ) + ".png"
        output_path = ensure_long_path(save_dir / file_name)
        fig.savefig(
            output_path,
            bbox_inches="tight",
            format="png",
            dpi=900 if high_quality else 100,
        )
        summary_df["output path"] = str(output_path)

    if show:
        plt.show()
    plt.close(fig)

    return summary_df


if __name__ == "__main__":
    master_device = "GPU"
    count_and_fingerprint_result = pd.read_pickle(
        MASTER_DATA / f"Tree_and_GP_{master_device}_count_and_fingerprint.pkl"
    )
    count_results = pd.read_pickle(
        MASTER_DATA / f"Tree_and_GP_COUNT_only_{master_device}.pkl"
    )
    permuted_fp_result = pd.read_pickle(
        MASTER_DATA / f"GP_random_fp_permutation_{master_device}.pkl"
    )
    
    Separate_datasets_save = (
        HERE
        / "result_analysis"
        / "absolute_metric"
        / "hybridization_comparison"
        / "separate_dataset"
    )


    plot_hybridization_performance_vs_data_number(
        df=count_and_fingerprint_result,
        metric="OOF_R2",
        model="GPytorchMAP",
        fp_kernels=[
            "RBF",
            "Matern32",
            "Matern52",
        ],
        count_kernels=["Matern32", "Matern52", "RBF"],
        show_all_targets=True,
        mixing_methods=[
            "sum",
            "product",
            # "(count:+)x(fp:x)",
            "(count:+)x(fp:+)",  # train on this
            "(count:x)+(fp:x)",
        ],
        y_label="R² (OOF)",
        fontsize=17,
        figsize=(11, 8),
        show=True,
        save_dir=(
            HERE
            / "result_analysis"
            / "absolute_metric"
            / "hybridization_comparison"
        ),
        file_name="R2OOF_GPytorchMAP_Bitwise_all_config_vs_data_number_all_targets.png",
    )


    # plot_hybridization_performance_vs_data_number(
    #     df=count_and_fingerprint_result,
    #     metric="OOF_R2",
    #     model="MGK",
    #     fp_kernels=["Graph"],
    #     count_kernels=["Matern32", "Matern52", "RBF"],
    #     show_all_targets=True,
    #     mixing_methods=[
    #         "sum",
    #         "product",
    #         "(count:+)x(Graph:x)",# train on this
    #         "(count:x)+(Graph:x)",
    #     ],
    #     y_label="R² (OOF)",
    #     fontsize=17,
    #     figsize=(11, 8),
    #     show=True,
    #     save_dir=(
    #         HERE
    #         / "result_analysis"
    #         / "absolute_metric"
    #         / "hybridization_comparison"
    #     ),
    #     file_name="R2OOF_MGK_all_config_vs_data_number_all_targets.png",
    # )


    plot_regular_vs_permuted_fp_performance_by_target(
        regular_df=count_and_fingerprint_result,
        permuted_df=permuted_fp_result,
        count_only_df=count_results,
        metric="OOF_R2",
        model="GPytorchMAP",
        fp_kernels=[
            "Matern32",
        ],
        count_kernels=["Matern32"],
        mixing_methods=[
            "sum",
            "product",
            "(count:+)x(fp:+)",
            "(count:x)+(fp:x)",
        ],
        y_label="R² (OOF)",
        fontsize=15,
        figsize=(6, 5),
        show=False,
        save_dir=Separate_datasets_save,
    )
