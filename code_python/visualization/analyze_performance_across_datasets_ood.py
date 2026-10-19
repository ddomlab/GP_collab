"""Compare model and hybridization performance on OOD or matched-IID cases.

Dataset sizes refer to the full dataset-target sizes in the regular analysis
metadata, rather than cluster test sizes or the size of the training split.
"""

from pathlib import Path
import re
from typing import Any, Optional
import warnings

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

from visualization_setting import (
    ensure_long_path,
    lighten_color,
    model_type_color,
    selection_values,
)
from visualize_gp import _add_datapoint_counts
from visualize_ood_gp import (
    EVALUATIONS,
    HYBRIDIZATION_METHOD_COLORS,
    _hybridization_input,
    _hybridization_metadata,
    _metric_name,
    _normalize_evaluations,
    _keep_common_complete_clusters,
    _ordered_metadata,
    _select_scores,
)


HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent.parent / "results"
MASTER_DATA = RESULTS / "master_ood_performance_data"
DEFAULT_SAVE_DIR = HERE / "ood_result_analysis"


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
    show_values: bool = False,
    show: bool = True,
    high_quality: bool = False,
    save_dir: Optional[Path] = DEFAULT_SAVE_DIR,
    file_name: Optional[str] = None,
    ood_and_iid: str = "OOD",
    complete_cases: bool = True,
) -> pd.DataFrame:
    """Compare hybridization methods by dataset size for OOD or IID only.

    Pass ``ood_and_iid="OOD"`` or ``"IID"``; no availability from the
    other split is required. By default only clusters with scores for every
    selected model/kernel/method configuration in that split are retained.

    Repeats are first averaged within each cluster and kernel configuration.
    Selected kernel pairs are then averaged within each cluster/method.
    With ``show_all_targets=True``, each dataset-target has a separate
    horizontal bar group; bars show mean +/- sample std across its clusters.

    With ``show_all_targets=False``, tasks with the same full dataset size
    share a vertical bar group. By default each dataset-target contributes its
    mean across clusters equally, and std is across those task means. Set
    ``dataset_as_experiment_points=False`` to instead pool cluster means
    within each size group. IID seeds never receive extra cluster weight.

    ``y_lim`` limits the score axis (horizontal in the per-target view).
    Returned columns include evaluation, metric, method, mean, std and count.
    Annotations and high-resolution saving are disabled by default.
    """
    if ood_and_iid is None:
        raise ValueError("ood_and_iid must be either 'OOD' or 'IID'.")
    evaluation = _normalize_evaluations(ood_and_iid)[0]
    metric = _metric_name(metric)
    if y_lim is not None:
        if len(y_lim) != 2 or not np.isfinite(y_lim).all() or not y_lim[0] < y_lim[1]:
            raise ValueError("y_lim must contain two finite values with lower < upper.")
        if log_y and y_lim[0] <= 0:
            raise ValueError("log_y=True requires a positive lower score limit.")

    scores, metadata = _select_scores(
        _hybridization_input(df),
        metric=metric,
        model=model,
        fp_kernels=fp_kernels,
        count_kernels=count_kernels,
        mixing_methods=mixing_methods,
        kernel_triples=kernel_triples,
        include_kernel_config=True,
        complete_cases=complete_cases,
        ood_and_iid=evaluation,
    )
    metadata = _hybridization_metadata(metadata, model, mixing_methods)
    scores = scores.merge(
        metadata[["configuration key", "hybridization method"]],
        on="configuration key",
        validate="many_to_one",
    )
    groups = ["paper", "target", "model", "mixing method", "hybridization method"]
    kernel_cases = scores.groupby(
        [*groups, "fp kernel", "count kernel", "cluster"],
        as_index=False,
        dropna=False,
    )["score"].mean()
    cases = kernel_cases.groupby(
        [*groups, "cluster"], as_index=False, dropna=False,
    )["score"].mean().rename(columns={"paper": "dataset"})
    cases = _add_datapoint_counts(cases)

    method_groups = ["model", "mixing method", "hybridization method"]
    if show_all_targets:
        values = cases
        summary_groups = ["dataset", "target", "n datapoints", *method_groups]
    else:
        values = cases
        if dataset_as_experiment_points:
            values = cases.groupby(
                ["dataset", "target", "n datapoints", *method_groups],
                as_index=False,
                dropna=False,
            )["score"].mean()
        summary_groups = ["n datapoints", *method_groups]

    summary = values.groupby(
        summary_groups, as_index=False, dropna=False,
    )["score"].agg(mean="mean", std="std", count="count")
    summary["std"] = summary["std"].fillna(0.0)
    summary["evaluation"] = evaluation
    summary["metric"] = metric
    # Retain the metric-specific summary column names used by regular analysis.
    summary = summary.rename(columns={
        "mean": f"{metric}_mean",
        "std": f"{metric}_std",
        "count": f"{metric}_count",
    })
    mean_col, std_col = f"{metric}_mean", f"{metric}_std"
    if log_y and summary[mean_col].le(0).any():
        raise ValueError("log_y=True requires all plotted means to be positive.")

    method_order = metadata["hybridization method"].drop_duplicates().tolist()
    palette = {
        row["hybridization method"]: HYBRIDIZATION_METHOD_COLORS.get(
            str(row["mixing method"]), "#808080",
        )
        for _, row in metadata.iterrows()
    }
    if show_all_targets:
        target_order = summary[["dataset", "target", "n datapoints"]].drop_duplicates()
        target_order = target_order.sort_values(
            ["n datapoints", "dataset", "target"], kind="stable",
        )
        # Use task identities for plotting; target display names need not be unique.
        task_index = pd.MultiIndex.from_frame(target_order[["dataset", "target"]])
        positions = np.arange(len(task_index), dtype=float)
        tick_labels = [
            f"{str(target).removeprefix('target_')}\n(n={size:g})"
            for target, size in zip(target_order["target"], target_order["n datapoints"])
        ]
        plot_figsize = (figsize[0], max(figsize[1], 0.8 * len(task_index) + 2))
        lookup_columns = ["dataset", "target"]
    else:
        sizes = sorted(summary["n datapoints"].unique())
        task_index = pd.Index(sizes, name="n datapoints")
        positions = np.arange(len(sizes), dtype=float)
        tick_labels = [f"{size:g}" for size in sizes]
        plot_figsize = figsize
        lookup_columns = ["n datapoints"]

    summary["_method_order"] = summary["hybridization method"].map(
        {method: i for i, method in enumerate(method_order)}
    )
    summary = summary.sort_values(
        ["n datapoints", "dataset", "target", "_method_order"]
        if show_all_targets else ["n datapoints", "_method_order"],
        kind="stable",
    ).drop(columns="_method_order").reset_index(drop=True)

    fig, ax = plt.subplots(figsize=plot_figsize)
    width = 0.72 / len(method_order)
    for index, method in enumerate(method_order):
        method_values = summary[summary["hybridization method"].eq(method)]
        values = method_values.set_index(lookup_columns).reindex(task_index)
        means = values[mean_col].to_numpy(dtype=float)
        stds = values[std_col].to_numpy(dtype=float)
        bar_positions = positions + (index - (len(method_order) - 1) / 2) * width
        options = dict(
            color=palette[method],
            capsize=3,
            error_kw={"elinewidth": 1.1, "capthick": 1.1},
        )
        if show_all_targets:
            ax.barh(bar_positions, means, height=width, xerr=stds, **options)
        else:
            ax.bar(bar_positions, means, width=width, yerr=stds, **options)

    score_label = y_label or f"Mean {metric.upper()} ({evaluation})"
    if show_all_targets:
        ax.set_yticks(positions)
        ax.set_yticklabels(tick_labels)
        ax.invert_yaxis()
        ax.set_xlabel(score_label, fontsize=fontsize, fontweight="bold")
        ax.set_ylabel(
            "Target (number of datapoints)" if x_label == "Number of datapoints" else x_label,
            fontsize=fontsize, fontweight="bold",
        )
        if log_y:
            ax.set_xscale("log")
        if y_lim is not None:
            ax.set_xlim(*y_lim)
    else:
        ax.set_xticks(positions)
        ax.set_xticklabels(
            tick_labels, rotation=x_tick_rotation,
            ha="right" if x_tick_rotation else "center",
        )
        ax.set_xlabel(x_label, fontsize=fontsize, fontweight="bold")
        ax.set_ylabel(score_label, fontsize=fontsize, fontweight="bold")
        if log_y:
            ax.set_yscale("log")
        if y_lim is not None:
            ax.set_ylim(*y_lim)
    if title:
        ax.set_title(title, fontsize=fontsize + 2)
    ax.tick_params(axis="both", labelsize=fontsize - 2)
    ax.legend(
        handles=[Patch(facecolor=palette[method], label=method) for method in method_order],
        loc="lower center", bbox_to_anchor=(0.5, 1.02),
        ncol=min(4, len(method_order)), frameon=False, fontsize=max(fontsize - 3, 7),
    )
    if show_values:
        # Clip optional labels to the axes, so tight saving cannot grow the canvas.
        for container in ax.containers:
            if not hasattr(container, "patches"):
                continue
            ax.bar_label(
                container, fmt="%.2f", padding=2,
                fontsize=max(fontsize - 6, 6), clip_on=True,
            )

    try:
        fig.tight_layout()
        if save_dir is not None:
            output_dir = ensure_long_path(Path(save_dir))
            output_dir.mkdir(parents=True, exist_ok=True)
            models = "_".join(str(value) for value in (selection_values(model) or ["all"]))
            output_name = file_name or (
                f"{models}_{metric}_hybridization_performance_vs_data_number.png"
            )
            # Always qualify custom filenames as well, preventing split overwrites.
            output_path = Path(output_name)
            if not output_path.stem.lower().endswith(f"_{evaluation.lower()}"):
                output_path = output_path.with_name(
                    f"{output_path.stem}_{evaluation.lower()}{output_path.suffix}"
                )
            fig.savefig(
                ensure_long_path(output_dir / output_path),
                bbox_inches="tight", dpi=900 if high_quality else 150,
            )
        if show:
            plt.show()
    finally:
        plt.close(fig)
    return summary


def _target_filename_component(value: Any) -> str:
    """Use readable filename components without hash suffixes."""
    text = str(value)
    readable = re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("._") or "unnamed"
    return readable[:70]


def plot_model_performance_by_target(
    df: pd.DataFrame,
    metric: str = "r2",
    model: Any = None,
    kernel_triples: Any = None,
    fp_kernels: Any = None,
    count_kernels: Any = None,
    mixing_methods: Any = None,
    include_kernel_config: bool = True,
    complete_cases: bool = True,
    figsize: tuple = (9, 5),
    fontsize: int = 12,
    x_label: str = "Model",
    y_label: Optional[str] = None,
    x_tick_rotation: int = 35,
    y_lim: Optional[tuple] = None,
    show: bool = False,
    high_quality: bool = False,
    save_dir: Optional[Path] = DEFAULT_SAVE_DIR / "absolute_metric" / "model_comparison",
) -> pd.DataFrame:
    """Save one paired OOD/IID model-comparison bar plot per dataset/target.

    IID repeats are averaged within each cluster. For both evaluations, bars
    show mean +/- sample std across cluster scores within the task. Each GP
    kernel configuration gets its own bar (e.g. GP-MAP Bitwise and SK), while
    kernel-free tree/GNN/MLP models survive kernel filtering automatically.

    By default each task uses clusters shared by both evaluations and all
    configurations available for that task. Configurations without both splits
    are reported and omitted for that task; they do not remove other tasks.
    A single cluster has std=0. No value annotations or gradients are drawn.

    OOD bars use model colors and adjacent IID bars use lighter matching colors.
    Filenames encode dataset, target, metric and ood_vs_iid. Use separate
    save directories for structure and continuous clustering. Returns one row
    per plotted task/configuration/evaluation with mean, std, counts and path.
    """
    metric = _metric_name(metric)
    if y_lim is not None:
        if len(y_lim) != 2 or not np.isfinite(y_lim).all() or not y_lim[0] < y_lim[1]:
            raise ValueError("y_lim must contain two finite values with lower < upper.")

    scores, metadata = _select_scores(
        df,
        metric=metric,
        model=model,
        fp_kernels=fp_kernels,
        count_kernels=count_kernels,
        mixing_methods=mixing_methods,
        kernel_triples=kernel_triples,
        include_kernel_config=include_kernel_config,
        complete_cases=False,  # Enforce completeness independently within each task.
        ood_and_iid=None,
    )
    metadata = _ordered_metadata(metadata, model)
    configuration_columns = [
        "configuration", "model", "fp kernel", "count kernel", "mixing method",
    ]
    summaries = []
    for (paper, target), task_scores in scores.groupby(
        ["paper", "target"], sort=False, dropna=False,
    ):
        coverage = task_scores.groupby("configuration key")["evaluation"].nunique()
        keys = set(coverage.index[coverage.eq(len(EVALUATIONS))])
        present = metadata["configuration key"].map(lambda key: key in keys)
        task_metadata = metadata[present]
        if not present.all():
            warnings.warn(
                f"{paper!r}, {target!r}: missing OOD or IID {metric} scores for "
                f"{metadata.loc[~present, 'configuration'].tolist()}; omitting those bars.",
                RuntimeWarning,
                stacklevel=2,
            )
        task_scores = task_scores[
            task_scores["configuration key"].map(lambda key: key in keys)
        ].copy()
        if task_scores.empty:
            continue
        if complete_cases:
            task_scores = _keep_common_complete_clusters(
                task_scores, evaluations=EVALUATIONS,
            )
        if task_scores.empty:
            warnings.warn(
                f"Skipping {paper!r}, {target!r}: no common OOD/IID clusters "
                "remain for the selected configurations.",
                RuntimeWarning,
                stacklevel=2,
            )
            continue

        cases = task_scores.groupby(
            ["evaluation", "cluster", *configuration_columns], as_index=False, dropna=False,
        )["score"].mean()
        summary = cases.groupby(
            ["evaluation", *configuration_columns], as_index=False, dropna=False,
        )["score"].agg(mean="mean", std="std", n_cases="count")
        summary["std"] = summary["std"].fillna(0.0)
        counts = task_scores.groupby(["evaluation", "configuration"]).size()
        summary = summary.merge(
            counts.rename("n_scores").reset_index(),
            on=["evaluation", "configuration"], validate="one_to_one",
        )
        summary["dataset"] = paper
        summary["target"] = target
        summary["metric"] = metric
        summary = _add_datapoint_counts(summary)
        order = task_metadata["configuration"].tolist()
        summary = summary.assign(
            _order=summary["configuration"].map({label: i for i, label in enumerate(order)}),
            _evaluation=summary["evaluation"].map({"OOD": 0, "IID": 1}),
        ).sort_values(["_order", "_evaluation"], kind="stable").drop(
            columns=["_order", "_evaluation"],
        ).reset_index(drop=True)
        summary["output path"] = None

        fig, ax = plt.subplots(figsize=figsize)
        try:
            positions = np.arange(len(order), dtype=float)
            width = 0.33
            for evaluation, offset in (("OOD", -width / 2), ("IID", width / 2)):
                values = summary[summary["evaluation"].eq(evaluation)].set_index(
                    "configuration"
                ).reindex(order)
                colors = [model_type_color(value) for value in values["model"]]
                if evaluation == "IID":
                    colors = [lighten_color(color) for color in colors]
                ax.bar(
                    positions + offset,
                    values["mean"].to_numpy(dtype=float),
                    width=width,
                    yerr=values["std"].to_numpy(dtype=float),
                    capsize=4,
                    color=colors,
                    edgecolor="none",
                    error_kw={"elinewidth": 1.2, "capthick": 1.2},
                )
            ax.legend(
                handles=[
                    Patch(facecolor="#666666", label="OOD"),
                    Patch(facecolor=lighten_color("#666666"), label="IID"),
                ],
                frameon=False,
                title="Evaluation",
            )
            ax.set_xticks(positions)
            ax.set_xticklabels(
                order, rotation=x_tick_rotation,
                ha="right" if x_tick_rotation else "center",
            )
            ax.set_xlabel(x_label, fontsize=fontsize, fontweight="bold")
            ax.set_ylabel(
                y_label or f"Mean {metric.upper()}",
                fontsize=fontsize, fontweight="bold",
            )
            ax.tick_params(axis="both", labelsize=fontsize - 2)
            if y_lim is not None:
                ax.set_ylim(*y_lim)
            fig.tight_layout()
            if save_dir is not None:
                output_dir = ensure_long_path(Path(save_dir))
                output_dir.mkdir(parents=True, exist_ok=True)
                filename = (
                    f"{_target_filename_component(paper)}_"
                    f"{_target_filename_component(target)}_"
                    f"{_target_filename_component(metric)}_ood_vs_iid.png"
                )
                output_path = ensure_long_path(output_dir / filename)
                fig.savefig(
                    output_path, bbox_inches="tight", dpi=900 if high_quality else 150,
                )
                summary["output path"] = str(output_path)
            if show:
                plt.show()
        finally:
            plt.close(fig)
        summaries.append(summary)
    if not summaries:
        raise ValueError("No dataset/target has usable cases for the model comparison.")
    return pd.concat(summaries, ignore_index=True)


if __name__ == "__main__":
    models = [
        "RF", "XGBR", "NGB", "GPytorchMAP", "GNN-GP", "GNN-GP-SN",
        "GNN-MLP", "GNN-MLP-Ens", "MLP-GP", "MLP-GP-SN",
    ]
    kernel_triples = [
        ("TanimotoMatern32", "Matern32", "(count:x)+(fp:x)"),
        ("Matern32", "Matern32", "(count:x)+(fp:x)"),
    ]
    requested_methods = [
        # "sum",
        "product",
        "(count:+)x(fp:+)",
        "(count:x)+(fp:x)",
    ]
    for ood_type in ("structure", "continuous"):
        result_df = pd.read_pickle(
            MASTER_DATA / f"Tree_GP_MLP_{ood_type}_cluster_OOD.pkl"
        )
        metrics = "cvpp_ama"
        plot_model_performance_by_target(
            df=result_df,
            metric=metrics,
            model=models,
            kernel_triples=kernel_triples,
            y_label="AMA",
            y_lim=(0, .5),
            fontsize=15,
            figsize=(9, 5),
            show=False,
            high_quality=False,
            save_dir=(
                DEFAULT_SAVE_DIR / "absolute_metric" / "model_comparison"
                / "separate_datasets" / ood_type/f"{metrics.capitalize()}"
            ),
        )

        # for evaluation in EVALUATIONS:
            # plot_hybridization_performance_vs_data_number(
            #     df=result_df,
            #     metric="r2",
            #     model="GPytorchMAP",
            #     fp_kernels="TanimotoMatern32",
            #     count_kernels="Matern32",
            #     mixing_methods=requested_methods,
            #     ood_and_iid=evaluation,
            #     show_all_targets=True,
            #     y_label=f"R²",
            #     y_lim=(0, 1) if evaluation == "IID" else (-2, 1),
            #     fontsize=15,
            #     figsize=(11, 8),
            #     show=False,
            #     high_quality=False,
            #     save_dir=(
            #         DEFAULT_SAVE_DIR / "absolute_metric" / "hybridization_comparison"
            #     ),
            #     file_name=f"R2_{ood_type}_GPytorchMAP_SK_vs_data_number.png",
            # )
