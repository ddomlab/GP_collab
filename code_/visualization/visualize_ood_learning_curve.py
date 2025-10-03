import json
from pathlib import Path
from typing import List, Optional, Dict, Tuple
import os 
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns
from visualization_setting import set_plot_style, save_img_path
from utils_uncertainty_calibration import compute_all_uncertainty_metrics
set_plot_style()

#--------------HELPER FUNCTIONS----------------#

def get_score(scores: Dict, cluster: str, score_metric: str) -> float:
    """
    Helper function to get the mean or std score for a given cluster and score type.
    Returns 0 if not available.
    """
    mean_key = f"test_{score_metric}_mean"
    std_key = f"test_{score_metric}_std"
    return (
        scores.get(cluster, {}).get("summary_stats", {}).get(mean_key, 0),
        scores.get(cluster, {}).get("summary_stats", {}).get(std_key, 0),
    )


def ensure_long_path(path):
        """Ensures Windows handles long paths by adding '\\?\' if needed."""
        path_str = str(path)
        if os.name == 'nt' and len(path_str) > 250:  
            return Path(f"\\\\?\\{path_str}")
        return path


def get_uncertainty_in_learning_curve(pred_file: Dict, method: Optional[str] = None) -> pd.DataFrame:
    results = []

    for cluster, ratios in pred_file.items():
        if not cluster.startswith("CO_"):
            continue

        cluster_size = ratios.get("training size", 0)
        y_true = np.array(ratios.get("y_true", []))

        for ratio_key, seeds in ratios.items():
            if ratio_key in {"y_true", "training size"}:
                continue

            train_ratio = float(ratio_key.replace("ratio_", ""))
            train_set_size = round(train_ratio * cluster_size)

            all_metrics = []

            for seed, predictions in seeds.items():
                y_pred = np.array(predictions.get("y_test_pred", []))
                y_uncertainty = np.array(predictions.get("y_test_uncertainty", []))

                metrics = compute_all_uncertainty_metrics(y_true, y_pred, y_uncertainty, method=method)
                all_metrics.append(metrics)

            # Aggregate
            if method:
                metric_values = np.array(all_metrics)
                results.append({
                    "Cluster": cluster,
                    "Train Set Size": train_set_size,
                    f"{method} mean": np.mean(metric_values),
                    f"{method} std": np.std(metric_values),
                })
            else:
                # Compute mean/std for each metric in the returned dictionary
                aggregated = {"Cluster": cluster, "Train Set Size": train_set_size}
                keys = all_metrics[0].keys()
                for k in keys:
                    vals = [m[k] for m in all_metrics]
                    aggregated[f"{k} mean"] = np.mean(vals)
                    aggregated[f"{k} std"] = np.std(vals)
                results.append(aggregated)

    return pd.DataFrame(results)


def process_learning_curve_scores(summary_scores,
                                   metric="rmse") -> Tuple[pd.DataFrame, pd.DataFrame]:
    data = []
    cluster_train_sizes = {}

    for cluster, ratios in summary_scores.items():
        if cluster.startswith("CO_") or cluster.startswith("IID_"):
            cluster_size = ratios.get("training size", 0)
            cluster_train_sizes[cluster] = set()
            for ratio, stats in ratios.items():
                if ratio == "training size":
                    continue
                train_ratio = float(ratio.replace("ratio_", ""))
                train_set_size = round(train_ratio * cluster_size, 1)
                cluster_train_sizes[cluster].add(train_set_size)

                metric_test = f"test_{metric}_mean"
                metric_train = f"train_{metric}_mean"
                std_test = f"test_{metric}_std"
                std_train = f"train_{metric}_std"

                if metric_test in stats["test_summary_stats"] and metric_train in stats["train_summary_stats"]:
                    data.append({
                        "Cluster": cluster,
                        "Train Ratio": train_ratio,
                        "Train Set Size": train_set_size,
                        "Score": stats["test_summary_stats"][metric_test],
                        "Std": stats["test_summary_stats"][std_test],
                        "Score Type": "Test"
                    })
                    data.append({
                        "Cluster": cluster,
                        "Train Ratio": train_ratio,
                        "Train Set Size": train_set_size,
                        "Score": stats["train_summary_stats"][metric_train],
                        "Std": stats["train_summary_stats"][std_train],
                        "Score Type": "Train"
                    })

    full_df = pd.DataFrame(data)
    # Determine common train sizes from CO clusters only
    co_only_sizes = {k: v for k, v in cluster_train_sizes.items() if k.startswith("CO_")}
    if co_only_sizes:
        common_train_sizes = set.intersection(*co_only_sizes.values())
    else:
        common_train_sizes = set()

    equal_train_size_df = full_df[full_df["Train Set Size"].isin(common_train_sizes)].copy()

    return full_df, equal_train_size_df


def process_learning_curve_distances(all_distances: Dict[str, Dict]) -> pd.DataFrame:
    """
    Convert raw learning curve distance results into a DataFrame with mean/std across seeds.
    
    Handles OOD (flat dict), IID (nested dict), and edge cases (direct value).
    """
    records = []
    for cluster_key, cluster_dict in all_distances.items():
        # skip metadata
        train_size = cluster_dict.get("training size", None)

        for ratio_key, seeds_dict in cluster_dict.items():
            if ratio_key == "training size":
                continue

            # Handle ratio key type
            if isinstance(ratio_key, str) and ratio_key.startswith("ratio_"):
                train_ratio = float(ratio_key.replace("ratio_", ""))
            else:
                train_ratio = float(ratio_key)

            values = []

            if isinstance(seeds_dict, dict):
                if cluster_key.startswith("IID_"):
                    # Nested dict: seed -> test_seed -> value
                    for _, test_seeds in seeds_dict.items():
                        if isinstance(test_seeds, dict):
                            values.extend(test_seeds.values())
                        else:  # already a value
                            values.append(test_seeds)
                else:  
                    # OOD: dict of {seed: value}
                    values.extend(seeds_dict.values())
            else:
                # Edge case: already a single value
                values.append(seeds_dict)

            if not values:
                continue

            mean_val = np.mean(values)
            std_val = np.std(values)
            effective_size = int(round(train_ratio * train_size))
            records.append({
                "Cluster": cluster_key,
                "Prefix": "CO" if cluster_key.startswith("CO_") else "IID",
                "Train Ratio": round(train_ratio,2),
                "Train Set Size": effective_size,
                "Mean": mean_val,
                "Std": std_val
            })
    
    return pd.DataFrame(records)



#--------------PLOT FUNCTIONS----------------#

def plot_bar_ood_iid(data: pd.DataFrame, ml_score_metric: str, 
                    saving_path, file_name: str,
                    figsize=(12, 7), text_size=14,
                    ncol=3, title:Optional[str]=None,
                    cluster_order:Optional[list]=None) -> None:
    
    def parse_cluster_info(cluster_label):
        if cluster_label.startswith("CO_"):
            return cluster_label.split("_")[1], "OOD"
        elif cluster_label.startswith("ID_") or cluster_label.startswith("IID_"):
            return cluster_label.split("_")[1], "IID"
        else:
            return cluster_label, "Unknown"

    parsed_info = data['Cluster'].apply(parse_cluster_info)
    data['BaseCluster'] = parsed_info.apply(lambda x: x[0])
    data['SplitType'] = parsed_info.apply(lambda x: x[1])

    # Use given order or fallback to unique values
    if cluster_order is None:
        clusters = data['BaseCluster'].unique()
    else:
        clusters = cluster_order

    models = data['Model'].unique()
    palette = sns.color_palette("Set2", n_colors=len(models))
    model_colors = {model: palette[i] for i, model in enumerate(models)}

    bar_width = 0.45
    n_models = len(models)
    group_width = n_models * 2 * bar_width + 0.5

    fig, ax = plt.subplots(figsize=figsize)
    cluster_ticks = []
    cluster_labels = []

    for i, cluster in enumerate(clusters):
        cluster_offset = i * group_width
        for j, model in enumerate(models):
            subset = data[(data["BaseCluster"] == cluster) & (data["Model"] == model)]
            if subset.empty:
                continue

            base_x = cluster_offset + j * 2 * bar_width
            color = model_colors[model]

            for _, row in subset.iterrows():
                if row["SplitType"] == "OOD":
                    ax.bar(base_x, row["Score"], bar_width, yerr=row["Std"],
                           error_kw={'elinewidth': 1, 'capthick': 1}, capsize=3,
                           color=color, alpha=1, label=f"{model} OOD")
                elif row["SplitType"] == "IID":
                    ax.bar(base_x + bar_width, row["Score"], bar_width, yerr=row["Std"],
                           error_kw={'elinewidth': 1, 'capthick': 1}, capsize=3,
                           color=color, alpha=.5, label=f"{model} IID")

        center_of_group = cluster_offset + (n_models * bar_width)
        cluster_ticks.append(center_of_group)
        cluster_labels.append(cluster)

    # De-duplicate legend
    handles, labels = ax.get_legend_handles_labels()
    seen = {l: h for h, l in zip(handles, labels)}
    ax.legend(seen.values(), seen.keys(), fontsize=text_size - 2,
              frameon=True, loc='upper left', ncol=ncol)

    # Y axis
    max_score = np.nanmax(data["Score"].values)
    ymax = np.ceil(max_score / .2) * .2
    ax.set_ylim(0, 1.2)
    ax.set_yticks(np.arange(0, 1.4 , .2))

    # X axis
    ax.set_xlabel("Cluster", fontsize=text_size, fontweight='bold')
    ax.set_xticks(cluster_ticks)
    ax.set_xticklabels(cluster_labels, fontsize=text_size, rotation=45)

    ax.set_ylabel(f'{ml_score_metric.upper()}', fontsize=text_size, fontweight='bold')
    ax.tick_params(axis='y', labelsize=text_size)

    plt.tight_layout()
    save_img_path(saving_path, f"{file_name}.png")
    plt.show()
    plt.close()


def plot_ood_learning_uncertainty(
    model_predictions: Dict[str, Dict],
    folder: Path = None,
    file_name: str = 'model_comparison',
    title: Optional[str] = None,
    fontsize: int = 18,
    uncertainty_methods: Tuple[str, str] = ("RUSC", "Cv")
) -> None:
    method1, method2 = uncertainty_methods
    color_1 = '#E74747'  # For method1
    color_2 = '#367CE2'  # For method2

    model_styles = {
        'RF': '-',     # solid
        'XGBR': '--',  # dashed
        'NGB': '-.'    # dash-dot (optional, you can use '--' if only 2 models)
    }

    model_markers = {'RF': 'o', 'XGBR': '^', 'NGB': 'D'}

    all_df_1, all_df_2 = [], []

    for model_name, prediction in model_predictions.items():
        df_1 = get_uncertainty_in_learning_curve(prediction, method=method1)
        df_1["Model"] = model_name
        all_df_1.append(df_1)

        df_2 = get_uncertainty_in_learning_curve(prediction, method=method2)
        df_2["Model"] = model_name
        all_df_2.append(df_2)

    df_1 = pd.concat(all_df_1, ignore_index=True)
    df_2 = pd.concat(all_df_2, ignore_index=True)

    co_clusters = sorted([c for c in df_1['Cluster'].unique() if c.startswith("CO_")])
    col_wrap = 4 if len(co_clusters) > 4 else len(co_clusters)

    g = sns.FacetGrid(df_1[df_1['Cluster'].isin(co_clusters)],
                      col="Cluster", col_wrap=col_wrap, height=4, sharey=False)
    twin_axes = []

    def plot_model_uncertainties(data, **kwargs):
        cluster = data["Cluster"].iloc[0]
        ax = plt.gca()
        ax2 = ax.twinx()
        twin_axes.append((ax, ax2))

        for model_name in data["Model"].unique():
            marker = model_markers.get(model_name, 'o')
            line_style = model_styles.get(model_name, '-')

            data_1 = df_1[(df_1["Cluster"] == cluster) & (df_1["Model"] == model_name)].sort_values("Train Set Size")
            data_2 = df_2[(df_2["Cluster"] == cluster) & (df_2["Model"] == model_name)].sort_values("Train Set Size")

            ax.plot(
                data_1["Train Set Size"],
                data_1[f"{method1} mean"],
                label=None,
                marker=marker,
                linewidth=2.5,
                linestyle=line_style,
                color=color_1,
                markersize=8
            )
            ax.fill_between(
                data_1["Train Set Size"],
                data_1[f"{method1} mean"] - data_1[f"{method1} std"],
                data_1[f"{method1} mean"] + data_1[f"{method1} std"],
                alpha=0.2,
                color=color_1,
            )

            ax2.plot(
                data_2["Train Set Size"],
                data_2[f"{method2} mean"],
                label=None,
                marker=marker,
                linewidth=2.5,
                linestyle=line_style,
                color=color_2,
                markersize=8
            )
            ax2.fill_between(
                data_2["Train Set Size"],
                data_2[f"{method2} mean"] - data_2[f"{method2} std"],
                data_2[f"{method2} mean"] + data_2[f"{method2} std"],
                alpha=0.2,
                color=color_2
            )

        xticks = np.arange(0, 251, 50)
        ax.set_xticks(xticks)
        ax.set_ylim(-1, 1)
        ax.set_yticks(np.arange(-1, 1.5, 0.5))
        ax2.set_ylim(0, .5)
        ax2.set_yticks(np.arange(0, .6, 0.1))

        ax.set_xlabel("Training Set Size", fontsize=fontsize, fontweight='bold')
        ax.tick_params(axis='both', labelsize=fontsize - 2)
        ax2.tick_params(axis='y', labelsize=fontsize - 2)

    g.map_dataframe(plot_model_uncertainties)

    for i, (ax, ax2) in enumerate(twin_axes):
        if i == 0:
            ax.set_ylabel(method1, fontsize=fontsize, color=color_1, fontweight='bold')
        else:
            ax.set_ylabel("")
            ax.set_yticklabels([])

        if i == len(twin_axes) - 1:
            ax2.set_ylabel(method2, fontsize=fontsize, color=color_2, fontweight='bold')
        else:
            ax2.set_ylabel("")
            ax2.set_yticklabels([])
        
        raw_title = g.col_names[i]  
        clean_title = raw_title.replace("CO_", "Cluster ")
        ax.set_title(clean_title, fontsize=fontsize + 4)

    # Legend for models (using linestyle instead of color)
    legend_lines = [
        Line2D([0], [0], color='black', marker=model_markers[model], linestyle=model_styles[model], markersize=8,label=model)
        for model in df_1["Model"].unique()
    ]

    if title:
        g.fig.suptitle(title, fontsize=fontsize + 6, fontweight='bold', y=1.12)

    plt.figlegend(
        handles=legend_lines,
        loc='upper center',
        bbox_to_anchor=(0.5, 1.12),
        ncol=len(legend_lines),
        fontsize=fontsize,
        frameon=True,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.98])

    if folder:
        save_img_path(folder, f"uncertainty comparison ({file_name}).png")

    plt.show()
    plt.close()



def plot_ood_learning_accuracy_only(
    all_scores: Dict[str, Dict],
    metric: str = "rmse",
    folder: Path = None,
    file_name: str = 'model_accuracy_comparison',
    title: Optional[str] = None,
    fontsize: int = 14
) -> None:
    # Custom colors
    model_colors = {
        'RF': "#D32A41",    # orange
        'XGBR': "#6190bb",  # green
        'NGB': "#7570b3",   # purple as fallback
    }

    # Line/marker styles
    domain_styles = {
        "CO_": "-",      # OOD solid
        "IID_": "--",    # IID dashed
    }
    domain_markers = {
        "CO_": 'o',      # OOD circle
        "IID_": '^',     # IID triangle
    }

    all_dfs = []
    for model_name, summary_scores in all_scores.items():
        score_df, _ = process_learning_curve_scores(summary_scores, metric)
        score_df["Model"] = model_name
        all_dfs.append(score_df)

    score_df = pd.concat(all_dfs, ignore_index=True)

    if score_df.empty:
        raise ValueError(f"No data found for metric '{metric}'.")

    # Keep only Test scores
    score_df = score_df[score_df["Score Type"] == "Test"]

    max_score = np.ceil(score_df["Score"].max() * 10) / 10
    score_yticks = np.arange(0, max_score + 0.4, 0.4)

    co_clusters = sorted(set(c.replace("CO_", "") for c in score_df['Cluster'].unique() if c.startswith("CO_")))
    col_wrap = 4 if len(co_clusters) > 4 else len(co_clusters)

    g = sns.FacetGrid(score_df[score_df['Cluster'].str.startswith("CO_")],
                      col="Cluster", col_wrap=col_wrap, height=4, sharey=False)

    def plot_accuracy(data, **kwargs):
        cluster_oid = data["Cluster"].iloc[0]
        cluster_idx = cluster_oid.replace("CO_", "")
        ax = plt.gca()

        for model in data["Model"].unique():
            for prefix in ["CO_", "IID_"]:  # OOD and IID
                cluster_name = f"{prefix}{cluster_idx}"
                sub_df = score_df[
                    (score_df["Cluster"] == cluster_name) &
                    (score_df["Model"] == model)
                ]
                if sub_df.empty:
                    continue

                sub_df = sub_df.sort_values("Train Set Size")
                color = model_colors.get(model, 'black')
                linestyle = domain_styles[prefix]
                marker = domain_markers[prefix]

                ax.plot(
                    sub_df["Train Set Size"],
                    sub_df["Score"],
                    label=None,
                    marker=marker,
                    linestyle=linestyle,
                    linewidth=2.5,
                    color=color,
                    markersize=8
                )
                ax.fill_between(
                    sub_df["Train Set Size"],
                    sub_df["Score"] - sub_df["Std"],
                    sub_df["Score"] + sub_df["Std"],
                    alpha=0.15,
                    color=color
                )

        ax.set_xticks(np.arange(0, 251, 50))
        ax.set_ylim(0, max_score)
        ax.set_yticks(score_yticks)
        ax.tick_params(axis='x', labelsize=fontsize-1)
        ax.tick_params(axis='y', labelsize=fontsize-1)
        ax.set_xlabel("Training Set Size", fontsize=fontsize + 1, fontweight='bold')

        for spine in ['top', 'right', 'bottom', 'left']:
            ax.spines[spine].set_visible(True)

    g.map_dataframe(plot_accuracy)

    for i, ax in enumerate(g.axes.flat):
        if i == 0:
            ax.set_ylabel(metric.upper(), fontsize=fontsize + 1, fontweight='bold')
        else:
            ax.set_ylabel("")
            ax.set_yticklabels([])

        raw_title = g.col_names[i]
        clean_title = raw_title.replace("CO_", "Cluster ")
        ax.set_title(clean_title, fontsize=fontsize + 4)

    # Legend: models + domain
    color_legend = [
        Line2D([0], [0], color=color, lw=3, label=model)
        for model, color in model_colors.items()
        if model in score_df["Model"].unique()
    ]
    domain_legend = [
        Line2D([0], [0], color='black', lw=2.5, linestyle=ls, marker=mk, markersize=8, label=label)
        for label, (ls, mk) in zip(["OOD Test", "IID Test"], [("-", "o"), ("--", "^")])
    ]

    if title:
        g.fig.suptitle(title, fontsize=fontsize + 10, fontweight='bold', y=1.12)

    plt.figlegend(
        handles=color_legend + domain_legend,
        loc='upper center',
        bbox_to_anchor=(0.5, 1.11),
        ncol=len(color_legend + domain_legend),
        fontsize=fontsize,
        frameon=True
    )

    plt.tight_layout(rect=[0, 0, 1, 0.98])

    if folder:
        save_img_path(folder, f"accuracy comparison ({file_name}).png")

    plt.show()
    plt.close()


def plot_ood_learning_distances(
    all_distances: Dict[str, Dict],
    folder: Optional[Path] = None,
    file_name: str = 'ood_distance_comparison',
    title: Optional[str] = None,
    fontsize: int = 14
) -> None:

    df = all_distances.copy()

    if df.empty:
        raise ValueError("No distance data found.")

    # Style
    set2_palette = sns.color_palette("Set2")
    domain_styles = {"CO": "-", "IID": "--"}  # OOD solid, IID dashed

    max_val = np.ceil(df["Mean"].max() * 10) / 10
    yticks = np.arange(0, 11, 2)

    # Extract clusters (only CO to anchor panels)
    co_clusters = sorted(c.replace("CO_", "") for c in df['Cluster'].unique() if c.startswith("CO_"))
    col_wrap = 4 if len(co_clusters) > 4 else len(co_clusters)

    g = sns.FacetGrid(df[df['Cluster'].str.startswith("CO_")],
                      col="Cluster", col_wrap=col_wrap, height=4, sharey=False)

    def plot_cluster(data, **kwargs):
        cluster_oid = data["Cluster"].iloc[0]
        cluster_idx = cluster_oid.replace("CO_", "")
        ax = plt.gca()

        for prefix in ["CO", "IID"]:
            cluster_name = f"{prefix}_{cluster_idx}"
            sub_df = df[df["Cluster"] == cluster_name]
            if sub_df.empty:
                continue

            sub_df = sub_df.sort_values("Train Set Size")

            linestyle = domain_styles[prefix]
            color = set2_palette[0] if prefix == "CO" else set2_palette[1]

            ax.plot(
                sub_df["Train Set Size"],
                sub_df["Mean"],
                label=None,
                marker="o",
                linestyle=linestyle,
                linewidth=2.5,
                color=color,
                markersize=8
            )
            ax.fill_between(
                sub_df["Train Set Size"],
                sub_df["Mean"] - sub_df["Std"],
                sub_df["Mean"] + sub_df["Std"],
                alpha=0.2,
                color=color
            )



        xticks = np.arange(0, 251, 50)
        ax.set_xticks(xticks)
        # ax.set_xticks(sorted(sub_df["Train Set Size"].unique()))
        ax.set_ylim(0, max_val)
        ax.set_yticks(yticks)
        ax.tick_params(axis='both', labelsize=fontsize+2)
        ax.set_xlabel("Train Set Size", fontsize=fontsize + 2, fontweight='bold')

        for spine in ['top', 'right', 'bottom', 'left']:
            ax.spines[spine].set_visible(True)

    g.map_dataframe(plot_cluster)

    # y-axis labels
    for i, ax in enumerate(g.axes.flat):
        if i == 0:
            ax.set_ylabel("Wasserstein Distance", fontsize=fontsize , fontweight='bold')
        else:
            ax.set_ylabel("")
            ax.set_yticklabels([])


        raw_title = g.col_names[i]  
        clean_title = raw_title.replace("CO_", "Cluster ")
        ax.set_title(clean_title, fontsize=fontsize + 3)
        # ax.set_title(g.col_names[i], fontsize=fontsize + 4)

    # Legend
    domain_legend = [
        Line2D([0], [0], color=set2_palette[0], lw=2.5, linestyle="-", marker="o", label="OOD"),
        Line2D([0], [0], color=set2_palette[1], lw=2.5, linestyle="--", marker="o", label="IID"),
    ]

    if title:
        g.fig.suptitle(title, fontsize=fontsize + 8, fontweight='bold', y=1.08)

    plt.figlegend(
        handles=domain_legend,
        loc='upper center',
        bbox_to_anchor=(0.5, 1.11),
        ncol=len(domain_legend),
        fontsize=fontsize,
        frameon=True
    )

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    if folder:
        save_img_path(folder, f"{file_name}.png")

    plt.show()
    plt.close()



# comparison_of_features_lc = {
#     '(Mordred-Xn-Mw-PDI-concentration-temperature-polymer dP-polymer dD-polymer dH-solvent dP-solvent dD-solvent dH)_RF_hypOFF_Standard_lc':'Xn + + Mordred ',
#     '(Mordred-Xn-Mw-PDI-concentration-temperature-polymer dP-polymer dD-polymer dH-solvent dP-solvent dD-solvent dH-light exposure-aging time-aging temperature-prep temperature-prep time)_RF_hypOFF_Standard_lc':'Mordred+continuous+aging',
#     '(Xn-Mw-PDI-concentration-temperature-polymer dP-polymer dD-polymer dH-solvent dP-solvent dD-solvent dH)_RF_hypOFF_Standard_lc':'continuous',
#     '(Xn-Mw-PDI-concentration-temperature-polymer dP-polymer dD-polymer dH-solvent dP-solvent dD-solvent dH-light exposure-aging time-aging temperature-prep temperature-prep time)_RF_hypOFF_Standard_lc':'continuous+aging',
# }

def get_comparison_of_features(model: str, suffix: str):
    return {
        # f'(concentration-temperature-solvent dP-solvent dD-solvent dH)_{model}{suffix}': 'solvent_properties + solvent_HSPs',
        # f'(concentration-temperature-solvent dP-solvent dD-solvent dH-light exposure-aging time-aging temperature-prep temperature-prep time)_{model}{suffix}': 'solvent_properties + solvent_HSPs + environmental.thermal history',
        # f'(light exposure-aging time-aging temperature-prep temperature-prep time)_{model}{suffix}': 'environmental.thermal history',
        # f'(Xn)_{model}{suffix}': 'Xn',
        # f'(Xn-Mw-PDI)_{model}{suffix}': 'Xn + polysize',
        # f'(Xn-Mw-PDI-polymer dP-polymer dD-polymer dH)_{model}{suffix}': 'Xn + polysize + polymer HSPs',
        f'(Xn-Mw-PDI-concentration-temperature-polymer dP-polymer dD-polymer dH-solvent dP-solvent dD-solvent dH-light exposure-aging time-aging temperature-prep temperature-prep time)_{model}{suffix}': 'combination of all',
        # f'(Mordred-Xn-Mw-PDI)_{model}{suffix}': 'Xn + polysize + Mordred',
        # f'(MACCS-Xn-Mw-PDI)_{model}{suffix}': 'Xn + polysize + MACCS',
        # f'(ECFP3.count.512-Xn-Mw-PDI)_{model}{suffix}': 'Xn + polysize + ECFP6.count.512',
    }



if __name__ == "__main__":
    HERE: Path = Path(__file__).resolve().parent
    results_path = HERE.parent.parent / 'results'/ 'OOD_target_log Rg (nm)'
    cluster_list = [
                    # 'HBD3 MACCS cluster',
                    # 'substructure cluster',
                    'KM4 ECFP6_Count_512bit cluster',	
                    'KM3 Mordred cluster',
                    # 'KM4 polymer_solvent HSP cluster',
                    'KM4 Mordred_Polysize cluster',
                    # 'Polymers cluster'
                    ]


    for cluster in cluster_list:

        # Plot ood_learning_uncertainty and score new
        all_model_predictions = {}
        all_model_scores = {}
        models = ['RF', 'XGBR']

        for model in models:
            scores_folder_path = results_path / cluster / 'scaler'
            score_file_lc = scores_folder_path / f'(Xn-Mw-PDI-concentration-temperature-polymer dP-polymer dD-polymer dH-solvent dP-solvent dD-solvent dH-light exposure-aging time-aging temperature-prep temperature-prep time)_{model}_hypOFF_Standard_lc_scores.json'
            predictions_file_lc = scores_folder_path / f'(Xn-Mw-PDI-concentration-temperature-polymer dP-polymer dD-polymer dH-solvent dP-solvent dD-solvent dH-light exposure-aging time-aging temperature-prep temperature-prep time)_{model}_hypOFF_Standard_lc_predictions.json'
            predictions_file_lc = ensure_long_path(predictions_file_lc)
            score_file_lc = ensure_long_path(score_file_lc)
            if not os.path.exists(score_file_lc) or not os.path.exists(predictions_file_lc):
                print(f"File not found: {predictions_file_lc}")
                continue  

            with open(score_file_lc, "r") as f:
                scores_lc = json.load(f)
            with open(predictions_file_lc, "r") as s:
                predictions_lc = json.load(s)

            all_model_predictions[model] = predictions_lc
            all_model_scores[model] = scores_lc
        saving_folder = results_path / cluster / 'comparison of models learning curve (score and uncertainty)'
        plot_ood_learning_uncertainty(
            model_predictions=all_model_predictions,
            folder=saving_folder,
            file_name=f"combination of all features",
            fontsize=20,
            uncertainty_methods=("RUSC", "AMA"),
            # title=f"Uncertainty Comparison for {file_discription}"
        )

        print("Saved model uncertainty comparison plots.")


        plot_ood_learning_accuracy_only(
            all_scores=all_model_scores,
            metric="rmse",
            folder=saving_folder,
            file_name="accuracy only combination of all features_RF",
            # title="Model Accuracy Comparison",
            fontsize=20
        )

        print("Saved model accuracy-only learning curve plot.")



            # Plot OOD vs IID barplot at the same training size for comparison of features
            
        #     _, ood_iid_eq_tarining_size_df = process_learning_curve_scores(scores_lc, metric="rmse")
        #     ood_iid_eq_tarining_size_df = ood_iid_eq_tarining_size_df[ood_iid_eq_tarining_size_df["Score Type"] == "Test"]
        #     ood_iid_eq_tarining_size_df['Model'] = file_discription
        #     all_score_eq_training_size.append(ood_iid_eq_tarining_size_df)
        # all_score_eq_training_size: pd.DataFrame = pd.concat(all_score_eq_training_size, ignore_index=True)
        # saving_folder = results_path / cluster / f'OOD-IID bar plot at equal training set (comparison of features)'
        # plot_bar_ood_iid(all_score_eq_training_size, 'rmse', 
        #                  saving_folder, file_name=f'rmse_{model}_comparison of feature',
        #                      text_size=16, figsize=(15, 8), ncol=3)

        # print("save OOD vs IID bar plot at equal training size for comparison of features")

        # plot OOD vs IID barplot at the same training size 
        # accuracy_metric = "rmse"
        # all_score_eq_training_size = []
        # for model in ['XGBR', 'RF']:

        #         scores_folder_path = results_path / cluster / 'scaler'
        #         scores_lc = scores_folder_path / f'(Xn-Mw-PDI-concentration-temperature-polymer dP-polymer dD-polymer dH-solvent dP-solvent dD-solvent dH-light exposure-aging time-aging temperature-prep temperature-prep time)_{model}_hypOFF_Standard_lc_scores.json'
        #         scores_lc = ensure_long_path(scores_lc)
        #         if not os.path.exists(scores_lc):
        #             print(f"File not found: {scores_lc}")
        #             continue  
        #         with open(scores_lc, "r") as f:
        #             scores = json.load(f)
        #         _, ood_iid_eq_tarining_size_df = process_learning_curve_scores(scores, metric=accuracy_metric)
        #         ood_iid_eq_tarining_size_df = ood_iid_eq_tarining_size_df[ood_iid_eq_tarining_size_df["Score Type"] == "Test"]
        #         ood_iid_eq_tarining_size_df['Model'] = model

        #         all_score_eq_training_size.append(ood_iid_eq_tarining_size_df)
        # all_score_eq_training_size:pd.DataFrame = pd.concat(all_score_eq_training_size, ignore_index=True)
        # saving_folder = scores_folder_path/  f'OOD-IID bar plot(equal training set)'
        
        # figsize = (16, 10) if cluster == 'Polymers cluster' else (8, 6)
        # plot_bar_ood_iid(all_score_eq_training_size, accuracy_metric, saving_folder, file_name=f'metric-{accuracy_metric}',
        #                     text_size=22,figsize=figsize, ncol=2)
        # print("save OOD vs IID bar plot at equal training size")






        









