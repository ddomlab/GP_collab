import json
from pathlib import Path
from typing import List, Optional, Any, Dict
import os 
import re

# import cmcrameri.cm as cmc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import krippendorff

from visualization_setting import set_plot_style, save_img_path, ensure_long_path

set_plot_style()

HERE = Path(__file__).resolve().parent
# DATASETS = HERE.parent.parent / "datasets" / "Validation datasets"
RESULTS = HERE.parent.parent / "results"


# 
#----------- Helper functions for heatmap visualization -----------#
def generate_annotations(num: float) -> str:
    """
    Args:
        num: Number to annotate

    Returns:
        String to annotate heatmap
    """
    if isinstance(num, float) and not np.isnan(num):
        num_txt: str = f"{round(num, 2)}"
    else:
        num_txt = "NaN"
    return num_txt


def plot_manual_heatmap(
    root_dir: Path,
    score_metric: str,
    score_to_show: pd.DataFrame,
    figsize: tuple[int, int],
    fig_title: str,
    x_title: str,
    y_title: str,
    fname: str,
    vmin: float = None,
    vmax: float = None,
    feature_order: list[str] = None,
    model_order: list[str] = None,
    num_ticks: int = 3,
    **kwargs,
) -> None:
    
    # def wrap_label(label: str, max_words_per_line: int = 2) -> str:
    #     words = label.split(" + ")
    #     return "\n".join(
    #         [" + ".join(words[i:i + max_words_per_line]) for i in range(0, len(words), max_words_per_line)]
    #     )

    # Pivot the DataFrame
    avg_scores = score_to_show.pivot(index='features', columns='model', values='score')
    annotations = score_to_show.pivot(index='features', columns='model', values='annotations')

    if feature_order is not None:
        avg_scores = avg_scores.reindex(feature_order)
        annotations = annotations.reindex(feature_order)

    if model_order is not None:
        avg_scores = avg_scores[model_order]
        annotations = annotations[model_order]

    fig, ax = plt.subplots(figsize=figsize)

    palette = "viridis" if score_metric in ["r", "r2"] else "viridis_r"
    custom_cmap = sns.color_palette(palette, as_cmap=True)
    custom_cmap.set_bad(color="lightgray")

    # Create heatmap without automatic colorbar
    hmap = sns.heatmap(
        avg_scores,
        annot=annotations,
        fmt="",
        cmap=custom_cmap,
        cbar=False,
        vmin=vmin,
        vmax=vmax,
        ax=ax,
        mask=avg_scores.isnull(),
        annot_kws={"fontsize": kwargs['fontsize']},
    )

    # Flip y-ticks to the right side
    ax.yaxis.tick_right()

    # Set tick labels
    ax.set_xticks(np.arange(len(avg_scores.columns)) + 0.5)
    ax.set_yticks(np.arange(len(avg_scores.index)) + 0.5)

    x_tick_labels = avg_scores.columns.tolist()
    # y_tick_labels = [wrap_label(label) for label in avg_scores.index]

    ax.set_xticklabels(x_tick_labels, rotation=45, ha="right", fontsize=kwargs['fontsize'])
    ax.set_yticklabels(avg_scores.index, rotation=0, ha="left", fontsize=kwargs['fontsize'])
    ax.tick_params(axis='y', pad=5)

    # Titles
    plt.title(fig_title, fontsize=kwargs['fontsize'] + 2, fontweight='bold')
    ax.set_xlabel(x_title, fontsize=kwargs['fontsize'] + 2, fontweight='bold')
    ax.set_ylabel(y_title, fontsize=kwargs['fontsize'] + 2, fontweight='bold')

    # Colorbar manually on the left
    var_titles = {"stdev": "Stdev"}
    score_txt = "$R^2$" if score_metric == "r2" else score_metric

    if vmin is None:
        vmin = np.nanmin(avg_scores.values)
    if vmax is None:
        vmax = np.nanmax(avg_scores.values)

    ticks = np.linspace(vmin, vmax, num_ticks)

    # Get the QuadMesh from heatmap and create colorbar manually
    im = hmap.get_children()[0]
    cbar = fig.colorbar(im, ax=ax, orientation='vertical', location='left', shrink=0.6, pad=0.02)
    cbar.set_ticks(np.round(ticks, 1))
    cbar.set_label(
        f"Average {score_txt.upper()} ± {var_titles.get('stdev', 'Stdev')}",
        rotation=90,
        labelpad=10,
        fontsize=kwargs['fontsize']+2,
        fontweight='bold',
    )
    cbar.ax.tick_params(labelsize=kwargs['fontsize'])

    visualization_folder_path = root_dir / "heatmap"
    # os.makedirs(visualization_folder_path, exist_ok=True)
    plt.tight_layout()
    save_img_path(visualization_folder_path, fname)
    plt.show()
    plt.close()


def _parse_property_string(file_name: str):

    feat_match = re.search(r"\((.*?)\)", file_name)
    if not feat_match:
        return None
    feature_str = feat_match.group(1)
    features = re.sub(r"[-]+", "+", feature_str)  # replace '-' with '+'


    # extract model name (e.g. RF, XGBR, sklearn-GPR, etc.)
    model_match = re.search(r"\)_(.*?)_", file_name)
    model = model_match.group(1) if model_match else None

    imp_match = re.search(rf"{model}_(\w+)_hypOFF", file_name)
    imputation = imp_match.group(1) if imp_match else None

    # Exclude false matches like 'hypOFF' itself
    if imputation and imputation.lower() not in ["hypoff"]:
        features = f"{features} ({imputation} imp)"

    return features, model


def get_results_from_file(
    file_path: Path,
    score_metric: str,
    # impute: bool = False,
) :
    """
    Args:
        root_dir: Root directory containing all results
        representation: Representation for which to get scores
        model: Model for which to get scores.
        score: Score to plot
        var: Variance to plot

    Returns:
        Average and variance of score
    """
    file_path = ensure_long_path(file_path)
    if not file_path.exists():
            print('not exists')
            features, model = None, None
            avg, std = np.nan, np.nan
        
    else:
        features, model = _parse_property_string(file_path.name)
         # Read JSON file
       
        with open(file_path, "r") as f:
            data = json.load(f)

        avg = data[f"{score_metric}_avg"]
        std = data[f"{score_metric}_stdev"]
        
        if score_metric in ["mae", "rmse"]:
            avg, std = abs(avg), abs(std)
        return features, model, avg, std



def _get_count_fp_comparison(
                            target_folder: Path,
                            scoring_metric: str,
                            selected_models: set[str],
                            features_to_draw: set[str] = None,
                            ) -> pd.DataFrame:

    scores_to_report: List = []
    pattern: str = "*_scores.json"

    # for value in comparison_value:
    # value_folder = os.path.join(target_folder, value)
    score_files = list(Path(target_folder).rglob(pattern))
    for score_path in score_files:
        if "generalizability" in score_path.name or "test" in score_path.parts or 'lc_scores' in score_path.name:
            continue

        feats, model, av, std = get_results_from_file(file_path=score_path, score_metric=scoring_metric)
        # Only keep selected features
        if features_to_draw is None:
            # draw all features
            pass
        elif feats not in features_to_draw:
            continue

        if model not in selected_models and selected_models is not None:
            continue

        anot = f"{np.round(av, 2)}\n±{np.round(std, 2)}"
        scores_to_report.append({
            "features": feats,
            "model": model,
            "score": np.round(av, 2),
            "annotations": anot
        })

    return pd.DataFrame(scores_to_report)


def creat_count_fp_heatmap(
                            target_dir:Path,
                            score_metric:str,
                            features_to_draw: List[str] = None,
                            models_to_draw: set[str] = {"RF","XGBR"},
                            ) -> None:
    scores_to_show:pd.DataFrame = _get_count_fp_comparison(
                                                            target_folder=target_dir,
                                                            scoring_metric=score_metric,
                                                            selected_models=models_to_draw,
                                                            features_to_draw=features_to_draw
                                                            )
    print(scores_to_show)
    fname= f"model vs features ({score_metric})"
    if score_metric == "r2":
        vmax= 1
        vmin= 0
        n_cbar_tick = 6 
    elif score_metric == "mae":
        vmax= .5
        vmin= 0.1
        n_cbar_tick = 5  
    elif score_metric == "rmse":
        vmax= 1
        vmin= 0
        n_cbar_tick = 6 
    plot_manual_heatmap(root_dir=target_dir/"comparison heatmap for polymer properties",
                        score_metric=score_metric,
                        score_to_show=scores_to_show,
                        figsize=(7, 7),
                        fig_title=f" \n ",
                        x_title="Models",
                        y_title="",
                        fname=fname,
                        vmin=vmin,
                        vmax=vmax,
                        # feature_order=,
                        # model_order=['RF','DT','MLR'],
                        num_ticks=n_cbar_tick,
                        fontsize=16,
                        )



def plot_average_feature_importances(scores_data: Dict[str, Any], save_loc: Path, file_extension: str, figsize: tuple=(10,7)) -> None:
    all_model_fi: List[Dict[str, float]] = []
    all_shap_fi: List[Dict[str, float]] = []

    # Collect model importances
    if isinstance(scores_data.get("feature_importance_model"), list):
        all_model_fi.extend(scores_data["feature_importance_model"])
    else:
        for value in scores_data.values():
            if isinstance(value, dict) and isinstance(value.get("feature_importance_model"), list):
                all_model_fi.extend(value["feature_importance_model"])

    # Collect SHAP importances
    if isinstance(scores_data.get("feature_importance_SHAP"), list):
        all_shap_fi.extend(scores_data["feature_importance_SHAP"])
    else:
        for value in scores_data.values():
            if isinstance(value, dict) and isinstance(value.get("feature_importance_SHAP"), list):
                all_shap_fi.extend(value["feature_importance_SHAP"])

    # Build DataFrames
    df_model = pd.DataFrame(all_model_fi) if all_model_fi else pd.DataFrame()
    df_shap = pd.DataFrame(all_shap_fi) if all_shap_fi else pd.DataFrame()
    if df_model.empty and df_shap.empty:
        raise ValueError("No feature importance data found in scores_data.")

    # Compute means and stds
    mean_model, std_model, mean_shap, std_shap = None, None, None, None
    if not df_model.empty:
        mean_model = df_model.mean()
        std_model = df_model.std()
    if not df_shap.empty:
        mean_shap = df_shap.mean()
        std_shap = df_shap.std()

    # Determine feature order (prefer model if available)
    if mean_model is not None:
        features = mean_model.sort_values(ascending=False).index
    else:
        features = mean_shap.sort_values(ascending=False).index

    # Align both importances
    if mean_model is not None:
        mean_model = mean_model.reindex(features, fill_value=0)
        std_model = std_model.reindex(features, fill_value=0)
    if mean_shap is not None:
        mean_shap = mean_shap.reindex(features, fill_value=0)
        std_shap = std_shap.reindex(features, fill_value=0)

    # ---- Select top 15 most important features ----
    if mean_model is not None and mean_shap is not None:
        combined_mean = (mean_model + mean_shap) / 2
    elif mean_model is not None:
        combined_mean = mean_model
    else:
        combined_mean = mean_shap

    top_features = combined_mean.sort_values(ascending=False).head(15).index

    # Restrict to top 15
    mean_model = mean_model[top_features] if mean_model is not None else None
    std_model = std_model[top_features] if std_model is not None else None
    mean_shap = mean_shap[top_features] if mean_shap is not None else None
    std_shap = std_shap[top_features] if std_shap is not None else None
    features = top_features
    # ------------------------------------------------------------

    # ----- NEW VIOLIN + SWARM PLOT SECTION -----
    fig, ax = plt.subplots(figsize=figsize)

    # Build long-format dataframe
    plot_data = []

    if not df_model.empty:
        for feat in features:
            for val in df_model[feat].dropna().values:
                plot_data.append([feat, val, "Model Importance"])

    if not df_shap.empty:
        for feat in features:
            for val in df_shap[feat].dropna().values:
                plot_data.append([feat, val, "SHAP Importance"])

    plot_df = pd.DataFrame(plot_data, columns=["Feature", "Value", "Type"])

    # Colors: red for model, blue for shap
    palette = {
        "Model Importance": "#BB3A5A",  # red
        "SHAP Importance": "#1f77b4",   # blue
    }

    # VIOLIN PLOTS (side-by-side)
    sns.violinplot(
        data=plot_df,
        x="Feature",
        y="Value",
        hue="Type",
        palette=palette,
        inner=None,
        linewidth=1.2,
        cut=0,
        dodge=True,
        ax=ax,
    )

    # SWARMPLOTS
    sns.swarmplot(
        data=plot_df,
        x="Feature",
        y="Value",
        hue="Type",
        dodge=True,
        palette=palette, 
        size=3,
        alpha=0.7,
        ax=ax,
    )

    # Fix duplicate legends (remove swarm's entries)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], frameon=False)

    ax.set_xticklabels(features, rotation=75, ha="right")
    ax.set_ylabel("Feature Importance")
    ax.set_xlabel("Top 15 Features")
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    plt.tight_layout()
    save_img_path(save_loc / "feature importance", f"feature_importance_top15_{file_extension}.png")
    # plt.show()
    plt.close()
    return df_model, df_shap







def krippendorff_alpha_by_feature(df, save_loc, file_extension, n_seeds=7, folds_per_seed=5, figsize=(12,6)):
    total_needed = n_seeds * folds_per_seed

    if len(df) < total_needed:
        raise ValueError(f"Not enough rows: need {total_needed}, but got {len(df)}")

    # Ensure we use exactly N rows (or you can shuffle before slicing)
    df_cut = df.iloc[:total_needed]

    alphas = {}

    # Loop over all features
    for feature in df_cut.columns:
        values = df_cut[feature].values

        # Reshape rows into (n_seeds × folds_per_seed)
        ratings = values.reshape(n_seeds, folds_per_seed)

        # Compute alpha
        try:
            alpha = krippendorff.alpha(
                reliability_data=ratings,
                level_of_measurement='interval'
            )
        except Exception:
            alpha = np.nan

        alphas[feature] = alpha

    # Convert to DataFrame
    alphas_df = pd.DataFrame({
        "Feature": list(alphas.keys()),
        "Alpha": list(alphas.values())
    }).sort_values("Alpha", ascending=False)

    # ---- Plot bar chart ----
    plt.figure(figsize=figsize)
    plt.bar(alphas_df["Feature"], alphas_df["Alpha"], color="#0b81a5")
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Krippendorff’s α")
    # plt.title("Krippendorff’s Alpha for Each Feature")
    plt.tight_layout()
    save_img_path(save_loc / "feature importance", f"feature_krippendorff_stability_{file_extension}.png")
    # plt.show()
    plt.close()
    return alphas_df





if __name__ == "__main__":

    PAPER = {
            "Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation": ["target_log Rg (nm)"],
            # "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices": ["target_calculated PCE (%)"],
            # "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery": ["target_log (Separation factor)","target_log (Total flux)"],
            # "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation": ["target_log (Separation factor)","target_log (Total flux)"],
            # "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning": [
            # "target_flux decline ratio (%)",
            # "target_flux recovery ratio (%)",
            # "target_irreversible fouling ratio(%)",
            # "target_organic compound removal (%)",
            # "target_reversible fouling ratio (%)",
            # r"target_water permeability (LMH\bar)",
            # ],
            }
    

    for paper_name, target_list in PAPER.items():
        for target in target_list:
    #         print(paper_name, target)
    #         creat_count_fp_heatmap(
    #                                 target_dir=RESULTS/paper_name/target,
    #                                 score_metric='r2',
    #                                 # comparison_value=['scaler', 'Trimer_scaler'],
    #                                 )

            paper_loc: Path = RESULTS / paper_name / target
            file_name = "(ECFP3_count_512-COUNT)_XGBR_hypOFF_Standard_Standard_scores"
            score_path = ensure_long_path(paper_loc / f"{file_name}.json")
            with open(score_path, "r") as f:
                scores = json.load(f)

            rf_imp, shap_imp = plot_average_feature_importances(scores_data=scores,
                                            save_loc=paper_loc,
                                            file_extension=file_name,
                                            figsize=(8,7.5)
                                            )
            feature_means = shap_imp.abs().mean()

            top15_features = feature_means.sort_values(ascending=False).head(15).index

            # 3. Filter the DataFrame to these 15 features
            df_top15 = shap_imp[top15_features]
            # print(df_top15)
            # plot_top15_feature_stability(
            #                     scores_data=scores,
            #                     # save_loc=paper_loc,
            #                     # file_extension=file_name,
            #                     # top_n=15,
            #                     # figsize=(8,6)
            #                     )
            krippendorff_alpha_by_feature(
                                        df=df_top15,             
                                        save_loc=paper_loc,
                                        file_extension=file_name,
                                        figsize=(9,6)
                                        )