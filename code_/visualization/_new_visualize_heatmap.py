import json
from pathlib import Path
from typing import List, Optional
import os 
import re
# import cmcrameri.cm as cmc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
# from matplotlib import rc
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

if __name__ == "__main__":

    PAPER = {
            # "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices": ["target_calculated PCE (%)"],
            "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery": ["target_log (Separation factor)","target_log (Total flux)"],
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
            print(paper_name, target)
            creat_count_fp_heatmap(
                                    target_dir=RESULTS/paper_name/target,
                                    score_metric='r2',
                                    # comparison_value=['scaler', 'Trimer_scaler'],
                                    )