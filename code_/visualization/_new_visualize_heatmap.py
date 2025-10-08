import json
from itertools import product
from pathlib import Path
from typing import List, Optional
import os 
# import cmcrameri.cm as cmc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
# from matplotlib import rc
from visualization_setting import set_plot_style, save_img_path


HERE: Path = Path(__file__).resolve().parent
RESULTS: Path = HERE.parent.parent/ 'results'



#----------- Helper functions for heatmap visualization -----------#
def ensure_long_path(path: Path) -> Path:
    """Ensures Windows handles long paths by adding '\\?\' if needed."""
    path_str = str(path)
    if os.name == 'nt' and len(path_str) > 250:
        return Path(f"\\\\?\\{path_str}")
    return path


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


def _create_heatmap(
    root_dir: Path,
    score: str,
    avg_scores: pd.DataFrame,
    annotations: pd.DataFrame,
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
    fontsize: int = 16,
    **kwargs,
) -> None:
    """
    Args:
        root_dir: Root directory containing all results
        score: Score to plot ("r", "r2", "mae", etc.)
        var: Variance type to annotate ("std", "sem", etc.)
        avg_scores: DataFrame of average scores
        annotations: DataFrame of annotations (e.g., standard deviations)
        figsize: Tuple of figure size (width, height)
        fig_title: Title for the plot
        x_title: Label for the x-axis
        y_title: Label for the y-axis
        fname: Filename to save the figure (without extension)
        vmin: Minimum value for color scaling
        vmax: Maximum value for color scaling
        feature_order: List specifying the desired order of features (columns)
        model_order: List specifying the desired order of ML models (rows)
    """

    # Reorder DataFrames if specific order is provided
    if feature_order is not None:
        avg_scores = avg_scores.loc[feature_order]
        annotations = annotations.loc[feature_order]
    if model_order is not None:
        avg_scores = avg_scores[model_order]
        annotations = annotations[model_order]

    # Create heatmap
    fig, ax = plt.subplots(figsize=figsize)
    palette: str = "viridis" if score in ["r", "r2"] else "viridis_r"
    custom_cmap = sns.color_palette(palette, as_cmap=True)
    custom_cmap.set_bad(color="lightgray")

    hmap = sns.heatmap(
        avg_scores,
        annot=annotations,
        fmt="",
        cmap=custom_cmap,
        cbar=True,
        vmin=vmin,
        vmax=vmax,
        ax=ax,
        mask=avg_scores.isnull(),
        annot_kws={"fontsize": fontsize},
    )

    # Set axis labels and tick labels
    ax.set_xticks(np.arange(len(avg_scores.columns)) + 0.5)
    ax.set_yticks(np.arange(len(avg_scores.index)) + 0.5)
    x_tick_labels: list[str] = avg_scores.columns.tolist()
    y_tick_labels: list[str] = avg_scores.index.tolist()

    ax.set_xticklabels(x_tick_labels, rotation=45, ha="right", fontsize=fontsize)
    ax.set_yticklabels(y_tick_labels, rotation=0, ha="right", fontsize=fontsize)

    # Set plot and axis titles
    plt.title(fig_title, fontsize=fontsize+2, fontweight='bold')
    ax.set_xlabel(x_title, fontsize=fontsize+2, fontweight='bold')
    ax.set_ylabel(y_title, fontsize=fontsize+2, fontweight='bold')

    # Set colorbar title and custom ticks
    score_txt: str = "$R^2$" if score == "r2" else score
    cbar = hmap.collections[0].colorbar

    if vmin is None:
        vmin = np.nanmin(avg_scores.values)
    if vmax is None:
        vmax = np.nanmax(avg_scores.values)

    num_ticks = num_ticks
    ticks = np.linspace(vmin, vmax, num_ticks)
    cbar.set_ticks(np.round(ticks,1))

    cbar.set_label(
        f"Average {score_txt.upper()} ± stdev",
        rotation=270,
        labelpad=20,
        fontsize=fontsize,
        fontweight='bold',
    )
    cbar.ax.tick_params(labelsize=fontsize)

    # Save the figure
    visualization_folder_path = root_dir / "heatmap"
    os.makedirs(visualization_folder_path, exist_ok=True)
    plt.tight_layout()
    plt.savefig(visualization_folder_path / f"{fname}.png", dpi=600)

    plt.show()
    plt.close()


def get_results_from_file(
    file_path: Path,
    score: str,
    peak_number:int=None,
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
        # for just scaler features
        if "scaler" == file_path.parent.name:
            model:str = file_path.name.split("_")[1] 
            # features:str = file_path.name.split("_")[0].replace("(", "", 1)[::-1].replace(")", "", 1)[::-1]
            # print(file_path.name)
            features = parse_property_string(file_path.name)
            # print(features)
        # for mixture of scaler and structural 
        elif "scaler" in file_path.parent.name and file_path.parent.name != "scaler":
            # features:str = file_path.name.split("_")[0].replace("(", "").replace(")", "")
            features = parse_property_string(file_path.name)
            # print(features)
            model:str = file_path.name.split("_")[1]

        # for structural only
        else:
            features:str = file_path.name.split("_")[0].replace("(", "").replace(")", "")
            model:str = file_path.name.split("_")[1]
            features=correct_ecfp_name(features)

       
        with open(file_path, "r") as f:
            data = json.load(f)

        avg = data[f"{score}_avg"]
        std = data[f"{score}_stdev"]
        

        avg = avg[peak_number] if isinstance(avg, list) else avg
        std = std[peak_number] if isinstance(std, list) else std
        if score in ["mae", "rmse"]:
            avg, std = abs(avg), abs(std)
        return features, model, avg, std
