import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path
from typing import Union

import matplotlib.pyplot as plt
import seaborn as sns

def set_plot_style(
    font_family="sans-serif",
    font_sans_serif="Arial",
    title_size=20,
    label_size=16,
    tick_size=14,
    legend_size=14,
    seaborn_style=None,
    seaborn_palette=None,
    ax=None,
    fig=None,
    x_title=None,
    y_title=None,
    fig_title=None,
    x_tick_labels=None,
    y_tick_labels=None,
    rotate_xticks=45,
    rotate_yticks=0
):
    """Configure global and optional local (Axes/Figure) font and style settings."""
    
    # Global settings
    sns.set_style(seaborn_style)  # e.g., "white", "darkgrid", etc.
    sns.set_palette(seaborn_palette)

    plt.rc("font", family=font_family)
    plt.rc("font", **{"sans-serif": [font_sans_serif]})

    plt.rc("axes", titlesize=title_size)
    plt.rc("axes", labelsize=label_size)
    plt.rc("xtick", labelsize=tick_size)
    plt.rc("ytick", labelsize=tick_size)
    plt.rc("legend", fontsize=legend_size)

    sns.set_context("notebook", rc={
        "axes.titlesize": title_size,
        "axes.labelsize": label_size,
        "xtick.labelsize": tick_size,
        "ytick.labelsize": tick_size,
        "legend.fontsize": legend_size,
    })


def save_img_path(folder_path:Union[str, Path], file_name:str)->None:
    folder_path = Path(folder_path)
    os.makedirs(folder_path, exist_ok=True)

    save_path = folder_path / file_name
    if os.name == 'nt':  # Only for Windows
        save_path = Path(f"\\\\?\\{save_path.resolve()}")

    plt.savefig(save_path, dpi=800, bbox_inches='tight')