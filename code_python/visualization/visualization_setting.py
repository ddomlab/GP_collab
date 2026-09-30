import os
from pathlib import Path
import re
from typing import Any, Mapping, Optional, Union

import matplotlib.pyplot as plt
from matplotlib.colors import to_hex, to_rgb
import numpy as np
import pandas as pd
import seaborn as sns


TREE_MODELS = {"RF", "XGBR", "NGB"}
GNN_MODELS = {"GNN", "GCN", "GAT", "GIN", "MPNN", "DMPNN"}
MODEL_TYPE_COLORS = {
    "tree": "#7093B9",
    "gp": "#E45756",
    "mgk": "#9C1E1E",
    "gnn": "#72B7B2",
}
MODEL_DISPLAY_NAMES = {
    "XGBR": "XGB",
    "GPytorchMAP": "GP-MAP",
    "GpyroHMC": "GP-HMC",
}
FEATURE_STABILITY_COLUMNS = [
    "lengthscale_kendalls_w",
    "feature_importance_MDI_kendalls_w",
    "feature_importance_SHAP_kendalls_w",
]
HIGHER_IS_BETTER_METRICS = {"r2", "oof_r2", "rusc", "oof_rusc"}

def set_plot_style(
    # font_family="sans-serif",
    # font_sans_serif="Arial",
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

    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
    # plt.rc("font", **{"sans-serif": [font_sans_serif]})

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


def ensure_long_path(path: Path) -> Path:
    """Ensures Windows handles long paths by adding '\\\\?\\' if needed."""
    path_str = str(path)
    if os.name == 'nt' and len(path_str) > 250 and not path_str.startswith('\\\\?\\'):
        return Path(f"\\\\?\\{path_str}")
    return path


def save_img_path(folder_path: Union[str, Path], file_name: str) -> None:
    """Creates the folder (with long path support) and saves a plot image safely."""
    folder_path = ensure_long_path(Path(folder_path))
    os.makedirs(folder_path, exist_ok=True)

    save_path = ensure_long_path(folder_path / file_name)
    plt.savefig(save_path, dpi=1200, bbox_inches='tight')


def selection_values(values: Any) -> Optional[list[Any]]:
    if values is None:
        return None
    if isinstance(values, str):
        return [values]
    return list(values)


def filter_selection(
    df: pd.DataFrame,
    column: str,
    values: Any,
    keep_missing: bool = False,
) -> pd.DataFrame:
    selected = selection_values(values)
    if selected is None:
        return df

    selected_lower = {str(value).lower() for value in selected}
    mask = df[column].astype(str).str.lower().isin(selected_lower)
    if keep_missing:
        mask |= df[column].isna()
    return df[mask].copy()


def display_model_name(model: Any) -> str:
    model = str(model)
    return MODEL_DISPLAY_NAMES.get(model, model)


def is_tree_model(model: Any) -> bool:
    return str(model).upper() in TREE_MODELS


def is_feature_stability_metric(metric: Any) -> bool:
    metric_key = re.sub(r"[\s\-]+", "_", str(metric).strip().lower())
    return metric_key in {
        "feature_stability",
        "feature_importance_stability",
        "feature_importance_kendalls_w",
        "feature_importance_stability_kendalls_w",
        "tree_feature_importance_stability",
    }


def metric_higher_is_better(metric: Any) -> bool:
    """Return the optimization direction used by the visualization scripts."""
    metric_text = str(metric).strip()
    return (
        metric_text.lower() in HIGHER_IS_BETTER_METRICS
        or metric_text in FEATURE_STABILITY_COLUMNS
        or is_feature_stability_metric(metric_text)
    )


def model_config_label(
    row: pd.Series,
    include_kernel_config: bool = True,
) -> str:
    raw_model = str(row["model"])
    model = display_model_name(raw_model)
    if not include_kernel_config or raw_model == "MGK":
        return model
    if pd.isna(row["fp kernel"]) and pd.isna(row["count kernel"]):
        return model

    suffix = "SK" if "tanimoto" in str(row["fp kernel"]).lower() else "Bitwise"
    return f"{model} ({suffix})"


def model_config_sort_key(
    row: pd.Series,
    model_order: Optional[Mapping[str, int]] = None,
) -> tuple[int, int, str]:
    model = str(row["model"])
    model_rank = len(model_order) if model_order is not None else 0
    if model_order is not None:
        model_rank = model_order.get(model.lower(), model_rank)

    fp_kernel = row["fp kernel"]
    if is_tree_model(model):
        family_rank = 0
    elif model == "MGK":
        family_rank = 3
    elif pd.notna(fp_kernel) and "tanimoto" not in str(fp_kernel).lower():
        family_rank = 1
    elif pd.notna(fp_kernel) and "tanimoto" in str(fp_kernel).lower():
        family_rank = 2
    else:
        family_rank = 4
    return (family_rank, model_rank, model)


def model_type_color(
    model: Any,
    palette: Optional[Mapping[str, str]] = None,
) -> str:
    colors = MODEL_TYPE_COLORS if palette is None else palette
    model_name = str(model)
    model_upper = model_name.upper()
    if is_tree_model(model_name):
        return colors["tree"]
    if model_name == "MGK":
        return colors["mgk"]
    if model_name in GNN_MODELS or any(token in model_upper for token in GNN_MODELS):
        return colors["gnn"]
    if "GP" in model_upper:
        return colors["gp"]
    return "#808080"


def lighten_color(color: str, amount: float = 0.45) -> str:
    """Blend a color toward white by ``amount`` in the interval [0, 1]."""
    if not 0 <= amount <= 1:
        raise ValueError("amount must be between 0 and 1.")
    rgb = np.asarray(to_rgb(color), dtype=float)
    return to_hex(rgb + (1.0 - rgb) * amount)
