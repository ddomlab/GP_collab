import json
from pathlib import Path
from typing import List, Optional, Any, Dict
import os 
import re
from collections import defaultdict
import matplotlib.colors as mcolors

# visualization imports
# import cmcrameri.cm as cmc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

#statistics imports
# import krippendorff
import pingouin as pg
from scipy.stats import kendalltau
#docs
from docx import Document
from docx.shared import Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH

# internal imports
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
    print('yes')

    score_files = list(Path(target_folder).glob(pattern))
    for score_path in score_files:
        if "heatmap" in score_path.parts:
            continue

        if (
            "generalizability" in score_path.name
            or "test" in score_path.parts
            or "lc_scores" in score_path.name
            or "dropped_nans" in score_path.name
            ):
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
                            figsize: tuple=(6,6)
                            ) -> None:
    scores_to_show:pd.DataFrame = _get_count_fp_comparison(
                                                            target_folder=target_dir,
                                                            scoring_metric=score_metric,
                                                            selected_models=models_to_draw,
                                                            features_to_draw=features_to_draw
                                                            )
                                                            
    print(scores_to_show)
    fname= f"model vs features_{score_metric}"
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
                        figsize=figsize,
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
    if isinstance(scores_data.get("feature_importance_MDI"), list):
        all_model_fi.extend(scores_data["feature_importance_MDI"])
    else:
        for value in scores_data.values():
            if isinstance(value, dict) and isinstance(value.get("feature_importance_MDI"), list):
                all_model_fi.extend(value["feature_importance_MDI"])

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
    # fig, ax = plt.subplots(figsize=figsize)

    # # Build long-format dataframe
    # plot_data = []

    # if not df_model.empty:
    #     for feat in features:
    #         for val in df_model[feat].dropna().values:
    #             plot_data.append([feat, val, "Model Importance"])

    # if not df_shap.empty:
    #     for feat in features:
    #         for val in df_shap[feat].dropna().values:
    #             plot_data.append([feat, val, "SHAP Importance"])

    # plot_df = pd.DataFrame(plot_data, columns=["Feature", "Value", "Type"])

    # Colors: red for model, blue for shap
    # palette = {
    #     "Model Importance": "#BB3A5A",  # red
    #     "SHAP Importance": "#1f77b4",   # blue
    # }

    # # VIOLIN PLOTS (side-by-side)
    # sns.violinplot(
    #     data=plot_df,
    #     x="Feature",
    #     y="Value",
    #     hue="Type",
    #     palette=palette,
    #     inner=None,
    #     linewidth=1.2,
    #     cut=0,
    #     dodge=True,
    #     ax=ax,
    # )

    # # SWARMPLOTS
    # sns.swarmplot(
    #     data=plot_df,
    #     x="Feature",
    #     y="Value",
    #     hue="Type",
    #     dodge=True,
    #     palette=palette, 
    #     size=3,
    #     alpha=0.7,
    #     ax=ax,
    # )

    # # Fix duplicate legends (remove swarm's entries)
    # handles, labels = ax.get_legend_handles_labels()
    # ax.legend(handles[:2], labels[:2], frameon=False)

    # ax.set_xticklabels(features, rotation=75, ha="right")
    # ax.set_ylabel("Feature Importance")
    # ax.set_xlabel("Top 15 Features")
    # ax.grid(axis="y", linestyle="--", alpha=0.4)

    # plt.tight_layout()
    # save_img_path(save_loc / "feature importance", f"feature_importance_top15_{file_extension}.png")
    # # plt.show()
    # plt.close()
    return df_model, df_shap



def get_lengthscale_stat(
        scores: dict,
        expert_rank=None,
        feature_stability=True,
        mean_std=True
    ) -> Dict:

    all_rows = []
    feat_samples = {}   # <-- collect raw samples for mean/std over seeds, folds, draws

    for key, value in scores.items():
        if isinstance(value, dict) and "test_lengthscale" in value:
            fold_list = value["test_lengthscale"]

            for fold_dict in fold_list:
                row_data = {}

                for feat_name, samples in fold_dict.items():
                    if samples is not None:
                        samples = np.asarray(samples)

                        # Row dataframe: mean across draws only for df_ls
                        row_data[feat_name] = np.mean(samples)

                        # For global stats: keep ALL draws across folds + seeds
                        if feat_name not in feat_samples:
                            feat_samples[feat_name] = []
                        feat_samples[feat_name].append(samples.flatten())
                    else:
                        row_data[feat_name] = np.nan

                all_rows.append(row_data)

    # main df like before (mean over draws per fold/seed)
    df_ls = pd.DataFrame(all_rows)

    # ---- NEW: mean + std over ALL seeds × folds × draws ----
    mean_std_stats = {}
    if mean_std or expert_rank is not None:
        for feat, arr_list in feat_samples.items():
            arr = np.concatenate(arr_list)
            mean_std_stats[feat] = {
                "mean": np.mean(arr),
                "std": np.std(arr, ddof=1)
            }
            df_mean_std = pd.DataFrame(mean_std_stats).T

    kendalls_w_result=  kendalls_w(df_ls, tie_corrected=False)["Kendall's W"] if feature_stability else None
        

    kendall_tau_result = None
    if expert_rank is not None:
        expert_series = expert_rank.iloc[0].copy()

        assert set(expert_series.index) == set(df_mean_std.index), \
            "Features in expert_rank and length scale mean must match exactly!"

        ls_rank = df_mean_std["mean"].rank(method="average", ascending=True)
        ls_rank_aligned = ls_rank.loc[expert_series.index]
        print("expert rank\n", expert_series)
        print("model rank (aligned)\n", ls_rank_aligned)
        print("real ls values\n", df_mean_std["mean"].loc[expert_series.index])
        kendall_tau_result = kendalltau(
            expert_series.values,
            ls_rank_aligned.values
        ).statistic

    stat_results = {
        "kendalls_w":kendalls_w_result,
        "kendall_tau": kendall_tau_result,
        "mean_std": mean_std_stats if mean_std else None,
    }
    # kendallTau_results = kendalltau(x=, y=)
    return stat_results


def kendalls_w(df_input, tie_corrected=True):
    """
    Kendall's W for agreement across raters (rows) on items (columns).
    If tie_corrected=True, applies the standard tie correction.
    """
    m = len(df_input)          # raters / folds
    n = len(df_input.columns)  # items / features

    # Rank each row
    ranks = df_input.rank(axis=1, ascending=False, method="average")

    # Sum of ranks for each feature
    R = ranks.sum(axis=0)
    R_bar = m * (n + 1) / 2.0

    # S term
    S = np.sum((R - R_bar) ** 2)

    if not tie_corrected:
        # Classic denominator
        S_max = (m**2 * (n**3 - n)) / 12.0
        W = S / S_max if S_max != 0 else np.nan
        chi_sq = m * (n - 1) * W
        df_chi = n - 1
        return {
            "m (Runs)": m,
            "n (Features)": n,
            "Kendall's W": W,
            "Chi-square": chi_sq,
            "Degrees of Freedom": df_chi,
            "Tie Corrected": False
        }

    # -------- Tie corrected version --------
    T_total = 0
    for _, row in df_input.iterrows():
        counts = row.value_counts()
        ties = counts[counts > 1]
        if len(ties) > 0:
            T_total += np.sum(ties**3 - ties)

    denom = (m**2 * (n**3 - n)) - (m * T_total)
    W = (12 * S) / denom if denom != 0 else np.nan

    chi_sq = m * (n - 1) * W
    df_chi = n - 1

    return {
        "m (Runs)": m,
        "n (Features)": n,
        "Kendall's W": W,
        "Chi-square": chi_sq,
        "Degrees of Freedom": df_chi,
        "Tie Corrected": True,
        "Total Tie Correction Term (Σ(t^3−t))": T_total
    }



def get_scores(data: Dict, metric: List[str]):
    special = {"rmse", "r2", "mae"}

    return {
        score: (
            f"{round(data[f'{score}_avg'], 2)} ± {round(data[f'{score}_stdev'], 2)}"
            if score.lower() in special
            else f"{int(round(data[score]))}"
        )
        for score in metric
    }


def create_word_table_table(rows_data, folder_path, file_name="results_gp_table.docx"):
    document = Document()

    rows = len(rows_data) + 3
    cols = 9   # <-- 9 columns now (Combination + 2 models × 3 metrics)
    table = document.add_table(rows=rows, cols=cols)
    table.style = "Table Grid"

    # ---------------- HEADER LEVEL 1 ----------------
    h = table.rows[0].cells
    h[0].text = "Combination"
    h[3].text = "GPyro HMC"
    h[6].text = "GPyTorch MAP"

    h[0].merge(h[2])   # Combination spans Count/FP/Mixing
    h[3].merge(h[5])   # GPyro HMC spans RMSE / R² / runtime
    h[6].merge(h[8])   # GPyTorch MAP spans RMSE / R² / runtime

    # ---------------- HEADER LEVEL 2 ----------------
    h2 = table.rows[1].cells
    h2[0].text = "Kernel"
    h2[3].text = "RMSE"
    h2[4].text = "R²"
    h2[5].text = "Run time (sec)"
    h2[6].text = "RMSE"
    h2[7].text = "R²"
    h2[8].text = "Run time (sec)"
    h2[0].merge(h2[1])

    # ---------------- HEADER LEVEL 3 ----------------
    h3 = table.rows[2].cells
    h3[0].text = "Count"
    h3[1].text = "FP"
    h3[2].text = "Mixing method"

    # ---- VERTICAL MERGES ----
    table.rows[1].cells[2].merge(table.rows[2].cells[2])
    for i in range(3, 9):
        table.rows[1].cells[i].merge(table.rows[2].cells[i])

    # ---------------- DATA ROWS ----------------
    r = 3
    for row in rows_data:
        cells = table.rows[r].cells
        for c, v in enumerate(row):
            cells[c].text = str(v)
        r += 1

    # ---------------- STYLING ----------------
    for row in table.rows:
        for cell in row.cells:
            for p in cell.paragraphs:
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                for run in p.runs:
                    run.font.size = Pt(10)

    os.makedirs(folder_path, exist_ok=True)
    document.save(folder_path / file_name)
    print(f"Saved Word table → {file_name}")



baseline_kernel = ["RBF", "Matern32", "Matern52"]
tanimoto_kernel = ["TanimotoMatern32", "TanimotoMatern52", "TanimotoRBF","Tanimoto",  ]
count_kernel = ["Matern32", "Matern52", "RBF"]
mixing_methods = ["sum", "product", "averageProduct"]

# PAPER = {
#         "Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation": ["target_log Rg (nm)"],
#         "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices": ["target_calculated PCE (%)"],
#         "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery": ["target_log (Separation factor)","target_log (Total flux)"],
#         "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation": ["target_log (Separation factor)","target_log (Total flux)"],
#         "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning": [
#         "target_flux decline ratio (%)",
#         "target_flux recovery ratio (%)",
#         "target_irreversible fouling ratio(%)",
#         "target_organic compound removal (%)",
#         "target_reversible fouling ratio (%)",
#         r"target_water permeability (LMH\bar)",
#         ],
#         }



PLS_Ranks = pd.DataFrame([{
                "Xn": 1,
                "Mw (g/mol)": 2,
                "Concentration (mg/ml)": 3,
                "Temperature SANS/SLS/DLS/SEC (K)": 4,
                "polymer dD": 5,
                "solvent dD": 6,
                "polymer dP": 7,
                "solvent dP": 8,
                "polymer dH": 9,
                "solvent dH": 10,
                "fp": 11,
                "PDI": 12
            }], 
            index=["expert_rank"]
            )

BMS_Ranks = pd.DataFrame([{
                "HOMO_D (eV)": 1,
                "LUMO_A (eV)": 2,
                "Eg_D (eV)": 3,
                "Eg_A (eV)": 4,
                "Ehl_D (eV)": 5,
                "Ehl_A (eV)": 6,
                "HOMO_A (eV)": 7,
                "LUMO_D (eV)": 8,
                "fp_Donor": 9,
                "fp_Acceptor": 10,
                "D:A ratio (m/m)": 11,
                "temperature of thermal annealing": 12,
                "solvent additive conc. (% v/v)": 13,
                "HTL energy level (eV)": 14,
                "ETL energy level (eV)": 15
            }],
            index=["expert_rank"]
            )

PAPER = {
    "Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation": {
        "target": ["target_log Rg (nm)"],
        "expert_impt": PLS_Ranks
    },
    "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices": {
        "target": ["target_calculated PCE (%)"],
        "expert_impt": BMS_Ranks
    },
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery": {
        "target": ["target_log (Separation factor)", "target_log (Total flux)"],
        "expert_impt": None
    },
    "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation": {
        "target": ["target_log (Separation factor)", "target_log (Total flux)"],
        "expert_impt": None
    },
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning": {
        "target": [
            "target_flux decline ratio (%)",
            "target_flux recovery ratio (%)",
            "target_irreversible fouling ratio(%)",
            "target_organic compound removal (%)",
            "target_reversible fouling ratio (%)",
            r"target_water permeability (LMH\bar)"
        ],
        # "expert_impt": None
    }
}


models = [
            "GpyroMCMC", 
            "GPytorchMAP"
            ]




tanimoto_kernel = ["TanimotoMatern32", "TanimotoMatern52", "TanimotoRBF", "Tanimoto"]
count_kernel    = ["Matern32", "Matern52", "RBF"]
mixing_methods  = ["sum", "product", "averageProduct"]
mixing_labels   = {"sum": "Sum", "product": "Product", "averageProduct": "Av(count)×FP"}


def load_r2(paper_loc, model, fp_k, count_k, mix_method):
    """Return r2_avg or None if the file is missing."""
    for suffix in ["", "_mean"]:
        tmpl = (
            f"(ECFP3_count_512-COUNT)_"
            f"({model}_{fp_k}-{count_k}_{mix_method})"
            f"{suffix}_hypOFF_Standard_Standard_scores"
        )
        p = ensure_long_path(paper_loc / f"{tmpl}.json")
        if p.exists():
            with open(p) as f:
                data = json.load(f)
            return data.get("r2_avg")
    return None

def build_heatmap_data(model):
    scores = defaultdict(list)

    for paper_name, paper_info in PAPER.items():
        for target in paper_info["target"]:
            paper_loc = RESULTS / paper_name / target
            for fp_k in tanimoto_kernel:
                for count_k in count_kernel:
                    for mix in mixing_methods:
                        r2 = load_r2(paper_loc, model, fp_k, count_k, mix)
                        if r2 is not None:
                            scores[(fp_k, count_k, mix)].append(r2)

    combos = [f"{fp}:\n{c}" for fp in tanimoto_kernel for c in count_kernel]
    rows = [mixing_labels[m] for m in mixing_methods]
    pairs = [(fp, c) for fp in tanimoto_kernel for c in count_kernel]

    mean_mat = np.full((len(mixing_methods), len(combos)), np.nan)
    std_mat  = np.full((len(mixing_methods), len(combos)), np.nan)
    annot_mat = np.full((len(mixing_methods), len(combos)), "", dtype=object)

    for j, (fp_k, count_k) in enumerate(pairs):
        for i, mix in enumerate(mixing_methods):
            vals = scores[(fp_k, count_k, mix)]
            if vals:
                avg = np.mean(vals)
                std = np.std(vals)   # or np.std(vals, ddof=1) for sample std if len(vals) > 1

                mean_mat[i, j] = avg
                std_mat[i, j] = std

                avg_txt = f"{avg:.2f}"
                std_txt = f"{std:.2f}"
                annot_mat[i, j] = f"{avg_txt}\n±{std_txt}"

    mean_df = pd.DataFrame(mean_mat, index=rows, columns=combos)
    std_df = pd.DataFrame(std_mat, index=rows, columns=combos)
    annot_df = pd.DataFrame(annot_mat, index=rows, columns=combos)

    return mean_df, std_df, annot_df


def plot_heatmap(model, save_dir: Path):
    df, std_df, annot_df = build_heatmap_data(model)

    df = df.dropna(axis=0, how="all").dropna(axis=1, how="all")
    annot_df = annot_df.loc[df.index, df.columns]

    fig, ax = plt.subplots(figsize=(10, 3.5))

    sns.heatmap(
        df,
        ax=ax,
        cmap="Reds",
        vmin=0.5,
        vmax=0.9,
        mask=~np.isfinite(df),
        annot=annot_df,
        fmt="",
        annot_kws={"fontsize": 8, "fontweight": "bold"},
        linewidths=0.5,
        linecolor="white",
        cbar_kws={
            # "label": "Average R² ± Stdev",
            # "pad": 0.02,
            "ticks": [0.5, 0.6, 0.7, 0.8, 0.9],
        },
    )

    cbar = ax.collections[0].colorbar
    cbar.set_label("Average R² ± Stdev", rotation=270, labelpad=20, fontsize=14)
    cbar.set_ticks([0.5, 0.6, 0.7, 0.8, 0.9])
    cbar.set_ticklabels(["0.5", "0.6", "0.7", "0.8", "0.9"])

    n_count = len(count_kernel)
    for g in range(1, len(tanimoto_kernel)):
        ax.axvline(g * n_count, color="white", linewidth=2)

    ax.set_xticklabels(ax.get_xticklabels(), fontsize=8, ha="right", rotation=30)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=10, rotation=0)
    ax.set_xlabel("FP Kernel : Count Kernel", fontsize=11, labelpad=10)
    ax.set_ylabel("Hybridization Configuration", fontsize=11)

    plt.tight_layout()
    save_img_path(save_dir, f"heatmap_r2_{model}.png")
    plt.close(fig)

if __name__ == "__main__":

    # model_stats = {}
    # for paper_name, paper_info in PAPER.items():
    #     for target in paper_info["target"]:
    #         print(paper_name, target)
            # creat_count_fp_heatmap(
            #                         target_dir=RESULTS/paper_name/target,
            #                         score_metric='rmse',
            #                         figsize=(7,4.5),
            #                         # comparison_value=['scaler', 'Trimer_scaler'],
            #                         )
            # feat_impt_expert = paper_info["expert_impt"]
            # for model in ["RF","XGBR", "NGB"]:    
                # paper_loc: Path = RESULTS / paper_name / target
            #     file_name = f"(ECFP3_count_512-COUNT)_{model}_hypOFF_Standard_Standard_chain1_scores"
            #     score_path = ensure_long_path(paper_loc / f"{file_name}.json")
            #     if not score_path.exists():
            #         file_name = f"(ECFP3_count_512-COUNT)_{model}_hypOFF_Standard_Standard_product_chain1_scores"
            #         score_path = ensure_long_path(paper_loc / f"{file_name}.json")
            #     if not score_path.exists():
            #         file_name = f"(ECFP3_count_512-COUNT)_{model}_mean_hypOFF_Standard_Standard_product_chain1_scores"
            #         score_path = ensure_long_path(paper_loc / f"{file_name}.json")
            #     # else:
            #     #     print("file not found:", score_path)
            #     with open(score_path, "r") as f:
            #         scores = json.load(f)
            #     MDI_imp, shap_imp = plot_average_feature_importances(scores_data=scores,
            #                                     save_loc=paper_loc,
            #                                     file_extension=file_name,
            #                                     figsize=(8,7.5)
            #                                     )
                # shap_feature_means = shap_imp.abs().mean()
                # df_top15_shap_features = shap_imp[shap_feature_means.sort_values(ascending=False).head(15).index]

                # mdi_feature_means = MDI_imp.mean()
                # df_top15_mdi_features = MDI_imp[mdi_feature_means.sort_values(ascending=False).head(15).index]

                # 3. Filter the DataFrame to these 15 features
                # plot_top15_feature_stability(
                #                     scores_data=scores,
                #                     # save_loc=paper_loc,
                #                     # file_extension=file_name,
                #                     # top_n=15,
                #                     # figsize=(8,6)
                #                     )
                # krippendorff_alpha_by_feature(
                #                             df=df_top15,             
                #                             save_loc=paper_loc,
                #                             file_extension=file_name,
                #                             figsize=(9,6)
                #                             )
                # ls_stats = get_lengthscale_stat(scores, expert_rank=feat_impt_expert, feature_stability=True, mean_std=False)
                # print(df_ls)
                # model_stats.setdefault(paper_name, {}).setdefault(target, {}).setdefault(model, {})["length scale"] = ls_stats
                # model_stats.setdefault(paper_name, {}).setdefault(target, {}).setdefault(model, {})["MDI"] = kendalls_w(MDI_imp)["Kendall's W"]
                # model_stats.setdefault(paper_name, {}).setdefault(target, {}).setdefault(model, {})["SHAP"] = kendalls_w(shap_imp)["Kendall's W"]

                # print(pg.friedman(df_top15))
                # print(MDI_imp)
    # with open(RESULTS / "model_stats" / "tree_model_stability.json", "w") as f:
    #     json.dump(model_stats, f, indent=2)

    # with open(RESULTS / "model_stats" / "model_stability_ls.json", "w") as f:
    #     json.dump(model_stats, f, indent=2)

            ### SCORE TABLE CREATION ###
                # rows_data = []
                # for count_k in count_kernel:
                #     for fp_k in baseline_kernel:
                #         for mix_method in mixing_methods:
                #             row = [count_k, fp_k, "Av(co)*FP" if mix_method=="averageProduct" else mix_method]
                #             for model in models:
                #                 score_file = None
                #                 file_template = f"(ECFP3_count_512-COUNT)_({model}_{fp_k}-{count_k}_{mix_method})_hypOFF_Standard_Standard_scores"
                #                 score_path = ensure_long_path(paper_loc / f"{file_template}.json")
                #                 if not score_path.exists():
                #                     file_template = f"(ECFP3_count_512-COUNT)_({model}_{fp_k}-{count_k}_{mix_method})_mean_hypOFF_Standard_Standard_scores"
                #                     score_path = ensure_long_path(paper_loc / f"{file_template}.json")
                #                     if not score_path.exists():
                #                         file_template = f"(ECFP3_count_512-COUNT)_({model}_{fp_k}-{count_k}_{mix_method})_hypOFF_Standard_Standard_scores"
                #                         score_path = ensure_long_path(paper_loc / f"{file_template}.json")
                #                         if not score_path.exists():
                #                             file_template = f"(ECFP3_count_512-COUNT)_({model}_{fp_k}-{count_k}_{mix_method})_mean_hypOFF_Standard_Standard_scores"
                #                             score_path = ensure_long_path(paper_loc / f"{file_template}.json")
                #                 if not score_path.exists():
                #                     print(f"❌ Missing score: {paper_name}\n{target}\nfile name: {file_template}")

                #                 else:
                #                     with open(score_path, "r") as f:
                #                         score_file = json.load(f)

                #                 if score_file is None:
                #                     rmse_value = ""
                #                     r2_value = ""
                #                     run_time_value = ""
                #                 else:
                #                     score_annot = get_scores(score_file, metric=['rmse','r2', "run_time_sec"])
                #                     rmse_value = score_annot["rmse"]
                #                     r2_value = score_annot["r2"]
                #                     run_time_value = score_annot["run_time_sec"]

                #                 row.extend([rmse_value, r2_value, run_time_value])
                #             rows_data.append(row)
                # create_word_table_table(rows_data, folder_path=paper_loc/"tabular results", file_name=f"kernel_combination_scores.docx")
                for model in models:
                    plot_heatmap(model, save_dir=HERE/"result_analysis")
                                