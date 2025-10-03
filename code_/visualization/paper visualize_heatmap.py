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
import re
# from matplotlib import rc
from visualization_setting import set_plot_style, save_img_path
from visualize_ood_learning_curve import ensure_long_path

HERE: Path = Path(__file__).resolve().parent
RESULTS: Path = HERE.parent.parent/ 'results'

target_list = [
    # 'target_Rg',
    # 'target_Rh' 
    # "target_multimodal Rh",
    # "target_multimodal Rh_without_log",
    # "target_multimodal Rh with padding",
    # "target_multimodal Rh (e-5 place holder)_with_Log",
    # "target_multimodal Rh (e-5 place holder)",
    # "target_log10 multimodal Rh (e-5 place holder)"
    # "target_Rh First Peak",
    # "target_Rh Second Peak",
    # "target_Rh Third Peak",
    # "target_Rh First Peak_with_Log",
    # "target_Rh Second Peak_with_Log",
    # "target_Rh Third Peak_with_Log",
    "target_log Rg (nm)"
    # 'target_Rh (1_1000 nm) (highest intensity)_LogFT'
    ]

transformer_list = [
    "Standard",
    # "Robust Scaler",
    # "transformerOFF"
                    ]

scores_list: list = [
                    "r2", 
                    "mae", 
                    "rmse"
                     ] 
var_titles: dict[str, str] = {"stdev": "Standard Deviation", "stderr": "Standard Error"}




def parse_property_string(prop_string):
    import re

    # Define the feature categories
    categories = {
        'Xn': ['Xn'],
        'polysize': ['Mw', 'PDI'],
        'solvent_properties': ['concentration', 'temperature'],
        'polymer_HSPs': ['polymer dP', 'polymer dD', 'polymer dH'],
        'solvent_HSPs': ['solvent dP', 'solvent dD', 'solvent dH'],
        'Ra': ['Ra'],
        'HSPs differences': [
            'abs(solvent dD - polymer dD)',
            'abs(solvent dP - polymer dP)',
            'abs(solvent dH - polymer dH)'
        ],
        'environmental.thermal history': [
            'light exposure', 'aging time', 'aging temperature',
            'prep temperature', 'prep time'
        ],
        'ECFP6.count.512': ['ECFP3.count.512'],
        'MACCS': ['MACCS'],
        'Mordred': ['Mordred'],
    }

    # ✅ Step 1: Split by `)_` and remove leading `(` to isolate feature block
    try:
        feature_block = prop_string.split(')_')[0][1:]
    except IndexError:
        return 'unknown format'

    # ✅ Step 2: Tokenize into features while preserving abs(...) blocks
    components = []
    i = 0
    while i < len(feature_block):
        if feature_block[i:i+4] == 'abs(':
            j = i + 4
            bracket_level = 1
            while j < len(feature_block) and bracket_level > 0:
                if feature_block[j] == '(':
                    bracket_level += 1
                elif feature_block[j] == ')':
                    bracket_level -= 1
                j += 1
            components.append('abs(' + feature_block[i+4:j-1].strip() + ')')
            i = j + 1 if j < len(feature_block) and feature_block[j] == '-' else j
        else:
            j = i
            while j < len(feature_block) and feature_block[j] != '-':
                j += 1
            token = feature_block[i:j].strip()
            if token:
                components.append(token)
            i = j + 1

    # ✅ Step 3: Match categories exactly
    present_categories = []
    for cat_name, feat_list in categories.items():
        if all(f in components for f in feat_list):
            present_categories.append(cat_name)

    return ' + '.join(present_categories) if present_categories else 'unknown format'

    
    
    

    


def get_total(prop_string: str) -> Optional[str]:
    """
    Maps full property string to a simplified label for structural_scaler features.
    """
    mapping = {
        'polysize + solvent_properties + polymer_HSPs + solvent_HSPs + Mordred':
            'Mordred + continuous',
        'polysize + solvent_properties + polymer_HSPs + solvent_HSPs + MACCS':
            'MACCS + continuous',
        'polysize + solvent_properties + polymer_HSPs + solvent_HSPs + ECFP6.count.512':
            'ECFP6.count.512 + continuous',

    }
    return mapping.get(prop_string)


def correct_ecfp_name(prop_string):
    if 'ECFP' in prop_string:
        parts = prop_string.split('.')
        ecfp_part = parts[0]  # e.g., 'ECFP3'
        ecfp_text = ecfp_part[:4]  # 'ECFP'
        ecfp_num = int(ecfp_part[4:])  # extract number
        corrected_ecfp = f"{ecfp_text}{ecfp_num * 2}"  # multiply by 2
        return '.'.join([corrected_ecfp] + parts[1:])
    return prop_string
       
    
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
        f"Average {score_txt.upper()} ± {var_titles.get('stdev')}",
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



def get_polymer_propeties_comparison(target_folder: Path,
                                     score: str,
                                     comparison_value: List[str],
                                     features_to_draw: set[str] = None,
                                    ) -> pd.DataFrame:
    scores_to_report: List = []
    pattern: str = "*_scores.json"

    selected_models: set = {
        "RF",
        "NGB",}

    for value in comparison_value:
        value_folder = os.path.join(target_folder, value)
        score_files = list(Path(value_folder).rglob(pattern))

        for score_path in score_files:
            if "generalizability" in score_path.name or "test" in score_path.name or 'lc_scores' in score_path.name:
                continue

            feats, model, av, std = get_results_from_file(file_path=score_path, score=score)
            # Only keep selected features
            if feats not in features_to_draw:
                continue

            if model not in selected_models:
                continue

            anot = f"{np.round(av, 2)}\n±{np.round(std, 2)}"
            scores_to_report.append({
                "features": feats,
                "model": model,
                "score": np.round(av, 2),
                "annotations": anot
            })

    return pd.DataFrame(scores_to_report)


def plot_manual_heatmap(
    root_dir: Path,
    score: str,
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
    
    def wrap_label(label: str, max_words_per_line: int = 2) -> str:
        words = label.split(" + ")
        return "\n".join(
            [" + ".join(words[i:i + max_words_per_line]) for i in range(0, len(words), max_words_per_line)]
        )

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

    palette = "viridis" if score in ["r", "r2"] else "viridis_r"
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
    y_tick_labels = [wrap_label(label) for label in avg_scores.index]

    ax.set_xticklabels(x_tick_labels, rotation=45, ha="right", fontsize=kwargs['fontsize'])
    ax.set_yticklabels(y_tick_labels, rotation=0, ha="left", fontsize=kwargs['fontsize'])
    ax.tick_params(axis='y', pad=5)

    # Titles
    plt.title(fig_title, fontsize=kwargs['fontsize'] + 2, fontweight='bold')
    ax.set_xlabel(x_title, fontsize=kwargs['fontsize'] + 2, fontweight='bold')
    ax.set_ylabel(y_title, fontsize=kwargs['fontsize'] + 2, fontweight='bold')

    # Colorbar manually on the left
    var_titles = {"stdev": "Stdev"}
    score_txt = "$R^2$" if score == "r2" else score

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

    plt.tight_layout()
    save_img_path(root_dir, f"{fname}.png")
    plt.show()
    plt.close()



def creat_result_df(target_dir: Path,
                    score: str,
                    data_type: str,
                    transformer_type: str,
                    regressor_model: Optional[str] = None,
                    peak_number: Optional[int] = None
) -> tuple[pd.DataFrame, pd.DataFrame, list]:
    
    avg_scores: pd.DataFrame = pd.DataFrame()
    std_scores: pd.DataFrame = pd.DataFrame()
    annotations: pd.DataFrame = pd.DataFrame()
    models = set()
    pattern: str = "*_scores.json"
    
    for representation in os.listdir(target_dir):
        score_files = []
        if data_type == 'structural':
            if 'test' not in representation and 'scaler' not in representation:
                filterd_rep = representation
                representation_dir = os.path.join(target_dir, representation)
                score_files = list(Path(representation_dir).rglob(pattern))

        elif data_type == 'scaler':
            if 'scaler' == representation:
                filterd_rep = representation
                representation_dir = os.path.join(target_dir, representation)
                score_files = list(Path(representation_dir).rglob(pattern))

        elif data_type == 'structural_scaler':
            if 'scaler' in representation and 'scaler' != representation:
                filterd_rep = representation
                representation_dir = os.path.join(target_dir, representation)
                score_files = list(Path(representation_dir).rglob(pattern))

        for file_path in score_files:
            if "generalizability" in file_path.name or "test" in file_path.name or 'lc_scores' in file_path.name:
                continue

            if data_type == 'structural':
                if regressor_model != file_path.name.split("_")[1]:
                    continue
                if transformer_type not in file_path.name:
                    continue

                feats, model, av, std = get_results_from_file(file_path=file_path, score=score, peak_number=peak_number)
                models.add(model)

                if data_type == 'structural':
                    if feats not in ['MACCS', 'Mordred', 'ECFP6.count.512']:
                        continue

                  

            elif data_type in ['scaler', 'structural_scaler']:  
                if transformer_type not in file_path.name:
                    continue

                feats, model, av, std = get_results_from_file(file_path=file_path, score=score, peak_number=peak_number)
                models.add(model)
                
                if data_type == 'scaler':
                    exclude_feats = ["unknown format", "polysize + solvent_properties + polymer_HSPs + solvent_HSPs + Ra"]
                    if feats in exclude_feats:
                        continue
                else:
                    simplified_feats = get_total(feats)
                    if simplified_feats is None:
                        continue
                    feats = simplified_feats

            if data_type in ['scaler', 'structural_scaler']:
                if data_type== 'structural_scaler':
                    x= model
                    y= feats
                else:
                    x = feats
                    y = model
                if feats not in avg_scores.columns:
                    avg_scores.loc[x, y] = av
                    std_scores.loc[x, y] = std
                else:
                    avg_scores.at[x, y] = av
                    std_scores.at[x, y] = std
            else:  # structural
                if feats not in avg_scores.columns:
                    avg_scores.loc[filterd_rep, feats] = av
                    std_scores.loc[filterd_rep, feats] = std
                else:
                    avg_scores.at[filterd_rep, feats] = av
                    std_scores.at[filterd_rep, feats] = std

    for x, y in product(avg_scores.columns.to_list(), avg_scores.index.to_list()):
        avg: float = avg_scores.loc[y, x]
        std: float = std_scores.loc[y, x]
        avg_txt: str = generate_annotations(avg)
        std_txt: str = generate_annotations(std)
        annotations.loc[y, x] = f"{avg_txt}\n±{std_txt}"

    avg_scores = avg_scores.astype(float)
    annotations = annotations.astype(str)
    return avg_scores, annotations, list(models)


def create_structural_result(target_dir:Path,
                             regressor_model:str,
                             target:str,
                             score:str,
                            #  var:str,
                             data_type:str,
                             transformer_type:str,
                             peak_num:int=None
                             ) -> None:
    ave, anot, model = creat_result_df(target_dir=target_dir,score=score,data_type=data_type
                                       ,regressor_model=regressor_model, transformer_type=transformer_type,
                                       peak_number=peak_num)
    
    model_in_title:str =  ",".join(model)
    score_txt: str = "$R^2$" if score == "r2" else score.upper()
    reg_name = f'{regressor_model} on peak {peak_num+1}' if peak_num else regressor_model
    fname= f"selected PolymerRepresentation vs Fingerprint trained by {reg_name} search heatmap_{score} score"
    if score == "r2":
        vmax= .9
        vmin= .1
        n_cbar_tick = 9
    elif score == "mae":
        vmax= .5
        vmin= 0.1
        n_cbar_tick = 5  
    elif score == "rmse":
        vmax= .6
        vmin= 0.2
        n_cbar_tick = 5
    _create_heatmap(root_dir=HERE,
                    score=score,
                    # var=var,
                    avg_scores=ave,
                    annotations=anot,
                    figsize=(9, 8),
                    fig_title=f"",
                    x_title="Molecular Representations",
                    y_title="Polymer Unit Representation",
                    fname=fname,
                    vmin=vmin,
                    vmax=vmax,
                    num_ticks=n_cbar_tick,
                    )


def create_structural_scaler_result(target_dir:Path,
                                    # regressor_model:str,
                                    target:str,
                                    score:str,
                                    # var:str,
                                    data_type:str,
                                    transformer_type:str,
                                    peak_num:int=None
                                    ) -> None:

    ave, anot, model = creat_result_df(target_dir=target_dir,score=score,data_type=data_type
                                       , transformer_type=transformer_type,
                                       peak_number=peak_num)
    model_in_title:str =  ",".join(model)
    score_txt: str = "$R^2$" if score == "r2" else score.upper()
    # fname = f'{fname} on peak {peak_num+1}' if peak_num else fname
    fname= f"selected Model vs all features search heatmap (structural+continuous)_{score} score"
    if score == "r2":
        vmax= .9
        vmin= .3
        n_cbar_tick = 4  
    elif score == "mae":
        vmax= .5
        vmin= 0.1
        n_cbar_tick = 5  
    elif score == "rmse":
        vmax= .6
        vmin= 0.2 
        n_cbar_tick = 5 
    _create_heatmap(root_dir=HERE,
                    score=score,
                    # var=var,
                    avg_scores=ave,
                    annotations=anot,
                    figsize=(10,7),
                    fig_title='',
                    x_title="Feature Space",
                    y_title="Regression Models",
                    fname=fname,
                    num_ticks=n_cbar_tick,
                    vmin=vmin,
                    vmax=vmax,
                    )

#    'XGBR','RF','NGB'"GPR.matern", "GPR.rbf" "GPR"
complex_models = ['NGB']


# for transformer in transformer_list:
#     for model in complex_models: 
#         for target_folder in target_list:
#             for i in scores_list:
#                 create_structural_scaler_result(target_dir=RESULTS/target_folder,regressor_model= model,target=f'{target_folder} with',
#                                                 score=i,var='stdev',data_type='structural_scaler', transformer_type=transformer)
                # create_structural_result(target_dir=RESULTS/target_folder,regressor_model= model,target=f'{target_folder} with',
                #                             score=i,data_type='structural', transformer_type=transformer)



def create_scaler_result(target_dir:Path,
                        score:str,
                        data_type:str,
                        transformer_type:str,
                        peak_num:int=None
                        )->None:

    ave, anot, model = creat_result_df(target_dir=target_dir,score=score,data_type=data_type,
                                       regressor_model=None,transformer_type=transformer_type,
                                       peak_number=peak_num)
    model_in_title:str =  ",".join(model)
    score_txt: str = "$R^2$" if score == "r2" else score.upper()
    # reg_name = f'{regressor_model} on {peak_num}' if peak_num else regressor_model
    fname= f"selected Regression Models vs numerical features search heatmap_{score}"
    fname = f'{fname} on peak {peak_num+1}' if peak_num else fname


    if score == "r2":
        vmax= .9
        vmin= .1
        n_cbar_tick = 9  
    elif score == "mae":
        vmax= .5
        vmin= 0.1
        n_cbar_tick = 5  
    elif score == "rmse":
        vmax= .6
        vmin= 0.2
        n_cbar_tick = 5

    _create_heatmap(root_dir=HERE,
                    score=score,
                    # var=var,
                    avg_scores=ave,
                    annotations=anot,
                    figsize=(10,9),
                    fig_title=f" \n ",
                    x_title="Regression Models",
                    y_title="Feature Space",
                    fname=fname,
                    vmin=vmin,
                    vmax=vmax,
                    feature_order=['Xn', 'polysize', 'Xn + polysize', 'solvent_properties','polymer_HSPs','solvent_HSPs', 'polymer_HSPs + solvent_HSPs', 'HSPs differences','Ra','environmental.thermal history', ],
                    model_order=['RF', 'ElasticNet','MLR'],
                    num_ticks=n_cbar_tick,
                    # vmin=0.4,
                    # vmax=0.6,
                    )



# for transformer in transformer_list:
#     for target_folder in target_list:
#         for i in scores_list:
            # create_structural_scaler_result(target_dir=RESULTS/target_folder,target=f'{target_folder} with',
            #                                     score=i,data_type='structural_scaler', transformer_type=transformer)
            # create_scaler_result(target_dir=RESULTS/target_folder,
            #                     score=i,data_type='scaler',transformer_type=transformer)
            



            
# selected_features: set = [
#     # "solvent_properties + solvent_HSPs",
#     "polysize + solvent_properties + solvent_HSPs",
#     "polysize + solvent_properties + polymer_HSPs + solvent_HSPs",
#     # "solvent_properties + polymer_HSPs + solvent_HSPs",
#     "polysize + solvent_properties + solvent_HSPs + MACCS",
#     "polysize + solvent_properties + solvent_HSPs + Mordred",
#     "polysize + solvent_properties + solvent_HSPs + ECFP6.count.512",
# ]

# def creat_polymer_properties_comparison(target_dir:Path,
#                                         score:str,
#                                         comparison_value:List[str],
#                                         ) -> None:
    # scores_to_show:pd.DataFrame = get_polymer_propeties_comparison(target_folder=target_dir,
    #                                                                score=score,
    #                                                                comparison_value=comparison_value,
    #                                                                features_to_draw=selected_features)
#     score_txt: str = "$R^2$" if score == "r2" else score.upper()
#     fname= f"model vs features in {score}"
#     plot_manual_heatmap(root_dir=target_dir/"comparison heatmap for polymer properties",
#                         score=score,
#                         score_to_show=scores_to_show,
#                         figsize=(9, 8),
#                         fig_title=f" \n ",
#                         x_title="Feature Space",
#                         y_title="Models",
#                         fname=fname,
#                         vmin=.2,
#                         vmax=.6,
#                         feature_order=selected_features,
#                         # model_order=['RF','DT','MLR'],
#                         num_ticks=5,
#                         )


# creat_polymer_properties_comparison(target_dir=RESULTS/'target_log Rg (nm)',
#                                     score='r2',
#                                     comparison_value=['scaler', 'Trimer_scaler'],
#                                     )



def get_aging_comparison(target_folder: Path,
                            score: str,
                            comparison_value: List[str],
                            features_to_draw: set[str] = None,
                            models_to_draw: set[str] = None,
                            special_namings: List[str] = None,
                        ) -> pd.DataFrame:
    scores_to_report: List = []
    pattern: str = "*_scores.json"

    for value in comparison_value:
        value_folder = os.path.join(target_folder, value)
        score_files = list(Path(value_folder).rglob(pattern))

        for score_path in score_files:
            if "generalizability" in score_path.name or "test" in score_path.name or 'lc_scores' in score_path.name:
                continue
            feats, model, av, std = get_results_from_file(file_path=score_path, score=score)

            # Only keep selected features
            # print(feats)
            if feats not in features_to_draw:
                continue
            if model not in models_to_draw:
                continue
            if special_namings:
                for special_name in special_namings:
                    if special_name in score_path.name:
                        feats = f'{feats} ({special_name})'
                    else:
                        feats = feats


            anot = f"{np.round(av, 2)}\n±{np.round(std, 2)}"
            scores_to_report.append({
                "features": feats,
                "model": model,
                "score": np.round(av, 2),
                "annotations": anot
            })
    # print(scores_to_report)
    return pd.DataFrame(scores_to_report) 


def creat_aging_comparison_heatmap(target_dir:Path,
                                        score_metrics:str,
                                        comparison_value:List[str],
                                        features_to_draw: List[str] = None,
                                        models_to_draw: set[str] = None,
                                        special_namings: List[str] = None,
                                        ) -> None:
    
    scores_to_show:pd.DataFrame = get_aging_comparison(target_folder=target_dir,
                                                        score=score_metrics,
                                                        comparison_value=comparison_value,
                                                        features_to_draw=features_to_draw,
                                                        models_to_draw=models_to_draw,
                                                        special_namings=special_namings,
                                                        )
    # score_txt: str = "$R^2$" if score_metrics == "r2" else score_metrics.upper()
    fname= f"model vs features in {score_metrics} for separate comparison (three criteria)"
    if score_metrics == "r2":
        vmax= .9
        vmin= .1
        n_cbar_tick = 9  
    elif score_metrics == "mae":
        vmax= .5
        vmin= 0.1
        n_cbar_tick = 5  
    elif score_metrics == "rmse":
        vmax= .6
        vmin= 0.2 
        n_cbar_tick = 4 
    plot_manual_heatmap(root_dir=target_dir/'comparison heatmap for polymer properties',
                        score=score_metrics,
                        score_to_show=scores_to_show,
                        figsize=(10,12),
                        fig_title=f" ",
                        x_title="Feature Space",
                        y_title="Models",
                        fname=fname,
                        vmin=vmin,
                        vmax=vmax,
                        feature_order=features_to_draw,
                        model_order=['RF','NGB','XGBR'],
                        num_ticks=5,
                        fontsize=16,
                        )


aging_features: List = [
    'Xn + polysize',
    'Xn + polysize + Mordred',
    'Xn + polysize + MACCS',
    'Xn + polysize + ECFP6.count.512',
    'Xn + polysize + polymer_HSPs',
    # 'environmental.thermal history',
    'solvent_properties + solvent_HSPs',
    'solvent_properties + solvent_HSPs + environmental.thermal history',
    'Xn + polysize + solvent_properties + polymer_HSPs + solvent_HSPs',
    'Xn + polysize + solvent_properties + polymer_HSPs + solvent_HSPs + environmental.thermal history',
    # 'Xn + polysize + solvent_properties + environmental.thermal history',

    # 'Xn + polysize + solvent_properties + polymer_HSPs + solvent_HSPs',
    # 'Xn + polysize + solvent_properties + polymer_HSPs + solvent_HSPs + Mordred',
]

creat_aging_comparison_heatmap(target_dir=RESULTS/'target_log Rg (nm)',
                                    score_metrics='r2',
                                    comparison_value=['scaler', 'Trimer_scaler'],
                                    features_to_draw=aging_features,
                                    models_to_draw={'RF','NGB','XGBR'},
                                    # special_namings=['aging_imputed']
                                    )