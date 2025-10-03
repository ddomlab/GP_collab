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
    var: str,
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
    if not file_path.exists():
        features, model = None, None
        avg, std = np.nan, np.nan
    else:
        # for just scaler features
        if "scaler" == file_path.parent.name:
            model:str = file_path.name.split("_")[1] 
            features:str = file_path.name.split("_")[0].replace("(", "", 1)[::-1].replace(")", "", 1)[::-1]
        # for mixture of scaler and structural 
        elif "scaler" in file_path.parent.name and file_path.parent.name != "scaler":
            features:str = file_path.name.split("_")[0].replace("(", "").replace(")", "")
            model:str = file_path.name.split("_")[1]
        # for structural only
        else:
            features:str = file_path.name.split("_")[0].replace("(", "").replace(")", "")
            model:str = file_path.name.split("_")[1]
            features=correct_ecfp_name(features)

       
        with open(file_path, "r") as f:
            data = json.load(f)


        avg = data[f"{score}_avg"]
        # print(avg)
        if var == "stdev":
            std = data[f"{score}_stdev"]
        elif var == "stderr":
            std = data[f"{score}_stderr"]
        else:
            raise ValueError(f"Unknown variance type: {var}")
        

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
    var: str,
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
    **kwargs,
) -> None:
    """
    Args:
        root_dir: Root directory containing all results
        score: Score to plot
        var: Variance to plot
        x_labels: Labels for x-axis
        y_labels: Labels for y-axis
        figsize: Figure size
        fig_title: Figure title
        x_title: X-axis title
        y_title: Y-axis title
        fname: Filename to save figure
    """

    if feature_order is not None:
        avg_scores = avg_scores[feature_order]
        annotations = annotations[feature_order]

    if model_order is not None:
        avg_scores = avg_scores.loc[model_order]
        annotations = annotations.loc[model_order]

    fig, ax = plt.subplots(figsize=figsize)
    palette: str = "viridis" if score in ["r", "r2"] else "viridis_r"
    custom_cmap = sns.color_palette(palette, as_cmap=True)

    # If r2: force vmin and prepare masked grey values
    if score == "r2":
        
        grey_mask = avg_scores < 0.0

    elif score == "rmse":
        grey_mask = avg_scores > 0.6
    elif score == "mae":
        grey_mask = avg_scores > .5
    else:
        grey_mask = pd.DataFrame(False, index=avg_scores.index, columns=avg_scores.columns)

    # Draw full heatmap first
    # custom_cmap.set_bad(color="lightgray")
    hmap = sns.heatmap(
        avg_scores.mask(grey_mask),  # Mask < 0 temporarily to color them grey
        annot=annotations,
        fmt="",
        cmap=custom_cmap,
        cbar=True,
        vmin=vmin,
        vmax=vmax,
        ax=ax,
        annot_kws={"fontsize": 16},
    )

    # Overlay grey boxes for < 0 values in r2
    # if grey_mask.any().any():
    #     for i in range(avg_scores.shape[0]):
    #         for j in range(avg_scores.shape[1]):
    #             if grey_mask.iloc[i, j]:
    #                 rect = plt.Rectangle((j, i), 1, 1, facecolor="lightgray", edgecolor="white", lw=0.5)
    #                 ax.add_patch(rect)
    #                 ax.text(
    #                     j + 0.5,
    #                     i + 0.5,
    #                     f"{annotations.iloc[i, j]}",
    #                     ha="center",
    #                     va="center",
    #                     fontsize=16,
    #                     color="black",
    #                     # fontweight='bold',
    #                 )

    # Set tick positions and labels
    # ax.set_xticks(np.arange(len(avg_scores.columns)) + 0.5)
    # ax.set_yticks(np.arange(len(avg_scores.index)) + 0.5)

    x_tick_labels: list[str] = [col for col in avg_scores.columns]
    y_tick_labels: list[str] = avg_scores.index.to_list()

    ax.set_xticklabels(x_tick_labels, rotation=45, ha="right", fontsize=16)
    ax.set_yticklabels(y_tick_labels, rotation=0, ha="right", fontsize=16)

    # Titles
    plt.title(fig_title, fontsize=22, fontweight='bold',pad=10)
    ax.set_xlabel(x_title, fontsize=20, fontweight='bold')
    ax.set_ylabel(y_title, fontsize=20, fontweight='bold')

    # Colorbar
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
        f"Average {score_txt.upper()} ± {var_titles[var]}",
        rotation=270,
        labelpad=20,
        fontsize=18,
        fontweight='bold',
    )

    cbar.ax.tick_params(labelsize=18)

    # Save
    visualization_folder_path = root_dir / "heatmap"
    os.makedirs(visualization_folder_path, exist_ok=True)
    plt.tight_layout()
    plt.savefig(visualization_folder_path / f"{fname}.png", dpi=600)
    plt.show()
    plt.close()




def creat_result_df(target_dir: Path,
                    score: str,
                    var: str,
                    data_type:str,
                    transformer_type:str,
                    regressor_model:Optional[str]=None,
                    peak_number:Optional[int]=None
) -> tuple[pd.DataFrame,pd.DataFrame]:
    
    avg_scores: pd.DataFrame = pd.DataFrame()
    std_scores: pd.DataFrame = pd.DataFrame()
    annotations: pd.DataFrame = pd.DataFrame()
    models = set()
    pattern: str = "*_scores.json"
    for representation in os.listdir(target_dir):
        score_files = []
        if data_type == 'structural':
            
            if 'test' not in representation  and 'scaler' not in representation:
                filterd_rep = representation
                representation_dir = os.path.join(target_dir, representation)
                score_files: list[Path] = list(Path(representation_dir).rglob(pattern))
        elif data_type == 'scaler':
        
            if  'scaler' == representation:
                filterd_rep = representation
                representation_dir = os.path.join(target_dir, representation)
                score_files: list[Path] = list(Path(representation_dir).rglob(pattern))
        
        elif data_type=='structural_scaler':
            
            if  'scaler' in representation and 'scaler' != representation:
                filterd_rep = representation
                representation_dir = os.path.join(target_dir, representation)
                score_files: list[Path] = list(Path(representation_dir).rglob(pattern))

        for file_path in score_files:
            # for structural and mix of structural-scaler
            if "generalizability" not in file_path.name and "test" not in file_path.name:

                if data_type=='structural_scaler' or data_type=='structural':
                    if regressor_model == file_path.name.split("_")[1]:
                        if transformer_type in file_path.name:
                            # print(transformer_type)
                            feats, model, av , std = get_results_from_file(file_path=file_path, score=score, var=var,peak_number=peak_number)                
                            models.add(model)
                        else:
                            continue
                    
                    else:
                        continue
                    # for just scaler 
                
                else:
                    if transformer_type in file_path.name:
                        print(transformer_type)

                        feats, model, av , std = get_results_from_file(file_path=file_path, score=score, var=var, peak_number=peak_number)  
                        models.add(model)
                    else:
                            continue
                # for scaler features only
                # print(type(av))

                if data_type=='scaler':
                    if feats not in avg_scores.columns:
                        
                        avg_scores.loc[model,feats] = av
                        std_scores.loc[model,feats] = std
                    else:
                        avg_scores.at[model,feats] = av
                        std_scores.at[model,feats] = std
                # for mixtures and fingerprints 
                else:
                    if feats not in avg_scores.columns:

                        avg_scores.loc[filterd_rep,feats] = av
                        std_scores.loc[filterd_rep,feats] = std
                    else:
                        avg_scores.at[filterd_rep,feats] = av
                        std_scores.at[filterd_rep,feats] = std
        


    for x, y in product(avg_scores.columns.to_list(), avg_scores.index.to_list()):
        avg: float = avg_scores.loc[y, x]
        std: float = std_scores.loc[y, x]
        avg_txt: str = generate_annotations(avg)
        std_txt: str = generate_annotations(std)
        annotations.loc[y, x] = f"{avg_txt}\n±{std_txt}"

    avg_scores = avg_scores.astype(float)
    annotations = annotations.astype(str)
    return avg_scores, annotations, list(models)


# def set_minmax_cbar(score_type):


#     return vmin, vmax

def create_structural_result(target_dir:Path,
                             regressor_model:str,
                             target:str,
                             score:str,
                             var:str,
                             data_type:str,
                             transformer_type:str,
                             peak_num:int=None
                             ) -> None:
    ave, anot, model = creat_result_df(target_dir=target_dir,score=score, var=var,data_type=data_type
                                       ,regressor_model=regressor_model, transformer_type=transformer_type,
                                       peak_number=peak_num)
    
    model_in_title:str =  ",".join(model)
    score_txt: str = "$R^2$" if score == "r2" else score.upper()
    reg_name = f'{regressor_model} on peak {peak_num+1}' if peak_num else regressor_model
    fname= f"PolymerRepresentation vs Fingerprint trained by {reg_name}  with {transformer_type} search heatmap_{score} score"
    # vmin, vmax = set_minmax_cbar(score)
   
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


    _create_heatmap(root_dir=target_dir,
                    score=score,
                    var=var,
                    avg_scores=ave,
                    annotations=anot,
                    figsize=(14, 8),
                    fig_title=f"Algorithm: {model_in_title}",
                    x_title="Fingerprint Representations",
                    y_title="Polymer Unit Representation",
                    vmax=vmax,
                    vmin=vmin,
                    num_ticks=n_cbar_tick,
                    fname=fname)



# def create_structural_scaler_result(target_dir:Path,
#                                     regressor_model:str,
#                                     target:str,
#                                     score:str,
#                                     var:str,
#                                     data_type:str,
#                                     transformer_type:str,
#                                     peak_num:int=None
#                                     ) -> None:

#     ave, anot, model = creat_result_df(target_dir=target_dir,score=score, var=var,data_type=data_type,
#                                        regressor_model=regressor_model, transformer_type=transformer_type,
#                                        peak_number=peak_num)
#     model_in_title:str =  ",".join(model)
#     score_txt: str = "$R^2$" if score == "r2" else score.upper()
#     reg_name = f'{regressor_model} on peak {peak_num+1}' if peak_num else regressor_model
#     fname= f"PolymerRepresentation vs (Fingerprint-numerical) trained by {reg_name}  with {transformer_type} search heatmap_{score} score"
#     _create_heatmap(root_dir=target_dir,
#                     score=score,
#                     var=var,
#                     avg_scores=ave,
#                     annotations=anot,
#                     figsize=(20, 16),
#                     fig_title=f"Average {score_txt} Scores for Fingerprint-numerical Predicting {target} using {model_in_title} model",
#                     x_title="Fingerprint-numerical Representations",
#                     y_title="Polymer Unit Representation",
#                     fname=fname
#                     )

# for fixed PU with different models
def create_structural_scaler_result(target_dir:Path,
                                    # regressor_model:str,
                                    target:str,
                                    score:str,
                                    var:str,
                                    data_type:str,
                                    transformer_type:str,
                                    peak_num:int=None
                                    ) -> None:

    ave, anot, model = creat_result_df(target_dir=target_dir,score=score, var=var,data_type=data_type
                                       , transformer_type=transformer_type,
                                       peak_number=peak_num)
    model_in_title:str =  ",".join(model)
    score_txt: str = "$R^2$" if score == "r2" else score.upper()
    # fname = f'{fname} on peak {peak_num+1}' if peak_num else fname
    fname= f"all PolymerRepresentation vs all features search heatmap_{score} score"
    print(ave)
    _create_heatmap(root_dir=target_dir,
                    score=score,
                    var=var,
                    avg_scores=ave,
                    annotations=anot,
                    figsize=(12,6),
                    fig_title=f"\n",
                    x_title="Feature Space",
                    y_title="Regression Models",
                    fname=fname,
                    num_ticks=5,
                    vmin=0.2,
                    vmax=0.6,
                    )

        
        
        
        
        
complex_models = ['XGBR','RF', 'NGB']

for transformer in transformer_list:
    for model in complex_models: 
        for target_folder in target_list:
            for i in scores_list:
#                 create_structural_scaler_result(target_dir=RESULTS/target_folder,regressor_model= model,target=f'{target_folder} with',
#                                                 score=i,var='stdev',data_type='structural_scaler', transformer_type=transformer)
                create_structural_result(target_dir=RESULTS/target_folder,regressor_model= model,target=f'{target_folder}',
                                            score=i,var='stdev',data_type='structural', transformer_type=transformer)

# for peak in [0,1,2]:
#     for transformer in transformer_list:
#         for model in complex_models: 
#             for target_folder in target_list:
#                 for i in scores_list:
#                     # create_structural_scaler_result(target_dir=RESULTS/target_folder,regressor_model= model,target=f'{target_folder} with',
#                     #                                 score=i,var='stdev',data_type='structural_scaler', transformer_type=transformer,peak_num=peak)
                    
#                     create_structural_result(target_dir=RESULTS/target_folder,regressor_model= model,target=f'{target_folder} with',
#                                                 score=i,var='stdev',data_type='structural', transformer_type=transformer,peak_num=peak)



def create_scaler_result(target_dir:Path,
                        score:str,
                        target:str,
                        var:str,
                        data_type:str,
                        transformer_type:str,
                        peak_num:int=None
                        )->None:

    ave, anot, model = creat_result_df(target_dir=target_dir,score=score, var=var,data_type=data_type,
                                       regressor_model=None,transformer_type=transformer_type,
                                       peak_number=peak_num)
    model_in_title:str =  ",".join(model)
    score_txt: str = "$R^2$" if score == "r2" else score.upper()
    # reg_name = f'{regressor_model} on {peak_num}' if peak_num else regressor_model
    fname= f"Regression Models vs numerical features with {transformer_type} search heatmap_{score}"
    fname = f'{fname} on peak {peak_num+1}' if peak_num else fname
    if score == "r2":
        vmax= .6
        vmin= .4
    elif score == "mae":
        vmax= .55
        vmin= 0.45
    
    elif score == "rmse":
        vmax= .65
        vmin= 0.50 
    
    _create_heatmap(root_dir=target_dir,
                    score=score,
                    var=var,
                    avg_scores=ave,
                    annotations=anot,
                    figsize=(24, 14),
                    fig_title=f"Average {score_txt} of {target.split()[1]} with continuous features",
                    x_title="Feature Space",
                    y_title="Regression Models",
                    fname=fname,
                    model_order=['RF', 'DT', 'MLR'],
                    num_ticks=3,
                    vmax=vmax,
                    vmin=vmin,
                    )


# simple_models = ['MLR','DT','RF']

# for transformer in transformer_list:
#     for target_folder in target_list:
#         for i in scores_list:
            # create_scaler_result(target_dir=RESULTS/target_folder,target=f'{target_folder}',
            #                     score=i,var='stdev',data_type='scaler',transformer_type=transformer)
            

            # create_structural_scaler_result(target_dir=RESULTS/target_folder,target=f'{target_folder}',
            #                                     score=i,var='stdev',data_type='structural_scaler', transformer_type=transformer)
            
# for peak in [0,1,2]:
#     for transformer in transformer_list:
#         for target_folder in target_list:
#             for i in scores_list:
#                 create_scaler_result(target_dir=RESULTS/target_folder,target=f'{target_folder} with',
#                                     score=i,var='stdev',data_type='scaler',transformer_type=transformer,peak_num=peak)



