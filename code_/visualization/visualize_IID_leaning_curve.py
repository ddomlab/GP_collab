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
from visualize_ood_learning_curve import ensure_long_path
from visualization_setting import set_plot_style

HERE: Path = Path(__file__).resolve().parent
RESULTS: Path = HERE.parent.parent/ 'results'

# target_dir: Path = RESULTS/'target_Lp'


# with open(filters, "r") as f:
#     FILTERS: dict = json.load(f)

# score_bounds: dict[str, int] = {"r": 1, "r2": 1, "mae": 7.5, "rmse": 7.5}
var_titles: dict[str, str] = {"stdev": "Standard Deviation", "stderr": "Standard Error"}
target_list = ['target_log Rg (nm)']
models = [
    'XGBR'
        #   ,'RF',
        # 'NGB'
          ]



set_plot_style()


# file = Path(r'C:\Users\sdehgha2\Desktop\PhD code\pls-dataset-project\PLS-Dataset\results\target_Lp\Trimer\(ECFP8)_count_1024_RF_scores.json')
# print(file.name)

def clean_array(score_list):
    """
    Extracts the first element from the list if it's nested, 
    otherwise returns the value as is.
    
    Parameters:
    score_list (list or float): The score or list containing the score.
    
    Returns:
    np.array: Cleaned array of scores.
    """
    return np.array([score[0] if isinstance(score, list) else score for score in score_list])

def get_learning_data(
    file_path: Path,
) -> tuple[float, float]:
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
        learning_df = None
    else:
        # for just scaler features
        if "scaler" == file_path.parent.name:
            model:str = file_path.name.split("_")[1] 
            features:str = file_path.name.split("_")[0].replace("(", "").replace(")", "")
        # for mixture of scaler and structural 
        elif "scaler" in file_path.parent.name and file_path.parent.name != "scaler":
            features:str = file_path.name.split("_")[0].replace("(", "").replace(")", "")
            model:str = file_path.name.split("_")[1]
        # for structural only
        else:
            features:str = file_path.name.split("_")[0].replace("(", "").replace(")", "")
            model:str = file_path.name.split("_")[1]

       
        with open(file_path, "r") as f:
            data = json.load(f)
        train_sizes = data['aggregated_results']["train_sizes_fraction"]
        train_scores_mean = clean_array(data['aggregated_results']["train_scores_mean"])
        train_scores_std = clean_array(data['aggregated_results']["train_scores_std"])
        test_scores_mean = clean_array(data['aggregated_results']["test_scores_mean"])
        test_scores_std = clean_array(data['aggregated_results']["test_scores_std"])
        learning_df = pd.DataFrame({
                "train_sizes_fraction": train_sizes,
                "train_scores_mean": train_scores_mean,
                "train_scores_std": train_scores_std,
                "test_scores_mean": test_scores_mean,
                "test_scores_std": test_scores_std
            })
    return features, model, learning_df



def _save_path(features, result_type, model, root_dir, polymer_representation,hyp_status):
    visualization_folder_path =  root_dir/result_type/polymer_representation
    
    os.makedirs(visualization_folder_path, exist_ok=True)
    fname = f"({model})"
    fname = f"{fname}_hypeOFF" if hyp_status==False else fname 
    saving_path = visualization_folder_path/ f"{fname}.png"
    return saving_path




def _create_learning_curve(
    root_dir: Path,
    features:str,
    model:str,
    hyp_status:Optional[str],
    poly_representation:str,
    df:pd.DataFrame,
    # x_labels: list[str],
    # y_labels: list[str],
    figsize: tuple[int, int],
    fig_title: str,
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

    
    plt.figure(figsize=figsize)

    # Plot training scores
    ax = sns.lineplot(x="train_sizes_fraction", y="train_scores_mean", data=df, label="Training Score", marker='o', markersize=10,linewidth=5,color="#FF9637")
    sns.lineplot(x="train_sizes_fraction", y="test_scores_mean", data=df, label="Test Score", marker='v', markersize=10,linewidth=5,color="#418BBF")


    plt.fill_between(df["train_sizes_fraction"],
                 df["train_scores_mean"] - df["train_scores_std"],
                 df["train_scores_mean"] + df["train_scores_std"],
                 color="#FFA657", alpha=0.2)

    plt.fill_between(df["train_sizes_fraction"],
                 df["test_scores_mean"] - df["test_scores_std"],
                 df["test_scores_mean"] + df["test_scores_std"],
                 color="#4F92C5", alpha=0.2)


    # plt.title(fig_title, fontsize=22, fontweight='bold')  # Plot title with size and bold
    ax.set_xlabel("Training set size", fontsize=24, fontweight='bold')  # X-axis label
    ax.set_ylabel("$R^2$", fontsize=24, fontweight='bold')  # Y-axis label

    # Customizing tick labels
    ax.tick_params(axis='both', which='major', labelsize=21)  # Set tick label size
    # plt.xticks(fontsize=20)  # Set x-tick labels bold
    # plt.yticks(fontsize=20)  # Set y-tick labels bold

    plt.legend(fontsize=21)

    print(hyp_status)
    # make folder and save the files
    saving_file_path =_save_path(features=features,
                                       result_type=f"learning curves of {model}",
                                       model=model, 
                                       root_dir=root_dir, 
                                       hyp_status=hyp_status,
                                       polymer_representation=poly_representation)

    plt.tight_layout()
    plt.savefig(saving_file_path, dpi=600)

    # Show the heatmap
    # plt.show()
    plt.close()



def save_learning_curve(target_dir: Path,
                    # data_type:str,
                    # regressor_model:Optional[str]=None,
) -> tuple[pd.DataFrame,pd.DataFrame]:
    

    pattern: str = "*lc_scores.json"
    for representation in os.listdir(target_dir):
        score_files = []
            
        representation_dir = os.path.join(target_dir, representation)
        score_files: list[Path] = list(Path(representation_dir).rglob(pattern))
        
        for file_path in score_files:
            # for structural and mix of structural-scaler
                poly_representation_name = file_path.parent.name
                print(poly_representation_name)
                feats, model, learning_score_data = get_learning_data(file_path=file_path)  
                hyper_param_status = True if "hypOFF" not in file_path.name else False  
                _create_learning_curve(target_dir,
                                       feats,
                                       model,
                                       hyper_param_status,
                                       poly_representation_name,
                                       learning_score_data,
                                       figsize=(8, 6),
                                       fig_title ="Learning Curve"
                                       )


# (Mordred-Mw-PDI-concentration-temperature-polymer dP-polymer dD-polymer dH-solvent dP-solvent dD-solvent dH)_NGB_hypOFF_Standard_lc_scores.json

# for target_folder in target_list:
#     save_learning_curve(target_dir=RESULTS/target_folder)



# manual
# f_name = "(Xn-Mw-PDI-concentration-temperature-polymer dP-polymer dD-polymer dH-solvent dP-solvent dD-solvent dH-light exposure-aging time-aging temperature-prep temperature-prep time)_XGBR_hypOFF_Standard_lc_scores.json"
# poly_representation_name = 'scaler'
# target_to_an = 'target_log Rg (nm)'
# lc_file:Path = RESULTS/target_to_an/poly_representation_name/f_name
# lc_file = ensure_long_path(lc_file)
# feats, model, learning_score_data = get_learning_data(file_path=lc_file)  
# print(feats, model)
# hyper_param_status = True if "hypOFF" not in lc_file.name else False  
# _create_learning_curve(
#                     RESULTS/target_to_an,
#                     feats,
#                     model,
#                     hyper_param_status,
#                     poly_representation_name,
#                     learning_score_data,
#                     figsize=(8, 6),
#                     fig_title ="Learning Curve"
#                     )




# def create_structural_result(target_dir:Path,
#                              regressor_model:str,
#                              target:str,
#                              score:str,
#                              var:str,
#                              data_type:str
#                              ) -> None:
#     ave, anot, model = creat_result_df(target_dir=target_dir,score=score, var=var,data_type=data_type
#                                        ,regressor_model=regressor_model)
#     model_in_title:str =  ",".join(model)
#     score_txt: str = "$R^2$" if score == "r2" else score.upper()
#     _create_heatmap(root_dir=target_dir,
#                     score=score,
#                     var=var,
#                     avg_scores=ave,
#                     annotations=anot,
#                     figsize=(12, 8),
#                     fig_title=f"Average {score_txt} Scores for Fingerprint Predicting {target} using {model_in_title} model(s)",
#                     x_title="Fingerprint Representations",
#                     y_title="Polymer Unit Representation",
#                     fname=f"PolymerRepresentation vs Fingerprint trained by {regressor_model} search heatmap_{score} score")


# scores_list: list = {"r", "r2", "mae", "rmse"}
# var_titles: dict[str, str] = {"stdev": "Standard Deviation", "stderr": "Standard Error"}
# for i in scores_list:
#     create_structural_result(target_dir=target_dir,target='Lp (nm) with filteration on concentation and Lp',score=i,var='stdev',data_type='structural')


# feat, model, av, std = get_results_from_file(file,score='r2', var='stdev')



# for model in models: 
#     for target_folder in target_list:
#         for i in scores_list:
#             create_structural_result(target_dir=RESULTS/target_folder,regressor_model= model,target=f'{target_folder} with',
#                                             score=i,var='stdev',data_type='structural')





# def create_structural_scaler_result(target_dir:Path,
#                                     regressor_model:str,
#                                     target:str,
#                                     score:str,
#                                     var:str,
#                                     data_type:str
#                                     ) -> None:

#     ave, anot, model = creat_result_df(target_dir=target_dir,score=score, var=var,data_type=data_type,
#                                        regressor_model=regressor_model)
#     model_in_title:str =  ",".join(model)
    
#     score_txt: str = "$R^2$" if score == "r2" else score.upper()
#     _create_heatmap(root_dir=target_dir,
#                     score=score,
#                     var=var,
#                     avg_scores=ave,
#                     annotations=anot,
#                     figsize=(20, 16),
#                     fig_title=f"Average {score_txt} Scores for Fingerprint-numerical Predicting {target} using {model_in_title} model",
#                     x_title="Fingerprint-numerical Representations",
#                     y_title="Polymer Unit Representation",
#                     fname=f"PolymerRepresentation vs (Fingerprint-numerical) trained by {regressor_model} search heatmap_{score} score")
    



# for model in models: 
#     for target_folder in target_list:
#         for i in scores_list:
#             create_structural_scaler_result(target_dir=RESULTS/target_folder,regressor_model= model,target=f'{target_folder} with',
#                                             score=i,var='stdev',data_type='structural_scaler')




# def create_scaler_result(target_dir:Path,
#                         score:str,
#                         target:str,
#                         var:str,
#                         data_type:str
#                         )->None:

#     ave, anot, model = creat_result_df(target_dir=target_dir,score=score, var=var,data_type=data_type,
#                                        regressor_model=None)
#     model_in_title:str =  ",".join(model)
#     score_txt: str = "$R^2$" if score == "r2" else score.upper()
#     _create_heatmap(root_dir=target_dir,
#                     score=score,
#                     var=var,
#                     avg_scores=ave,
#                     annotations=anot,
#                     figsize=(20, 10),
#                     fig_title=f"Average {score_txt} Scores for numerical Predicting {target} using {model_in_title} model",
#                     x_title="numerical Representations",
#                     y_title="Regression Models",
#                     fname=f"Regression Models vs numerical features search heatmap_{score}")
    

# for target_folder in target_list:
#     for i in scores_list:
#         create_scaler_result(target_dir=RESULTS/target_folder,target=f'{target_folder} with',
#                              score=i,var='stdev',data_type='scaler')