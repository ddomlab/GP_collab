import json
from math import ceil
from pathlib import Path
import os 
import numpy as np
# import cmcrameri.cm as cmc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import List, Optional
from matplotlib import rc
from visualize_ood_learning_curve import ensure_long_path


HERE: Path = Path(__file__).resolve().parent
DATASETS: Path = HERE.parent.parent / "datasets"
RESULTS: Path = HERE.parent.parent/ 'results'
target_list = ['target_Rg']

from visualization_setting import set_plot_style, save_img_path

set_plot_style()


def get_file_info(file_path:Path)-> tuple[str, str]:
    if not file_path.exists():
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    return file_path.stem , file_path.parent.name

def get_prediction_plot(results_directory: Path, ground_truth_file: Path, target:str)->None:
    w_data = pd.read_pickle(ground_truth_file)
    true_values: pd.DataFrame = w_data[target].dropna()
    # true_values: pd.Series = pd.read_csv(ground_truth_file)["calculated PCE (%)"]
    
    pattern: str = "*_predictions.csv"
    for representation in os.listdir(results_directory):
        score_files = []
            
        representation_dir = os.path.join(results_directory, representation)
        score_files: list[Path] = list(Path(representation_dir).rglob(pattern))
        
        for pred_path in score_files:
            file_name, polymer_representation = get_file_info(pred_path)
            print(file_name)
            r2_avg, r2_stderr = get_scores(pred_path)
            draw_predictions_plot(true_values,
                                    pred_path,
                                    r2_avg,
                                    r2_stderr,
                                    results_directory,
                                    polymer_representation,
                                    file_name)



def get_scores(dir: Path,peak_number:Optional[int],score_metrics:str) -> tuple[float, float]:
    score_path = dir.with_name(dir.stem.replace("_prediction", "_score") + ".json")
    with open(score_path, "r") as f:
        scores: dict = json.load(f)
        
    avg = scores[f"{score_metrics}_avg"]
    std = scores[f"{score_metrics}_stdev"]
    avg = avg[peak_number] if isinstance(avg, list) else avg
    std = std[peak_number] if isinstance(std, list) else std
    r2_avg = round(avg, 2)
    r2_stderr = round(std, 2)
    return r2_avg, r2_stderr


def draw_predictions_plot(
                        #    true_values: pd.Series,
                        target:str,
                        predictions: Path,
                        # r2_avg: float,
                        # r2_stderr: float,
                        root_dir:Path,
                        poly_representation_name:str,
                        file_name,
                           ) -> None:
    # Load the data from CSV files
    
    predicted_values = pd.read_csv(predictions)
    seeds = predicted_values.drop(target, axis=1).columns 

    # There are 7 columns in predicted_values, each corresponding to a different seed
    # Create a Series consisting of the ground truth values repeated 7 times
    true_values_ext = pd.concat([predicted_values[target]] * len(seeds), ignore_index=True)
    # Create a Series consisting of the predicted values, with the column names as the index
    predicted_values_ext = pd.concat([predicted_values[col] for col in seeds], axis=0, ignore_index=True)

    ext_comb_df = pd.concat([true_values_ext, predicted_values_ext], axis=1)
    print(ext_comb_df)
    combined_data = np.log10(ext_comb_df+1)

    # print(ext_comb_df)
    # Create the hex-binned plot with value distributions for all y-axis columns

    g = sns.jointplot(data=combined_data, x=target, y=0,
                      kind="hex",
                    #   cmap="viridis",
                      # joint_kws={"gridsize": 50, "cmap": "Blues"},
                    #   joint_kws={"gridsize": (150,45)},
                      marginal_kws={"bins": 25},
                      )
    ax_max = ceil(max(combined_data.max()))
    # print(ax_max)
    g.ax_joint.plot([0, ax_max], [0, ax_max], ls="--", c=".3")
    # g.ax_joint.annotate(f"$R^2$ = {r2_avg} ± {r2_stderr}",
    #                     # xy=(0.1, 0.9), xycoords='axes fraction',
    #                     # ha='left', va='center',
    #                     # bbox={'boxstyle': 'round', 'fc': 'powderblue', 'ec': 'navy'}
    #                     )
    #  kwargs: linewidth=0.2, edgecolor='white',  mincnt=1
    # plt.text(0.95, 0.05, f"$R^2$ = {r2_avg} ± {r2_stderr}",
    #          horizontalalignment='right',
    #          verticalalignment='bottom',
    #          transform=g.ax_joint.transAxes,
    #          )
    # g.plot_marginals(sns.kdeplot, color="blue")
    # Set plot limits to (0, 15) for both axes
    g.set_axis_labels(f"log True {target} ", f"log Predicted {target}")
    g.ax_joint.set_xlim(0, ax_max)
    g.ax_joint.set_ylim(0, ax_max)
    # plt.tight_layout()
    # g.ax_joint.set_xscale("log")
    # g.ax_joint.set_yscale("log")
    # plt.show()


    visualization_folder_path =  root_dir/"parity plot"/poly_representation_name
    os.makedirs(visualization_folder_path, exist_ok=True)
    saving_path = visualization_folder_path/ f"{file_name}.png"

    # output: Path = predictions.parent / f"{predictions.stem}_plot.png"
    # plt.savefig(output, dpi=600)
    try:
        g.savefig(saving_path, dpi=600)
        print(f"Saved {saving_path}")
    except Exception as e:
        print(f"Failed to save {saving_path} due to {e}")

    # plt.show() (if needed)
    plt.close()


def draw_single_parity_plot(
                        target:str,
                        prediction_path: Path,
                        # r2_avg: float,
                        # r2_stderr: float,
                        root_dir:Path,
                        poly_representation_name:str,
                        file_name:str,
                        log_scale:bool,
                        x_y_target_name:str,
                        peak_num:Optional[int]=None,
                        score_metrics:str='r2',
                           ) -> None:
    # Load the data from CSV files
    r2_avg, r2_stderr = get_scores(prediction_path,peak_number=peak_num,score_metrics=score_metrics)
    predicted_values = pd.read_csv(prediction_path)
    seeds = predicted_values.drop(target, axis=1).columns 
    if peak_num:
        file_name = f"{file_name}_{peak_num}.png"  
        predicted_values = pd.DataFrame({
        col: predicted_values[col].values.reshape(-1, 3).tolist()  # Reshape each column
        for col in predicted_values.columns
        })
        predicted_values[target] = predicted_values[target].apply(lambda x: eval(x)[peak_num] if isinstance(x, str) else x[peak_num])

        predicted_values_ext = pd.concat([predicted_values[col].apply(lambda x: eval(x)[peak_num] if isinstance(x, str) else x[peak_num]) for col in seeds], axis=0, ignore_index=True)
    else:
        predicted_values_ext = pd.concat([predicted_values[col] for col in seeds], axis=0, ignore_index=True)
    true_values_ext = pd.concat([predicted_values[target]] * len(seeds), ignore_index=True)
    combined_data = pd.concat([true_values_ext, predicted_values_ext], axis=1)
    r_preason = combined_data.corr().iloc[0, 1]
    r_spearman = combined_data.corr(method='spearman').iloc[0, 1]
    if log_scale:
        combined_data = np.log10(combined_data)
    
    g = sns.jointplot(data=combined_data, x=target, y=0,
                      kind="hex",
                    #   cmap="Purples",
                      joint_kws={"gridsize": 19, "cmap": "Blues"},
                    #   joint_kws={"gridsize": (40,30)},
                      marginal_kws={"bins": 25},
                      )
    ax_max = ceil(max(combined_data.max()))
    ax_min = ceil(min(combined_data.min()))
    g.ax_joint.plot([0, ax_max], [0, ax_max], ls="--", c=".3")
    score_sign = '$R^2$' if score_metrics == 'r2' else '$RMSE$'
    g.ax_joint.annotate(f"{score_sign} = {r2_avg} ± {r2_stderr}",
                        xy=(0.1, 0.9), xycoords='axes fraction',
                        ha='left', va='center',
                        fontsize=16.5,
                        # bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray")
                        )

    # Set plot limits to (0, 15) for both axes
    xlabel = f"log True {x_y_target_name}" if log_scale else f"True {x_y_target_name}"
    ylabel = f"log Predicted {x_y_target_name}" if log_scale else f"Predicted {x_y_target_name}"

    g.ax_joint.set_xlabel(xlabel, fontdict={'fontsize': 18, 'fontweight': 'bold'})
    g.ax_joint.set_ylabel(ylabel, fontdict={'fontsize': 18, 'fontweight': 'bold'})
    g.ax_joint.set_xlim(ax_min, ax_max)
    g.ax_joint.set_ylim(ax_min, ax_max)
    # plt.tight_layout()
    # g.ax_joint.set_xscale("log")
    # g.ax_joint.set_yscale("log")


    visualization_folder_path =  root_dir/"parity plot"/poly_representation_name
    os.makedirs(visualization_folder_path, exist_ok=True)

    saving_path = visualization_folder_path/ file_name
    print(saving_path)
    try:
        save_img_path(visualization_folder_path, file_name)
    except Exception as e:
        print(f"Failed to save {saving_path} due to {e}")
    plt.show()
    plt.close()





if __name__ == "__main__":


    # training_df_dir: Path = DATASETS/ "training_dataset"/ "dataset_wo_block_cp_(fp-hsp)_added_additive_dropped.pkl"
    # for target_folder in target_list:
    #     get_prediction_plot(RESULTS/target_folder, training_df_dir,"Rg1 (nm)")


    target_to_an = 'target_log Rg (nm)'
    file_n = "(Xn-Mw-PDI-concentration-temperature-polymer dP-polymer dD-polymer dH-solvent dP-solvent dD-solvent dH-light exposure-aging time-aging temperature-prep temperature-prep time)_XGBR_Standard_predictions.csv"
    poly_representation_name = 'scaler'
    truth_val_file:Path = RESULTS/target_to_an/poly_representation_name/file_n
    truth_val_file = ensure_long_path(truth_val_file)
    print(truth_val_file.exists())
    fname,_ = get_file_info(truth_val_file)
    draw_single_parity_plot(
        target='log Rg (nm)',
        x_y_target_name='log Rg (nm)',
        prediction_path=truth_val_file,
        # r2_avg: float,
        # r2_stderr: float,
        root_dir=RESULTS/target_to_an,
        poly_representation_name=poly_representation_name,
        # peak_num=0,
        file_name=fname,
        log_scale=False,
        score_metrics='r2'
    )