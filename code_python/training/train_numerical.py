import pandas as pd
from pathlib import Path

from filter_data import _get_dataset_features, get_structural_info
from training_utils_gp_special import train_regressor
from typing import Callable, Optional, Union, Dict, Tuple
import numpy as np
from data_handling import save_results
from utils import parse_arguments
from train_structure_numerical import main_structural_numerical
# import sys

# sys.modules.setdefault("numpy._core",         np.core)
# sys.modules.setdefault("numpy._core.numeric", np.core.numeric)
# sys.modules.setdefault("numpy._core.multiarray", np.core.multiarray)
# sys.modules.setdefault("numpy._core.umath",   np.core.umath)



## add feature importance using gini or shap
HERE = Path(__file__).resolve().parent
DATASETS = HERE.parent.parent / "datasets" / "Validation datasets"
RESULTS = HERE.parent.parent / "results"

DEBUGGING = False

if __name__ == "__main__":
    if DEBUGGING is False:
        args = parse_arguments()
        PAPER = args.paper
        dataset_name =args.dataset
        w_data, feats, all_targets, _ = _get_dataset_features(DATASETS, PAPER, dataset_name)

        for targ in all_targets:
            main_structural_numerical(
                dataset=w_data,
                output_dir_name=PAPER,
                representation=None,
                polymer_unit=None,
                regressor_type=args.regressor_type, #"GPytorchMAP", #MGK-sklearn "GpyroHMC"
                target_features=[targ],  
                feat_transformer="Standard",
                target_transformer="Standard",
                numerical_feats=feats,
                hyperparameter_optimization=False,
                kernel_type={
                    "count": args.K_count #"Matern32"   #args.K_count
                    },
                kernel_mixing_method=args.Kernel_mixing_method,
                use_cuda=True, #True for GPU, False for CPU
                )


    # else:
    #     ##### DEBUGGING #####
    #     PAPER = "Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation"
    #     dataset_name ="Rg data with clusters aging imputed"
    #     w_data, feats, all_targets, polymer_unit = _get_dataset_features(DATASETS, PAPER, dataset_name)

    #     for targ in all_targets:
    #         # optuna_save_dir = RESULTS/PAPER/f"target_{targ}"/ "MGK_hyperprameters"/f"Graph_Matern32_{args.kernel_feature_mode}"
    #         main_structural_numerical(
    #             dataset=w_data,
    #             representation="ECFP", #"ECFP" # MG #SSK
    #             radius=3,
    #             vector="count",
    #             regressor_type="RF", #"GPytorchMAP", #MGK-sklearn "GpyroHMC"
    #             polymer_unit=polymer_unit,
    #             target_features=[targ],  
    #             feat_transformer="Standard",
    #             target_transformer="Standard",
    #             numerical_feats=feats,
    #             hyperparameter_optimization=False,
    #             # kernel_type={
    #             #     "fp": "TanimotoRBF", #Graph #"TanimotoRBF" #args.K_fp
    #             #     "count": "RBF" #"Matern32"   #args.K_count
    #             #     # "fp": "Graph", 
    #             #     # "count": "Matern32"
    #             #     },
    #             # kernel_mixing_method="product",
    #             # use_cuda=False, #True for GPU, False for CPU
    #             # kernel_feature_mode = args.kernel_feature_mode, #joint, #per_feature for MGK
    #             clustering_method="substructure cluster",
    #             )
            
