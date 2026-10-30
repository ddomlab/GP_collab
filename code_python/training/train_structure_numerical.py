import pandas as pd
from pathlib import Path

from filter_data import _get_dataset_features, get_structural_info
from training_utils_gp_special import train_regressor
from typing import Callable, Optional, Union, Dict, Tuple
import numpy as np
from data_handling import save_results
from utils import parse_arguments
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


def main_structural_numerical(
    dataset: pd.DataFrame,
    regressor_type: str,
    target_features: list[str],
    feat_transformer: str,
    target_transformer:str,
    representation:str,
    polymer_unit: str,
    output_dir_name: str,
    numerical_feats: Optional[list[str]]=None,
    hyperparameter_optimization: bool=False,
    columns_to_impute: Optional[list[str]]=None,
    special_impute: Optional[str]=None,
    radius:int=None,
    vector:str=None,
    n_bits:int=None,
    imputer:Optional[str]=None,
    **kwargs,
) -> None:
    
    structural_features, unroll_single_feat = get_structural_info(
        representation,
        polymer_unit,
        radius,
        vector,
        n_bits,
    )
    scores, predictions  =   train_regressor(
                                            dataset=dataset,
                                            features_impute=columns_to_impute,
                                            special_impute=special_impute,
                                            unroll=unroll_single_feat,
                                            structural_features=structural_features,
                                            numerical_feats=numerical_feats,
                                            target_features=target_features,
                                            regressor_type=regressor_type,
                                            target_transformer=target_transformer,
                                            feat_transformer=feat_transformer,
                                            hyperparameter_optimization=hyperparameter_optimization,
                                            imputer=imputer,
                                            Test=DEBUGGING,
                                            **kwargs,
                                            )
  
    save_results(scores,
                predictions=predictions,
                representation=representation,
                target_features=target_features,
                regressor_type=regressor_type,
                numerical_feats=numerical_feats,
                radius= radius,
                vector =vector,
                n_bits=n_bits,
                hypop=hyperparameter_optimization,
                transform_type=feat_transformer,
                target_transformer=target_transformer,
                imputer=imputer,
                output_dir_name= output_dir_name,
                TEST=DEBUGGING,
                # special_folder_name='HDF',
                special_file_name=(
                    "SHAP" if kwargs.get("return_gp_shap", False) else None
                ),
                **kwargs,
                )


if __name__ == "__main__":
    if DEBUGGING is False:
        args = parse_arguments()
        PAPER = args.paper
        dataset_name = args.dataset
        model_name = args.regressor_type
        w_data, feats, all_targets, polymer_unit = _get_dataset_features(
            DATASETS,
            PAPER,
            dataset_name,
        )

        is_mgk = model_name == "MGK"
        is_kernel_mix_gp = model_name in {"GPytorchMAP", "GpyroHMC"}

        if is_mgk:
            representation_options = {
                "representation": "MG",
                "radius": None,
                "vector": None,
            }
        else:
            selected_representation = args.representation or "ECFP"
            representation_options = {
                "representation": selected_representation,
                "radius": (args.radius or 3)
                if selected_representation in {"ECFP", "HDF"} else None,
                "vector": args.vector
                if selected_representation == "ECFP" else None,
                "n_bits": (
                    64 if args.n_bits is None else args.n_bits
                ) if selected_representation == "HDF" else None,
            }

        model_options = {}
        if is_mgk:
            model_options = {
                "kernel_type": {
                    "fp": "Graph",
                    "count": args.K_count or "Matern32",
                },
                "kernel_mixing_method": (
                    args.Kernel_mixing_method or "product"
                ),
                "kernel_feature_mode": (
                    args.kernel_feature_mode or "per_feature"
                ),
                "use_cuda": True,
            }
        elif is_kernel_mix_gp:
            model_options = {
                "kernel_type": {
                    "fp": args.K_fp or "TanimotoRBF",
                    "count": args.K_count or "Matern32",
                },
                "kernel_mixing_method": (
                    args.Kernel_mixing_method or "product"
                ),
                "use_cuda": True,
                "return_gp_shap": args.return_gp_shap,
            }

        targets_to_train = all_targets
        if args.target_features:
            if args.target_features not in all_targets:
                raise ValueError(
                    f"Unknown target {args.target_features!r}. "
                    f"Available targets: {all_targets}"
                )
            targets_to_train = [args.target_features]

        for targ in targets_to_train:
            main_structural_numerical(
                dataset=w_data,
                regressor_type=model_name,
                polymer_unit=polymer_unit,
                target_features=[targ],
                feat_transformer="Standard",
                target_transformer="Standard",
                numerical_feats=feats,
                hyperparameter_optimization=False,
                output_dir_name=PAPER,
                permute_features=False,
                clustering_method=args.clustering_method,
                **representation_options,
                **model_options,
                )


    else:
        ##### DEBUGGING #####
        PAPER = "Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation"
        dataset_name ="Rg data with clusters aging imputed"
        w_data, feats, all_targets, polymer_unit = _get_dataset_features(DATASETS, PAPER, dataset_name)

        for targ in all_targets:
            # optuna_save_dir = RESULTS/PAPER/f"target_{targ}"/ "MGK_hyperprameters"/f"Graph_Matern32_{args.kernel_feature_mode}"
            main_structural_numerical(
                dataset=w_data,
                representation="HDF", #"ECFP" # MG #SSK
                n_bits=128,
                radius=3,
                # radius=3,
                # vector="count",
                regressor_type="RF", #"GPytorchMAP", #MGK-sklearn "GpyroHMC"
                polymer_unit=polymer_unit,
                target_features=[targ],  
                feat_transformer="Standard",
                target_transformer="Standard",
                numerical_feats=feats,
                hyperparameter_optimization=False,
                output_dir_name=PAPER,

                # kernel_type={
                #     "fp": "TanimotoRBF", #Graph #"TanimotoRBF" #args.K_fp
                #     "count": "RBF" #"Matern32"   #args.K_count
                #     # "fp": "Graph", 
                #     # "count": "Matern32"
                #     },
                # kernel_mixing_method="product",
                # use_cuda=False, #True for GPU, False for CPU
                # kernel_feature_mode = args.kernel_feature_mode, #joint, #per_feature for MGK
                # clustering_method="substructure cluster", #continuous_cluster  #structure_cluster
                )
            
