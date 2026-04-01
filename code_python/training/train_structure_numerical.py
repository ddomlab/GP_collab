import pandas as pd
from pathlib import Path

import torch
from filter_data import _get_dataset_features, get_structural_info
from training_utils_gp_special import train_regressor
from typing import Callable, Optional, Union, Dict, Tuple
import numpy as np
import sys
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

TEST = False


def main_structural_numerical(
    dataset: pd.DataFrame,
    regressor_type: str,
    target_features: list[str],
    feat_transformer: str,
    target_transformer:str,
    representation:str,
    polymer_unit: str,
    numerical_feats: Optional[list[str]]=None,
    hyperparameter_optimization: bool=False,
    columns_to_impute: Optional[list[str]]=None,
    special_impute: Optional[str]=None,
    radius:int=None,
    vector:str=None,
    imputer:Optional[str]=None,
    **kwargs,
) -> None:

    structural_features, unroll_single_feat = get_structural_info(representation, polymer_unit, radius, vector)
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
                                            Test=TEST,
                                            **kwargs,
                                            )
  
    save_results(scores,
                predictions=predictions,
                representation=unroll_single_feat["representation"],
                target_features=target_features,
                regressor_type=regressor_type,
                numerical_feats=numerical_feats,
                radius= radius,
                vector =vector,
                hypop=hyperparameter_optimization,
                transform_type=feat_transformer,
                target_transformer=target_transformer,
                imputer=imputer,
                output_dir_name= PAPER,
                TEST=TEST,
                # special_folder_name='hp_RF_differences',
                special_file_name='noPriors',
                **kwargs,
                )


if __name__ == "__main__":
    args = parse_arguments()
    PAPER = args.paper
    # # "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation"
    dataset_name = args.dataset
    # # "cleaned_dataset_pervaporation_membranes_wang"
    # #non_imputed_dropped_nan_Rg_data
    w_data, feats, all_targets, polymer_unit = _get_dataset_features(DATASETS, PAPER, dataset_name)

    for targ in all_targets:
                main_structural_numerical(
                    dataset=w_data,
                    representation="MG",
                    # radius=3,
                    # vector="count",
                    regressor_type="MGK",
                    # GPytorchMixMCMC
                    # GPMixMCMC
                    polymer_unit=polymer_unit,
                    target_features=[targ],  
                    feat_transformer=None,
                    target_transformer='Standard',
                    numerical_feats=feats,
                    hyperparameter_optimization=False,
                    kernel_type={
                        "fp": "Graph",
                        "count": "RBF"
                        # "fp":args.K_fp,
                        # "count":args.K_count
                                },
                    kernel_mixing_method=args.Kernel_mixing_method,
                    # imputer="mean",
                    # columns_to_impute=['P_MW','surface tension (mN/m)','pore maker molecular weight (Da)','organic compound size (Da)','solubility parameter (MPa1/2)',]
                    )

 
        # 'log (Total flux)',
        # 'log (Separation factor)'

    # print(w_data['log (Total flux)'].describe())
    # 'water permeability (LMH/bar)',
    # 'organic compound removal (%)',
    # 'flux decline ratio (%)',
    # 'flux recovery ratio (%) ',
    # 'reversible fouling ratio (%)',
    # 'irreversible fouling ratio(%)',
    # import torch
    # train_x = torch.linspace(0, 1, 4)
    # print(train_x)

    # PAPER = "Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation"
    # # "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation"
    # dataset_name = "Rg data with clusters aging imputed"
    # # "cleaned_dataset_pervaporation_membranes_wang"
    # #non_imputed_dropped_nan_Rg_data
    # w_data, feats, all_targets, polymer_unit = _get_dataset_features(DATASETS, PAPER, dataset_name)
    # for targ in all_targets:
    #     # GPMixMCMC
    #     # GPytorchMCMC
    #     # GPytorchMAP

    #     # for model in ["RF", "XGBR", "NGB"]:
    #         # print(w_data["Monomer SMILES"])
    #         main_structural_numerical(
    #             dataset=w_data,
    #             representation="SMILES",
    #             # radius=3,
    #             # vector="count",
    #             regressor_type="GPytorchMAP",
    #             # GPytorchMixMCMC
    #             # GPMixMCMC
    #             polymer_unit=polymer_unit,
    #             target_features=[targ],  
    #             feat_transformer='Standard',
    #             target_transformer='Standard',
    #             numerical_feats=feats,
    #             hyperparameter_optimization=False,
    #             kernel_type={"fp":"SSK",
    #                          "count":"RBF"
    #                          },
    #             kernel_mixing_method="product",
    #             # imputer="mean",
    #             # columns_to_impute=['P_MW','surface tension (mN/m)','pore maker molecular weight (Da)','organic compound size (Da)','solubility parameter (MPa1/2)',]
    #             )





