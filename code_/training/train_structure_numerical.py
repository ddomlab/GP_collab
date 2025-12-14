import pandas as pd
from pathlib import Path
from filter_data import _get_dataset_features, get_structural_info
from training_utils_gp_special import train_regressor
from typing import Callable, Optional, Union, Dict, Tuple
import numpy as np
import sys
from data_handling import save_results
from utils import parse_arguments

import sys

sys.modules.setdefault("numpy._core",         np.core)
sys.modules.setdefault("numpy._core.numeric", np.core.numeric)
sys.modules.setdefault("numpy._core.multiarray", np.core.multiarray)
sys.modules.setdefault("numpy._core.umath",   np.core.umath)



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
) -> None:

    structural_features, unroll_single_feat = get_structural_info(representation,polymer_unit,radius,vector)
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
                                            )
  
    save_results(scores,
                predictions=predictions,
                representation=representation,
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
                special_file_name='sum_chain1',
                )



if __name__ == "__main__":
    PAPER = "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices"
    #non_imputed_dropped_nan_Rg_data
    w_data,feats, all_targets = _get_dataset_features(DATASETS, PAPER, "Beyond molecular structure_seifrid_imputed")

    # args = parse_arguments()

    # main_structural_numerical(
    #     dataset=w_data,
    #     representation=args.representation,
    #     radius=args.radius,
    #     vector=args.vector,
    #     polymer_unit=args.polymer_unit,
    #     regressor_type=args.regressor_type,
    #     kernel=args.kernel,
    #     target_features=[args.target_features],  
    #     feat_transformer='Standard',
    #     second_transformer=None,
    #     hyperparameter_optimization=True,
    #     columns_to_impute=args.columns_to_impute,  
    #     special_impute=args.special_impute,
    #     numerical_feats=args.numerical_feats,  
    #     imputer=args.imputer,
    #     cutoff=None,  
    # )
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
    for targ in all_targets:
        # GPMixMCMC
        # for model in ["RF", "XGBR", "NGB"]:
            
            main_structural_numerical(
                dataset=w_data,
                representation="ECFP",
                radius=3,
                vector="count",
                regressor_type="GPMixMCMC",
                # GPytorchMixMCMC
                # GPMixMCMC
                # kernel="matern32_j_rbf_mix",
                polymer_unit=["Donor", "Acceptor"],
                target_features=[targ],  
                feat_transformer='Standard',
                target_transformer='Standard',
                numerical_feats=feats,
                hyperparameter_optimization=False,
                # imputer="mean",
                # columns_to_impute=['P_MW','surface tension (mN/m)','pore maker molecular weight (Da)','organic compound size (Da)','solubility parameter (MPa1/2)',]
                )





