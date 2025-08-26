import pandas as pd
from pathlib import Path
from training_utils import train_regressor
from all_factories import radius_to_bits,cutoffs
from typing import Callable, Optional, Union, Dict, Tuple
import numpy as np
import sys
sys.path.append("../cleaning")
from argparse import ArgumentParser
from data_handling import save_results
from train_structure_numerical import parse_arguments

HERE: Path = Path(__file__).resolve().parent
DATASETS: Path = HERE.parent.parent / "datasets"
RESULTS = Path = HERE.parent.parent / "results"

training_df_dir: Path = DATASETS/ "training_dataset"/ "Rg data with clusters.pkl"
w_data = pd.read_pickle(training_df_dir)

TEST = False

def main_numerical_only(
    dataset: pd.DataFrame,
    regressor_type: str,
    target_features: list[str],
    feat_transformer: str,
    hyperparameter_optimization: bool,
    numerical_feats: Optional[list[str]],
    kernel:str=None,
    cutoff:Optional[str]=None,
    target_transformer:str=None,
    classification:bool=False,
) -> None:


    scores, predictions,ls  = train_regressor(
                                            dataset=dataset,
                                            structural_features=None,
                                            unroll=None,
                                            numerical_feats=numerical_feats,
                                            target_features=target_features,
                                            regressor_type=regressor_type,
                                            kernel=kernel,
                                            cutoff=cutoff,
                                            hyperparameter_optimization=hyperparameter_optimization,
                                            feat_transformer=feat_transformer,
                                            target_transformer=target_transformer,
                                            Test=TEST,
                                            classification=classification,
                                            )
    
    save_results(scores,
                predictions=predictions,
                feat_length_scale=ls,
                representation= None,
                pu_type= None,
                target_features=target_features,
                regressor_type=regressor_type,
                kernel=kernel,
                numerical_feats=numerical_feats,
                cutoff=cutoffs,
                TEST=TEST,
                hypop=hyperparameter_optimization,
                transform_type=feat_transformer,
                target_transformer=target_transformer,
                classification=classification,
                # special_folder_name='hp_RF_differences'
                )


    # columns_to_impute: list[str] = ["PDI","Temperature SANS/SLS/DLS/SEC (K)","Concentration (mg/ml)"]
    # special_column: str = "Mw (g/mol)"
    # numerical_feats: list[str] = ["Mn (g/mol)", "Mw (g/mol)", "PDI", "Temperature SANS/SLS/DLS/SEC (K)","Concentration (mg/ml)"]
    # imputer = "mean"
    # transform_type= "Standard"
    # target_features= ['Lp (nm)']
    


if __name__ == "__main__":
    # if TEST==False:

    #     args = parse_arguments()
    #     main_numerical_only(
    #         dataset=w_data,
    #         regressor_type=args.regressor_type,
    #         kernel=args.kernel,
    #         target_features=[args.target_features],  
    #         transform_type='Standard',
    #         hyperparameter_optimization=True,
    #         columns_to_impute=args.columns_to_impute,  
    #         special_impute=args.special_impute,
    #         numerical_feats=args.numerical_feats,  
    #         imputer=args.imputer,
    #         cutoff=None,  
    #         second_transformer=None,
    #         classification=False
    #     )
    # else:
        main_numerical_only(
            dataset=w_data,
            regressor_type="sklearn-GPR",
            kernel= "a_matern",
            target_features=['log Rg (nm)'],
            feat_transformer='Standard',
            target_transformer='Standard', 
            hyperparameter_optimization=False,
            numerical_feats=[
                            'Xn', 'Mw (g/mol)', 'PDI', "Concentration (mg/ml)", 
                            "Temperature SANS/SLS/DLS/SEC (K)", 
                            "polymer dP", "polymer dD", "polymer dH",
                            'solvent dP', 'solvent dD', 'solvent dH'
                             ],
            )

    # columns_to_impute: list[str] = ["PDI","Temperature SANS/SLS/DLS/SEC (K)","Concentration (mg/ml)"]
    # special_column: str = "Mw (g/mol)"
    # numerical_feats: list[str] = ["Mn (g/mol)", "Mw (g/mol)", "PDI", "Temperature SANS/SLS/DLS/SEC (K)","Concentration (mg/ml)"]

# "intensity weighted average over log(Rh (nm))"



