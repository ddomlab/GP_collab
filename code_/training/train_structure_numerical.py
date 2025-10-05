import pandas as pd
from pathlib import Path
from filter_data import _get_dataset
from training_utils_gp_special import train_regressor
from all_factories import radius_to_bits,cutoffs
from typing import Callable, Optional, Union, Dict, Tuple
import numpy as np
import sys
from argparse import ArgumentParser
from data_handling import save_results

## add feature importance using gini or shap
HERE = Path(__file__).resolve().parent
DATASETS = HERE.parent.parent / "datasets" / "Validation datasets"
PAPER = "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices"
RESULTS = HERE.parent.parent / "results"
w_data = _get_dataset(DATASETS / PAPER, "Beyond molecular structure_seifrid_imputed")

TEST = False


def get_appropriate_features(PAPER:str):
     
    if PAPER == "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices":
        return [
                "HOMO_D (eV)", "LUMO_D (eV)", "Eg_D (eV)", "Ehl_D (eV)",
                "HOMO_A (eV)", "LUMO_A (eV)", "Eg_A (eV)", "Ehl_A (eV)",
                "D:A ratio (m/m)", "solvent additive conc. (% v/v)", 
                "temperature of thermal annealing",
                "HTL energy level (eV)", "ETL energy level (eV)"
                ]

    
    if PAPER == "Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation":
        return [
                'Xn', 'Mw (g/mol)', 'PDI', "Concentration (mg/ml)", 
                "Temperature SANS/SLS/DLS/SEC (K)", 
                "polymer dP", "polymer dD", "polymer dH",
                'solvent dP', 'solvent dD', 'solvent dH'
                ],
    
    else:
        print("No predefined features for this dataset. Please specify numerical features manually.")
        sys.exit(1)


def get_structural_info(fp:str,poly_unit_name:list[str],radius:int=None,vector:str=None)->Tuple:
       
        if fp == "Mordred" or fp == "MACCS":
            fp_features = [f"{unit}_{fp}" for unit in poly_unit_name]
            unrolling_featurs = {"representation": fp,
                                "unit_name":poly_unit_name,
                                "col_names": fp_features}
            return fp_features, unrolling_featurs
        
        if fp == "ECFP":
            n_bits = radius_to_bits[radius]
            fp_features = [
                        f"{unit}_{fp}{2 * radius}_{vector}_{n_bits}bits"
                        for unit in poly_unit_name
                        ]
            unrolling_featurs = {
                                "representation": fp,
                                "radius": radius,
                                "n_bits": n_bits,
                                "vector_type": vector,
                                "unit_name": poly_unit_name,
                                "col_names": fp_features,
                                }
            return fp_features, unrolling_featurs
        else:
              return None, None


def main_structural_numerical(
    dataset: pd.DataFrame,
    regressor_type: str,
    target_features: list[str],
    feat_transformer: str,
    target_transformer:str,
    representation:str,
    polymer_unit: str,
    numerical_feats: Optional[list[str]],
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
                output_dir_name= PAPER,
                TEST=TEST,

                # special_folder_name='hp_RF_differences',
                # special_file_name='v2_(max_feat_all_leaf_smaller)',
                )



def parse_arguments():
    parser = ArgumentParser(description="Process some data for numerical-only regression.")
    
    parser.add_argument(
        '--target_features',
        # choices=['Lp (nm)', 'Rg1 (nm)', 'Rh (IW avg log)'],  
        required=True,
        help="Specify a single target for the analysis."
    )

    parser.add_argument(
        '--regressor_type', 
        type=str, 
        choices=['RF', 'DT', 'MLR', 'SVR', 'XGBR','KNN', 'GPR', 'NGB', 'sklearn-GPR', 'MLP'], 
        required=True, 
        help="Regressor type required"
    )

    parser.add_argument(
        '--numerical_feats',
        type=str,
        choices=['Xn','Mn (g/mol)', 'Mw (g/mol)', 'PDI', 'Temperature SANS/SLS/DLS/SEC (K)',
                  'Concentration (mg/ml)','solvent dP',	'polymer dP',	'solvent dD',	'polymer dD',	'solvent dH',	'polymer dH', 'Ra',
                  "abs(solvent dD - polymer dD)", "abs(solvent dP - polymer dP)", "abs(solvent dH - polymer dH)"],

        nargs='+',  # Allows multiple choices
        required=None,
        help="Numerical features: choose"
    )
    
    parser.add_argument(
        '--columns_to_impute',
        type=str,
        choices=['Xn','Mn (g/mol)', 'Mw (g/mol)', 'PDI', 'Temperature SANS/SLS/DLS/SEC (K)',
                  'Concentration (mg/ml)','solvent dP',	'polymer dP',	'solvent dD',	'polymer dD',	'solvent dH',	'polymer dH', 'Ra'],

        nargs='*',  # This allows 0 or more values
        default=None,  
        help="imputation features: choose"
    )

    parser.add_argument(
        '--imputer',
        choices=['mean', 'median', 'most_frequent',"distance KNN", None],  
        nargs='?',  # This allows the argument to be optional
        default=None,  
        help="Specify the imputation strategy or leave it as None."
    )

    parser.add_argument(
        '--special_impute',
        choices=['Mw (g/mol)', None],  
        nargs='?',  # This allows the argument to be optional
        default=None,  # Set the default value to None
        help="Specify the imputation strategy or leave it as None."
    )

    parser.add_argument(
        "--feat_transformer", 
        type=str, 
        choices=["Standard", "Robust Scaler"], 
        default= "Standard", 
        help="transform type required"
    )

    parser.add_argument(
        "--kernel", 
        type=str,
        default=None,
        help='kernel for GP is optinal'
    )

    parser.add_argument(
        '--representation', 
        type=str, 
        choices=['ECFP', 'MACCS', 'Mordred'], 
        required=None, 
        help="Fingerprint required"
    )

    parser.add_argument(
        '--polymer_unit', 
        type=str, 
        choices=['Monomer', 'Dimer', 'Trimer', 'RRU Monomer', 'RRU Dimer', 'RRU Trimer'], 
        required=None, 
        help="Fingerprint required"
    )

    parser.add_argument(
        '--radius',
        type=int,
        choices=[3, 4, 5, 6],
        nargs='?',  # This allows the argument to be optional
        default=None,  # Set the default value to None
        help='Radius for ECFP'
    )

    parser.add_argument(
        '--vector',
        type=str,
        choices=['count', 'binary'],
        nargs='?',  # This allows the argument to be optional
        default='count',  # Set the default value to None
        help='Type of vector (default: count)'
    )

    parser.add_argument(
        '--clustering_method',
        type=str,
        nargs='?',
        help='Type of clustering method'
    )


    return parser.parse_args()

if __name__ == "__main__":
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

    main_structural_numerical(
        dataset=w_data,
        representation="ECFP",
        radius=3,
        vector="count",
        regressor_type="RF",
        # kernel="matern32_j_rbf_mix",
        polymer_unit=["Donor","Acceptor"],
        target_features=['calculated PCE (%)'],  
        feat_transformer='Standard',
        target_transformer='Standard',
        numerical_feats=get_appropriate_features(PAPER),
        hyperparameter_optimization=False,
    )





