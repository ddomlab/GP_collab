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

HERE: Path = Path(__file__).resolve().parent
DATASETS: Path = HERE.parent.parent / "datasets"
RESULTS = Path = HERE.parent.parent / "results"

training_df_dir: Path = DATASETS/ "training_dataset"/"Rg data with clusters.pkl"
w_data = pd.read_pickle(training_df_dir)

TEST = False

def get_structural_info(fp:str,poly_unit:str,radius:int=None,vector:str=None)->Tuple:
       
        if fp == "Mordred" or fp == "MACCS":
            fp_features: list[str] = [f"{poly_unit}_{fp}"]
            unrolling_featurs = {"representation": fp,
                                "oligomer_representation":poly_unit,
                                "col_names": fp_features}
            return fp_features, unrolling_featurs
        
        if fp == "ECFP":
            n_bits = radius_to_bits[radius]
            fp_features: list[str] = [
            f"{poly_unit}_{fp}{2 * radius}_{vector}_{n_bits}bits"]
            unrolling_featurs = {
                                "representation": fp,
                                "radius": radius,
                                "n_bits": n_bits,
                                "vector_type": vector,
                                "oligomer_representation": poly_unit,
                                "col_names": fp_features,
                            }
            return fp_features, unrolling_featurs
        else:
              return None, None


def main_structural_numerical(
    dataset: pd.DataFrame,
    regressor_type: str,
    target_features: list[str],
    transform_type: str,
    representation:str,
    oligomer_representation: str,
    columns_to_impute: Optional[list[str]]=None,
    special_impute: Optional[str]=None,
    numerical_feats: Optional[list[str]]=None,
    hyperparameter_optimization: bool=True,
    radius:int=None,
    vector:str=None,
    kernel:str = None,
    imputer:Optional[str]=None,
    cutoff:Optional[str]=None,
    second_transformer:str=None,
) -> None:
    
    structural_features, unroll_single_feat = get_structural_info(representation,oligomer_representation,radius,vector)
    scores, predictions,ls  = train_regressor(
                                                    dataset=dataset,
                                                    features_impute=columns_to_impute,
                                                    special_impute=special_impute,
                                                    representation=representation,
                                                    structural_features=structural_features,
                                                    unroll=unroll_single_feat,
                                                    numerical_feats=numerical_feats,
                                                    target_features=target_features,
                                                    regressor_type=regressor_type,
                                                    kernel=kernel,
                                                    transform_type=transform_type,
                                                    second_transformer=second_transformer,
                                                    cutoff=cutoff,
                                                    hyperparameter_optimization=hyperparameter_optimization,
                                                    imputer=imputer,
                                                    Test=TEST,
                                                )
    save_results(scores,
                predictions=predictions,
                imputer=imputer,
                feat_length_scale=ls,
                representation= representation,
                pu_type= oligomer_representation,
                target_features=target_features,
                regressor_type=regressor_type,
                numerical_feats=numerical_feats,
                radius= radius,
                vector =vector,
                kernel=kernel,
                cutoff=cutoff,
                TEST=TEST,
                hypop=hyperparameter_optimization,
                transform_type=transform_type,
                second_transformer=second_transformer,
                # special_folder_name='hp_RF_differences',
                # special_file_name='v2_(max_feat_all_leaf_smaller)',
                )



def parse_arguments():
    parser = ArgumentParser(description="Process some data for numerical-only regression.")
    
    # Argument for regressor_type
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
        "--transform_type", 
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
        '--oligomer_representation', 
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
    #     oligomer_representation = args.oligomer_representation,
    #     regressor_type=args.regressor_type,
    #     kernel=args.kernel,
    #     target_features=[args.target_features],  
    #     transform_type='Standard',
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
        regressor_type="sklearn-GPR",
        kernel="j_rbf",
        oligomer_representation="Monomer",
        target_features=['log Rg (nm)'],  
        transform_type="Robust Scaler",
        second_transformer=None,
        # columns_to_impute=["PDI", "Temperature SANS/SLS/DLS/SEC (K)", "Concentration (mg/ml)"],
        # special_impute="Mw (g/mol)",
        # numerical_feats=['Mw (g/mol)', 'PDI', 'Concentration (mg/ml)', 'Temperature SANS/SLS/DLS/SEC (K)', 'solvent dP', 'solvent dD', 'solvent dH'],
        # imputer='mean',
        hyperparameter_optimization=False,
        # oligomer_representation="Trimer",
    )





