import pandas as pd
from pathlib import Path
from typing import Callable, Optional, Union, Dict, Tuple
from learning_curve_utils import get_generalizability_predictions
from all_factories import radius_to_bits,cutoffs
import sys
import json
import numpy as np
sys.path.append("../cleaning")
from argparse import ArgumentParser
from data_handling import save_results
from train_structure_numerical import get_structural_info,parse_arguments

HERE: Path = Path(__file__).resolve().parent
DATASETS: Path = HERE.parent.parent / "datasets"
RESULTS = Path = HERE.parent.parent / "results"

training_df_dir: Path = DATASETS/ "training_dataset"/ "Rg data with clusters aging imputed.pkl"
w_data = pd.read_pickle(training_df_dir)
# ['Monomer', 'Dimer', 'Trimer', 'RRU Monomer', 'RRU Dimer', 'RRU Trimer']
TEST=False




def main_structural_numerical(
    dataset: pd.DataFrame,
    regressor_type: str,
    target_features: list[str],
    transform_type: str,
    representation:str=None,
    oligomer_representation: str=None,
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

    structural_features, unroll_single_feat= get_structural_info(representation,oligomer_representation,radius,vector)

    learning_score  = get_generalizability_predictions(
                                                    dataset=dataset,
                                                    features_impute=columns_to_impute,
                                                    special_impute=special_impute,
                                                    representation=representation,
                                                    structural_features=structural_features,
                                                    unroll=unroll_single_feat,
                                                    numerical_feats=numerical_feats,
                                                    target_features=target_features,
                                                    regressor_type=regressor_type,
                                                    transform_type=transform_type,
                                                    cutoff=cutoffs,
                                                    hyperparameter_optimization=hyperparameter_optimization,
                                                    imputer=imputer,
                                                    second_transformer=second_transformer,
                                                    Test=TEST,
                                                )
    _=save_results(
                scores=learning_score,
                imputer=imputer,
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
                target_transformer=second_transformer,
                learning_curve=True
                )



if __name__ == '__main__':
    if TEST==False:
        args = parse_arguments()

        main_structural_numerical(
        dataset=w_data,
        target_features=[args.target_features],  
        representation=args.representation,
        radius=args.radius,
        vector=args.vector,
        oligomer_representation = args.oligomer_representation,
        regressor_type=args.regressor_type,
        transform_type='Standard',
        second_transformer=None,
        hyperparameter_optimization=False,
        columns_to_impute=args.columns_to_impute,  
        special_impute=args.special_impute,
        numerical_feats=args.numerical_feats,  
        imputer=args.imputer,
        cutoff=None,  
        )

    else:
        main_structural_numerical(
            dataset=w_data,
            # representation="Mordred",
            # radius=3,
            # vector="count",
            regressor_type="NGB",
            target_features=['log Rg (nm)'],  
            transform_type='Standard',
            second_transformer=None,
            # columns_to_impute=["PDI", "Temperature SANS/SLS/DLS/SEC (K)", "Concentration (mg/ml)"],
            # special_impute="Mw (g/mol)",
            numerical_feats=['Mw (g/mol)', 'PDI', 'Concentration (mg/ml)', 'Temperature SANS/SLS/DLS/SEC (K)', 'solvent dP', 'solvent dD', 'solvent dH'],
            # imputer='mean',
            hyperparameter_optimization=False,
            # oligomer_representation="Trimer",
        )






