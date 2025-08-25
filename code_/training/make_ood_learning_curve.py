import pandas as pd
from pathlib import Path
from get_ood_split_learning_curve import train_ood_learning_curve
from train_structure_numerical import get_structural_info
from typing import Callable, Optional, Union, Dict, Tuple
import numpy as np
import sys
import os

from argparse import ArgumentParser
from data_handling import save_results
from train_structure_numerical import parse_arguments
HERE: Path = Path(__file__).resolve().parent
DATASETS: Path = HERE.parent.parent / "datasets"
RESULTS = Path = HERE.parent.parent / "results"

training_df_dir: Path = DATASETS/ "training_dataset"/"Rg data with clusters.pkl"
w_data = pd.read_pickle(training_df_dir)

# clusters = 'KM4 ECFP6_Count_512bit cluster'	
# 'KM3 Mordred cluster'	
# 'substructure cluster'
# 'EG-Ionic-Based Cluster'
# 'KM5 polymer_solvent HSP and polysize cluster'	
# 'KM4 polymer_solvent HSP cluster'
# 'KM4 Mordred_Polysize cluster'


TEST = False


def main_structural_numerical(
    dataset: pd.DataFrame,
    regressor_type: str,
    target_features: list[str],
    transform_type: str,
    numerical_feats: Optional[list[str]],
    representation:str=None,
    oligomer_representation: str=None,
    radius:int=None,
    vector:str=None,
    second_transformer:str=None,
    clustering_method:str=None,
) -> None:
    
    structural_features, unroll_single_feat = get_structural_info(representation,oligomer_representation,radius,vector)
    scores, predictions  = train_ood_learning_curve(
                                                    dataset=dataset,
                                                    target_features=target_features,
                                                    representation=representation,
                                                    structural_features=structural_features,
                                                    unroll=unroll_single_feat,
                                                    numerical_feats=numerical_feats,
                                                    regressor_type=regressor_type,
                                                    transform_type=transform_type,
                                                    second_transformer=second_transformer,
                                                    clustering_method=clustering_method,
                                                    Test=TEST,
                                                )

    _ = save_results(
                    scores,
                    predictions=predictions,
                    representation= representation,
                    pu_type= oligomer_representation,
                    target_features=target_features,
                    regressor_type=regressor_type,
                    numerical_feats=numerical_feats,
                    radius= radius,
                    vector =vector,
                    TEST=TEST,
                    hypop=False,
                    transform_type=transform_type,
                    target_transformer=second_transformer,
                    clustering_method=clustering_method,
                    learning_curve=True,
                    # special_folder_name='mean_aggregated',
                    # special_file_name='v1_(max_feat_sqrt)'
                    )
    


if __name__ == "__main__":
    if TEST==False:
        args = parse_arguments()

        main_structural_numerical(
            dataset=w_data,
            representation=args.representation,
            radius=args.radius,
            vector=args.vector,
            oligomer_representation = args.oligomer_representation,
            regressor_type=args.regressor_type,
            target_features=[args.target_features],  
            transform_type='Standard',
            second_transformer=None,
            numerical_feats=args.numerical_feats, 
            clustering_method=args.clustering_method,
        )


    else:
            main_structural_numerical(
            dataset=w_data,
            representation="Mordred",
            # radius=3,
            # vector="count",
            regressor_type="RF",
            target_features=['log Rg (nm)'],  
            transform_type='Standard',
            second_transformer=None,
            numerical_feats=['Xn', 'Mw (g/mol)', 'PDI', 'Concentration (mg/ml)', 'Temperature SANS/SLS/DLS/SEC (K)', 'polymer dP', 'polymer dD' , 'polymer dH', 'solvent dP', 'solvent dD', 'solvent dH'],
            oligomer_representation="Trimer",
            clustering_method='canonical_name'
            )
        
            # main_structural_numerical(
            # dataset=w_data,
            # representation=None,
            # # radius=3,
            # # vector="count",
            # regressor_type="RF",
            # target_features=['log Rg (nm)'],  
            # transform_type='Standard',
            # second_transformer=None,
            # numerical_feats=['Mw (g/mol)', 'PDI', 'Concentration (mg/ml)', 'Temperature SANS/SLS/DLS/SEC (K)', 'polymer dP', 'polymer dD' , 'polymer dH', 'solvent dP', 'solvent dD', 'solvent dH'],
            # oligomer_representation=None,
            # clustering_method='substructure cluster'
            # )