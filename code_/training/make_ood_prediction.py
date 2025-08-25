import pandas as pd
from pathlib import Path
from get_ood_split import train_regressor
from train_structure_numerical import get_structural_info
from all_factories import radius_to_bits,cutoffs
from typing import Callable, Optional, Union, Dict, Tuple
import numpy as np
import sys
import os

# sys.path.append("../cleaning")
# sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../visualization")))
# from visualize_ood_scores import plot_splits_scores, plot_splits_parity
from argparse import ArgumentParser
from data_handling import save_results
from all_factories import generate_acronym_string
from train_structure_numerical import parse_arguments
HERE: Path = Path(__file__).resolve().parent
DATASETS: Path = HERE.parent.parent / "datasets"
RESULTS = Path = HERE.parent.parent / "results"

training_df_dir: Path = DATASETS/ "training_dataset"/"Rg data with clusters.pkl"
w_data = pd.read_pickle(training_df_dir)

# clusters = 'KM4 ECFP6_Count_512bit cluster'	
# 'HBD3 MACCS cluster'  do this as new
# 'KM3 Mordred cluster'	
# 'substructure cluster'
# 'KM5 polymer_solvent HSP and polysize cluster'	redo this
# 'KM4 polymer_solvent HSP and polysize cluster'    do this as new one
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
    hyperparameter_optimization: bool=True,
    radius:int=None,
    vector:str=None,
    kernel:str = None,
    cutoff:Optional[str]=None,
    second_transformer:str=None,
    clustering_method:str=None,
) -> None:
    
    structural_features, unroll_single_feat = get_structural_info(representation,oligomer_representation,radius,vector)
    scores, predictions, cluster_y_truth  = train_regressor(
                                                    dataset=dataset,
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
                                                    clustering_method=clustering_method,
                                                    Test=TEST,
                                                )
    

    _ = save_results(scores,
                predictions=predictions,
                representation= representation,
                ground_truth=cluster_y_truth,
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
                clustering_method=clustering_method,
                # special_folder_name='mean_aggregated',
                )
    #TODO: Plot the results
    # scores_criteria: list= ['mad', 'mae', 'rmse',
    #                         'r2', 'ystd', 'pearson_r', 
    #                         'pearson_p_value', 'spearman_r',
    #                         'spearman_p_value', 'kendall_r', 'kendall_p_value']
    # suffix = f"{regressor_type}_{representation}" if representation else f'{regressor_type}'
    # suffix = f"{suffix}_{transform_type}" if transform_type else suffix
    # feats_abbv = generate_acronym_string(numerical_feats) if numerical_feats else None
    # suffix = f"{suffix}_{feats_abbv}" if feats_abbv else suffix
    # score_plot_folder = saving_folder/ f'comparitive cluster scores ({suffix})'
    # plot_splits_scores(scores=scores, scores_criteria=scores_criteria, folder=score_plot_folder)

    # print("-"*30
    #       ,"\nPlotted Comparitive Cluster Scores!")
    

    # parity_folder = saving_folder/ f'parity plot ({suffix})'
    # plot_splits_parity(predicted_values=predictions,
    #                     ground_truth=cluster_y_truth,
    #                     score=scores,
    #                     folder=parity_folder)
    
    # print("_"*30,
    #       "\nPlotted Parity Plots!")


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
            kernel=args.kernel,
            target_features=[args.target_features],  
            transform_type='Standard',
            second_transformer=None,
            hyperparameter_optimization=True,
            numerical_feats=args.numerical_feats, 
            clustering_method=args.clustering_method,
        )

# # TODO: Update the dataset to include the clustering methods
    else:
        main_structural_numerical(
            dataset=w_data,
            representation="MACCS",
            # radius=3,
            # vector="count",
            regressor_type="MLP",
            target_features=['log Rg (nm)'],  
            transform_type='Standard',
            second_transformer=None,
            numerical_feats=['Xn', 'Mw (g/mol)', 'PDI', 'Concentration (mg/ml)', 'Temperature SANS/SLS/DLS/SEC (K)', 'polymer dP', 'polymer dD' , 'polymer dH', 'solvent dP', 'solvent dD', 'solvent dH'],
            hyperparameter_optimization=False,
            oligomer_representation="Trimer",
            clustering_method='substructure cluster'
        )