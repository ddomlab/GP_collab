import pandas as pd
from pathlib import Path
from training_utils import train_regressor
from all_factories import radius_to_bits,cutoffs
import sys
# import numpy as np
from typing import Callable, Optional, Union, Dict, Tuple
sys.path.append("../cleaning")
# from clean_dataset import open_json
from argparse import ArgumentParser
from data_handling import save_results
from train_structure_numerical import get_structural_info
from utils import parse_arguments

HERE: Path = Path(__file__).resolve().parent
DATASETS: Path = HERE.parent.parent / "datasets"
RESULTS = Path = HERE.parent.parent / "results"

training_df_dir: Path = DATASETS/ "training_dataset"/"Rg data with clusters aging imputed.pkl"
w_data = pd.read_pickle(training_df_dir)


TEST=False


def main_structural(
    dataset: pd.DataFrame,
    regressor_type: str,
    target_features: list[str],
    transform_type: str,
    representation:str,
    oligomer_representation: str,
    hyperparameter_optimization: bool=True,
    radius:int=None,
    vector:str=None,
    kernel:str = None,
    cutoff:Optional[str]=None,
    second_transformer:str=None,
) -> None:
    
    structural_features, unroll_single_feat = get_structural_info(representation,oligomer_representation,radius,vector)
    scores, predictions,data_shapes  = train_regressor(
                                                    dataset=dataset,
                                                    representation=representation,
                                                    structural_features=structural_features,
                                                    unroll=unroll_single_feat,
                                                    features_impute=None,
                                                    special_impute=None,
                                                    numerical_feats=None,
                                                    imputer=None,
                                                    target_features=target_features,
                                                    regressor_type=regressor_type,
                                                    kernel=kernel,
                                                    feat_transformer=transform_type,
                                                    target_transformer=second_transformer,
                                                    cutoff=cutoff,
                                                    hyperparameter_optimization=hyperparameter_optimization,
                                                    Test=TEST,
                                                )
    save_results(scores,
                predictions=predictions,
                df_shapes=data_shapes,
                representation= representation,
                pu_type= oligomer_representation,
                target_features=target_features,
                regressor_type=regressor_type,
                radius= radius,
                vector =vector,
                kernel=kernel,
                cutoff=cutoff,
                TEST=TEST,
                hypop=hyperparameter_optimization,
                transform_type=transform_type,
                target_transformer=second_transformer,
                )





if __name__ == "__main__":
    args = parse_arguments()

    main_structural(
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
                )











# if __name__ == '__main__':
    # main()
    # main_structural(
    #     dataset=w_data,
    #     representation="MACCS",
    #     # radius=3,
    #     # vector="count",
    #     regressor_type="RF",
    #     target_features=['log First Peak (e-5 place holder)'],  
    #     transform_type=None,
    #     second_transformer=None,
    #     hyperparameter_optimization=True,
    #     oligomer_representation="Monomer",
    # )

