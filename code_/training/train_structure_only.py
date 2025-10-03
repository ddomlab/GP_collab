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
        choices=['RF', 'DT', 'MLR', 'SVR', 'XGBR','KNN', 'GPR', 'NGB', 'sklearn-GPR'], 
        required=True, 
        help="Regressor type required"
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
        required=True, 
        help="Fingerprint required"
    )

    parser.add_argument(
        '--oligomer_representation', 
        type=str, 
        choices=['Monomer', 'Dimer', 'Trimer', 'RRU Monomer', 'RRU Dimer', 'RRU Trimer'], 
        required=True, 
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


    return parser.parse_args()

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

