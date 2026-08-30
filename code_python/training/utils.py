import glob
import numpy as np
import pandas as pd
import random
from pyparsing import Union
import scipy.stats
from argparse import ArgumentParser
# import torch

from sklearn.decomposition import PCA 
from sklearn.model_selection import KFold, train_test_split
from sklearn.metrics import mean_squared_error as mse
from sklearn.metrics import r2_score
import torch



def read_split_files(dataset, data_dir = 'data'):
    fnames = glob.glob(f'{data_dir}/{dataset}_splits*.npz')
    splits = []
    for f in fnames:
        splits.append(np.load(f))
        
    return splits


def parse_arguments():
    parser = ArgumentParser(description="Process some data for numerical-only regression.")
    
    parser.add_argument(
        '--target_features',
        # choices=['Lp (nm)', 'Rg1 (nm)', 'Rh (IW avg log)'],  
        required=False,
        help="Specify a single target for the analysis."
    )

    parser.add_argument(
        '--cuda',
        # choices=['Lp (nm)', 'Rg1 (nm)', 'Rh (IW avg log)'],  
        required=False,
        type=bool,
        help="Specify whether to use CUDA (GPU) or not (CPU)."
    )

    parser.add_argument(
        '--paper',
        # choices=['Lp (nm)', 'Rg1 (nm)', 'Rh (IW avg log)'],  
        required=False,
        help="Specify a single target for the analysis."
    )

    parser.add_argument(
        '--dataset',
        # choices=['Lp (nm)', 'Rg1 (nm)', 'Rh (IW avg log)'],  
        required=False,
        help="Specify dataset name"
    )

    parser.add_argument(
        '--K_fp',
        # choices=['Lp (nm)', 'Rg1 (nm)', 'Rh (IW avg log)'],  
        required=False,
        help="Specify fingerprint kernel"
    )

    parser.add_argument(
        '--K_count',
        # choices=['Lp (nm)', 'Rg1 (nm)', 'Rh (IW avg log)'],  
        required=False,
        help="Specify count kernel"
    )

    parser.add_argument(
        '--Kernel_mixing_method',
        # choices=['Lp (nm)', 'Rg1 (nm)', 'Rh (IW avg log)'],  
        required=False,
        help="Specify a hybridization method for the analysis."
    )

    parser.add_argument(
        '--kernel_feature_mode',
        required=False,
        help="Specify a kernel feature mode (joint(ARD)/per_feature) for the analysis."
    )

    parser.add_argument(
        '--regressor_type', 
        type=str, 
        required=False, 
        help="Regressor type required"
    )

    parser.add_argument(
        '--numerical_feats',
        type=str,
        choices=['Xn','Mn (g/mol)', 'Mw (g/mol)', 'PDI', 'Temperature SANS/SLS/DLS/SEC (K)',
                  'Concentration (mg/ml)','solvent dP',	'polymer dP',	'solvent dD',	'polymer dD',	'solvent dH',	'polymer dH', 'Ra',
                  "abs(solvent dD - polymer dD)", "abs(solvent dP - polymer dP)", "abs(solvent dH - polymer dH)"],

        nargs='+',  # Allows multiple choices
        required=False,
        help="Numerical features: choose"
    )
    
    parser.add_argument(
        '--columns_to_impute',
        type=str,
        choices=['Xn','Mn (g/mol)', 'Mw (g/mol)', 'PDI', 'Temperature SANS/SLS/DLS/SEC (K)',
                  'Concentration (mg/ml)','solvent dP',	'polymer dP',	'solvent dD',	'polymer dD',	'solvent dH',	'polymer dH', 'Ra'],

        nargs='*',  # This allows 0 or more values
        default=False,  
        help="imputation features: choose"
    )

    parser.add_argument(
        '--imputer',
        choices=['mean', 'median', 'most_frequent',"distance KNN", None],  
        nargs='?',  # This allows the argument to be optional
        default=False,  
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



def split_for_training(
    data: Union[pd.DataFrame, np.ndarray, pd.Series, torch.Tensor],
    indices: np.ndarray
) -> Union[pd.DataFrame, np.ndarray, pd.Series, torch.Tensor]:

    # ensure numpy index array
    indices = np.asarray(indices)
    # pandas DataFrame
    if isinstance(data, pd.DataFrame):
        return data.iloc[indices].copy()
    # pandas Series
    elif isinstance(data, pd.Series):
        return data.iloc[indices].copy()

    # numpy array
    elif isinstance(data, np.ndarray):
        return data[indices]
    # torch tensor
    elif isinstance(data, torch.Tensor):
        return data[indices]
    else:
        raise ValueError(
            f"Unsupported data type {type(data)}. "
            "Must be DataFrame, Series, numpy array, or torch Tensor."
        )