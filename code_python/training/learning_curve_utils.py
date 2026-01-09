
from pathlib import Path
from typing import Callable, Optional, Union, Dict, Tuple

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer, TransformedTargetRegressor
from sklearn.model_selection import KFold
from sklearn.pipeline import Pipeline

import matplotlib.pyplot as plt

from filter_data import filter_dataset
from all_factories import (
                            regressor_factory,
                            transforms)

from training_utils import (
     _optimize_hyperparams,
     split_for_training,
     get_target_transformer
    )
from all_factories import optimized_models
from imputation_normalization import preprocessing_workflow
from scoring import (
    process_learning_score,
    get_incremental_split
    )


def set_globals(Test: bool=False) -> None:
    global SEEDS, N_FOLDS, BO_ITER
    if not Test:
        SEEDS = [6, 13, 42, 69, 420, 1234567890, 473129]
        N_FOLDS = 5
        BO_ITER = 42
    else:
        SEEDS = [42,13]
        N_FOLDS = 2
        BO_ITER = 1


def get_generalizability_predictions(
    dataset: pd.DataFrame,
    features_impute: Optional[list[str]],
    special_impute: Optional[str],
    representation: Optional[str],
    structural_features: Optional[list[str]],
    numerical_feats: Optional[list[str]],
    unroll: Union[dict[str, str], list[dict[str, str]], None],
    regressor_type: str,
    target_features: str,
    transform_type: str,
    second_transformer:str=None,
    hyperparameter_optimization: bool=True,
    imputer: Optional[str] = None,
    cutoff:Dict[str, Tuple[Optional[float], Optional[float]]]=None,
    Test:bool=False,
    # output_dir_name: str = "results",
    ) -> None:
        """
        you should change the name here for prepare
        """
            #seed scores and seed prediction
        set_globals(Test)
        learning_score = _prepare_data(
                                        dataset=dataset,
                                        features_impute= features_impute,
                                        special_impute= special_impute,
                                        representation=representation,
                                        structural_features=structural_features,
                                        unroll=unroll,
                                        numerical_feats = numerical_feats,
                                        target_features=target_features,
                                        regressor_type=regressor_type,
                                        transform_type=transform_type,
                                        second_transformer=second_transformer,
                                        imputer=imputer,
                                        cutoff=cutoff,
                                        hyperparameter_optimization=hyperparameter_optimization,
                                        )
        
        
        scores = process_learning_score(learning_score)
        return scores
        



def _prepare_data(
    dataset: pd.DataFrame,
    target_features: str,
    regressor_type: str,
    features_impute: Optional[list[str]]=None,
    special_impute: Optional[str]=None,
    representation: Optional[str]=None,
    structural_features: Optional[list[str]]=None,
    numerical_feats: Optional[list[str]]=None,
    unroll: Union[dict, list, None] = None,
    transform_type: str = "Standard",
    hyperparameter_optimization: bool = True,
    imputer: Optional[str] = None,
    cutoff: Dict[str, Tuple[Optional[float], Optional[float]]]=None,
    second_transformer:str=None,
    **kwargs,
    ) -> tuple[dict[int, dict[str, float]], pd.DataFrame]:


    """
    here you should change the names
    """



    X, y, unrolled_feats, _ = filter_dataset(
                                        raw_dataset=dataset,
                                        structure_feats=structural_features,
                                        scalar_feats=numerical_feats,
                                        target_feats=target_features,
                                        cutoff=cutoff,
                                        dropna = True,
                                        unroll=unroll,
                                        )

    # Pipline workflow here and preprocessor
    preprocessor: Pipeline = preprocessing_workflow(imputer=imputer,
                                                    feat_to_impute=features_impute,
                                                    representation = representation,
                                                    numerical_feat=numerical_feats,
                                                    structural_feat = unrolled_feats,
                                                    special_column=special_impute,
                                                    scaler=transform_type)
    
    preprocessor.set_output(transform="pandas")
    learning_score = run_leaning(
                                    X,
                                    y,
                                    preprocessor=preprocessor,
                                    regressor_type=regressor_type,
                                    transform_type=transform_type,
                                    hyperparameter_optimization=hyperparameter_optimization,
                                    second_transformer=second_transformer,
                                    **kwargs,
                                    )

    return learning_score

def run_leaning(
    X, y, preprocessor: Union[ColumnTransformer, Pipeline], regressor_type: str,
    transform_type: str, second_transformer:str=None, hyperparameter_optimization: bool = True,
    **kwargs,
    ) -> dict[int, dict[str, float]]:

    seed_learning_curve_scores:dict[int, dict] = {}

    for seed in SEEDS:
        print(f' REGRESSOR {regressor_type}, \tSEED:, {seed}')
        cv_outer = KFold(n_splits=N_FOLDS, shuffle=True, random_state=seed)
        y_transform = get_target_transformer(transform_type,second_transformer)

        if hyperparameter_optimization:
            y_transform_regressor = TransformedTargetRegressor(
                        regressor=regressor_factory[regressor_type],
                        transformer=y_transform,
                )
            
        else:
            model = optimized_models(regressor_type)
            y_transform_regressor = TransformedTargetRegressor(
                        regressor=model,
                        transformer=y_transform,
                )
        new_preprocessor = 'passthrough' if len(preprocessor.steps) == 0 else preprocessor
        regressor :Pipeline= Pipeline(steps=[
                    ("preprocessor", new_preprocessor),
                    ("regressor", y_transform_regressor),
                        ])
        
        regressor.set_output(transform="pandas")
        steps=0.1
        train_ratio:np.ndarray = np.linspace(0.1, 1, int(0.9 / steps))
        if hyperparameter_optimization:
                best_estimator, regressor_params = _optimize_hyperparams(
                    X,
                    y,
                    cv_outer=cv_outer,
                    seed=seed,
                    regressor_type=regressor_type,
                    regressor=regressor,
                )
                train_sizes, train_scores, test_scores = get_incremental_split(best_estimator,
                                                                                X,
                                                                                y,
                                                                                cv_outer,
                                                                                train_ratio=train_ratio,
                                                                                random_state=seed)


        else:
                train_sizes, train_scores, test_scores = get_incremental_split(regressor,
                                                                                X,
                                                                                y,
                                                                                cv_outer,
                                                                                train_ratio=train_ratio,
                                                                                random_state=seed)


        seed_learning_curve_scores[seed] = {
                    "train_sizes": train_sizes,   # 1D array of training sizes used
                    "train_sizes_fraction": train_sizes/len(X),
                    "train_scores": train_scores,  # 2D array of training scores
                    "test_scores": test_scores,  # 2D array of validation (cross-validation) scores
                    "best_params": regressor_params if hyperparameter_optimization else "Default"
                }
    return seed_learning_curve_scores






