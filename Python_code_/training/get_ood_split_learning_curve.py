from pathlib import Path
from typing import Callable, Optional, Union, Dict, Tuple
from sklearn.pipeline import Pipeline
from sklearn.compose import ColumnTransformer, TransformedTargetRegressor
from joblib import Parallel, delayed
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold, StratifiedKFold, train_test_split
from skopt import BayesSearchCV
from data_handling import remove_unserializable_keys, save_results
from filter_data import filter_dataset
from all_factories import (
                            transforms,
                            get_regressor_search_space)
from collections import Counter

from imputation_normalization import preprocessing_workflow
from training_utils import get_target_transformer, split_for_training,set_globals
from scoring import (
                    train_and_predict_ood,
                    process_ood_learning_curve_score,
                    get_incremental_split
                )
from all_factories import optimized_models

from get_ood_split import (StratifiedKFoldWithLabels,
                            get_loco_splits)

import warnings
from sklearn.exceptions import DataConversionWarning
warnings.filterwarnings(action='ignore', category=DataConversionWarning)

def set_globals(Test: bool=False) -> None:
    global N_FOLDS, BO_ITER
    if not Test:
        # SEEDS = [6, 13, 42, 69, 420, 1234567890, 473129]
        N_FOLDS = 5
        BO_ITER = 42
        MAINRASIOS = []
    else:
        # SEEDS = [42,13]
        N_FOLDS = 2
        BO_ITER = 1



def train_ood_learning_curve(
                            dataset: pd.DataFrame,
                            target_features: str,
                            representation: Optional[str],
                            structural_features: Optional[list[str]],
                            numerical_feats: Optional[list[str]],
                            unroll: Union[dict[str, str], list[dict[str, str]], None],
                            regressor_type: str,
                            transform_type: str=None,
                            second_transformer:str=None,
                            clustering_method:str= None,
                            Test:bool=False,
                            ):
                          

    set_globals(Test)
    score, predictions = prepare_data(
                                    dataset=dataset,
                                    target_features=target_features,
                                    regressor_type=regressor_type,
                                    representation=representation,
                                    structural_features=structural_features,
                                    numerical_feats=numerical_feats,
                                    unroll=unroll,
                                    transform_type=transform_type,
                                    second_transformer=second_transformer,
                                    clustering_method=clustering_method,
                                    )
    score = process_ood_learning_curve_score(score)
    return score, predictions

def prepare_data(    
                dataset: pd.DataFrame,
                target_features: str,
                regressor_type: str,
                representation: Optional[str]=None,
                structural_features: Optional[list[str]]=None,
                numerical_feats: Optional[list[str]]=None,
                unroll: Union[dict, list, None] = None,
                transform_type: str = None,
                second_transformer:str = None,
                clustering_method:str= None,
                ):

    X, y, unrolled_feats, cluster_labels, _ = filter_dataset(
                                                    raw_dataset=dataset,
                                                    structure_feats=structural_features,
                                                    scalar_feats=numerical_feats,
                                                    target_feats=target_features,
                                                    dropna = True,
                                                    unroll=unroll,
                                                    cluster_type=clustering_method,
                                                    )

    preprocessor: Pipeline = preprocessing_workflow(imputer=None,
                                                    feat_to_impute=None,
                                                    representation = representation,
                                                    numerical_feat=numerical_feats,
                                                    structural_feat = unrolled_feats,
                                                    special_column=None,
                                                    scaler=transform_type
                                                    )
    preprocessor.set_output(transform="pandas")

    return run_ood_learning_curve(X, y,
                                    cluster_labels,
                                    regressor_type,
                                    transform_type,
                                    second_transformer,
                                    preprocessor
                                    )


def run_single_ood(
    seed: int,
    train_ratio: float,
    X_train_val,
    y_train_val,
    X_test,
    y_test,
    cluster_train_val_labels_OOD: np.ndarray,
    model_name: str,
    transform_type: str,
    second_transformer: str,
    preprocessor: Optional[Pipeline],
) -> Tuple[int, Tuple]:
    if train_ratio == 1:
        X_train_OOD, y_train_OOD = X_train_val, y_train_val
    else:
        try:
            X_train_OOD, _, y_train_OOD, _ = train_test_split(
                X_train_val, y_train_val, train_size=train_ratio,
                random_state=seed, stratify=cluster_train_val_labels_OOD, shuffle=True,
            )
        
        except ValueError:
             X_train_OOD, _, y_train_OOD, _ = train_test_split(
                X_train_val, y_train_val, train_size=train_ratio,
                random_state=seed, shuffle=True,
            )
    return seed, fit_and_eval(
        X_train_OOD, y_train_OOD, X_test, y_test,
        model_name, transform_type, second_transformer, preprocessor
    )


def run_single_iid(
    seed: int,
    test_seed: int,
    train_ratio: float,
    X_full,
    y_full,
    y_test_len: int,
    model_name: str,
    transform_type: str,
    second_transformer: str,
    preprocessor: Optional[Pipeline],
) -> Tuple[int, int, Tuple]:
    X_train_val_IID, X_test_IID, y_train_val_IID, y_test_IID = train_test_split(
        X_full, y_full, test_size=y_test_len,
        random_state=test_seed, shuffle=True,
    )
    if train_ratio == 1:
        X_train_IID, y_train_IID = X_train_val_IID, y_train_val_IID
    else:
        X_train_IID, _, y_train_IID, _ = train_test_split(
            X_train_val_IID, y_train_val_IID,
            train_size=train_ratio, random_state=seed, shuffle=True
        )
    return seed, test_seed, fit_and_eval(
        X_train_IID, y_train_IID, X_test_IID, y_test_IID,
        model_name, transform_type, second_transformer, preprocessor
    )


def fit_and_eval(
    X_train, y_train, X_test, y_test,
    model_name, transform_type, second_transformer, preprocessor
):
    regressor = build_regressor(
        algorithm=model_name,
        first_tran=transform_type,
        second_tran=second_transformer,
        preprocessor=preprocessor,
    )
    unc_preprocessor = preprocessor
    test_scores, train_scores, y_pred, y_unc = train_and_predict_ood(
        regressor, X_train, y_train, X_test, y_test,
        return_train_pred=True, algorithm=model_name,
        manual_preprocessor=unc_preprocessor,
    )
    return train_scores, test_scores, y_pred.flatten(), (
        None if y_unc is None else y_unc.flatten()
    )

def run_ood_learning_curve(
    X,
    y,
    cluster_labels: Union[np.ndarray, dict[str, np.ndarray]],
    model_name: str,
    transform_type: str = None,
    second_transformer: str = None,
    preprocessor: Optional[Pipeline] = None,
    n_jobs: int = -1,
) -> Tuple[dict, dict]:

    loco_split_idx = get_loco_splits(cluster_labels)
    learning_curve_scores, learning_curve_predictions = {}, {}

    train_sizes = [len(tv_idxs) for tv_idxs, _ in loco_split_idx.values()]
    min_train_size = min(train_sizes)
    for cluster, (tv_idx, test_idx) in loco_split_idx.items():
        if cluster == 'Polar':
            cluster_tv_labels_OOD = split_for_training(cluster_labels['Side Chain Cluster'], tv_idx)
        elif cluster in ['Fluorene', 'PPV', 'Thiophene']:
            cluster_tv_labels_OOD = split_for_training(cluster_labels['substructure cluster'], tv_idx)
        else:
            # cluster_tv_labels_OOD = split_for_training(cluster_labels, tv_idx)
            raw_labels = split_for_training(cluster_labels, tv_idx)
            cluster_tv_labels_OOD = merge_small_clusters(raw_labels, min_count=2)
            # name, counts = np.unique(cluster_tv_labels_OOD, return_counts=True)
        print(cluster)
            

        X_tv_OOD = split_for_training(X, tv_idx)
        y_tv_OOD = split_for_training(y, tv_idx)
        X_test_OOD = split_for_training(X, test_idx)
        y_test_OOD = split_for_training(y, test_idx)

        learning_curve_predictions.setdefault(f'CO_{cluster}', {})['y_true'] = y_test_OOD.flatten()
        train_ratios = [
            .1, 0.3, 0.5, 0.7,
                         .9,1]
        min_ratio = min_train_size / len(X_tv_OOD)
        if min_ratio not in train_ratios:
            train_ratios.append(min_ratio)

        for train_ratio in train_ratios:
            seeds = random_state_generator(train_ratio)
            print("-------- OOD Parallel --------")
            ood_results = Parallel(n_jobs=n_jobs)(
                delayed(run_single_ood)(
                    seed, train_ratio, X_tv_OOD, y_tv_OOD,
                    X_test_OOD, y_test_OOD, cluster_tv_labels_OOD,
                    model_name, transform_type, second_transformer, preprocessor
                ) for seed in seeds
            )

            for seed, (train_scores, test_scores, y_pred, y_unc) in ood_results:
                learning_curve_scores.setdefault(f'CO_{cluster}', {}).setdefault(f'ratio_{train_ratio}', {})[f'seed_{seed}'] = (train_scores, test_scores)
                learning_curve_predictions.setdefault(f'CO_{cluster}', {}).setdefault(f'ratio_{train_ratio}', {})[f'seed_{seed}'] = {
                    'y_test_pred': y_pred,
                    'y_test_uncertainty': y_unc,
                }
            
            print("-------- IID Parallel --------")
            y_test_len = len(y_test_OOD)
            test_ratio = y_test_len / len(y)
            test_seeds = random_state_generator(test_ratio)
            iid_results = Parallel(n_jobs=n_jobs)(
                delayed(run_single_iid)(
                    seed, test_seed, train_ratio, X, y,
                    y_test_len, model_name, transform_type,
                    second_transformer, preprocessor
                ) for seed in seeds for test_seed in test_seeds
            )

            for seed, test_seed, (train_scores, test_scores, y_pred, y_unc) in iid_results:
                scores = learning_curve_scores.setdefault(f'IID_{cluster}', {}).setdefault(f'ratio_{train_ratio}', {}).setdefault(f'seed_{seed}', {})
                preds = learning_curve_predictions.setdefault(f'IID_{cluster}', {}).setdefault(f'ratio_{train_ratio}', {}).setdefault(f'seed_{seed}', {})
                scores[f'test_set_seed_{test_seed}'] = (train_scores, test_scores)
                preds[f'test_set_seed_{test_seed}'] = {
                    'y_test_pred': y_pred,
                    'y_test_uncertainty': y_unc,
                }

        for prefix in ['CO', 'IID']:
            learning_curve_scores[f'{prefix}_{cluster}']['training size'] = len(X_tv_OOD)
            learning_curve_predictions[f'{prefix}_{cluster}']['training size'] = len(X_tv_OOD)

    return learning_curve_scores, learning_curve_predictions



def build_regressor(algorithm,first_tran,second_tran,preprocessor):
    y_transform = get_target_transformer(first_tran, second_tran)
    model = optimized_models(algorithm)

    y_model = TransformedTargetRegressor(
                        regressor=model,
                        transformer=y_transform,
                        )
    new_preprocessor = 'passthrough' if len(preprocessor.steps) == 0 else preprocessor
    regressor :Pipeline= Pipeline(steps=[
                        ("preprocessor", new_preprocessor),
                        ("regressor", y_model),
                        ])
    regressor.set_output(transform="pandas")
    return regressor


def random_state_generator(train_ratio: float) -> np.ndarray:
    if train_ratio >0.7:
        random_state_list = np.arange(5)
    elif train_ratio >0.5:
        random_state_list = np.arange(7)
    elif train_ratio >0.3:
        random_state_list = np.arange(12)
    elif train_ratio >=0.1:
        random_state_list = np.arange(30)
    else:
        random_state_list = np.arange(40)

    return random_state_list



def merge_small_clusters(labels, min_count=2):
    label_counts = Counter(labels)
    labels = np.array(labels)
    merged_labels = np.array([
        label if label_counts[label] >= min_count else 'Other'
        for label in labels
    ])
    return merged_labels



