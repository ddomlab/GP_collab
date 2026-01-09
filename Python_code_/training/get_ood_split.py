from pathlib import Path
from typing import Callable, Optional, Union, Dict, Tuple
from sklearn.pipeline import Pipeline
from sklearn.compose import ColumnTransformer, TransformedTargetRegressor
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold, StratifiedKFold
from skopt import BayesSearchCV
from data_handling import remove_unserializable_keys, save_results
from filter_data import filter_dataset
from all_factories import (
                            regressor_factory,
                            optimized_models,
                            # regressor_search_space,
                            transforms,
                            construct_kernel,
                            get_regressor_search_space)

from imputation_normalization import preprocessing_workflow
from training_utils import get_target_transformer, split_for_training,_optimize_hyperparams,set_globals
from scoring import (
                    cross_validate_regressor,
                    process_scores,
                    train_and_predict_ood,
                    process_ood_scores
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


class StratifiedKFoldWithLabels(StratifiedKFold):
    def __init__(self, labels, n_splits=5, shuffle=True, random_state=None):
        super().__init__(n_splits=n_splits, shuffle=shuffle, random_state=random_state)
        self.labels = labels

    def split(self, X, y=None, groups=None):
        # Use the labels instead of the target variable (y)
        return super().split(X, self.labels)




def train_regressor(
    dataset: pd.DataFrame,
    representation: Optional[str],
    structural_features: Optional[list[str]],
    numerical_feats: Optional[list[str]],
    unroll: Union[dict[str, str], list[dict[str, str]], None],
    regressor_type: str,
    target_features: str,
    features_impute: Optional[list[str]]=None,
    special_impute: Optional[str]=None,
    transform_type: str=None,
    second_transformer:str=None,
    hyperparameter_optimization: bool=True,
    imputer: Optional[str] = None,
    cutoff:Dict[str, Tuple[Optional[float], Optional[float]]]=None,
    kernel: Optional[str] = None,
    Test:bool=False,
    clustering_method:str= None,
    classification:bool=False,
    ) -> None:
        """
        you should change the name here for prepare
        """
            #seed scores and seed prediction
        set_globals(Test)
        scores, predictions,cluster_y_ture = _prepare_data(
                                            dataset=dataset,
                                            representation=representation,
                                            structural_features=structural_features,
                                            unroll=unroll,
                                            numerical_feats = numerical_feats,
                                            target_features=target_features,
                                            regressor_type=regressor_type,
                                            transform_type=transform_type,
                                            second_transformer=second_transformer,
                                            cutoff=cutoff,
                                            hyperparameter_optimization=hyperparameter_optimization,
                                            kernel=kernel,
                                            clustering_method=clustering_method,
                                            )
        scores = process_ood_scores(scores)
        return scores, predictions, cluster_y_ture



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
    transform_type: str = None,
    second_transformer:str = None,
    hyperparameter_optimization: bool = True,
    imputer: Optional[str] = None,
    cutoff: Dict[str, Tuple[Optional[float], Optional[float]]]=None,
    clustering_method:str= None,
    kernel:Optional[str] = None,
    **kwargs,
    ) -> tuple[dict[int, dict[str, float]], pd.DataFrame]:


    """
    here you should change the names
    """


    X, y, unrolled_feats, cluster_labels, X_y_shape = filter_dataset(
                                                    raw_dataset=dataset,
                                                    structure_feats=structural_features,
                                                    scalar_feats=numerical_feats,
                                                    target_feats=target_features,
                                                    cutoff=cutoff,
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
    # preprocessor = transformer?
    score, predication, cluster_y_ture  = run_loco_cv(
                                            X,
                                            y,
                                            preprocessor=preprocessor,
                                            second_transformer=second_transformer,
                                            regressor_type=regressor_type,
                                            transform_type=transform_type,
                                            hyperparameter_optimization=hyperparameter_optimization,
                                            kernel=kernel,
                                            cluster_labels=cluster_labels,
                                            **kwargs,
                                            )
    score['overall data shape'] = X_y_shape
    return score, predication, cluster_y_ture




def run_loco_cv(X, y, 
                preprocessor: Union[ColumnTransformer, Pipeline], 
                transform_type: Optional[str],
                second_transformer:Optional[str],
                regressor_type: str,
                cluster_labels: Union[np.ndarray, dict[str, np.ndarray]],
                hyperparameter_optimization: bool = True,
                kernel:Optional[str] = None,
                **kwargs,
                ) -> tuple[dict[int, dict[str, float]], pd.DataFrame]:

    cluster_scores: dict[int, dict[int, dict[str, float]]] = {}
    cluster_predictions: dict[int, dict[int, np.ndarray]] = {}
    cluster_y_test: dict[int, np.ndarray] = {}

    search_space = get_regressor_search_space(regressor_type,kernel)
    kernel = construct_kernel(regressor_type, kernel)

    cluster_y_test[f'ID_y_true'] = y.flatten()
    loco_split_idx:Dict[int,tuple[np.ndarray]] = get_loco_splits(cluster_labels)
    for cluster, (tv_idx,test_idx) in loco_split_idx.items():
        # print(cluster)
        if cluster=='Polar':
            cluster_tv_labels = split_for_training(cluster_labels['Side Chain Cluster'],tv_idx)
        elif cluster in ['Fluorene', 'PPV', 'Thiophene']:
            cluster_tv_labels = split_for_training(cluster_labels['substructure cluster'],tv_idx)
        else:
            cluster_tv_labels = split_for_training(cluster_labels,tv_idx)
        X_tv, y_tv = split_for_training(X, tv_idx), split_for_training(y,tv_idx)
        X_test, y_test = split_for_training(X, test_idx), split_for_training(y, test_idx)

        #OOD section
        cluster_scores[f'CO_{cluster}'] = {}
        cluster_predictions[f'CO_{cluster}'] = {}
        cluster_y_test[f'CO_{cluster}']=y_test.flatten()
        cluster_scores[f'CO_{cluster}']['cluster size (%)'] = round(len(y_test)/len(y)*100)

        ##IID section
        cluster_scores[f'ID_{cluster}'] = {}
        cluster_predictions[f'ID_{cluster}'] = {}
        cluster_scores[f'ID_{cluster}']['cluster size (%)'] = round((1/round(len(y)/len(y_test)))*100)




        print("\n\n")
        print("-"*50,
            "\nOOD TEST ON", cluster
        )
        
        for seed in SEEDS:
            y_transform = get_target_transformer(transform_type,second_transformer)
            skop_scoring = "neg_root_mean_squared_error"

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

            ID_n_split = round(len(y)/len(y_test))
            IID_cv_baseline = KFold(n_splits=ID_n_split, shuffle=True, random_state=seed)
            uncertainty_preprocessr: Pipeline = Pipeline(steps=[
                    ("preprocessor", new_preprocessor),
                ]) 
            if hyperparameter_optimization:
                OOD_best_estimator, OOD_regressor_params = optimize_ood_hp(
                    X_tv,
                    y_tv,
                    seed=seed,
                    regressor_type=regressor_type,
                    search_space=search_space,
                    regressor=regressor,
                    scoring=skop_scoring,
                    cluster_lables=cluster_tv_labels,
                )


                
                OOD_scores, OOD_predictions, uncertenty_predictions = train_and_predict_ood(OOD_best_estimator, X_tv, y_tv, X_test,
                                                                        y_test, return_train_pred=False, algorithm=regressor_type, manual_preprocessor=uncertainty_preprocessr) 
                OOD_scores["best_params"] = OOD_regressor_params

                ID_cv_outer  =  KFold(n_splits=N_FOLDS, shuffle=True, random_state=seed)
                ID_cv_in = KFold(n_splits=N_FOLDS, shuffle=True, random_state=seed)
                best_estimator, regressor_params = _optimize_hyperparams(
                    X,
                    y,
                    cv_outer=ID_cv_outer,
                    cv_in=ID_cv_in,
                    seed=seed,
                    n_iter=BO_ITER,
                    regressor_type=regressor_type,
                    search_space=search_space,
                    regressor=regressor,
                    scoring=skop_scoring,
                    classification=False
                )

                ID_scores, ID_predictions = cross_validate_regressor(
                    best_estimator, X, y, IID_cv_baseline, classification=False
                )
                ID_scores["best_params"] = regressor_params

            else:
                OOD_scores, OOD_predictions, uncertenty_predictions = train_and_predict_ood(regressor, X_tv, y_tv, X_test,
                                                                        y_test,return_train_pred=False,
                                                                        algorithm=regressor_type,
                                                                        manual_preprocessor=uncertainty_preprocessr)
                
                ID_scores, ID_predictions = cross_validate_regressor(
                                            regressor, X, y, IID_cv_baseline, classification=False
                                            )


            cluster_scores[f'CO_{cluster}'][seed] = OOD_scores
            cluster_predictions[f'CO_{cluster}'][seed] = {
                'y_test_prediction':OOD_predictions.flatten(),
                'y_test_uncertainty': uncertenty_predictions.flatten() if uncertenty_predictions is not None else None,
                }
            # print(cluster_predictions)                                  
            cluster_scores[f'ID_{cluster}'][seed] = ID_scores
            cluster_predictions[f'ID_{cluster}'][seed] = ID_predictions.flatten()

    return cluster_scores, cluster_predictions, cluster_y_test


def _get_splits(cluster_labels: np.ndarray) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    """Helper function to generate LOO splits from a cluster label array."""
    splits = {}
    cluster_names, counts = np.unique(cluster_labels, return_counts=True)
    
    if len(cluster_names) > 2:
        for cluster, count in zip(cluster_names, counts):
            if count <3:
                continue

            mask = cluster_labels == cluster
            test_idxs = np.where(mask)[0]
            tv_idxs = np.where(~mask)[0]
            splits[cluster] = (tv_idxs, test_idxs)
    else:
        rare_cluster = cluster_names[np.argmin(counts)]  # Identify the rarest cluster
        mask = cluster_labels == rare_cluster
        test_idxs = np.where(mask)[0]
        tv_idxs = np.where(~mask)[0]
        splits[rare_cluster] = (tv_idxs, test_idxs)
    
    return splits


def get_loco_splits(cluster_type: Union[np.ndarray, dict[str, np.ndarray]]) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    """
    Generates Leave-One-Cluster-Out (LOCO) splits.
    
    Args:
        cluster_type (Union[np.ndarray, dict[str, np.ndarray]]): 
            - If `np.ndarray`, it is treated as a cluster label array.
            - If `dict`, it should contain cluster labels under keys 'substructure cluster' and 'side chain'.
    
    Returns:
        dict[int, tuple[np.ndarray, np.ndarray]]: 
            A dictionary where keys are cluster labels, and values are tuples of (train/val indices, test indices).
    """
    if isinstance(cluster_type, np.ndarray):
        return _get_splits(cluster_type)

    if isinstance(cluster_type, dict):
        if "substructure cluster" not in cluster_type or "Side Chain Cluster" not in cluster_type:
            raise ValueError("Expected dictionary keys: 'substructure cluster' and 'Side Chain Cluster'.")

        splits = _get_splits(cluster_type["substructure cluster"])

        # Find rarest side chain cluster and create a split
        side_chain_labels = cluster_type["Side Chain Cluster"]
        unique_labels, counts = np.unique(side_chain_labels, return_counts=True)
        
        if len(unique_labels) > 0:
            rare_cluster = unique_labels[np.argmin(counts)]
            mask = side_chain_labels == rare_cluster
            test_idxs = np.where(mask)[0]
            tv_idxs = np.where(~mask)[0]
            splits[rare_cluster] = (tv_idxs, test_idxs)

        return splits

    raise TypeError("Input must be either an np.ndarray or a dictionary with cluster labels.")


def optimize_ood_hp(
                    x_train_val,
                    y_train_val,
                    regressor,
                    regressor_type,
                    search_space,
                    scoring,
                    cluster_lables:np.ndarray,
                    seed)-> tuple[Pipeline, dict]:
    
    estimators: list[BayesSearchCV] = []
    print("\n\n")
    print("-"*50,
        "\nOPTIMIZING HYPERPARAMETERS FOR REGRESSOR", regressor_type, "\tSEED:", seed
    )
    n_cluster = len(np.unique(cluster_lables))
    if n_cluster>1:
        # print('yes it is stratified')
        cv = StratifiedKFoldWithLabels(n_splits=N_FOLDS, labels=cluster_lables,shuffle=True, random_state=seed)
        for train_index, val_index in cv.split(x_train_val, cluster_lables):
            X_t = split_for_training(x_train_val,train_index)
            y_t = split_for_training(y_train_val,train_index)
            cluster_label_t = split_for_training(cluster_lables,train_index)
            cv_inner = StratifiedKFoldWithLabels(n_splits=N_FOLDS, labels=cluster_label_t,shuffle=True, random_state=seed)

            bayes = BayesSearchCV(
            regressor,
            search_space,
            n_iter=BO_ITER,
            cv=cv_inner,
            n_jobs=-1,
            random_state=seed,
            refit=True,
            scoring=scoring,
            return_train_score=True,
            )

            bayes.fit(X_t, y_t)
            # print(f"\n\nBest parameters: {bayes.best_params_}\n\n")
            estimators.append(bayes)
    else: 
        cv = KFold(n_splits=N_FOLDS, shuffle=True, random_state=seed)
        for train_index, val_index in cv.split(x_train_val, y_train_val):
            X_t = split_for_training(x_train_val,train_index)
            y_t = split_for_training(y_train_val,train_index)
            cv_inner = KFold(n_splits=N_FOLDS, shuffle=True, random_state=seed)
            bayes = BayesSearchCV(
                regressor,
                search_space,
                n_iter=BO_ITER,
                cv=cv_inner,
                n_jobs=-1,
                random_state=seed,
                refit=True,
                scoring=scoring,
                return_train_score=True,
                )
            
            bayes.fit(X_t, y_t)
            # print(f"\n\nBest parameters: {bayes.best_params_}\n\n")
            estimators.append(bayes)

    best_idx: int = np.argmax([est.best_score_ for est in estimators])
    best_estimator: Pipeline = estimators[best_idx].best_estimator_
    try:
        regressor_params: dict = best_estimator.named_steps.regressor.get_params()
        regressor_params = remove_unserializable_keys(regressor_params)
    except:
        regressor_params = {"bad params": "couldn't get them"}

    return best_estimator, regressor_params



