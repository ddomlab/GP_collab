import json
from pathlib import Path
import time
from typing import Callable, Optional, Union, Dict, Tuple

import numpy as np
import pandas as pd
from collections import defaultdict

from sklearn.compose import ColumnTransformer, TransformedTargetRegressor
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from skopt import BayesSearchCV
from sklearn.preprocessing import FunctionTransformer
from sklearn.multioutput import MultiOutputRegressor,MultiOutputClassifier
# from optuna.integration import OptunaSearchCV
from scipy.stats import invgamma

# from iterstrat.ml_stratifiers import MultilabelStratifiedKFold
from data_handling import remove_unserializable_keys, save_results
from filter_data import filter_dataset
from all_factories import (
                            regressor_factory,
                            # regressor_search_space,
                            transforms,
                            get_regressor_search_space,
                            _get_gp_kernel)
from all_factories import optimized_models

from imputation_normalization import preprocessing_workflow
from scoring import (
    cross_validate_regressor,
    process_scores,
    _average_ls,
    gp_cross_validate_regressor
)

from utils import split_for_training


HERE: Path = Path(__file__).resolve().parent


def set_globals(Test: bool=False) -> None:
    global SEEDS, N_FOLDS, BO_ITER
    if not Test:
        SEEDS = [6, 13, 42]
        N_FOLDS = 5
        BO_ITER = 42
    else:
        SEEDS = [42,13]
        N_FOLDS = 2
        BO_ITER = 1



def train_regressor(
    dataset: pd.DataFrame,
    features_impute: Optional[list[str]],
    special_impute: Optional[str],
    structural_features: Optional[list[str]],
    numerical_feats: Optional[list[str]],
    unroll: Union[dict[str, str], list[dict[str, str]], None],
    regressor_type: str,
    target_features: str,
    feat_transformer: str=None,
    target_transformer:str=None,
    hyperparameter_optimization: bool=True,
    cutoff:Dict[str, Tuple[Optional[float], Optional[float]]]=None,
    imputer: Optional[str] = None,
    Test:bool=False,
    **keyword,
    ) -> None:
        """
        you should change the name here for prepare
        """
            #seed scores and seed prediction
        set_globals(Test)
        start = time.time()
        scores, predictions = _prepare_data(
                                            dataset=dataset,
                                            features_impute= features_impute,
                                            special_impute= special_impute,
                                            structural_features=structural_features,
                                            unroll=unroll,
                                            numerical_feats = numerical_feats,
                                            target_features=target_features,
                                            regressor_type=regressor_type,
                                            transform_type=feat_transformer,
                                            second_transformer=target_transformer,
                                            imputer=imputer,
                                            cutoff=cutoff,
                                            hyperparameter_optimization=hyperparameter_optimization,
                                            **keyword
                                            )
        scores = process_scores(scores)
        end = time.time()
        scores["run_time_sec"] = np.round((end - start)/len(SEEDS), 3)
        # if ls is not None:
        #     ls = _average_ls(ls)
  
        return scores, predictions
        

def create_feature_groups(
    unrolled_feats: Optional[list[str]], 
    unroll_info: Union[dict, list, None], 
    numerical_feats: Optional[list[str]]
) -> dict[str, list[str]]:
    feat_group = {}

    if unrolled_feats and unroll_info:
        units = []
        if isinstance(unroll_info, dict):
            units = unroll_info.get("unit_name", [])
        elif isinstance(unroll_info, list):
            for d in unroll_info:
                u = d.get("unit_name", [])
                if isinstance(u, list):
                    units.extend(u)
                else:
                    units.append(u)

        for unit in units:
            unit_cols = [col for col in unrolled_feats if col.startswith(f"{unit}_")]
            if unit_cols:
                feat_group[f'fp_{unit}'] = unit_cols

    if numerical_feats:
        feat_group['count'] = numerical_feats

    return feat_group


def _prepare_data(
    dataset: pd.DataFrame,
    target_features: str,
    regressor_type: str,
    features_impute: Optional[list[str]]=None,
    special_impute: Optional[str]=None,
    structural_features: Optional[list[str]]=None,
    numerical_feats: Optional[list[str]]=None,
    unroll: Union[dict, list, None] = None,
    transform_type: str = "Standard",
    second_transformer:str = None,
    hyperparameter_optimization: bool = True,
    imputer: Optional[str] = None,
    cutoff: Dict[str, Tuple[Optional[float], Optional[float]]]=None,
    # kernel_type: Optional[str]=None,
    # kernel_mixing_method: Optional[str]=None,
    **kwargs,
    ) -> tuple[dict[int, dict[str, float]], pd.DataFrame]:


    """
    here you should change the names
    """



    X, y, unrolled_feats = filter_dataset(
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
                                                    numerical_feat=numerical_feats,
                                                    structural_feat=unrolled_feats,
                                                    special_column=special_impute,
                                                    scaler=transform_type
                                                    )
    


    preprocessor.set_output(transform="pandas")
    # if isinstance(regressor_type, dict):
    feat_group = create_feature_groups(unrolled_feats, unroll, numerical_feats)
    # else:
    #     feat_idx = None
   
    score,predication = run(
                            X,
                            y,
                            features_group=feat_group,
                            preprocessor=preprocessor,
                            second_transformer=second_transformer,
                            regressor_type=regressor_type,
                            hyperparameter_optimization=hyperparameter_optimization,
                            **kwargs,
                            )
    # print(X_y_shape)
    y_frame = pd.DataFrame(y.flatten(),columns=target_features)
    combined_prediction_ground_truth = pd.concat([predication, y_frame], axis=1)

    return score, combined_prediction_ground_truth


def run(
    X, y, 
    features_group: Optional[dict[str, list[int]]],
    preprocessor: Union[ColumnTransformer, Pipeline], 
    second_transformer:str, 
    regressor_type: str,
    hyperparameter_optimization: bool = True,
    **kwargs,
    ) -> tuple[dict[int, dict[str, float]], pd.DataFrame]:

    seed_scores: dict[int, dict[str, float]] = {}
    seed_predictions: dict[int, np.ndarray] = {}
    seed_indices: dict[int, np.ndarray] = {}

    for seed in SEEDS:

        print(f"Running seed: {seed}")
        cv_outer = get_default_kfold_splitter(n_splits=N_FOLDS,random_state=seed)
        y_transform = get_target_transformer(second_transformer)

        if hyperparameter_optimization:
            search_space = get_regressor_search_space(regressor_type)
            skop_scoring = "neg_root_mean_squared_error"

            if y.shape[1] > 1:
                y_transform_regressor = TransformedTargetRegressor(
                regressor = MultiOutputRegressor(
                estimator= regressor_factory[regressor_type]
                ),
                transformer=y_transform,
                )
                
                search_space = {
                f"regressor__regressor__estimator__{key.split('__')[-1]}": value
                for key, value in search_space.items()
                    }
            else:
                y_transform_regressor = TransformedTargetRegressor(
                        regressor= regressor_factory[regressor_type],
                        transformer=y_transform,
                )

            new_preprocessor = 'passthrough' if len(preprocessor.steps) == 0 else preprocessor
            regressor :Pipeline= Pipeline(steps=[
                        ("preprocessor", new_preprocessor),
                        ("regressor", y_transform_regressor),
                            ])

            regressor.set_output(transform="pandas")
            cv_in = get_default_kfold_splitter(n_splits=N_FOLDS,random_state=seed)
            best_estimator, regressor_params = _optimize_hyperparams(
                                                                    X,
                                                                    y,
                                                                    cv_outer=cv_outer,
                                                                    cv_in=cv_in,
                                                                    n_iter=BO_ITER,
                                                                    seed=seed,
                                                                    regressor_type=regressor_type,
                                                                    search_space=search_space,
                                                                    regressor=regressor,
                                                                    scoring=skop_scoring,
                                                                    )
            
            scores, predictions = cross_validate_regressor(
                                    best_estimator, X, y, cv_outer
                                    )
            scores["best_params"] = regressor_params


        else:
            print("No hyperparameter optimization")
            if isinstance(regressor_type, dict):
                kernel = _get_gp_kernel(regressor_type, idx=features_group)
                model = optimized_models(regressor_type['model'], kernel=kernel)
            else:
                model = optimized_models(
                                        regressor_type,
                                        feat_group=features_group,
                                        **kwargs
                                        )

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
            # y = y.ravel()
            y = y.flatten()
            # return_importance = False if "GP" in regressor_type else True
            if "gp" in regressor_type.lower():
                scores, predictions = gp_cross_validate_regressor(regressor, X, y, cv_outer, return_ls=True)
            else:
                scores, predictions = cross_validate_regressor(regressor, X, y, cv_outer,
                                                                custom=True,
                                                                early_stopping=True,
                                                                return_estimator=False,
                                                                return_feature_importances=True,
                                                                
                                                                )
        seed_scores[seed] = scores.copy()
        seed_scores[seed].pop("estimator", None)
        # seed_indices[seed] = indices
        # length_scale_fitted_model = regressor.named_steps["regressor"].regressor.get_params()["estimator"].kernel_.length_scale
        seed_predictions[seed] = predictions.flatten()

        seed_predictions: pd.DataFrame = pd.DataFrame.from_dict(
                        seed_predictions, orient="columns")
    # seed_indices_full_dropped_nans.json
    
    # save_folder = Path(r"C:\Users\sdehgha2\Desktop\phd-code\gp_collab\GP_collab\results\Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation\target_log Rg (nm)")
    # with open(save_folder / "seed_indices_full_dropped_nans.json", "w") as f:
    #     json.dump(seed_indices, f, cls=NumpyArrayEncoder, indent=2)

    return seed_scores, seed_predictions



import torch
def _to_torch_tensor(data):
    """
    Converts feature data (X) from DataFrame/NumPy to a float tensor
    and moves it to CUDA if self.use_cuda is True.
    """
    if isinstance(data, pd.DataFrame):
        data = data.to_numpy(dtype=np.float32)
    elif isinstance(data, np.ndarray):
        data = data.astype(np.float32)
    
    # Convert to tensor
    tensor = torch.tensor(data, dtype=torch.float)
    
    # Move to device
    return tensor



class NumpyArrayEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, tuple):
            return list(obj)
        else:
            return super(NumpyArrayEncoder, self).default(obj)


def _optimize_hyperparams(
    X, y, 
    cv_outer, 
    cv_in,seed: int,
    n_iter:int,
    regressor_type:str,
    search_space:dict,
    regressor: Pipeline,
    scoring:Union[str,Callable]) -> tuple:

    # Splitting for outer cross-validation loop
    estimators: list[BayesSearchCV] = []
    for train_index, _ in cv_outer.split(X, y):

        X_train = split_for_training(X, train_index)
        y_train = split_for_training(y, train_index)
        # print(X_train)
        # Splitting for inner hyperparameter optimization loop
        # cv_inner = KFold(n_splits=N_FOLDS, shuffle=True, random_state=seed)
        print("\n\n")
        print(
            "OPTIMIZING HYPERPARAMETERS FOR REGRESSOR", regressor_type, "\tSEED:", seed
        )
        # Bayesian hyperparameter optimization
        bayes = BayesSearchCV(
            regressor,
            search_space,
            n_iter=n_iter,
            cv=cv_in,
            n_jobs=-1,
            random_state=seed,
            refit=True,
            scoring=scoring,
            return_train_score=True,
        )
        bayes.fit(X_train, y_train)

        print(f"\n\nBest parameters: {bayes.best_params_}\n\n")
        estimators.append(bayes)

    # Extract the best estimator from hyperparameter optimization
    best_idx: int = np.argmax([est.best_score_ for est in estimators])
    best_estimator: Pipeline = estimators[best_idx].best_estimator_
    try:
        regressor_params: dict = best_estimator.named_steps.regressor.get_params()
        regressor_params = remove_unserializable_keys(regressor_params)
    except:
        regressor_params = {"bad params": "couldn't get them"}

    return best_estimator, regressor_params



def _pd_to_np(data):
    if isinstance(data, pd.DataFrame):
        return data.values
    elif isinstance(data, np.ndarray):
        return data
    else:
        raise ValueError("Data must be either a pandas DataFrame or a numpy array.")

# def custom_function(x):
#     x = x.astype(float)  
#     x[:, 2] = np.log10(x[:, 2] + 1e-6)  
#     return x




def get_target_transformer(y_transformer:str) -> Pipeline:
    return Pipeline(steps=[
            ("y scaler", transforms[y_transformer])  # StandardScaler to standardize the target
            ])



def get_default_kfold_splitter(n_splits: int,random_state:int):
        
    return KFold(n_splits=n_splits, shuffle=True, random_state=random_state)

