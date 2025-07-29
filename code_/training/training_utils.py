from pathlib import Path
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
from sklearn.gaussian_process.kernels import RBF, Matern, RationalQuadratic, ConstantKernel as C
from scipy.stats import invgamma

# from iterstrat.ml_stratifiers import MultilabelStratifiedKFold
from data_handling import remove_unserializable_keys, save_results
from filter_data import filter_dataset
from all_factories import (
                            regressor_factory,
                            # regressor_search_space,
                            transforms,
                            construct_kernel,
                            get_regressor_search_space)
from all_factories import optimized_models

from imputation_normalization import preprocessing_workflow
from scoring import (
    cross_validate_regressor,
    process_scores,
    _average_ls
)

from _custom_kernel import AdditiveRBF, JaccardRBF

HERE: Path = Path(__file__).resolve().parent


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



def train_regressor(
    dataset: pd.DataFrame,
    features_impute: Optional[list[str]],
    special_impute: Optional[str],
    representation: Optional[str],
    structural_features: Optional[list[str]],
    numerical_feats: Optional[list[str]],
    unroll: Union[dict[str, str], list[dict[str, str]], None],
    regressor_type: str,
    target_features: str,
    transform_type: str=None,
    second_transformer:str=None,
    hyperparameter_optimization: bool=True,
    imputer: Optional[str] = None,
    cutoff:Dict[str, Tuple[Optional[float], Optional[float]]]=None,
    kernel: Optional[str] = None,
    Test:bool=False,
    classification:bool=False,
    ) -> None:
        """
        you should change the name here for prepare
        """
            #seed scores and seed prediction
        set_globals(Test)
        scores, predictions, ls = _prepare_data(
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
                                                    kernel=kernel,
                                                    classification=classification,
                                                    )
        scores = process_scores(scores,classification)
        ls =_average_ls(ls)
  
        return scores, predictions, ls
        



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
    second_transformer:str = None,
    hyperparameter_optimization: bool = True,
    imputer: Optional[str] = None,
    cutoff: Dict[str, Tuple[Optional[float], Optional[float]]]=None,
    classification:bool=False,
    kernel:Optional[str] = None,
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
    score,predication, ls = run(
                                X,
                                y,
                                preprocessor=preprocessor,
                                second_transformer=second_transformer,
                                regressor_type=regressor_type,
                                transform_type=transform_type,
                                hyperparameter_optimization=hyperparameter_optimization,
                                kernel=kernel,
                                classification=classification,
                                **kwargs,
                                )
    # print(X_y_shape)
    y_frame = pd.DataFrame(y.flatten(),columns=target_features)
    combined_prediction_ground_truth = pd.concat([predication, y_frame], axis=1)

    return score, combined_prediction_ground_truth, ls

def run(
    X, y, preprocessor: Union[ColumnTransformer, Pipeline], classification:bool,second_transformer:str, regressor_type: str,
    transform_type: str, hyperparameter_optimization: bool = True,
    kernel:Optional[str] = None,**kwargs,
    ) -> tuple[dict[int, dict[str, float]], pd.DataFrame]:

    seed_scores: dict[int, dict[str, float]] = {}
    seed_predictions: dict[int, np.ndarray] = {}
    length_scale_fitted_model = defaultdict(dict)
    # if 'Rh (IW avg log)' in target_features:
    #     y = np.log10(y)

    
    for seed in SEEDS:
        cv_outer = get_default_kfold_splitter(n_splits=N_FOLDS,classification=classification,random_state=seed)
        #   cv_outer = KFold(n_splits=N_FOLDS, shuffle=True, random_state=seed)
        y_transform = get_target_transformer(transform_type,second_transformer)
        if hyperparameter_optimization:
            search_space = get_regressor_search_space(regressor_type,kernel)
            kernel = construct_kernel(regressor_type, kernel)
            if classification:
                skop_scoring = "f1"

                if y.shape[1] > 1:
                    y_transform_regressor = MultiOutputClassifier(regressor_factory[regressor_type],n_jobs=-1)
                    search_space = {
                    f"regressor__estimator__{key.split('__')[-1]}": value
                    for key, value in search_space.items()
                        }
                    
                else:
                    y_transform_regressor = regressor_factory[regressor_type](kernel=kernel) if kernel!=None else regressor_factory[regressor_type]
                                
                

            else:
                skop_scoring = "neg_root_mean_squared_error"
                if y.shape[1] > 1:
                    y_transform_regressor = TransformedTargetRegressor(
                    regressor = MultiOutputRegressor(
                    estimator=regressor_factory[regressor_type](kernel=kernel) if kernel is not None
                    else regressor_factory[regressor_type]
                    ),
                    transformer=y_transform,
                        )
                    
                    search_space = {
                    f"regressor__regressor__estimator__{key.split('__')[-1]}": value
                    for key, value in search_space.items()
                        }
                else:
                    y_transform_regressor = TransformedTargetRegressor(
                            regressor=regressor_factory[regressor_type](kernel=kernel) if kernel!=None
                                    else regressor_factory[regressor_type],
                            transformer=y_transform,
                    )
            new_preprocessor = 'passthrough' if len(preprocessor.steps) == 0 else preprocessor
            regressor :Pipeline= Pipeline(steps=[
                        ("preprocessor", new_preprocessor),
                        ("regressor", y_transform_regressor),
                            ])

            regressor.set_output(transform="pandas")
            cv_in = get_default_kfold_splitter(n_splits=N_FOLDS,classification=classification,random_state=seed)
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
                classification=classification
            )
            scores, predictions = cross_validate_regressor(
                best_estimator, X, y, cv_outer, classification=classification
            )
            scores["best_params"] = regressor_params


        else:
            if regressor_type == "sklearn-GPR":
            
                length_scale_init = np.full(X.shape[1], .5) if X.shape[1] < 500 else 1.5
                length_scale_bounds_inv_gamma = (invgamma.ppf(0.05, a=5, scale=5), invgamma.ppf(0.95, a=5, scale=5))
                if kernel =='rbf':
                    my_kernel = RBF(length_scale=length_scale_init,
                                    #  length_scale_bounds=length_scale_bounds_inv_gamma
                                     )
                elif kernel == 'matern':
                    my_kernel = Matern(length_scale=length_scale_init, nu=0.1,
                                        # length_scale_bounds=length_scale_bounds_inv_gamma
                                        )

                elif kernel == 'j_rbf':
                    my_kernel = JaccardRBF(length_scale=length_scale_init,
                                        # length_scale_bounds=length_scale_bounds_inv_gamma
                                        )
                elif kernel == 'a_rbf':
                    my_kernel = AdditiveRBF(length_scale=length_scale_init)


                else:
                    raise ValueError(f"kernel required, unsupported")
                
            model = optimized_models(regressor_type,kernel=my_kernel)
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
                
            scores, predictions = cross_validate_regressor(regressor, X, y, cv_outer)
        # print(scores)
        seed_scores[seed] = scores.copy()
        seed_scores[seed].pop("estimator", None)
        # length_scale_fitted_model = regressor.named_steps["regressor"].regressor.get_params()["estimator"].kernel_.length_scale
        seed_predictions[seed] = predictions.flatten()

        seed_predictions: pd.DataFrame = pd.DataFrame.from_dict(
                        seed_predictions, orient="columns")
        if regressor_type == "sklearn-GPR":
            if X.shape[1] < 500:
                for i, est in enumerate(scores["estimator"]):
                    feature_names = est.named_steps["preprocessor"].get_feature_names_out()
                    length_scales = est.named_steps["regressor"].regressor_ .kernel_.length_scale
                    length_scale_fitted_model[seed][i] = dict(zip(feature_names, length_scales))
            else:
                length_scale_fitted_model = None
        else:
            length_scale_fitted_model = None
    return seed_scores, seed_predictions, length_scale_fitted_model


def _optimize_hyperparams(
    X, y, cv_outer, cv_in,seed: int,n_iter:int,  regressor_type:str, search_space:dict, regressor: Pipeline, classification:bool,
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


def split_for_training(
    data: Union[pd.DataFrame, np.ndarray,pd.Series], indices: np.ndarray
) -> Union[pd.DataFrame, np.ndarray, pd.Series]:
    if isinstance(data, pd.DataFrame):
        split_data = data.iloc[indices]
    elif isinstance(data, np.ndarray):
        split_data = data[indices]
    elif isinstance(data, pd.Series):
        split_data = data.iloc[indices]
    else:
        raise ValueError("Data must be either a pandas DataFrame, Series, or a numpy array.")
    return split_data


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




def get_target_transformer(transformer:str,extra_transformer:str) -> Pipeline:
    
    if extra_transformer:
        return Pipeline(steps=[
            ("extra y transform", transforms[extra_transformer]), 
            ("y scaler", transforms[transformer])  
            ])
    else:
        return Pipeline(steps=[
            ("y scaler", transforms[transformer])  # StandardScaler to standardize the target
            ])



def get_default_kfold_splitter(n_splits: int, classification:bool, random_state:int):
    if classification:
        return StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    else:
        return KFold(n_splits=n_splits, shuffle=True, random_state=random_state)
