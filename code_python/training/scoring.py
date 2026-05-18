from itertools import product
from typing import Callable, Union, Dict, List
from xml.parsers.expat import model
from joblib import Parallel, delayed
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, kendalltau
from sklearn.pipeline import Pipeline
from sklearn.metrics._scorer import r2_scorer
from sklearn.model_selection import (
    cross_val_predict,
    cross_validate,
    train_test_split,
    learning_curve
    )

# from _validation import multioutput_cross_validate 
from utils import split_for_training
import shap
from sklearn.base import clone
from sklearn.metrics import (
    make_scorer,
    mean_absolute_error,
    mean_squared_error,
    root_mean_squared_error,
    r2_score,
    roc_auc_score,
    accuracy_score,
    f1_score,
    recall_score,
    precision_score,
)

from utils_uncertainty_calibration import (
    compute_ece,
    compute_cdf_ama,
    compute_cvpp_ama,
    gaussian_nll,
    compute_RUSC,
    )


import inspect
import copy


def pearson(y_true: pd.Series, y_pred: np.ndarray) -> float:
    if isinstance(y_true, pd.Series) or isinstance(y_true, pd.DataFrame):
        y_true = y_true.to_numpy()
    y_true = y_true.flatten()
    y_pred = y_pred.flatten()
    r = pearsonr(y_true, y_pred)[0]
    return r
r_scorer = make_scorer(pearson, greater_is_better=True)

def rmse_score(y_test: pd.Series, y_pred: pd.Series) -> float:
    """
    Calculate the root mean squared error.

    Args:
        y_test: Test targets.
        y_pred: Predicted targets.

    Returns:
        Root mean squared error.
    """
    return root_mean_squared_error(y_test, y_pred)


def flatten_array(arr):
    if isinstance(arr, (pd.Series, pd.DataFrame)):
        arr = arr.to_numpy()
    return arr.flatten()

def pearson_r(y_true, y_pred):
    y_true, y_pred = flatten_array(y_true), flatten_array(y_pred)
    return pearsonr(y_true, y_pred)[0]

def pearson_p(y_true, y_pred):
    y_true, y_pred = flatten_array(y_true), flatten_array(y_pred)
    return pearsonr(y_true, y_pred)[1]

def spearman_r(y_true, y_pred):
    y_true, y_pred = flatten_array(y_true), flatten_array(y_pred)
    return spearmanr(y_true, y_pred)[0]

def spearman_p(y_true, y_pred):
    y_true, y_pred = flatten_array(y_true), flatten_array(y_pred)
    return spearmanr(y_true, y_pred)[1]

def kendall_r(y_true, y_pred):
    y_true, y_pred = flatten_array(y_true), flatten_array(y_pred)
    return kendalltau(y_true, y_pred)[0]

def kendall_p(y_true, y_pred):
    y_true, y_pred = flatten_array(y_true), flatten_array(y_pred)
    return kendalltau(y_true, y_pred)[1]



# multiouput regression scorers and classification scorers
def r2_scorer_multi(y_true, y_pred):
    return r2_score(y_true, y_pred, multioutput="raw_values")

def rmse_scorer_multi(y_true, y_pred):
    return root_mean_squared_error(y_true, y_pred, multioutput="raw_values")

def mae_scorer_multi(y_true, y_pred):
    return mean_absolute_error(y_true, y_pred, multioutput="raw_values")


def accuracy_scorer_multi(y_true, y_pred):
    return accuracy_score(y_true, y_pred)


def f1_scorer_multi(y_true, y_pred):
    return f1_score(y_true, y_pred, average=None) 

def roc_auc_scorer_multi(y_true, y_pred):
    return roc_auc_score(y_true, y_pred, average=None)

def recall_scorer_multi(y_true, y_pred):
    return recall_score(y_true, y_pred, average=None)

def precision_scorer_multi(y_true, y_pred):
    return precision_score(y_true, y_pred, average=None)


r_scorer = make_scorer(pearson_r, greater_is_better=True)
spearman_scorer = make_scorer(spearman_r, greater_is_better=True)
kendall_scorer = make_scorer(kendall_r, greater_is_better=True)
rmse_scorer = make_scorer(root_mean_squared_error, greater_is_better=False)
mae_scorer = make_scorer(mean_absolute_error, greater_is_better=False)
r2_scorer = make_scorer(r2_score, greater_is_better=True)


# P-value scorers (not typically used for optimization but useful for logging)
pearson_p_scorer = make_scorer(pearson_p, greater_is_better=False)
spearman_p_scorer = make_scorer(spearman_p, greater_is_better=False)
kendall_p_scorer = make_scorer(kendall_p, greater_is_better=False)



f1_scorer = make_scorer(f1_score, greater_is_better=True)
roc_auc_scorer = make_scorer(roc_auc_score, greater_is_better=True)
recall_scorer = make_scorer(recall_score, greater_is_better=True)
precision_scorer = make_scorer(precision_score, greater_is_better=True)


def process_scores(
    scores: dict[int, dict[str, float]],
) -> dict[Union[int, str], dict[str, float]]:
    
    first_key = list(scores.keys())[0]
    
    score_types: list[str] = [
        key for key in scores[first_key].keys() 
        if key.startswith("test_") and "lengthscale" not in key
    ]

    arr = np.array(scores[first_key]["test_r2"])
    
    if arr.ndim > 1 and arr.shape[1] > 1:
        avg_r2 = np.round(np.mean(np.vstack([arr for seed in scores.values() for arr in seed["test_r2"]]), axis=0), 3)
        stdev_r2 = np.round(np.std(np.vstack([arr for seed in scores.values() for arr in seed["test_r2"]]), axis=0), 3)
        print("Average scores:\t", f"r2: {avg_r2}±{stdev_r2}")
        
        avgs: list[float] = [
            np.mean(np.vstack([arr for seed in scores.values() for arr in seed[score]]), axis=0) for score in score_types
        ]
        stdevs: list[float] = [
            np.std(np.vstack([arr for seed in scores.values() for arr in seed[score]]), axis=0) for score in score_types
        ]
    else:
        avg_rmse = round(np.mean([seed["test_rmse"] for seed in scores.values()]), 2)
        stdev_rmse = round(np.std([seed["test_rmse"] for seed in scores.values()]), 2)
        avg_r2 = round(np.mean([seed["test_r2"] for seed in scores.values()]), 2)
        stdev_r2 = round(np.std([seed["test_r2"] for seed in scores.values()]), 2)
        print("Average scores:\t",
            f"rmse: {abs(avg_rmse)}±{stdev_rmse}\t",
            f"r2: {avg_r2}±{stdev_r2}")

        avgs: list[float] = [
            np.mean([seed[score] for seed in scores.values()]) for score in score_types
        ]
        stdevs: list[float] = [
            np.std([seed[score] for seed in scores.values()]) for score in score_types
        ]

    clean_score_types: list[str] = [score.replace("test_", "") for score in score_types]
    for score, avg, stdev in zip(clean_score_types, avgs, stdevs):
        scores[f"{score}_avg"] = abs(avg) if score in ["rmse", "mae"] else avg
        scores[f"{score}_stdev"] = stdev
    
    if arr.ndim > 1 and arr.shape[1] > 1:
        for score in clean_score_types:
            scores[f"{score}_avg_aggregate"] = np.mean(scores[f"{score}_avg"])
            scores[f"{score}_stdev_aggregate"] = np.mean(scores[f"{score}_stdev"])

    return scores


def _average_ls(ls_data: Dict) -> None:
    feature_values = defaultdict(list)
    for _, folds in ls_data.items():
        for _, features in folds.items():
            for feature, value in features.items():
                feature_values[feature].append(value)

    stats = {
        feature: {
            "mean": float(np.mean(values)),
            "std": float(np.std(values, ddof=1))
        }
        for feature, values in feature_values.items()
    }

    ls_data["aggregated_ls"] = stats

    return ls_data


def compute_summary_stats(metrics: Dict[str, list[float]]) -> Dict[str, float]:
    """
    Helper function to compute the mean and standard deviation of the test/train metrics.
    """
    summary_stats = {}
    
    for metric, values in metrics.items():
        values = np.array(values, dtype=np.float64)  # Convert to NumPy array for handling NaNs
        valid_values = values[~np.isnan(values)]  # Remove NaNs before computing stats
        
        if valid_values.size > 0:
            # Use absolute mean for RMSE and MAE
            summary_stats[f"{metric}_mean"] = abs(np.mean(valid_values)) if metric in ["test_rmse", "train_rmse", "test_mae", "train_mae"] else np.mean(valid_values)
            summary_stats[f"{metric}_std"] = np.std(valid_values)
        else:
            summary_stats[f"{metric}_mean"] = np.nan
            summary_stats[f"{metric}_std"] = np.nan

    return summary_stats


def process_ood_scores(
    scores: Dict[int, Dict[str, Union[float, list]]]
) -> Dict[Union[int, str], Dict[str, float]]:
    for cluster, seed_data in scores.items():
        if cluster.startswith("CO_") or cluster.startswith("ID_"):
            test_metrics = {}

            # Collect test scores across seeds
            for seed, metrics in seed_data.items():
                if isinstance(seed, int):
                    for key, value in metrics.items():
                        if key.startswith("test_") and not isinstance(value, dict):
                            if isinstance(value, list):
                                test_metrics.setdefault(key, []).extend(value)
                            else:
                                test_metrics.setdefault(key, []).append(value)

            # Compute summary statistics for the cluster
            summary_stats = compute_summary_stats(test_metrics)
            
            # Store the results back in the scores dictionary
            scores[cluster]['summary_stats'] = summary_stats

    return scores


from collections import defaultdict

def process_ood_learning_curve_score(scores: dict) -> dict:
    """
    Process learning curve results for both CO_ and IID_ clusters:
    - CO_: aggregate across seeds
    - IID_: aggregate across seeds and test_set_seeds
    """
    for cluster, score_dict in scores.items():
        if cluster.startswith("CO_"):
            for train_ratio, seeds in score_dict.items():
                if train_ratio == "training size":
                    continue

                test_metrics = defaultdict(list)
                train_metrics = defaultdict(list)

                for seed, (train_results, test_results) in seeds.items():
                    for metric, value in test_results.items():
                        value = value.item() if isinstance(value, np.ndarray) else value
                        test_metrics[metric].append(value)
                    for metric, value in train_results.items():
                        value = value.item() if isinstance(value, np.ndarray) else value
                        train_metrics[metric].append(value)

                scores[cluster][train_ratio].update({
                    "test_summary_stats": compute_summary_stats(test_metrics),
                    "train_summary_stats": compute_summary_stats(train_metrics),
                })

        elif cluster.startswith("IID_"):
            for train_ratio, seeds_dict in score_dict.items():
                if train_ratio == "training size":
                    continue

                test_metrics = defaultdict(list)
                train_metrics = defaultdict(list)

                for seed, test_seed_dict in seeds_dict.items():
                    for test_seed, (train_results, test_results) in test_seed_dict.items():
                        for metric, value in test_results.items():
                            value = value.item() if isinstance(value, np.ndarray) else value
                            test_metrics[metric].append(value)
                        for metric, value in train_results.items():
                            value = value.item() if isinstance(value, np.ndarray) else value
                            train_metrics[metric].append(value)

                scores[cluster][train_ratio].update({
                    "test_summary_stats": compute_summary_stats(test_metrics),
                    "train_summary_stats": compute_summary_stats(train_metrics),
                })

    return scores



def process_learning_score(score: dict[int, dict[str, np.ndarray]]):
     # Initialize arrays for aggregation
    train_scores_mean = None
    train_scores_std = None
    test_scores_mean = None
    test_scores_std = None

    # Loop over seeds and accumulate results
    for _, results in score.items():
        if train_scores_mean is None:
            # Initialize mean and std with the first seed's results
            train_scores_mean = results["train_scores"].mean(axis=1, keepdims=True)
            train_scores_std = results["train_scores"].std(axis=1, keepdims=True)
            test_scores_mean = results["test_scores"].mean(axis=1, keepdims=True)
            test_scores_std = results["test_scores"].std(axis=1, keepdims=True)
        else:
            # Accumulate the means and stds
            train_scores_mean += results["train_scores"].mean(axis=1, keepdims=True)
            train_scores_std += results["train_scores"].std(axis=1, keepdims=True)
            test_scores_mean += results["test_scores"].mean(axis=1, keepdims=True)
            test_scores_std += results["test_scores"].std(axis=1, keepdims=True)

    # Calculate the averages over the number of seeds
    num_seeds = len(score)
    score['aggregated_results']= {
        "train_sizes": results["train_sizes"],
        "train_sizes_fraction": results["train_sizes_fraction"],
        "train_scores_mean": train_scores_mean / num_seeds,
        "train_scores_std": train_scores_std / num_seeds,
        "test_scores_mean": test_scores_mean / num_seeds,
        "test_scores_std": test_scores_std / num_seeds,
    }
    return score



def get_feature_importances_from_cv(score: dict, X: np.ndarray | None = None) -> Dict[str, List[Dict[str, float]]]:
    all_MDI_importances = []
    all_shap_importances = []

    first_est = score["estimator"][0]
    preprocessor = first_est.named_steps["preprocessor"]
    feature_names = preprocessor.get_feature_names_out()
    for est in score["estimator"]:
        model = est.named_steps["regressor"].regressor_

        # --- Model-based importances ---
        raw_fi = model.feature_importances_

        fi_model = raw_fi[0] if model.__class__.__name__ == "NGBRegressor" else raw_fi
        all_MDI_importances.append(dict(zip(feature_names, fi_model)))

        # --- SHAP importances ---
        if X is not None:
            X_transformed = preprocessor.transform(X)
            model_output = 0 if model.__class__.__name__ == "NGBRegressor" else "raw"
            try:
                explainer = shap.Explainer(model, X_transformed)
            except Exception:
                explainer = shap.TreeExplainer(model, model_output=model_output)
            shap_values = explainer(X_transformed, check_additivity=False)
            fi_shap = np.abs(shap_values.values).mean(axis=0)
            all_shap_importances.append(dict(zip(feature_names, fi_shap)))

    score["feature_importance_MDI"] = all_MDI_importances
    if all_shap_importances:
        score["feature_importance_SHAP"] = all_shap_importances




def _gp_fit_predict_score(estimator, model_type, X, y, train_idx, test_idx, scoring, return_ls: bool, UQ: bool):
    """
    Runs inside a parallel worker:
    - clone estimator
    - fit on train
    - predict on test
    - compute scores with provided scorers
    """
    if "mgk" in model_type.lower():
        est = copy.deepcopy(estimator)
    else:
        est = clone(estimator)

    # Use safe row selector for all data types
    X_train = split_for_training(X, train_idx)
    X_test  = split_for_training(X, test_idx)
    y_train = split_for_training(y, train_idx)
    y_test  = split_for_training(y, test_idx)

    # Fit model (with Pyro MCMC)
    est.fit(X_train, y_train)

    # Predict
    y_result = est.predict(X_test, return_std=UQ)
    results = {}
    if return_ls: 
        results["lengthscale"] = est.named_steps["regressor"]._get_lengthscale()
    # Compute scoring: scoring[name] is a scorer from make_scorer

    y_test = y_test.flatten()
    for name, scorer in scoring.items():
        results[name] = scorer(y_test, y_result["y_pred"])
    if UQ:
        UQ_scorers = {
            "ece": compute_ece,
            "RUSC": compute_RUSC,
            "cdf_ama": compute_cdf_ama,
            "cvpp_ama": compute_cvpp_ama,
            "nll": gaussian_nll,
        }
        for name, uq_scorer in UQ_scorers.items():
            results[name] = float(uq_scorer(y_test, y_result["y_pred"], y_result["y_std"]))

    return test_idx, y_result, results


def gp_cross_validate(
    estimator,
    model_type:str,
    X,
    y,
    cv,
    scoring,
    UQ,
    return_ls,
    n_jobs,
):
    """
    Returns:
        scores: Dictionary of lists containing fold-wise metrics.
        predictions: Dictionary containing out-of-fold predictions:
            - predictions["y_pred"]: predicted means, shape (n_samples,)
            - predictions["y_std"]: predicted standard deviations, shape (n_samples,)
              only available when returned by the model.
    """
    # if "mgk" in model_type.lower():
    #     parallel_results = Parallel(n_jobs=n_jobs, verbose=0, require="sharedmem")(
    #     delayed(_mgk_fit_and_score_fold)(estimator, X, y, train_idx, test_idx, scoring, UQ)
    #     for train_idx, test_idx in cv.split(X, y)
    #     )
    # else:
    parallel_results = Parallel(n_jobs=n_jobs, verbose=0, require="sharedmem")(
            delayed(_gp_fit_predict_score)(
                estimator, model_type, X, y, train_idx, test_idx, scoring, return_ls, UQ
            )
            for train_idx, test_idx in cv.split(X, y)
        )

    scores = defaultdict(list)
    n_samples = len(y)

    predictions = {
        "y_pred": np.full(n_samples, np.nan),
        "y_std": np.full(n_samples, np.nan),
    }

    for test_idx, y_result, fold_scores in parallel_results:
        predictions["y_pred"][test_idx] = y_result["y_pred"]

        if "y_std" in y_result:
            predictions["y_std"][test_idx] = y_result["y_std"]

        for key, val in fold_scores.items():
            scores[f"test_{key}"].append(val)

    return scores, predictions


def gp_cross_validate_regressor(
    regressor, model_type:str, X, y, cv, n_jobs=1, return_ls: bool = False, UQ: bool = False
    ) -> tuple[dict[str, float], np.ndarray]:


        scorers = {
            "rmse": root_mean_squared_error,
            "mae": mean_absolute_error,
            "r2": r2_score,
            }
        
        score, predictions = gp_cross_validate(
            regressor,
            model_type,
            X,
            y,
            cv=cv,
            scoring=scorers,
            n_jobs=n_jobs,
            return_ls=return_ls,
            UQ=UQ
            )
        return score, predictions



def _custom_fit_predict_score(
        estimator, X, y,
        train_idx, 
        test_idx, 
        scoring: Dict,
        return_estimator: bool,
        return_feature_importances: bool,
        early_stopping:bool
        ):
    """
    - clone estimator
    - split train into train / eval
    - fit with early stopping (XGB or NGBoost)
    - predict on test
    - compute scores
    """

    model = clone(estimator)

    X_train= split_for_training(X, train_idx)
    y_train = split_for_training(y, train_idx)
    X_test = split_for_training(X, test_idx)
    y_test = split_for_training(y, test_idx)

    # Internal validation split (ONLY from training fold)
    fit_kwargs = {}
    if  early_stopping:
        print("Using early stopping during training...")
        X_train, X_eval, y_train, y_eval = train_test_split(
            X_train,
            y_train,
            test_size=0.2,
            random_state=42,
        )

        # ---- Early stopping dispatch ----
        preprocessor = clone(model.named_steps["preprocessor"])
        reg = model.named_steps["regressor"]
        inner_model = reg.regressor
        yt  = clone(reg.transformer) 
        y_train_2d = np.asarray(y_train).reshape(-1, 1)
        y_eval_2d  = np.asarray(y_eval).reshape(-1, 1)
        yt.fit(y_train_2d)
        # y_train_tgt = yt.transform()
        y_eval_tgt  = yt.transform(y_eval_2d)
        preprocessor.fit(X_train)

        X_train_t = preprocessor.transform(X_train)
        X_eval_t  = preprocessor.transform(X_eval)
        X_test_t  = preprocessor.transform(X_test)



        fit_sig = inspect.signature(inner_model.fit).parameters

        # XGBoost / LightGBM style
        if "eval_set" in fit_sig:
            fit_kwargs["eval_set"] = [(X_eval_t, y_eval_tgt)]
            fit_kwargs["verbose"] = False

        elif "X_val" in fit_sig and "Y_val" in fit_sig:
            fit_kwargs["X_val"] = X_eval_t
            fit_kwargs["Y_val"] = y_eval_tgt
        else:
            raise ValueError("Early stopping requested but the model does not support eval_set or X_val/Y_val.")

        reg.fit(X_train_t, y_train, **fit_kwargs)
        # print(f"best score, {reg.regressor_.best_score}")
        # best_iteration = reg.regressor_.best_iteration
        # print(f"Best iteration: {best_iteration}")
        y_pred = reg.predict(X_test_t)
        scores = {
            f"test_{name}": scorer(y_test, y_pred)
            for name, scorer in scoring.items()
        }

        if return_feature_importances:
            # MDI_importances = []
            # shap_importances = []

            model_inner = reg.regressor_   # ✅ FIX
            feature_names = preprocessor.get_feature_names_out()

            raw_fi = model_inner.feature_importances_
            feat_imp = raw_fi[0] if model_inner.__class__.__name__ == "NGBRegressor" else raw_fi
            # MDI_importances.append(dict(zip(feature_names, feat_imp)))
            scores["feature_importance_MDI"] = dict(zip(feature_names, feat_imp))

            model_output = 0 if model_inner.__class__.__name__ == "NGBRegressor" else "raw"
            explainer = shap.TreeExplainer(model_inner, model_output=model_output)

            shap_values = explainer(X_train_t)
            fi_shap = np.abs(shap_values.values).mean(axis=0)
            # shap_importances.append(dict(zip(feature_names, fi_shap)))

            scores["feature_importance_SHAP"] = dict(zip(feature_names, fi_shap))
        if return_estimator:
            scores["estimator"] = {
                "preprocessor": preprocessor,
                "regressor": reg,
            }
    else:
        model.fit(X_train, y_train, **fit_kwargs)
        y_pred = model.predict(X_test)
        scores = {
            f"test_{name}": scorer(y_test, y_pred)
            for name, scorer in scoring.items()
        }
    
        if return_feature_importances:
            
            MDI_importances = []
            shap_importances = []
            preprocessor = model.named_steps["preprocessor"]
            feature_names = preprocessor.get_feature_names_out()
            model_inner = model.named_steps["regressor"].regressor_
            raw_fi = model_inner.feature_importances_
            feat_imp = raw_fi[0] if model_inner.__class__.__name__ == "NGBRegressor" else raw_fi
            # MDI_importances.append(dict(zip(feature_names, feat_imp)))
            scores["feature_importance_MDI"] = dict(zip(feature_names, feat_imp)) 
            model_output = 0 if model_inner.__class__.__name__ == "NGBRegressor" else "raw"
            explainer = shap.TreeExplainer(model_inner, model_output=model_output)
            x_t_transformed = preprocessor.transform(X_train)
            shap_values = explainer(x_t_transformed)
            fi_shap = np.abs(shap_values.values).mean(axis=0)
            # shap_importances.append(dict(zip(feature_names, fi_shap)))
            scores["feature_importance_SHAP"] = dict(zip(feature_names, fi_shap))
        if return_estimator:
            scores["estimator"] = model
    
    return test_idx, y_pred, scores


def custom_cross_validate(
    estimator,
    X,
    y,
    cv,
    scoring,
    n_jobs=-1,
    return_estimator: bool = False,
    return_feature_importances: bool = False,
    early_stopping:bool = False
):
    
    parallel_results = Parallel(n_jobs=n_jobs, verbose=0, pre_dispatch="all")(
        delayed(_custom_fit_predict_score)(
            estimator, X, y, train_idx, test_idx, scoring, return_estimator, return_feature_importances, early_stopping
        )
        for train_idx, test_idx in cv.split(X, y)
    )

    scores = defaultdict(list)
    n_samples = len(y)
    predictions = np.full(n_samples, np.nan)
        
    for test_idx, y_pred, fold_scores in parallel_results:
        predictions[test_idx] = y_pred

        for key, val in fold_scores.items():
                scores[key].append(val)
                
    return scores, predictions




def cross_validate_regressor(
    regressor, X, y, cv,
    custom: bool = False,
    early_stopping: bool = False,
    return_estimator: bool = False,
    return_feature_importances: bool = False,
    return_indices: bool = False,
    n_jobs: int = -1
    ) -> tuple[dict[str, float], np.ndarray]:

        # MULTIOUPUT 
        # if y.ndim > 1 and y.shape[1] > 1:
  
        #     scorers = {
        #             "r2": r2_scorer_multi,
        #             "rmse": rmse_scorer_multi,
        #             "mae": mae_scorer_multi
        #             }
        

        #     score =  multioutput_cross_validate(
        #         estimator= regressor,
        #         X=X,
        #         y= y,
        #         cv=cv,
        #         scorers=scorers,
        #         n_jobs=-1,
        #         verbose=0
        #         )

        # # SINGLE OUTPUT
        # else:
        

        if custom   :
            scorers = {
            "rmse": root_mean_squared_error,
            "mae": mean_absolute_error,
            "r2": r2_score,
            }
            score, predictions = custom_cross_validate(
                regressor,
                X,
                y,
                cv=cv,
                scoring=scorers,
                n_jobs=n_jobs,
                return_estimator=return_estimator,
                return_feature_importances=return_feature_importances,
                early_stopping=early_stopping
            )
        else:
            scorers = {
            "rmse": rmse_scorer,
            "mae": mae_scorer,
            "r2": r2_scorer,
            }

            score: dict[str, float] = cross_validate(
                regressor,
                X,
                y,
                cv=cv,
                scoring=scorers,
                return_estimator=True,
                n_jobs=n_jobs,
                return_indices=return_indices,
                )
            predictions: np.ndarray = cross_val_predict(
                regressor,
                X,
                y,
                cv=cv,
                n_jobs=n_jobs,
            )
            if return_feature_importances:
                get_feature_importances_from_cv(score, X=X)

            if return_indices:
                indexes = score.pop("indices")
                return score, predictions, indexes
        return score, predictions



def get_incremental_split(
        regressor_params, X, y, cv,
        random_state:int,
        train_ratio:np.ndarray,
        scoring:str='r2',
    ) -> tuple:
     
    training_sizes, training_scores, testing_scores = learning_curve(
                                                        regressor_params,
                                                        X,
                                                        y,
                                                        cv=cv,
                                                        n_jobs=-1,
                                                        train_sizes=train_ratio,
                                                        scoring=scoring,
                                                        shuffle=True,
                                                        random_state=random_state
                                                        )

 
    return training_sizes, training_scores, testing_scores


class PredictionUncertainty:
    def __init__(self, fitted_model):
        self.fitted_model = fitted_model

    def pred_dist(self, X) -> np.ndarray:
        X_array = X.values if hasattr(X, "values") else X
        all_preds = Parallel(n_jobs=-1)(
            delayed(tree.predict)(X_array) for tree in self.fitted_model.estimators_
        )
        return np.std(all_preds, axis=0)
    


def train_and_predict_ood(reg, X_train_val, y_train_val, X_test, y_test, 
                            algorithm:str, return_train_pred:bool=False,
                            manual_preprocessor:Pipeline=None) -> tuple:
    
    reg.fit(X_train_val, y_train_val)
    manual_preprocessor.fit(X_train_val)
    x_test_scaled = manual_preprocessor.transform(X_test)
    # print(x_test_scaled)
    if algorithm == 'NGB':
        y_test_predist = reg.named_steps['regressor'].regressor_.pred_dist(x_test_scaled)
        y_test_pred_uncertainty = np.array(np.sqrt(y_test_predist.var)).reshape(y_test.shape)
        y_test_pred = reg.predict(X_test)

    elif algorithm == 'RF':
        uncertainty_estimator = PredictionUncertainty(reg.named_steps['regressor'].regressor_)
        y_test_pred_uncertainty = uncertainty_estimator.pred_dist(x_test_scaled).reshape(y_test.shape)
        y_test_pred = reg.predict(X_test)
    else:
        y_test_pred = reg.predict(X_test)
        y_test_pred_uncertainty = None

    
    test_scores = get_prediction_scores(y_test, y_test_pred,'test')
    if return_train_pred:
        y_train_pred = reg.predict(X_train_val)
        train_scores = get_prediction_scores(y_train_val, y_train_pred,'train')
        return test_scores, train_scores, y_test_pred, y_test_pred_uncertainty
    
    return test_scores, y_test_pred, y_test_pred_uncertainty


def get_prediction_scores(y_true, y_pred, score_set:str='test'):
    return {
        f"{score_set}_mad": np.abs(y_pred - y_pred.mean()).mean(),
        f"{score_set}_ystd": y_pred.std(),
        f"{score_set}_mae": mean_absolute_error(y_true, y_pred),
        f"{score_set}_rmse": root_mean_squared_error(y_true, y_pred),
        f"{score_set}_r2": r2_score(y_true, y_pred),
        # f"{score_set}_pearson_r": pearsonr(y_test, y_pred)[0],
        # f"{score_set}_pearson_p_value": pearsonr(y_test, y_pred)[1],
        # f"{score_set}_spearman_r": spearmanr(y_test, y_pred)[0],
        # f"{score_set}_spearman_p_value": spearmanr(y_test, y_pred)[1],
        # f"{score_set}_kendall_r": kendalltau(y_test, y_pred)[0],
        # f"{score_set}_kendall_p_value": kendalltau(y_test, y_pred)[1],
    }


