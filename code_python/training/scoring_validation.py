from typing import Callable, Union, Dict, List, Optional
from joblib import Parallel, delayed
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, kendalltau
from sklearn.pipeline import Pipeline
from sklearn.metrics._scorer import r2_scorer
from sklearn.model_selection import cross_validate
from collections import defaultdict
from utils import split_for_training
from sklearn.base import clone
from sklearn.metrics import (
    make_scorer,
    mean_absolute_error,
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
    compute_Cv,
    compute_sharpness
    )

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
        if (
            key.startswith("test_")
            and "lengthscale" not in key
            and "feature_importance" not in key
        )
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


def _fit_predict_score(
    estimator,
    model_type,
    X,
    y,
    train_idx,
    test_idx,
    scoring,
    return_ls: bool,
    UQ: bool,
    return_estimator: bool = False,
    return_tree_importances: bool = False,
):
    """
    Shared GP/tree worker that runs inside a parallel worker:
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

    est.fit(X_train, y_train)

    # All project GP/tree wrappers return the same prediction dictionary.
    y_result = est.predict(X_test, return_std=UQ)
    if not isinstance(y_result, dict) or "y_pred" not in y_result:
        raise TypeError(
            f"{type(est).__name__}.predict must return a dictionary containing "
            "'y_pred' and 'y_std'."
        )

    y_result["y_pred"] = np.asarray(y_result["y_pred"]).ravel()
    if y_result.get("y_std") is not None:
        y_result["y_std"] = np.asarray(y_result["y_std"]).ravel()

    results = {}
    fitted_regressor = (
        est.named_steps["regressor"]
        if isinstance(est, Pipeline)
        else est
    )

    if return_ls:
        results["lengthscale"] = fitted_regressor._get_lengthscale()

    if return_tree_importances:
        results["feature_importance_MDI"] = fitted_regressor._get_MDI()
        results["feature_importance_SHAP"] = fitted_regressor._get_SHAP()

    if return_estimator:
        results["estimator"] = est

    y_test = np.asarray(y_test).ravel()
    for name, scorer in scoring.items():
        results[name] = scorer(y_test, y_result["y_pred"])

    if UQ and y_result.get("y_std") is not None:
        UQ_scorers = {
            "ece": compute_ece,
            "RUSC": compute_RUSC,
            "cdf_ama": compute_cdf_ama,
            "cvpp_ama": compute_cvpp_ama,
            "nll": gaussian_nll,
            "Cv": compute_Cv,
            "sharpness": compute_sharpness
        }
        for name, uq_scorer in UQ_scorers.items():
            if name in {"Cv", "sharpness"}:
                results[name] = float(uq_scorer(y_result["y_std"]))
            else:
                results[name] = float(uq_scorer(y_test, y_result["y_pred"], y_result["y_std"]))

    return test_idx, y_result, results


def _fit_predict_score_iid(
    estimator,
    model_type,
    X,
    y,
    reference_test_idx,
    fold_idx,
    group_name,
    base_seed,
    scoring,
    return_ls: bool,
    UQ: bool,
    return_estimator: bool = False,
    return_tree_importances: bool = False,
):
    n_samples = len(y)
    test_size = len(reference_test_idx)

    if test_size <= 0 or test_size >= n_samples:
        raise ValueError(
            f"Invalid IID test size {test_size} for dataset "
            f"containing {n_samples} samples."
        )

    # Different RNG sequence for each group
    iid_seed = base_seed + fold_idx * 10_000

    rng = np.random.default_rng(iid_seed)
    permutation = rng.permutation(n_samples)

    iid_test_idx = np.sort(permutation[:test_size])
    iid_train_idx = np.sort(permutation[test_size:])

    test_idx, y_result, results = _fit_predict_score(
        estimator=estimator,
        model_type=model_type,
        X=X,
        y=y,
        train_idx=iid_train_idx,
        test_idx=iid_test_idx,
        scoring=scoring,
        return_ls=return_ls,
        UQ=UQ,
        return_estimator=return_estimator,
        return_tree_importances=return_tree_importances,
    )

    return group_name, test_idx, y_result, results


def cross_validate(
    estimator,
    model_type,
    X,
    y,
    cv,
    scoring,
    UQ,
    n_jobs,
    return_ls,
    return_estimator,
    return_tree_importances,
    cluster_group=None,
    cluster_validation_mode="ood",
):
    """
    Cross-validation supporting:

        cluster_validation_mode="ood"
            Leave-one-group-out only.

        cluster_validation_mode="iid"
            Repeated random IID splits matched to each group's test-set size.

        cluster_validation_mode="both"
            Run both OOD and matched IID validation.

    If cluster_group is None, standard CV behavior is used.
    """

    if cluster_validation_mode not in {"ood", "iid", "both"}:
        raise ValueError(
            "validation_mode must be one of: "
            "'ood', 'iid', or 'both'."
        )

    grouped_cv = cluster_group is not None
    n_samples = len(y)

    # ---------------------------------------------------------
    # Parallel configuration
    # ---------------------------------------------------------

    parallel_kwargs = {
        "n_jobs": n_jobs,
        "verbose": 0,
    }

    if "gp" in model_type.lower() or "mgk" in model_type.lower():
        parallel_kwargs["require"] = "sharedmem"
    else:
        parallel_kwargs["pre_dispatch"] = "all"

    # =========================================================
    # STANDARD NON-GROUPED CV
    # =========================================================

    if not grouped_cv:

        splits = cv.split(X, y)

        parallel_results = Parallel(**parallel_kwargs)(
            delayed(_fit_predict_score)(
                estimator,
                model_type,
                X,
                y,
                train_idx,
                test_idx,
                scoring,
                return_ls,
                UQ,
                return_estimator,
                return_tree_importances,
            )
            for train_idx, test_idx in splits
        )

        scores = defaultdict(list)

        predictions = {
            "y_pred": np.full(n_samples, np.nan),
            "y_std": np.full(n_samples, np.nan),
        }

        for test_idx, y_result, fold_scores in parallel_results:

            predictions["y_pred"][test_idx] = y_result["y_pred"]

            if y_result.get("y_std") is not None:
                predictions["y_std"][test_idx] = y_result["y_std"]

            for key, val in fold_scores.items():

                score_key = (
                    key
                    if key == "estimator"
                    else f"test_{key}"
                )

                scores[score_key].append(val)

    # =========================================================
    # GROUP-BASED VALIDATION
    # =========================================================

    else:

        cluster_group_array = np.asarray(cluster_group)

        raw_splits = list(
            cv.split(
                X,
                y,
                groups=cluster_group
            )
        )

        # -----------------------------------------------------
        # Build fold metadata once
        # -----------------------------------------------------

        folds = []

        for fold_idx, (train_idx, test_idx) in enumerate(raw_splits):

            unique_groups = np.unique(
                cluster_group_array[test_idx]
            )

            if len(unique_groups) != 1:
                raise ValueError(
                    "Expected exactly one held-out group per fold, "
                    f"but found {len(unique_groups)}: {unique_groups}"
                )

            group_name = unique_groups[0]

            if isinstance(group_name, np.generic):
                group_name = group_name.item()

            folds.append(
                {
                    "fold_idx": fold_idx,
                    "group_name": group_name,
                    "train_idx": train_idx,
                    "test_idx": test_idx,
                }
            )

        scores = {}
        predictions = {}

        # =========================================================
        # OOD VALIDATION
        # =====================================================

        if cluster_validation_mode in {"ood", "both"}:

            ood_results = Parallel(**parallel_kwargs)(
                delayed(_fit_predict_score)(
                    estimator,
                    model_type,
                    X,
                    y,
                    fold["train_idx"],
                    fold["test_idx"],
                    scoring,
                    return_ls,
                    UQ,
                    return_estimator,
                    return_tree_importances,
                )
                for fold in folds
            )

            ood_scores = defaultdict(dict)

            ood_predictions = {
                "y_pred": {},
                "y_std": {},
                "test_idx": {},
            }

            for fold, (
                test_idx,
                y_result,
                fold_scores,
            ) in zip(folds, ood_results):

                group_name = fold["group_name"]

                # Store predictions
                ood_predictions["y_pred"][group_name] = (
                    y_result["y_pred"]
                )

                ood_predictions["y_std"][group_name] = (
                    y_result["y_std"]
                    if y_result.get("y_std") is not None
                    else None
                )

                ood_predictions["test_idx"][group_name] = (
                    np.asarray(test_idx)
                )

                # Store scores
                for key, val in fold_scores.items():

                    score_key = (
                        key
                        if key == "estimator"
                        else f"test_{key}"
                    )

                    ood_scores[score_key][group_name] = val

            scores["ood"] = ood_scores
            predictions["ood"] = ood_predictions

        # =====================================================
        # MATCHED IID VALIDATION
        # =====================================================

        if cluster_validation_mode in {"iid", "both"}:

            IID_SEEDS = [17, 29, 43, 71, 97]

            iid_parallel_results = Parallel(**parallel_kwargs)(
                delayed(_fit_predict_score_iid)(
                    estimator=estimator,
                    model_type=model_type,
                    X=X,
                    y=y,
                    reference_test_idx=fold["test_idx"],
                    fold_idx=fold["fold_idx"],
                    group_name=fold["group_name"],
                    base_seed=base_seed,
                    scoring=scoring,
                    return_ls=return_ls,
                    UQ=UQ,
                    return_estimator=return_estimator,
                    return_tree_importances=return_tree_importances,
                )
                for fold in folds
                for base_seed in IID_SEEDS
            )

            iid_scores = defaultdict(dict)

            iid_predictions = {
                "y_pred": {},
                "y_std": {},
                "test_idx": {},
            }

            for (
                group_name,
                test_idx,
                y_result,
                fold_scores,
            ) in iid_parallel_results:

                if group_name not in iid_predictions["y_pred"]:
                    iid_predictions["y_pred"][group_name] = []
                    iid_predictions["y_std"][group_name] = []
                    iid_predictions["test_idx"][group_name] = []

                # Store predictions
                iid_predictions["y_pred"][group_name].append(
                    y_result["y_pred"]
                )

                iid_predictions["y_std"][group_name].append(
                    y_result["y_std"]
                    if y_result.get("y_std") is not None
                    else None
                )

                iid_predictions["test_idx"][group_name].append(
                    np.asarray(test_idx)
                )

                # Store scores
                for key, val in fold_scores.items():

                    score_key = (
                        key
                        if key == "estimator"
                        else f"test_{key}"
                    )

                    if group_name not in iid_scores[score_key]:
                        iid_scores[score_key][group_name] = []

                    iid_scores[score_key][group_name].append(val)
            scores["iid"] = iid_scores
            predictions["iid"] = iid_predictions

    return scores, predictions


def cross_validate_regressor(
    regressor, 
    model_type:str, 
    X, y,
    cv, 
    UQ:bool=False,
    return_ls:bool=False,
    return_estimator:bool=False,
    return_tree_importances:bool=False,
    n_jobs:int=1,     
    cluster_group: Optional[str]=None,
    cluster_validation_mode: str="ood"
    ) -> tuple[dict[str, float], dict[str, np.ndarray]]:


        scorers = {
            "rmse": root_mean_squared_error,
            "mae": mean_absolute_error,
            "r2": r2_score,
            }
        
        score, predictions = cross_validate(
            regressor,
            model_type,
            X,
            y,
            cv=cv,
            scoring=scorers,
            n_jobs=n_jobs,
            return_ls=return_ls,
            UQ=UQ,
            return_estimator=return_estimator,
            return_tree_importances=return_tree_importances,
            cluster_group=cluster_group,
            cluster_validation_mode=cluster_validation_mode
            )
        return score, predictions
