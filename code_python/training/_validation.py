from sklearn.utils.parallel import Parallel, delayed
import numbers
import time
# import warnings
# from collections import Counter
# from contextlib import suppress
# from functools import partial
# from numbers import Real
# from traceback import format_exc
from sklearn.base import clone

import numpy as np
import scipy.sparse as sp
from joblib import logger

from sklearn.base import clone, is_classifier
from sklearn.exceptions import FitFailedWarning, UnsetMetadataPassedError
# from ..metrics import check_scoring, get_scorer_names
# from sklearn.metrics._scorer import _MultimetricScorer
from sklearn.utils import Bunch, _safe_indexing, check_random_state, indexable
# from ..utils._array_api import device, get_namespace
# from sklearn.utils._param_validation import (
#     HasMethods,
#     Integral,
#     Interval,
#     StrOptions,
#     validate_params,
# )
from sklearn.utils.metadata_routing import (
    MetadataRouter,
    MethodMapping,
    _routing_enabled,
    process_routing,
)
from sklearn.utils.metaestimators import _safe_split
from sklearn.model_selection._validation import (_check_method_params, 
                                      _check_params_groups_deprecation,
                                      )
from sklearn.model_selection._split import check_cv
import numbers




def multioutput_cross_validate(estimator, X, y,
                               scorers,
                               cv,
                               n_jobs=-1,
                               verbose=0,
                               pre_dispatch="2*n_jobs",
                               return_estimator=False,
                               return_train_score=False,
                               fit_params=None,
                               params=None,  
                               groups=None,  
                               ):
    """
    Custom cross-validation function supporting multiple custom scorers.

    Parameters:
    - estimator: The regressor or classifier to be evaluated.
    - X: Features matrix.
    - y: Target matrix (multioutput allowed).
    - scorers: Dictionary of scoring functions (e.g., {"r2": r2_scorer_func, "rmse": rmse_scorer_func}).
    - cv: Cross-validation strategy (e.g., KFold instance).
    - n_jobs: Number of jobs for parallel processing.
    - verbose: Verbosity level.
    - pre_dispatch: Controls the number of jobs dispatched during parallel execution.

    Returns:
    - results: Dictionary containing scores for each fold and each scorer.
    """
    params = _check_params_groups_deprecation(fit_params, params, groups)
    cv = check_cv(cv, y, classifier=is_classifier(estimator))

    X, y = indexable(X, y)
    if _routing_enabled():
        # For estimators, a MetadataRouter is created in get_metadata_routing
        # methods. For these router methods, we create the router to use
        # `process_routing` on it.
        router = (
            MetadataRouter(owner="cross_validate")
            .add(
                splitter=cv,
                method_mapping=MethodMapping().add(caller="fit", callee="split"),
            )
            .add(
                estimator=estimator,
                # TODO(SLEP6): also pass metadata to the predict method for
                # scoring?
                method_mapping=MethodMapping().add(caller="fit", callee="fit"),
            )
            .add(
                scorer=scorers,
                method_mapping=MethodMapping().add(caller="fit", callee="score"),
            )
        )
        try:
            routed_params = process_routing(router, "fit", **params)
        except UnsetMetadataPassedError as e:
            # The default exception would mention `fit` since in the above
            # `process_routing` code, we pass `fit` as the caller. However,
            # the user is not calling `fit` directly, so we change the message
            # to make it more suitable for this case.
            unrequested_params = sorted(e.unrequested_params)
            raise UnsetMetadataPassedError(
                message=(
                    f"{unrequested_params} are passed to cross validation but are not"
                    " explicitly set as requested or not requested for cross_validate's"
                    f" estimator: {estimator.__class__.__name__}. Call"
                    " `.set_fit_request({{metadata}}=True)` on the estimator for"
                    f" each metadata in {unrequested_params} that you"
                    " want to use and `metadata=False` for not using it. See the"
                    " Metadata Routing User guide"
                    " <https://scikit-learn.org/stable/metadata_routing.html> for more"
                    " information."
                ),
                unrequested_params=e.unrequested_params,
                routed_params=e.routed_params,
            )
    else:
        routed_params = Bunch()
        routed_params.splitter = Bunch(split={"groups": groups})
        routed_params.estimator = Bunch(fit=params)
        routed_params.scorer = Bunch(score={})

    # print('problem')
    indices = cv.split(X, y, **routed_params.splitter.split)
    # print(indices)
    # indices = list(cv.split(X))
    # print(indices)
    
    parallel = Parallel(n_jobs=n_jobs, verbose=verbose, pre_dispatch=pre_dispatch)
    results = parallel(
        delayed(_fit_and_score)(
            clone(estimator),
            X,
            y,
            scorers,
            train,
            test)
            for train, test in indices
    )
    # print(results)
    results = _aggregate_score_dicts(results)
    res = {}
    res["fit_time"] = results["fit_time"]
    res["score_time"] = results["score_time"]
    if return_estimator:
        res["estimator"] = results["estimator"]
    test_scores_dict = _normalize_score_results(results["test_scores"])
    if return_train_score:
        train_scores_dict = _normalize_score_results(results["train_scores"])
    # print(test_scores_dict)
    for name in test_scores_dict:
      res[f"test_{name}"] = test_scores_dict[name]
      # res["test_{}"] = test_scores_dict[name]

      if return_train_score:
          key = f"train_{name}"
          res[key] = train_scores_dict[name]
        
    return res


def _fit_and_score(estimator,
                  X,
                  y,
                  scorer,
                  train,
                  test,
                  parameters=None,
                  fit_params=None,
                  return_times=True,
                  return_parameters=False,
                  return_estimator=False,
                  return_train_score=False,
):
        # Clone the estimator for isolation
        fit_params = fit_params if fit_params is not None else {}
        fit_params = _check_method_params(X, params=fit_params, indices=train)
        if parameters is not None:
        # here we clone the parameters, since sometimes the parameters
        # themselves might be estimators, e.g. when we search over different
        # estimators in a pipeline.
        # ref: https://github.com/scikit-learn/scikit-learn/pull/26786
          estimator = estimator.set_params(**clone(parameters, safe=False))

        start_time = time.time()
        result = {}
        
        # Split data
        X_train, y_train = _safe_split(estimator, X, y, train)
        X_test, y_test = _safe_split(estimator, X, y, test, train)
        
        # Fit the model
        estimator.fit(X_train, y_train)
        fit_time = time.time() - start_time
        score_time = 0.0
        # Predict
        test_scores = _score(estimator, X_test,y_test, scorer)
        score_time = time.time() - start_time - fit_time
        if return_train_score:
            train_scores = _score(estimator, X_train, y_train, scorer)
        


        result["test_scores"] = test_scores
        if return_train_score:
          result["train_scores"] = train_scores
        if return_times:
            result["fit_time"] = fit_time
            result["score_time"] = score_time
        if return_parameters:
            result["parameters"] = parameters
        if return_estimator:
            result["estimator"] = estimator
        return result


def _score(estimator, X_test, y_test, score_space):
        y_pred = estimator.predict(X_test)
        
        # Compute scores for each scorer
        fold_scores = {}
        for name, sc in score_space.items():
            fold_scores[name] = sc(y_test, y_pred)
        
        return fold_scores


def _aggregate_score_dicts(scores):
    return {
        key: (
            np.asarray([score[key] for score in scores])
            if isinstance(scores[0][key], numbers.Number)
            else [score[key] for score in scores]
        )
        for key in scores[0]
    }


def _normalize_score_results(scores, scaler_score_key="score"):
    """Creates a scoring dictionary based on the type of `scores`"""
    if isinstance(scores[0], dict):
        # multimetric scoring
        return _aggregate_score_dicts(scores)
    # scaler
    return {scaler_score_key: scores}





############## test
# from sklearn.ensemble import RandomForestRegressor
# from sklearn.multioutput import MultiOutputRegressor
# from sklearn.datasets import make_regression
# from sklearn.model_selection import KFold
# from sklearn.metrics import r2_score, mean_squared_error,root_mean_squared_error
# import numpy as np
# import time

# # Define custom scorers
# def r2_scorer_func(y_true, y_pred):
#     return r2_score(y_true, y_pred, multioutput="raw_values")

# def rmse_scorer_func(y_true, y_pred):
#     return mean_squared_error(y_true, y_pred, multioutput="raw_values")

# def mse_scorer_func(y_true, y_pred):
#     return root_mean_squared_error(y_true, y_pred, multioutput="raw_values") 
# scorers = {
#     "r2": r2_scorer_func,
#     "rmse": mse_scorer_func,
#     "mae": mse_scorer_func
# }

# # Generate synthetic data
# X, y = make_regression(n_samples=100, n_features=10, n_targets=3, noise=0.1, random_state=25)

# # Define the regressor and cross-validation strategy
# regressor = MultiOutputRegressor(RandomForestRegressor(random_state=42), n_jobs=-1)
# cv = KFold(n_splits=5, shuffle=True, random_state=42)

# # Perform custom cross-validation
# score = multioutput_cross_validate(
#     estimator=regressor,
#     X=X,
#     y=y,
#     scorers=scorers,
#     cv=cv,
#     n_jobs=-1,
#     verbose=0
# )

# print(score)
# for i in ["r2","rmse"]:
#   score[f'{i}_average'] =  np.array(score[f'test_{i}']).mean(axis=0)

# print(score)
