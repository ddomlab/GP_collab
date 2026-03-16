#!/usr/bin/env python
# -*- coding: utf-8 -*-
from typing import List, Literal
import warnings
import numpy as np
import scipy
from sklearn.metrics import (
    mean_squared_error,
    root_mean_squared_error,
    mean_absolute_error,
    r2_score,
    accuracy_score,
    balanced_accuracy_score,
    roc_auc_score,
    precision_score,
    recall_score,
    f1_score,
    matthews_corrcoef
)


Metric = Literal['roc_auc', 'accuracy', 'balanced_accuracy', 'precision', 'recall', 'f1_score', 'mcc',
                 'rmse', 'mae', 'mse', 'r2', 'max', 'spearman', 'kendall', 'pearson']
AVAILABLE_METRICS_REGRESSION = ['rmse', 'mae', 'mse', 'r2', 'max', 'spearman', 'pearson', 'kendall']
AVAILABLE_METRICS_BINARY = ['roc_auc', 'accuracy', 'balanced_accuracy', 'precision', 'recall', 'f1_score', 'mcc']
AVAILABLE_METRICS_MULTICLASS = ['roc_auc', 'accuracy', 'balanced_accuracy', 'precision', 'recall', 'f1_score', 'mcc']


def metric_regression(y: List[float], y_pred: List[float], metric: Metric) -> float:
    """
    Calculate performance metrics for regression tasks.
    
    Args:
        y (List[float]): Reference list of target values.
        y_pred (List[float]): List of predicted values.
        metric (Metric): Metric to calculate.
    
    Returns:
        float: Calculated metric value.
    """
    if metric == 'rmse':
        return root_mean_squared_error(y, y_pred)
    elif metric == 'mae':
        return mean_absolute_error(y, y_pred)
    elif metric == 'mse':
        return mean_squared_error(y, y_pred)
    elif metric == 'r2':
        return r2_score(y, y_pred)
    elif metric == 'max':
        return np.max(abs(y - y_pred))
    elif metric == 'spearman':
        return scipy.stats.spearmanr(y, y_pred)[0]
    elif metric == 'pearson':
        return scipy.stats.pearsonr(y, y_pred)[0]
    elif metric == 'kendall':
        return scipy.stats.kendalltau(y, y_pred)[0]
    else:
        raise RuntimeError(f'Unsupported metrics {metric}')


def metric_binary(y: List[int], y_pred: List[float], metric: Metric) -> float:
    """
    Calculate performance metrics for binary classification tasks.
    
    Args:
        y (List[int]): Reference list of target values.
        y_pred (List[float]): List of predicted probabilities.
        metric (Metric): Metric to calculate.
    
    Returns:
        float: Calculated metric value.
    """
    sety = set(y)
    if sety == {0} or sety == {1}:
        warnings.warn('Only one class present in target values.')
    elif sety != {0, 1}:
        raise ValueError('Target values must be 0 or 1 for binary classification tasks.')
    y_pred_c = p2c(y, y_pred)
    if metric == 'roc_auc':
        return roc_auc_score(y, y_pred)
    elif metric == 'accuracy':
        return accuracy_score(y, y_pred_c)
    elif metric == 'balanced_accuracy':
        return balanced_accuracy_score(y, y_pred_c)
    elif metric == 'precision':
        return precision_score(y, y_pred_c)
    elif metric == 'recall':
        return recall_score(y, y_pred_c)
    elif metric == 'f1_score':
        return f1_score(y, y_pred_c)
    elif metric == 'mcc':
        return matthews_corrcoef(y, y_pred_c)
    else:
        raise ValueError(f'Unsupported metrics {metric}')


def p2c(y: List[float], y_pred: List[float]) -> List[float]:
    """
    Maps continuous predictions to nearest valid discrete values from reference set.
    
    Args:
        y (List[float]): Reference list containing valid target values. Used to 
            determine the set of allowed output values.
        y_pred (List[float]): List of predicted values to be mapped. Can contain 
            any continuous values.
    
    Returns:
        List[float]: Predicted values mapped to nearest values from reference set.
    
    Example:
        >>> y = [1.0, 2.0, 1.0, 3.0]
        >>> y_pred = [1.1, 2.3, 0.9, 2.7]
        >>> p2v(y, y_pred)
        [1.0, 2.0, 1.0, 3.0]
    """
    avail_values = list(set(y))
    y_pred = [min(avail_values, key=lambda x: abs(x - v)) for v in y_pred]
    return y_pred


def metric_multiclass(y: List[int], y_pred: List[float], metric: Metric) -> float:
    raise NotImplementedError('Multiclass classification metrics are not implemented yet.')
