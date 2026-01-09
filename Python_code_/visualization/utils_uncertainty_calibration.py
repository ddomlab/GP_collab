import numpy as np
from scipy.integrate import simpson
from typing import Callable, Optional, Tuple, OrderedDict, Dict, Union
from scipy.stats import pearsonr, spearmanr, norm
from sklearn.metrics import root_mean_squared_error


def compute_cvpp(y_true: np.ndarray, 
                y_pred: np.ndarray,
                y_err: np.ndarray,
                step: float = 0.05
                ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes the CVPP diagram values: observed coverage C(q) at various nominal confidence levels q.

    Parameters:
        y_true (np.ndarray): True target values.
        y_pred (np.ndarray): Predicted means.
        y_err (np.ndarray): Predicted standard deviations.
        step (float): Step size for confidence levels (default: 0.05).

    Returns:
        Tuple[np.ndarray, np.ndarray]: (qs, Cqs), where
            - qs: array of nominal confidence levels
            - Cqs: corresponding empirical coverage values
    """
    qs = np.arange(0, 1.0 + step, step)
    Cqs = np.empty_like(qs)

    for i, q in enumerate(qs):
        z = norm.ppf((1.0 + q) / 2.0)
        standardized_error  = np.abs((y_pred - y_true) / y_err)
        Cqs[i] = np.mean(standardized_error  < z)

    return qs, Cqs


def compute_ama(y_true: np.ndarray, 
                y_pred: np.ndarray,
                y_err: np.ndarray,
                step: float = 0.05) -> float:
    
    """
    Computes the Absolute Miscalibration Area (AMA) using CVPP.

    Returns:
        float: The AMA score.
    """
    qs, Cqs = compute_cvpp(y_true, y_pred, y_err, step=step)
    ama = simpson(np.abs(Cqs - qs), qs)
    return ama


def compute_residual_error_cal(y_true: np.ndarray, 
                                y_pred: np.ndarray,
                                y_err: np.ndarray,
                                ):

    res = np.abs(y_true-y_pred)
    correlation = spearmanr(res, y_err)[0]
    return correlation


def gaussian_nll(
                y_true,
                y_pred,
                y_err,
                eps:float=1e-6, 
                reduce='mean'
                ):
    """
    Computes Negative Log-Likelihood (NLL) using scipy.stats.norm.logpdf.

    Parameters:
    -----------
    y_true : np.ndarray
        True values, shape (n_samples,)
    y_pred_mean : np.ndarray
        Predicted means, shape (n_samples,)
    y_pred_std : np.ndarray
        Predicted standard deviations, shape (n_samples,)
    reduce : str
        One of 'mean', 'sum', or 'none' to control the reduction of NLL

    Returns:
    --------
    nll : float or np.ndarray
        The computed NLL value(s)
    """
    eps = 1e-6 
    y_err = np.clip(y_err, eps, None)

    log_probs = norm.logpdf(y_true, loc=y_pred, scale=y_err)
    nll = -log_probs  

    if reduce == 'mean':
        return np.mean(nll)
    elif reduce == 'sum':
        return np.sum(nll)
    elif reduce == 'none':
        return nll
    else:
        raise ValueError("reduce must be one of: 'mean', 'sum', or 'none'")



def compute_cv(stdevs):
    """
    Calculates the coefficient of variation (Cv) of predicted uncertainties.

    Parameters:
    -----------
    stdevs : np.ndarray
        Array of predicted standard deviations (uncertainty values), shape (n_samples,)

    Returns:
    --------
    cv : float
        Coefficient of variation of the uncertainty estimates
    """
    stdevs = np.asarray(stdevs)
    eps = 1e-8  # prevent division by zero
    mean_std = np.mean(stdevs)
    std_std = np.std(stdevs, ddof=1)  # unbiased estimator
    cv = std_std / (mean_std + eps)
    return cv


def compute_all_uncertainty_metrics(
    y_true: np.ndarray, 
    y_pred: np.ndarray, 
    y_err: np.ndarray, 
    step: float = 0.05, 
    method: Optional[str] = None
) -> Union[float, Dict[str, float]]:
    
    def sharpness(err): return np.sqrt(np.mean(err**2))

    metrics = {
        'NLL': lambda: gaussian_nll(y_true, y_pred, y_err, reduce='mean'),
        'Sharpness': lambda: sharpness(y_err),
        'Cv': lambda: compute_cv(y_err),
        'RUSC': lambda: compute_residual_error_cal(y_true, y_pred, y_err),
        'AMA': lambda: compute_ama(y_true, y_pred, y_err, step=step)
    }

    if method:
        return metrics[method]()
    
    return {name: func() for name, func in metrics.items()}
