import numpy as np
from scipy.integrate import simpson
from typing import Callable, Optional, Tuple, OrderedDict, Dict, Union
from scipy.stats import pearsonr, spearmanr, norm
from sklearn.metrics import root_mean_squared_error




def compute_ece(y_true, y_pred, y_std, n_bins=10):
    # Binned Expected Calibration Error (ECE)
    cdf_vals = norm.cdf(y_true, loc=y_pred, scale=y_std)
    ece = 0.0
    for i in range(1, n_bins + 1):
        p_expected = i / n_bins
        p_observed = np.mean(cdf_vals <= p_expected)
        ece += np.abs(p_expected - p_observed)
    return ece / n_bins


def compute_pit_values(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    y_std: np.ndarray,
    eps: float = 1e-8,
) -> np.ndarray:
    y_true = np.asarray(y_true, dtype=float).ravel()
    y_pred = np.asarray(y_pred, dtype=float).ravel()
    y_std = np.asarray(y_std, dtype=float).ravel()
    z_scores = (y_true - y_pred) / (y_std + eps)
    return norm.cdf(z_scores)


def compute_pit_cdf_curve(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    y_std: np.ndarray,
    eps: float = 1e-8,
) -> Tuple[np.ndarray, np.ndarray]:

    cdf_vals = compute_pit_values(y_true, y_pred, y_std, eps=eps)
    cdf_vals = np.asarray(cdf_vals, dtype=float).ravel()

    sorted_pit = np.sort(cdf_vals)
    empirical_cdf = np.arange(1, len(sorted_pit) + 1) / len(sorted_pit)

    return sorted_pit, empirical_cdf


def compute_cdf_ama(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    y_std: np.ndarray,
    eps: float = 1e-8,
) -> float:
    sorted_pit, empirical_cdf = compute_pit_cdf_curve(
        y_true, y_pred, y_std, eps=eps
    )

    return float(simpson(np.abs(empirical_cdf - sorted_pit), x=sorted_pit))



def compute_cvpp(y_true: np.ndarray, 
                y_pred: np.ndarray,
                y_std: np.ndarray,
                step: float = 0.05
                ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes the CVPP diagram values: observed coverage C(q) at various nominal confidence levels q.

    Parameters:
        y_true (np.ndarray): True target values.
        y_pred (np.ndarray): Predicted means.
        y_std (np.ndarray): Predicted standard deviations.
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
        standardized_error  = np.abs((y_pred - y_true) / y_std)
        Cqs[i] = np.mean(standardized_error  < z)

    return qs, Cqs


def compute_cvpp_ama(y_true: np.ndarray, 
                y_pred: np.ndarray,
                y_std: np.ndarray,
                step: float = 0.05) -> float:
    
    """
    Computes the Absolute Miscalibration Area (AMA) using CVPP.

    Returns:
        float: The AMA score.
    """
    qs, Cqs = compute_cvpp(y_test, y_pred, y_std, step=step)
    ama = simpson(np.abs(Cqs - qs), qs)
    return ama


def compute_residual_error_cal(y_true: np.ndarray, 
                                y_pred: np.ndarray,
                                y_std: np.ndarray,
                                ):

    res = np.abs(y_true-y_pred)
    correlation = spearmanr(res, y_std)[0]
    return correlation


def gaussian_nll(
                y_true,
                y_pred,
                y_std,
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
    y_std = np.clip(y_std, eps, None)

    log_probs = norm.logpdf(y_true, loc=y_pred, scale=y_std)
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
