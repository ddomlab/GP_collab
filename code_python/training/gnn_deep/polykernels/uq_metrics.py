"""Regression accuracy and predictive-uncertainty-calibration metrics.

All functions take 1-D arrays for a *single* target column. Multi-output problems
(e.g. dataset3's 6 membrane targets, dataset4's 2 pervaporation targets) call these
once per target dimension in the training script — nothing here aggregates across
columns, since pooling differently-scaled targets into one calibration number would be
close to meaningless.

The `gaussian_nll`/`compute_ece`/`compute_ama`/`compute_rusc`/`compute_cv`/
`compute_sharpness` metrics require a predictive std and are only meaningful for the
GP-headed models; the plain point-estimate model only gets `rmse_score`/`mae_score`/
`r2_score`.
"""

import numpy as np
from scipy import stats
from scipy.integrate import simpson
from sklearn.metrics import r2_score as _sklearn_r2_score

# Floor on predicted std, so a collapsed/degenerate posterior cannot divide by zero.
_STD_FLOOR = 1e-8


def _safe_std(std: np.ndarray) -> np.ndarray:
    """Clamp predicted std away from zero."""
    return np.maximum(np.asarray(std, dtype=float), _STD_FLOOR)


def rmse_score(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Root mean squared error."""
    return float(np.sqrt(np.mean((np.asarray(y_true) - np.asarray(y_pred)) ** 2)))


def mae_score(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Mean absolute error."""
    return float(np.mean(np.abs(np.asarray(y_true) - np.asarray(y_pred))))


def r2_score(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Coefficient of determination."""
    return float(_sklearn_r2_score(y_true, y_pred))


def spearman_score(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Spearman rank correlation between truth and prediction.

    A rank-accuracy counterpart to `r2_score`, for reporting alongside collaborators who
    score with it. Distinct from `compute_rusc`, which correlates |residual| against the
    predicted std and so measures whether uncertainty tracks error rather than accuracy.

    Returns NaN when either input is constant, since the correlation is undefined there.
    """
    return float(stats.spearmanr(np.asarray(y_true), np.asarray(y_pred)).statistic)


def gaussian_nll(y_true: np.ndarray, mean: np.ndarray, std: np.ndarray) -> float:
    """Gaussian negative log-likelihood of y_true under N(mean, std^2)."""
    var = _safe_std(std) ** 2
    return float(
        np.mean(0.5 * np.log(2 * np.pi * var) + (np.asarray(y_true) - mean) ** 2 / (2 * var))
    )


def compute_cvpp(
    y_true: np.ndarray, mean: np.ndarray, std: np.ndarray, n_points: int = 99
) -> tuple[np.ndarray, np.ndarray]:
    """Coverage-vs-confidence (CVPP) curve.

    For each nominal central-interval confidence level q, computes the empirical
    fraction of samples actually falling inside the predicted N(mean, std^2) interval
    at that confidence level. A perfectly calibrated model has empirical == nominal.

    Returns:
        qs: nominal confidence levels, linspace(0, 1, n_points + 2)[1:-1].
        empirical_coverage: matching empirical coverage at each q.
    """
    qs = np.linspace(0.0, 1.0, n_points + 2)[1:-1]
    abs_z = np.abs((np.asarray(y_true) - mean) / _safe_std(std))
    z_thresholds = stats.norm.ppf((1.0 + qs) / 2.0)
    empirical_coverage = np.array([np.mean(abs_z <= z) for z in z_thresholds])
    return qs, empirical_coverage


def compute_ama(y_true: np.ndarray, mean: np.ndarray, std: np.ndarray) -> float:
    """Absolute Miscalibration Area: integral of |empirical - nominal| coverage via Simpson's rule."""
    qs, empirical_coverage = compute_cvpp(y_true, mean, std)
    return float(simpson(np.abs(empirical_coverage - qs), x=qs))


def compute_ece(
    y_true: np.ndarray, mean: np.ndarray, std: np.ndarray, n_bins: int = 10
) -> float:
    """Binned Expected Calibration Error: mean |empirical - nominal| coverage across n_bins levels."""
    qs, empirical_coverage = compute_cvpp(y_true, mean, std, n_points=n_bins)
    return float(np.mean(np.abs(empirical_coverage - qs)))


def compute_rusc(y_true: np.ndarray, mean: np.ndarray, std: np.ndarray) -> float:
    """Spearman correlation between |residual| and predicted std.

    Returns NaN when the predicted std is constant, since the rank correlation is then
    undefined -- that happens for a degenerate/untrained posterior and is reported as
    NaN rather than silently coerced to 0.
    """
    abs_residual = np.abs(np.asarray(y_true) - mean)
    std = np.asarray(std, dtype=float)
    if np.ptp(std) == 0.0 or np.ptp(abs_residual) == 0.0:
        return float("nan")
    correlation, _ = stats.spearmanr(abs_residual, std)
    return float(correlation)


def compute_cv(std: np.ndarray) -> float:
    """Coefficient of variation (std/mean) of the predicted stds themselves."""
    std = _safe_std(std)
    return float(np.std(std) / np.mean(std))


def compute_sharpness(std: np.ndarray) -> float:
    """Mean predicted std -- lower is sharper (more confident) on average."""
    return float(np.mean(np.asarray(std, dtype=float)))
