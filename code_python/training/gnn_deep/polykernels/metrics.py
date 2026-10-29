"""Evaluation metrics for multilabel classification."""

from collections.abc import Callable, Sequence

import numpy as np
import torch
import torchmetrics.classification as tmc

from polykernels.utils import cast_random_state, round3


def _prepare_thresholded_probs(
    y_prob: torch.Tensor,
    threshold: float | Sequence[float],
) -> tuple[torch.Tensor, float]:
    if not isinstance(threshold, (int, float)):
        thresh_tensor = torch.as_tensor(threshold, dtype=y_prob.dtype, device=y_prob.device)
        return (y_prob >= thresh_tensor).to(y_prob.dtype), 0.5
    return y_prob, threshold


def auroc_score(y_true: torch.Tensor, y_prob: torch.Tensor) -> float:
    """Compute the macro-averaged Multilabel Area Under the Receiver Operating Characteristic (AUROC)."""
    n = y_true.shape[1]
    val = (
        tmc.MultilabelAUROC(num_labels=n, average="macro")
        .to(y_true.device)(y_prob, y_true.long())
        .item()
    )
    return round3(val)


def auprc_score(y_true: torch.Tensor, y_prob: torch.Tensor) -> float:
    """Compute the macro-averaged Multilabel Area Under the Precision-Recall Curve (AUPRC)."""
    n = y_true.shape[1]
    val = (
        tmc.MultilabelAveragePrecision(num_labels=n, average="macro")
        .to(y_true.device)(y_prob, y_true.long())
        .item()
    )
    return round3(val)


def f1_score(
    y_true: torch.Tensor,
    y_prob: torch.Tensor,
    threshold: float | Sequence[float] = 0.5,
) -> float:
    """Compute the macro-averaged Multilabel F1 Score with class-specific or uniform thresholds."""
    n = y_true.shape[1]
    y_prob_thresholded, eval_threshold = _prepare_thresholded_probs(y_prob, threshold)
    val = (
        tmc.MultilabelF1Score(num_labels=n, average="macro", threshold=eval_threshold)
        .to(y_true.device)(y_prob_thresholded, y_true.long())
        .item()
    )
    return round3(val)


def precision_score(
    y_true: torch.Tensor,
    y_prob: torch.Tensor,
    threshold: float | Sequence[float] = 0.5,
) -> float:
    """Compute the macro-averaged Multilabel Precision score."""
    n = y_true.shape[1]
    y_prob_thresholded, eval_threshold = _prepare_thresholded_probs(y_prob, threshold)
    val = (
        tmc.MultilabelPrecision(num_labels=n, average="macro", threshold=eval_threshold)
        .to(y_true.device)(y_prob_thresholded, y_true.long())
        .item()
    )
    return round3(val)


def recall_score(
    y_true: torch.Tensor,
    y_prob: torch.Tensor,
    threshold: float | Sequence[float] = 0.5,
) -> float:
    """Compute the macro-averaged Multilabel Recall score."""
    n = y_true.shape[1]
    y_prob_thresholded, eval_threshold = _prepare_thresholded_probs(y_prob, threshold)
    val = (
        tmc.MultilabelRecall(num_labels=n, average="macro", threshold=eval_threshold)
        .to(y_true.device)(y_prob_thresholded, y_true.long())
        .item()
    )
    return round3(val)


def evaluate(
    y_true: torch.Tensor,
    y_prob: torch.Tensor,
    threshold: float | Sequence[float] = 0.5,
    **_kwargs,
) -> dict[str, float]:
    """Calculate and return a dictionary of standard multilabel classification metrics (AUROC, AUPRC, F1, Precision, Recall)."""
    return {
        "auroc": auroc_score(y_true, y_prob),
        "auprc": auprc_score(y_true, y_prob),
        "f1": f1_score(y_true, y_prob, threshold=threshold),
        "precision": precision_score(y_true, y_prob, threshold=threshold),
        "recall": recall_score(y_true, y_prob, threshold=threshold),
    }


def bootstrap_scores(
    y_true: torch.Tensor,
    y_prob: torch.Tensor,
    metric_fn: Callable[[torch.Tensor, torch.Tensor], float],
    n_bootstraps: int = 1000,
    random_state: int | np.random.RandomState | None = None,
    *,
    drop_zero: bool = False,
    drop_nan: bool = True,
) -> np.ndarray:
    """Generate bootstrap scores for a metric, optionally dropping zeros/NaNs."""
    rng = cast_random_state(random_state)
    n_samples = y_true.shape[0]

    boot_metrics = []
    for _ in range(n_bootstraps):
        boot_idx = rng.randint(0, n_samples, size=n_samples)
        val = metric_fn(y_true[boot_idx], y_prob[boot_idx])
        boot_metrics.append(val)

    scores = np.asarray(boot_metrics)
    if drop_nan:
        scores = scores[~np.isnan(scores)]
    if drop_zero:
        scores = scores[scores != 0.0]
    return scores


def bootstrap_ci(
    y_true: torch.Tensor,
    y_prob: torch.Tensor,
    metric_fn: Callable[[torch.Tensor, torch.Tensor], float],
    n_bootstraps: int = 1000,
    ci_level: float = 0.95,
    random_state: int | np.random.RandomState | None = None,
    *,
    drop_zero: bool = False,
    drop_nan: bool = True,
) -> tuple[float, float, float]:
    """Compute the bootstrap mean and percentile confidence interval for a metric."""
    boot_metrics = bootstrap_scores(
        y_true=y_true,
        y_prob=y_prob,
        metric_fn=metric_fn,
        n_bootstraps=n_bootstraps,
        random_state=random_state,
        drop_zero=drop_zero,
        drop_nan=drop_nan,
    )
    if len(boot_metrics) == 0:
        return 0.0, 0.0, 0.0

    alpha = 1.0 - ci_level
    lower_pct = 100 * (alpha / 2.0)
    upper_pct = 100 * (1.0 - alpha / 2.0)

    boot_mean = round3(np.mean(boot_metrics))
    ci_lower = round3(np.percentile(boot_metrics, lower_pct))
    ci_upper = round3(np.percentile(boot_metrics, upper_pct))

    return boot_mean, ci_lower, ci_upper



