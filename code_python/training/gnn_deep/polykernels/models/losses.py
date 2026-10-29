"""Loss function definitions and dispatcher for multi-label classification and regression."""

from enum import StrEnum

import draccus
import torch
import torch.nn as nn
import torch.nn.functional as F  # noqa: N812


class ClassificationLossType(StrEnum):
    """Multi-label classification loss functions."""

    BCE = "bce"
    BCE_SMOOTHED = "bce_smoothed"
    FOCAL = "focal"
    ASL = "asymmetric"


class RegressionLossType(StrEnum):
    """Regression loss functions."""

    MSE = "mse"
    HUBER = "huber"


class LossType(StrEnum):
    """All supported loss functions, accepted by `get_loss_fn`.

    Note:
        Hyperparameter searches must sample from `ClassificationLossType` or
        `RegressionLossType` rather than this union, so a classification study never
        samples a regression loss (or vice versa).
    """

    BCE = "bce"
    BCE_SMOOTHED = "bce_smoothed"
    FOCAL = "focal"
    ASL = "asymmetric"
    MSE = "mse"
    HUBER = "huber"


draccus.encode.register(ClassificationLossType, lambda x: x.value)
draccus.encode.register(RegressionLossType, lambda x: x.value)
draccus.encode.register(LossType, lambda x: x.value)


class LabelSmoothedBCEWithLogitsLoss(nn.Module):
    """BCEWithLogitsLoss with label smoothing."""

    def __init__(self, smoothing: float = 0.1, reduction: str = "mean") -> None:
        super().__init__()
        self.smoothing = smoothing
        self.reduction = reduction

    def forward(self, logits: torch.Tensor, targets: torch.Tensor) -> torch.Tensor:
        smoothed_targets = targets * (1.0 - self.smoothing) + 0.5 * self.smoothing
        return F.binary_cross_entropy_with_logits(
            logits, smoothed_targets, reduction=self.reduction
        )


class FocalLossWithLogits(nn.Module):
    """Focal Loss for multi-label classification using logits.

    Reference:
        Focal Loss for Dense Object Detection (https://arxiv.org/abs/1708.02002)
    """

    def __init__(
        self, alpha: float = 0.25, gamma: float = 2.0, reduction: str = "mean"
    ) -> None:
        super().__init__()
        self.alpha = alpha
        self.gamma = gamma
        self.reduction = reduction

    def forward(self, logits: torch.Tensor, targets: torch.Tensor) -> torch.Tensor:
        bce_loss = F.binary_cross_entropy_with_logits(logits, targets, reduction="none")
        probs = torch.sigmoid(logits)
        p_t = probs * targets + (1 - probs) * (1 - targets)
        focal_weight = (1 - p_t) ** self.gamma

        if self.alpha >= 0:
            alpha_t = self.alpha * targets + (1 - self.alpha) * (1 - targets)
            focal_loss = alpha_t * focal_weight * bce_loss
        else:
            focal_loss = focal_weight * bce_loss

        if self.reduction == "mean":
            return focal_loss.mean()
        if self.reduction == "sum":
            return focal_loss.sum()
        return focal_loss


class AsymmetricLossWithLogits(nn.Module):
    """Asymmetric Loss for multi-label classification using logits.

    Reference:
        Asymmetric Loss for Multi-Label Classification (https://arxiv.org/abs/2009.14119)
    """

    def __init__(
        self,
        gamma_pos: float = 0.0,
        gamma_neg: float = 2.0,
        margin: float = 0.05,
        reduction: str = "mean",
    ) -> None:
        super().__init__()
        self.gamma_pos = gamma_pos
        self.gamma_neg = gamma_neg
        self.margin = margin
        self.reduction = reduction

    def forward(self, logits: torch.Tensor, targets: torch.Tensor) -> torch.Tensor:
        probs = torch.sigmoid(logits)

        # Positive loss (stable using logsigmoid)
        bce_pos = -targets * F.logsigmoid(logits)
        if self.gamma_pos > 0:
            bce_pos *= (1 - probs) ** self.gamma_pos

        # Negative loss with margin clipping
        if self.margin > 0:
            probs_neg = (probs - self.margin).clamp(min=0.0)
        else:
            probs_neg = probs

        bce_neg = -(1 - targets) * torch.log((1 - probs_neg).clamp(min=1e-8))
        if self.gamma_neg > 0:
            bce_neg *= probs_neg ** self.gamma_neg

        asl_loss = bce_pos + bce_neg

        if self.reduction == "mean":
            return asl_loss.mean()
        if self.reduction == "sum":
            return asl_loss.sum()
        return asl_loss


def get_loss_fn(
    loss_type: LossType | ClassificationLossType | RegressionLossType | str,
) -> nn.Module:
    """Dispatcher to construct loss function by type."""
    loss_type = LossType(loss_type)
    if loss_type == LossType.BCE:
        return nn.BCEWithLogitsLoss()
    if loss_type == LossType.BCE_SMOOTHED:
        return LabelSmoothedBCEWithLogitsLoss()
    if loss_type == LossType.FOCAL:
        return FocalLossWithLogits()
    if loss_type == LossType.ASL:
        return AsymmetricLossWithLogits()
    if loss_type == LossType.MSE:
        return nn.MSELoss()
    if loss_type == LossType.HUBER:
        return nn.HuberLoss()
    raise ValueError(f"Unsupported loss type: {loss_type}")
