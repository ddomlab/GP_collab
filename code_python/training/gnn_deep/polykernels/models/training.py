from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from safetensors.torch import save_file
import torch
import torch.nn as nn
import torch.optim.lr_scheduler as lr_scheduler
from torch.optim.swa_utils import AveragedModel, get_ema_multi_avg_fn
from torch.utils.data import DataLoader
from tqdm import tqdm

from polykernels.metrics import auroc_score
from polykernels.models.losses import LossType


@dataclass
class OptimizerConfig:
    """Optimizer and training-loop hyperparameters."""

    type: str = "adamw"
    lr: float = 0.005
    weight_decay: float = 6e-5
    epochs: int = 200
    batch_size: int = 256
    warmup: int = 20
    loss_type: LossType = LossType.BCE


def resolve_ema_decay(total_steps: int, requested_decay: float = 0.999) -> float:
    """Cap the EMA decay so the averaged weights actually converge within `total_steps`.

    A decay of d leaves d**total_steps of the weight on the random initialization. On
    small datasets (few batches per epoch) the requested 0.999 never washes that out --
    e.g. 200 steps leaves 82% of the init -- so the evaluated EMA weights would be
    mostly noise. Capping at 1 - 1/(total_steps/10) keeps the averaging horizon at
    roughly a tenth of training, which decays the init contribution to ~e^-10.
    """
    if total_steps <= 10:
        return 0.0
    return min(requested_decay, 1.0 - 10.0 / total_steps)


def unwrap_model(model: nn.Module) -> nn.Module:
    """Unwrap AveragedModel (EMA) and compiled model (OptimizedModule) to the base model."""
    if hasattr(model, "module"):
        model = model.module
    if hasattr(model, "_orig_mod"):
        model = model._orig_mod
    return model


class Checkpoint:
    """Manages retrospective checkpointing of EMA weights."""

    def __init__(self, filepath: str | Path, mode: str = "max", min_delta: float = 0.0) -> None:
        self.filepath = Path(filepath)
        self.mode = mode
        self.min_delta = min_delta
        self.best_score = float("-inf") if mode == "max" else float("inf")
        self.best_epoch = 0

    def step(self, score: float, epoch: int, ema_model: nn.Module) -> bool:
        """Saves EMA model state if the validation score improves by min_delta."""
        if self.mode == "max":
            improved = score > self.best_score + self.min_delta
        else:
            improved = score < self.best_score - self.min_delta

        if improved:
            self.best_score = score
            self.best_epoch = epoch
            self.filepath.parent.mkdir(parents=True, exist_ok=True)
            save_file(unwrap_model(ema_model).state_dict(), self.filepath)
            return True
        return False


def make_inferences(
    model: nn.Module,
    loader: DataLoader,
    device: torch.device,
    loss_fn: nn.Module,
) -> tuple[torch.Tensor, torch.Tensor, float]:
    """Collect model output raw logits, ground truth labels, and average loss over the loader."""
    model.eval()
    all_preds = []
    all_trues = []
    total_loss = 0.0
    with torch.no_grad():
        for x, y in loader:
            x, y = x.to(device), y.to(device)
            logits = model(x)
            total_loss += loss_fn(logits, y).item()
            all_preds.append(logits)
            all_trues.append(y)

    y_logits = torch.cat(all_preds, dim=0)
    y_true = torch.cat(all_trues, dim=0)
    return y_logits, y_true, total_loss / len(loader)


def make_regression_inferences(
    model: nn.Module,
    loader: DataLoader,
    device: torch.device,
    loss_fn: nn.Module,
) -> tuple[torch.Tensor, torch.Tensor, float]:
    """Collect raw regression predictions, targets, and average loss over the loader."""
    return make_inferences(model, loader, device, loss_fn)


def setup_optimizer_and_scheduler(
    model: nn.Module,
    config: OptimizerConfig,
    steps_per_epoch: int,
) -> tuple[torch.optim.Optimizer, lr_scheduler.OneCycleLR]:
    """Configure AdamW optimizer with grouped weight decay and get the OneCycleLR scheduler."""
    params = [
        p
        for p in model.parameters()
        if not isinstance(p, torch.nn.parameter.UninitializedParameter)
    ]
    decay_params = [p for p in params if p.dim() >= 2]
    no_decay_params = [p for p in params if p.dim() < 2]

    optimizer = torch.optim.AdamW(
        [
            {"params": decay_params, "weight_decay": config.weight_decay},
            {"params": no_decay_params, "weight_decay": 0.0},
        ],
        lr=config.lr,
    )
    scheduler = lr_scheduler.OneCycleLR(
        optimizer,
        max_lr=config.lr,
        steps_per_epoch=steps_per_epoch,
        epochs=config.epochs,
    )
    return optimizer, scheduler


def train_epoch(
    model: nn.Module,
    loader: DataLoader,
    optimizer: torch.optim.Optimizer,
    scheduler: lr_scheduler.LRScheduler,
    ema_model: AveragedModel,
    device: torch.device,
    compute_loss: Callable[[torch.Tensor, torch.Tensor], torch.Tensor],
) -> float:
    """Train one epoch, updating the EMA after every step. Returns the mean loss."""
    model.train()
    total_loss = 0.0
    for x, y in loader:
        x, y = x.to(device), y.to(device)
        optimizer.zero_grad()
        loss = compute_loss(x, y)
        loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
        optimizer.step()
        scheduler.step()
        total_loss += loss.item()
        ema_model.update_parameters(model)
    return total_loss / len(loader)


def train_model(
    model: nn.Module,
    train_loader: DataLoader,
    val_loader: DataLoader | None,
    optimizer: torch.optim.Optimizer,
    loss_fn: nn.Module,
    checkpoint: Checkpoint | None,
    epochs: int,
    device: torch.device,
    scheduler: lr_scheduler.LRScheduler,
    ema_decay: float = 0.999,
    epoch_callback: Callable[[int, float], None] | None = None,
) -> list[dict[str, float | int | str]]:
    """Train the model using train_loader, validate and log progress, and handle checkpointing."""
    ema_model = AveragedModel(
        model, multi_avg_fn=get_ema_multi_avg_fn(ema_decay)
    ).to(device)

    history = []
    pbar = tqdm(range(1, epochs + 1), desc="Training")
    for epoch in pbar:
        train_loss = train_epoch(
            model, train_loader, optimizer, scheduler, ema_model, device,
            lambda x, y: loss_fn(model(x), y),
        )
        current_lr = optimizer.param_groups[0]["lr"]

        pred_train_logits, y_train, _ = make_inferences(
            ema_model, train_loader, device, loss_fn
        )
        pred_train = torch.sigmoid(pred_train_logits)
        train_score = auroc_score(y_train, pred_train)

        history.append(
            {
                "epoch": epoch,
                "split": "train",
                "loss": train_loss,
                "auroc": train_score,
                "lr": current_lr,
            }
        )

        if val_loader is not None:
            pred_val_logits, y_val, val_loss = make_inferences(
                ema_model, val_loader, device, loss_fn
            )
            pred_val = torch.sigmoid(pred_val_logits)
            val_score = auroc_score(y_val, pred_val)

            history.append(
                {
                    "epoch": epoch,
                    "split": "val",
                    "loss": val_loss,
                    "auroc": val_score,
                    "lr": current_lr,
                }
            )

            if checkpoint is not None:
                checkpoint.step(val_score, epoch, ema_model)

            if epoch_callback is not None:
                epoch_callback(epoch, val_score)

        pbar.set_postfix(history[-1])

    # Unwrap AveragedModel (EMA) and compiled model (OptimizedModule)
    state_dict = unwrap_model(ema_model).state_dict()

    # Load into the unwrapped base model to avoid key mismatches
    unwrap_model(model).load_state_dict(state_dict)

    return history


def train_regression_model(
    model: nn.Module,
    train_loader: DataLoader,
    val_loader: DataLoader | None,
    optimizer: torch.optim.Optimizer,
    loss_fn: nn.Module,
    checkpoint: Checkpoint | None,
    epochs: int,
    device: torch.device,
    scheduler: lr_scheduler.LRScheduler,
    ema_decay: float = 0.999,
    epoch_callback: Callable[[int, float], None] | None = None,
    *,
    log_train_metrics: bool = False,
) -> list[dict[str, float | int | str]]:
    """Same shape as `train_model`, but scores on RMSE (minimized) instead of AUROC (maximized).

    Args:
        log_train_metrics: Also evaluate the EMA model on the training set each epoch to
            record train RMSE. Off by default -- it costs a full extra forward pass per
            epoch, which is wasted work when only the test-fold metrics are reported.
    """
    if checkpoint is not None and checkpoint.mode != "min":
        raise ValueError(
            f"train_regression_model reports RMSE (lower is better) but checkpoint.mode "
            f"is '{checkpoint.mode}'; pass Checkpoint(..., mode='min')."
        )

    ema_decay = resolve_ema_decay(epochs * len(train_loader), ema_decay)
    ema_model = AveragedModel(
        model, multi_avg_fn=get_ema_multi_avg_fn(ema_decay)
    ).to(device)

    history = []
    pbar = tqdm(range(1, epochs + 1), desc="Training")
    for epoch in pbar:
        train_loss = train_epoch(
            model, train_loader, optimizer, scheduler, ema_model, device,
            lambda x, y: loss_fn(model(x), y),
        )
        current_lr = optimizer.param_groups[0]["lr"]

        train_entry: dict[str, float | int | str] = {
            "epoch": epoch,
            "split": "train",
            "loss": train_loss,
            "lr": current_lr,
        }
        if log_train_metrics:
            pred_train, y_train, _ = make_regression_inferences(
                ema_model, train_loader, device, loss_fn
            )
            train_entry["rmse"] = float(
                torch.sqrt(torch.mean((pred_train - y_train) ** 2))
            )
        history.append(train_entry)

        if val_loader is not None:
            pred_val, y_val, val_loss = make_regression_inferences(
                ema_model, val_loader, device, loss_fn
            )
            val_rmse = float(torch.sqrt(torch.mean((pred_val - y_val) ** 2)))

            history.append(
                {
                    "epoch": epoch,
                    "split": "val",
                    "loss": val_loss,
                    "rmse": val_rmse,
                    "lr": current_lr,
                }
            )

            if checkpoint is not None:
                checkpoint.step(val_rmse, epoch, ema_model)

            if epoch_callback is not None:
                epoch_callback(epoch, val_rmse)

        pbar.set_postfix(history[-1])

    state_dict = unwrap_model(ema_model).state_dict()
    unwrap_model(model).load_state_dict(state_dict)

    return history
