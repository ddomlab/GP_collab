"""Gaussian-Process output heads for deep-kernel-learning regression models.

Wraps the existing `GNN`/`MLP` feature extractors (unchanged) with a per-target
stochastic variational GP (SVGP) output layer, giving a native predictive mean + std
for uncertainty quantification.
"""

from collections.abc import Callable

import gpytorch
import torch
import torch.nn as nn
import torch.optim.lr_scheduler as lr_scheduler
from torch.optim.swa_utils import AveragedModel, get_ema_multi_avg_fn
from torch.utils.data import DataLoader
import torch_geometric as pyg
from tqdm import tqdm

from polykernels.models.gnn import GNN
from polykernels.models.training import (
    Checkpoint,
    resolve_ema_decay,
    train_epoch,
    unwrap_model,
)


def cholesky_jitter() -> gpytorch.settings.cholesky_jitter:
    """Raise gpytorch's Cholesky jitter ceiling from 1e-6 to 1e-3 for float32.

    Inducing points sit on the thin feature manifold and drift during training, so K_zz
    can become near-singular; at the default ceiling that raises NotPSDError.
    """
    return gpytorch.settings.cholesky_jitter(float_value=1e-3, double_value=1e-6)


class GPOutputLayer(gpytorch.models.ApproximateGP):
    """Single-output stochastic variational GP with learnable inducing points.

    The `torch.randn` inducing points are a placeholder and must be re-seeded with
    `MultiTargetGPHead.init_inducing_from_features` before training (the pipeline does
    this in `warmup_lazy_modules`). ConstantMean beats LinearMean once they are seeded
    (dataset5_rg/mlp_gp: R2 0.84 vs 0.62).
    """

    def __init__(self, num_dim: int, num_inducing_points: int = 64) -> None:
        variational_distribution = gpytorch.variational.CholeskyVariationalDistribution(
            num_inducing_points
        )
        variational_strategy = gpytorch.variational.VariationalStrategy(
            self,
            torch.randn(num_inducing_points, num_dim),
            variational_distribution,
            learn_inducing_locations=True,
        )
        super().__init__(variational_strategy)
        self.mean_module = gpytorch.means.ConstantMean()
        self.covar_module = gpytorch.kernels.ScaleKernel(gpytorch.kernels.RBFKernel())

    def forward(self, x: torch.Tensor) -> gpytorch.distributions.MultivariateNormal:
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)


class MultiTargetGPHead(nn.Module):
    """One independent GP + Gaussian likelihood per target dimension, over shared features."""

    def __init__(
        self, feature_dim: int, output_dim: int, num_inducing_points: int = 256
    ) -> None:
        super().__init__()
        self.gps = nn.ModuleList(
            [
                GPOutputLayer(feature_dim, num_inducing_points)
                for _ in range(output_dim)
            ]
        )
        self.likelihoods = nn.ModuleList(
            [gpytorch.likelihoods.GaussianLikelihood() for _ in range(output_dim)]
        )
        self._mlls: list[gpytorch.mlls.VariationalELBO] | None = None
        self._mll_num_data: int | None = None

    @torch.no_grad()
    def init_inducing_from_features(self, features: torch.Tensor) -> None:
        """Put each GP's inducing points *and* its kernel lengthscale on the data's scale.

        Either default leaves K_xz ~ 0 for every input, so the GP falls back to its
        prior (constant mean, the same variance for every sample) however well the
        features fit. Random inducing points sit far from the thin region `embed()` maps
        data into, and the default lengthscale 0.693 is far below the embeddings'
        pairwise distances (~3-8), where the RBF kernel and its gradient are zero. On
        dataset5_rg/mlp_gp the median-distance lengthscale lifts R2 from -0.003 to 0.83.

        Points too close together make K_zz singular (NotPSDError) instead, so they are
        deduplicated, spread out by farthest-point traversal and jittered.
        """
        features = features.detach()
        unique = torch.unique(features, dim=0)
        scale = float(features.std()) or 1.0
        lengthscale = self._median_pairwise_distance(features)
        for gp in self.gps:
            inducing = gp.variational_strategy.inducing_points
            n = inducing.shape[0]
            sample = self._spread_out_sample(unique, n, scale)
            inducing.data.copy_(sample.to(inducing.dtype))
            gp.covar_module.base_kernel.lengthscale = lengthscale

    @staticmethod
    def _median_pairwise_distance(features: torch.Tensor, max_rows: int = 1024) -> float:
        """Median heuristic: the distance at which the kernel should be ~exp(-1/2).

        Subsampled because the full pairwise distance matrix is O(n**2) and this only
        needs to be right to within a factor of two. Falls back to the default lengthscale
        when there are too few distinct points to estimate from.
        """
        if features.shape[0] > max_rows:
            idx = torch.randperm(features.shape[0], device=features.device)[:max_rows]
            features = features[idx]
        if features.shape[0] < 2:
            return 0.693
        iu = torch.triu_indices(features.shape[0], features.shape[0], offset=1)
        median = torch.cdist(features, features)[iu[0], iu[1]].median().item()
        return median if median > 0 else 0.693

    @staticmethod
    def _spread_out_sample(
        pool: torch.Tensor, n: int, scale: float
    ) -> torch.Tensor:
        """Pick `n` maximally-separated rows from `pool`, padding with jitter if needed."""
        if pool.shape[0] > n:
            # Farthest-point traversal: start random, then repeatedly take the point
            # furthest from everything chosen so far.
            chosen = [int(torch.randint(pool.shape[0], (1,)).item())]
            dist = torch.cdist(pool, pool[chosen]).squeeze(-1)
            for _ in range(n - 1):
                nxt = int(torch.argmax(dist).item())
                chosen.append(nxt)
                dist = torch.minimum(dist, torch.cdist(pool, pool[[nxt]]).squeeze(-1))
            sample = pool[chosen]
        else:
            reps = -(-n // pool.shape[0])  # ceil division
            sample = pool.repeat(reps, 1)[:n]
        # Jitter proportional to the feature scale keeps every point distinct without
        # dragging it off the manifold.
        return sample + 1e-3 * scale * torch.randn_like(sample)

    def marginal_log_likelihoods(
        self, num_data: int
    ) -> list[gpytorch.mlls.VariationalELBO]:
        """Per-target VariationalELBO objectives, built once and cached per `num_data`."""
        if self._mlls is None or self._mll_num_data != num_data:
            self._mlls = [
                gpytorch.mlls.VariationalELBO(likelihood, gp, num_data=num_data)
                for gp, likelihood in zip(self.gps, self.likelihoods, strict=True)
            ]
            self._mll_num_data = num_data
        return self._mlls

    def forward(
        self, features: torch.Tensor
    ) -> list[gpytorch.distributions.MultivariateNormal]:
        return [gp(features) for gp in self.gps]

    def posterior_at(
        self, features: torch.Tensor
    ) -> list[tuple[torch.Tensor, torch.Tensor]]:
        """Per-target latent posterior (mean, FULL covariance) over a batch of features.

        `predict` keeps only the marginal std, which is all the metrics need. GP-SHAP needs
        the joint covariance across the batch, since the explanation covariance propagates
        through it. Returns one `(mean [n], covariance [n, n])` pair per target.
        """
        out = []
        with cholesky_jitter():
            for gp in self.gps:
                latent = gp(features)
                out.append((latent.mean, latent.covariance_matrix))
        return out

    def predict(self, features: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        """Returns (mean, std), each [n, output_dim], marginalized through each target's likelihood."""
        means, stds = [], []
        for gp, likelihood in zip(self.gps, self.likelihoods, strict=True):
            observed = likelihood(gp(features))
            means.append(observed.mean)
            stds.append(observed.stddev)
        return torch.stack(means, dim=-1), torch.stack(stds, dim=-1)


class GNNWithGPHead(nn.Module):
    """GNN feature extractor (unchanged, via `GNN.embed`) feeding a per-target GP head."""

    def __init__(
        self,
        gnn: GNN,
        feature_dim: int,
        output_dim: int,
        num_inducing_points: int = 256,
    ) -> None:
        super().__init__()
        self.gnn = gnn
        self.gp_head = MultiTargetGPHead(feature_dim, output_dim, num_inducing_points)

    def forward(
        self, data: pyg.data.Data
    ) -> list[gpytorch.distributions.MultivariateNormal]:
        return self.gp_head(self.gnn.embed(data))

    def predict(self, data: pyg.data.Data) -> tuple[torch.Tensor, torch.Tensor]:
        return self.gp_head.predict(self.gnn.embed(data))


class MLPWithGPHead(nn.Module):
    """Tabular MLP feature extractor (over x_cols only, no graph, via `MLP.embed`) feeding a per-target GP head."""

    def __init__(
        self,
        mlp: nn.Module,
        feature_dim: int,
        output_dim: int,
        num_inducing_points: int = 256,
    ) -> None:
        super().__init__()
        self.mlp = mlp
        self.gp_head = MultiTargetGPHead(feature_dim, output_dim, num_inducing_points)

    def forward(
        self, x: torch.Tensor
    ) -> list[gpytorch.distributions.MultivariateNormal]:
        return self.gp_head(self.mlp.embed(x))

    def predict(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        return self.gp_head.predict(self.mlp.embed(x))


def gp_elbo_loss(
    gp_head: MultiTargetGPHead,
    latents: list[gpytorch.distributions.MultivariateNormal],
    y: torch.Tensor,
    num_data: int,
) -> torch.Tensor:
    """Sum of per-target negative variational ELBO, each on that target's own (standardized) scale.

    Summing per-target terms after each has been computed on a comparably-scaled
    (standardized) target is what keeps the multi-output loss from being dominated by
    whichever target happens to have the largest raw variance.
    """
    loss = torch.zeros((), device=y.device, dtype=y.dtype)
    mlls = gp_head.marginal_log_likelihoods(num_data)
    for i, (mll, latent) in enumerate(zip(mlls, latents, strict=True)):
        loss = loss - mll(latent, y[:, i])
    return loss


def setup_gp_optimizer_and_scheduler(
    model: nn.Module,
    feature_lr: float,
    gp_lr: float,
    weight_decay: float,
    epochs: int,
    steps_per_epoch: int,
) -> tuple[torch.optim.Optimizer, lr_scheduler.OneCycleLR]:
    """AdamW with separate param groups/learning rates for the feature extractor vs GP head.

    Variational GP parameters (inducing locations, kernel hyperparameters) typically need
    a different learning rate than the upstream feature extractor to converge well. The
    feature extractor keeps the same decay/no-decay split as
    `polykernels.models.training.setup_optimizer_and_scheduler`, so norms and biases are
    not weight-decayed. GP parameters are never decayed.
    """
    base = unwrap_model(model)
    feature_extractor = base.gnn if hasattr(base, "gnn") else base.mlp
    # Skip still-uninitialized LazyLinear parameters: calling .dim() on an
    # UninitializedParameter raises. They materialize on the first forward pass, which
    # `warmup_lazy_modules` triggers before training. Mirrors the same guard in
    # `polykernels.models.training.setup_optimizer_and_scheduler`.
    def _initialized(p: torch.nn.Parameter) -> bool:
        return not isinstance(p, torch.nn.parameter.UninitializedParameter)

    # Frozen weights (a spectral-norm encoder's fixed input projection) are left out.
    trainable = [p for p in feature_extractor.parameters() if _initialized(p) and p.requires_grad]
    decay = [p for p in trainable if p.dim() >= 2]
    no_decay = [p for p in trainable if p.dim() < 2]

    optimizer = torch.optim.AdamW(
        [
            {"params": decay, "lr": feature_lr, "weight_decay": weight_decay},
            {"params": no_decay, "lr": feature_lr, "weight_decay": 0.0},
            {"params": base.gp_head.parameters(), "lr": gp_lr, "weight_decay": 0.0},
        ]
    )
    scheduler = lr_scheduler.OneCycleLR(
        optimizer,
        max_lr=[feature_lr, feature_lr, gp_lr],
        steps_per_epoch=steps_per_epoch,
        epochs=epochs,
    )
    return optimizer, scheduler


def make_gp_predictions(
    model: nn.Module,
    loader: DataLoader,
    device: torch.device,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Collect (mean, std, y_true) over a loader using the model's `.predict()`."""
    model.eval()
    means, stds, trues = [], [], []
    with torch.no_grad(), cholesky_jitter():
        for x, y in loader:
            x, y = x.to(device), y.to(device)
            mean, std = unwrap_model(model).predict(x)
            means.append(mean)
            stds.append(std)
            trues.append(y)
    return torch.cat(means, dim=0), torch.cat(stds, dim=0), torch.cat(trues, dim=0)


def gp_validation_loss(
    model: nn.Module, loader: DataLoader, device: torch.device, num_data: int
) -> float:
    """Average per-batch negative ELBO over a loader, for checkpointing/pruning during training."""
    model.eval()
    total_loss = 0.0
    with torch.no_grad(), cholesky_jitter():
        for x, y in loader:
            x, y = x.to(device), y.to(device)
            latents = model(x)
            total_loss += gp_elbo_loss(
                unwrap_model(model).gp_head, latents, y, num_data
            ).item()
    return total_loss / len(loader)


def train_gp_model(
    model: nn.Module,
    train_loader: DataLoader,
    val_loader: DataLoader | None,
    optimizer: torch.optim.Optimizer,
    checkpoint: Checkpoint | None,
    epochs: int,
    device: torch.device,
    scheduler: lr_scheduler.LRScheduler,
    num_train_samples: int,
    ema_decay: float = 0.999,
    epoch_callback: Callable[[int, float], None] | None = None,
) -> list[dict[str, float | int | str]]:
    """Train a GP-headed regression model (GNNWithGPHead/MLPWithGPHead) via variational ELBO.

    Mirrors `polykernels.models.training.train_model`'s shape (EMA, OneCycleLR,
    grad clipping, retrospective checkpointing), but scores on validation negative ELBO
    (minimized) instead of AUROC (maximized), since GP training optimizes an ELBO
    rather than a plain point-estimate loss.

    Note:
        `AveragedModel` deepcopies `model`, and gpytorch's variational strategy caches
        non-leaf tensors during a grad-enabled forward pass, which are not deepcopy-able.
        Any forward pass before this point must therefore run under `torch.no_grad()`
        (see `warmup_lazy_modules` in scripts/train_gnn_regression.py).
    """
    if checkpoint is not None and checkpoint.mode != "min":
        raise ValueError(
            f"train_gp_model reports negative ELBO (lower is better) but checkpoint.mode "
            f"is '{checkpoint.mode}'; pass Checkpoint(..., mode='min')."
        )

    ema_decay = resolve_ema_decay(epochs * len(train_loader), ema_decay)
    ema_model = AveragedModel(
        model, multi_avg_fn=get_ema_multi_avg_fn(ema_decay)
    ).to(device)

    def compute_loss(x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
        with cholesky_jitter():
            latents = model(x)
            gp_head = unwrap_model(model).gp_head
            return gp_elbo_loss(gp_head, latents, y, num_train_samples)

    history = []
    pbar = tqdm(range(1, epochs + 1), desc="Training (GP)")
    for epoch in pbar:
        train_loss = train_epoch(
            model, train_loader, optimizer, scheduler, ema_model, device, compute_loss
        )

        current_lr = optimizer.param_groups[0]["lr"]
        history.append(
            {"epoch": epoch, "split": "train", "loss": train_loss, "lr": current_lr}
        )

        if val_loader is not None:
            val_loss = gp_validation_loss(
                ema_model, val_loader, device, num_train_samples
            )
            history.append(
                {"epoch": epoch, "split": "val", "loss": val_loss, "lr": current_lr}
            )

            if checkpoint is not None:
                checkpoint.step(val_loss, epoch, ema_model)

            if epoch_callback is not None:
                epoch_callback(epoch, val_loss)

        pbar.set_postfix(history[-1])

    state_dict = unwrap_model(ema_model).state_dict()
    unwrap_model(model).load_state_dict(state_dict)

    return history
