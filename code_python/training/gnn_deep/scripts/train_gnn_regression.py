"""Train and evaluate regression models on the new polymer/membrane datasets.

Runs 3-seed x 5-fold IID cross-validation and, for problems with `cluster_cols`,
leave-one-cluster-out OOD splits + size-matched IID baselines, using hyperparameters
already tuned by `scripts/bayesopt_gnn_regression.py` (or the defaults).

OOD only, one held-out cluster of the collaborator's grouping, in their protocol (see
scripts/list_ood_tasks.py for the keys):
    python scripts/train_gnn_regression.py --problem dataset3_membrane --model_type gnn_gp \
        --split_types '[loco]' --ood_clusters structure_cluster \
        --only_cluster "structure_cluster=1" --seeds '[42]' \
        --results_dir results/gnn_regression_ood_collab/dataset3_membrane_gnn_gp/structure_cluster_1

Example:
    conda run -n polykernels python scripts/train_gnn_regression.py \
        --problem dataset3_membrane --model_type gnn_gp \
        --hparams_path results/gnn_regression_bayesopt/dataset3_membrane_gnn_gp/best_params.yaml

Test run:
    conda run -n polykernels python scripts/train_gnn_regression.py \
        --problem dataset3_membrane --model_type gnn_mlp --hparams_path <path> \
        --seeds '[6]' --n_folds 2 --results_dir results/test_gnn_regression
"""

from dataclasses import dataclass, field
from pathlib import Path
import time
from typing import Literal

import draccus
import numpy as np
import pandas as pd
from safetensors.torch import save_file
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Subset, TensorDataset
import torch_geometric.data
import torch_geometric.loader
import yaml

from polykernels import paths_and_constants, regression_data, splits, uq_metrics
from polykernels.datasets import (
    GraphDataset,
    overwrite_dataset_globals,
    prepare_graph_dataset,
)
import polykernels.models as models
from polykernels.models import ensemble
from polykernels.models.gp import (
    GNNWithGPHead,
    MLPWithGPHead,
    make_gp_predictions,
    setup_gp_optimizer_and_scheduler,
    train_gp_model,
)
from polykernels.molgraph import GraphType, get_template
from polykernels.utils import get_logger, set_seed

_P = paths_and_constants.Paths()
logger = get_logger()

ModelType = Literal["gnn_mlp", "gnn_gp", "mlp_gp", "gnn_gpsn", "mlp_gpsn", "gnn_mlpens"]
# `*_gpsn` are the GP models with a spectral-normalized, distance-preserving encoder
# (SNGP / DUE): the same hyperparameters as `gnn_gp` / `mlp_gp`, plus this bound on every
# weight's spectral norm -- the `norm_multiplier` of SNGP's reference implementation.
SPECTRAL_NORM_BOUND = 0.9
# `gnn_mlpens` is a deep ensemble of `gnn_mlp`s with a variance head (see
# polykernels.models.ensemble), giving the point-prediction model a predictive std. Like
# every model type it has exactly one underscore: run directories named
# <problem>_<model_type> are parsed by splitting off the last two `_`-separated tokens.
GRAPH_MODELS = ("gnn_mlp", "gnn_gp", "gnn_gpsn", "gnn_mlpens")
TABULAR_MODELS = ("mlp_gp", "mlp_gpsn")
SEEDS = [6, 13, 42]


@dataclass
class TrainConfig:
    """Master configuration for a regression CV run: which problem, model, and tuned hyperparameters."""

    problem: str = "dataset3_membrane"
    model_type: ModelType = "gnn_mlp"
    graph: GraphType = GraphType.CUSTOM
    # Empty (the default) trains from `get_default_hparams`; pass a path to override with
    # a tuned config. Note that a best_params.yaml written before the GNNConfig rename
    # (`global_type` -> `global_update`) will no longer load.
    hparams_path: str = ""  # noqa: E501 -- see resolve_hparams for why this is str, not Path | None
    results_dir: Path = _P.results_dir / "gnn_regression"
    seeds: list[int] = field(default_factory=lambda: list(SEEDS))
    n_folds: int = 5
    # Which job families to run: "iid" (seeds x n_folds random CV) and/or "loco" (one
    # held-out cluster per job, each with its size-matched IID baselines). "loco" is a
    # no-op for problems without cluster_cols. Passing only "loco" adds OOD results next
    # to an existing IID sweep without repeating it.
    split_types: list[str] = field(default_factory=lambda: ["iid", "loco"])
    # Restrict "loco" to one cluster key as printed by scripts/list_ood_tasks.py (e.g.
    # "Donor=PBDB-T"), so a SLURM array can run one held-out cluster per task.
    only_cluster: str = ""
    # Size-matched IID baseline folds per (cluster, seed); 0 runs the LOCO folds alone.
    matched_baseline_folds: int = 5
    # One of the collaborator's groupings, "structure_cluster" or "continuous_cluster"
    # (imported by scripts/import_collab_clusters.py). Set, "loco" follows their protocol
    # instead of the spec's cluster_cols: every cluster held out once, each matched by
    # five random splits of its size (see make_cv_jobs). They evaluate OOD with seed 42.
    ood_clusters: str = ""
    # Save every job's trained model (IID, loco and matched-IID folds alike) to
    # <results_dir>/models/<tag>/, named like its predictions file, together with the split
    # and standardization it used, in the layout the attribution scripts read
    # (model.safetensors, hparams.yaml, standardization.npz). The attribution scripts can
    # then explain the exact models behind the CV metrics. Saving happens after
    # evaluation, so the metrics are identical with or without it.
    save_models: bool = False



def get_default_hparams(model_type: ModelType) -> dict:
    """Default hyperparameter fallback matching pom-explore optimization setup."""
    if model_type in ("gnn_gpsn", "mlp_gpsn"):
        hparams = get_default_hparams(model_type.removesuffix("sn"))
        hparams["model"]["spectral_norm_bound"] = SPECTRAL_NORM_BOUND
        return hparams
    if model_type == "gnn_mlpens":
        # The gnn_mlp configuration unchanged, so the two differ only in the variance head,
        # the loss (beta-NLL instead of `loss_type`'s MSE) and the number of networks.
        hparams = get_default_hparams("gnn_mlp")
        del hparams["gp"]
        del hparams["optimizer"]["loss_type"]
        hparams["ensemble"] = {"n_members": ensemble.N_MEMBERS, "beta": ensemble.BETA}
        return hparams
    if model_type in ("gnn_mlp", "gnn_gp"):
        return {
            "model": {
                "node_dim": 128,
                "edge_dim": 32,
                "global_dim": 64,
                "emb_dim": 128,
                "global_update": "pna",
                "depth": 4,
                "dropout": 0.05,
                "num_att_heads": 4,
                # False so the tabular x_cols reach the model as graph globals. Under True
                # the encoder synthesizes its global state from node embeddings and never
                # reads `u`, making the process conditions dead inputs -- verified by
                # scrambling `u` and watching the output move by exactly 0.0.
                "seed_globals": False,
                "weight_tied": True,
            },
            "gp": {
                "num_inducing_points": 256,
                # Required by setup_gp_optimizer_and_scheduler. 0.05 is the median of the
                # values all 12 completed Optuna searches converged to (range 0.018-0.099);
                # the GP's variational parameters need a markedly higher LR than the
                # feature extractor's 0.005.
                "gp_lr": 0.05,
            },
            "optimizer": {
                "batch_size": 256,
                "lr": 0.005,
                "epochs": 200,
                "weight_decay": 6e-5,
                "loss_type": "mse",
            },
        }
    elif model_type == "mlp_gp":
        return {
            "model": {
                "hidden_dim": 128,
                "num_layers": 3,
                "dropout": 0.1,
            },
            "gp": {
                "num_inducing_points": 256,
                # Required by setup_gp_optimizer_and_scheduler. 0.05 is the median of the
                # values all 12 completed Optuna searches converged to (range 0.018-0.099);
                # the GP's variational parameters need a markedly higher LR than the
                # feature extractor's 0.005.
                "gp_lr": 0.05,
            },
            "optimizer": {
                "batch_size": 256,
                "lr": 0.005,
                "epochs": 200,
                "weight_decay": 6e-5,
            },
        }
    raise ValueError(f"Unknown model_type: {model_type}")


def resolve_hparams(config: "TrainConfig") -> dict:
    """Hyperparameters for this run: defaults unless an explicit override path is given.

    Kept as one helper so the CV loop and the run-record written to hparams.yaml can never
    disagree about which configuration was actually trained.

    `TrainConfig.hparams_path` is a `str` (empty meaning "unset") rather than the more
    natural `Path | None` because draccus turns each field's annotation into an argparse
    `type=` callable, and neither `Path | None` nor `Optional[Path]` is callable -- either
    one raises `TypeError` while building the parser, before main() runs at all.
    """
    if not config.hparams_path:
        logger.info(f"Using default hyperparameters for {config.model_type}.")
        return get_default_hparams(config.model_type)
    path = Path(config.hparams_path)
    if not path.exists():
        raise FileNotFoundError(
            f"--hparams_path {path} does not exist. Omit the flag to "
            f"train from get_default_hparams()."
        )
    logger.info(f"Loading hyperparameter overrides from {path}")
    with open(path) as f:
        return yaml.safe_load(f)





def standardize_columns(
    a: np.ndarray, train_idx: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-column zero-mean/unit-variance standardization, fit on `train_idx` only.

    Fitting on the training rows alone keeps the test fold out of the statistics. Used
    for both targets (so a multi-output loss is scale-balanced across targets) and
    tabular features (whose raw scales span orders of magnitude -- e.g. dataset3 mixes
    molecular weight ~1e5 with pressure ~1e0).
    """
    mean = a[train_idx].mean(axis=0)
    std = a[train_idx].std(axis=0)
    std = np.where(std < 1e-8, 1.0, std)
    return (a - mean) / std, mean, std


# Kept as a named alias: targets need the returned mean/std to un-standardize predictions.
standardize_targets = standardize_columns


def build_graph_dataset(
    smiles: np.ndarray,
    x_features: np.ndarray,
    y: np.ndarray,
    *,
    graph_type: GraphType = GraphType.CUSTOM,
    use_3d: bool = True,
) -> GraphDataset:
    """Featurize SMILES into graphs once using the specified graph template."""
    if graph_type == GraphType.CUSTOM:
        template = get_template(GraphType.CUSTOM, global_dim=x_features.shape[1])
    else:
        template = get_template(graph_type)
    return prepare_graph_dataset(smiles, y, graph_type=template, use_3d=use_3d)



def make_fold_loaders(
    model_type: ModelType,
    graph_dataset: GraphDataset | None,
    x_std: np.ndarray,
    y_std: np.ndarray,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    batch_size: int,
) -> tuple[DataLoader, DataLoader]:
    """Builds (train, test) loaders for one CV fold from fold-standardized features/targets."""
    if model_type in TABULAR_MODELS:
        x_t = torch.from_numpy(x_std).float()
        y_t = torch.from_numpy(y_std).float()
        train_ds = TensorDataset(x_t[train_idx], y_t[train_idx])
        test_ds = TensorDataset(x_t[test_idx], y_t[test_idx])
        return (
            DataLoader(train_ds, batch_size=batch_size, shuffle=True),
            DataLoader(test_ds, batch_size=batch_size, shuffle=False),
        )

    assert graph_dataset is not None
    fold_dataset = GraphDataset(
        graphs=graph_dataset.graphs, labels=torch.from_numpy(y_std).float()
    )
    if graph_dataset.graphs[0].u.shape[1] > 0:
        overwrite_dataset_globals(fold_dataset, np.ascontiguousarray(x_std))
    return (
        torch_geometric.loader.DataLoader(
            Subset(fold_dataset, train_idx), batch_size=batch_size, shuffle=True
        ),
        torch_geometric.loader.DataLoader(
            Subset(fold_dataset, test_idx), batch_size=batch_size, shuffle=False
        ),
    )



def build_model(
    model_type: ModelType, hparams: dict, output_dim: int, device: torch.device
) -> nn.Module:
    """Constructs the model for `model_type` from a tuned hyperparameter dict.

    The GP variants consume the backbone's `.embed()` output, so their wrapped model's
    own prediction head is unused -- it is sized to `output_dim` anyway (rather than
    left at the classification default of 138) to avoid carrying thousands of dead
    parameters through the optimizer, weight decay, and checkpoints.
    """
    if model_type == "gnn_mlp":
        gnn_cfg = models.GNNConfig(**hparams["model"], output_dim=output_dim)
        return models.GNN(gnn_cfg).to(device)
    if model_type == "gnn_mlpens":
        # Each member's head is twice as wide: a mean and a raw variance per target.
        gnn_cfg = models.GNNConfig(**hparams["model"], output_dim=2 * output_dim)
        members = [models.GNN(gnn_cfg) for _ in range(hparams["ensemble"]["n_members"])]
        return models.GNNEnsemble(members).to(device)
    # The spectral-norm variants differ only through `spectral_norm_bound` in their hparams.
    if model_type in ("gnn_gp", "gnn_gpsn"):
        gnn_cfg = models.GNNConfig(**hparams["model"], output_dim=output_dim)
        gnn = models.GNN(gnn_cfg)
        return GNNWithGPHead(
            gnn,
            feature_dim=gnn.latent_dim,
            output_dim=output_dim,
            num_inducing_points=hparams["gp"]["num_inducing_points"],
        ).to(device)
    if model_type in TABULAR_MODELS:
        mlp_cfg = models.MLPConfig(**hparams["model"], output_dim=output_dim)
        mlp = models.build_mlp(mlp_cfg)
        return MLPWithGPHead(
            mlp,
            feature_dim=mlp_cfg.hidden_dim,
            output_dim=output_dim,
            num_inducing_points=hparams["gp"]["num_inducing_points"],
        ).to(device)
    raise ValueError(f"Unknown model_type: {model_type}")


def warmup_lazy_modules(model: nn.Module, train_loader: DataLoader, device: torch.device) -> None:
    """Initialize LazyLinear layers so the optimizer can read their parameter shapes.

    For GP-headed models this also seeds the inducing points from real embeddings, which
    has to happen after the lazy layers exist but before training starts.

    Must stay under `no_grad`: gpytorch's variational strategy caches non-leaf tensors
    during a grad-enabled pass, which would then break `AveragedModel`'s deepcopy.
    """
    model.train()
    x, _ = next(iter(train_loader))
    x = x.to(device)
    with torch.no_grad():
        model(x)
        if hasattr(model, "gp_head"):
            _warmup_gp_head(model, x, train_loader, device)

    still_lazy = [
        n
        for n, p in model.named_parameters()
        if isinstance(p, torch.nn.parameter.UninitializedParameter)
    ]
    if still_lazy:
        raise RuntimeError(
            f"Parameters still uninitialized after warmup: {still_lazy}. The warmup "
            f"forward pass does not reach them; EMA averaging will fail on them."
        )


def _warmup_gp_head(
    model: nn.Module,
    x: torch.Tensor | torch_geometric.data.Data,
    train_loader: DataLoader,
    device: torch.device,
) -> None:
    """Materialize the unused backbone head and seed inducing points from real embeddings."""
    extractor = model.gnn if hasattr(model, "gnn") else model.mlp
    # A GP model's forward stops at `.embed()`, so the backbone's own prediction head is
    # never reached and its LazyLinear parameters stay uninitialized. They are excluded
    # from the optimizer, but `AveragedModel.update_parameters` still calls `.detach()`
    # on every parameter and raises on an uninitialized one -- so run the full backbone
    # once to materialize them.
    extractor(x)
    # A spectral-norm bound can only be registered on materialized weights, and must be in
    # place before the inducing points and lengthscale are fitted to these embeddings.
    if getattr(extractor, "spectral_norm_bound", 0.0) > 0:
        extractor.apply_spectral_norm()

    # Gather enough embeddings to seed every inducing point without duplicates:
    # num_inducing_points is searched up to 256 while a batch may hold only 128.
    needed = max(
        gp.variational_strategy.inducing_points.shape[0] for gp in model.gp_head.gps
    )
    chunks, total = [], 0
    for xb, _ in train_loader:
        emb = extractor.embed(xb.to(device))
        chunks.append(emb)
        total += emb.shape[0]
        if total >= needed:
            break
    model.gp_head.init_inducing_from_features(torch.cat(chunks, dim=0))


def train_one_fold(
    model_type: ModelType,
    model: nn.Module,
    hparams: dict,
    train_loader: DataLoader,
    val_loader: DataLoader | None,
    device: torch.device,
    num_train_samples: int,
    checkpoint: models.Checkpoint | None = None,
    epoch_callback=None,
) -> list[dict]:
    """Trains one model for the tuned number of epochs, dispatching by model type."""
    warmup_lazy_modules(model, train_loader, device)
    epochs = hparams["optimizer"]["epochs"]
    optimizer_config = models.OptimizerConfig(
        lr=hparams["optimizer"]["lr"],
        weight_decay=hparams["optimizer"]["weight_decay"],
        epochs=epochs,
        batch_size=hparams["optimizer"]["batch_size"],
    )
    if model_type == "gnn_mlp":
        optimizer, scheduler = models.setup_optimizer_and_scheduler(
            model, optimizer_config, steps_per_epoch=len(train_loader)
        )
        loss_fn = models.get_loss_fn(hparams["optimizer"]["loss_type"])
        return models.train_regression_model(
            model,
            train_loader,
            val_loader,
            optimizer,
            loss_fn,
            checkpoint,
            epochs,
            device,
            scheduler,
            epoch_callback=epoch_callback,
        )
    if model_type == "gnn_mlpens":
        # train_regression_model scores a validation loader as `pred - y`, which a member's
        # [N, 2T] output would broadcast into a wrong RMSE rather than fail on.
        if val_loader is not None or checkpoint is not None or epoch_callback is not None:
            raise ValueError(
                "gnn_mlpens trains for the fixed epoch budget only; it takes no validation "
                "loader, checkpoint or epoch callback."
            )
        loss_fn = models.BetaNLLLoss(hparams["ensemble"]["beta"])
        history = []
        # One member after another, each with its own optimizer, schedule and EMA, all
        # drawing on the RNG state set for this job: initializations and batch orders
        # differ between members but are reproducible.
        for i, member in enumerate(model.members):
            optimizer, scheduler = models.setup_optimizer_and_scheduler(
                member, optimizer_config, steps_per_epoch=len(train_loader)
            )
            member_history = models.train_regression_model(
                member, train_loader, None, optimizer, loss_fn, None, epochs, device, scheduler
            )
            history.extend({**record, "member": i} for record in member_history)
        return history

    optimizer, scheduler = setup_gp_optimizer_and_scheduler(
        model,
        feature_lr=hparams["optimizer"]["lr"],
        gp_lr=hparams["gp"]["gp_lr"],
        weight_decay=hparams["optimizer"]["weight_decay"],
        epochs=epochs,
        steps_per_epoch=len(train_loader),
    )
    return train_gp_model(
        model,
        train_loader,
        val_loader,
        optimizer,
        checkpoint,
        epochs,
        device,
        scheduler,
        num_train_samples,
        epoch_callback=epoch_callback,
    )


def evaluate_fold(
    model_type: ModelType,
    model: nn.Module,
    test_loader: DataLoader,
    device: torch.device,
    target_mean: np.ndarray,
    target_std: np.ndarray,
    target_names: list[str],
    loss_fn: nn.Module | None = None,
) -> tuple[list[dict], dict[str, np.ndarray]]:
    """Evaluates a trained fold, un-standardizes predictions, and returns per-target metrics.

    Returns:
        rows: one metrics dict per target dimension (never pooled across targets).
        predictions: {"y_true", "mean", "std"} arrays in the targets' original units,
            saved per fold so the reported metrics stay auditable after the run.
    """
    if model_type == "gnn_mlp":
        if loss_fn is None:
            raise ValueError("loss_fn is required to evaluate the gnn_mlp model type.")
        pred_std, y_true_std, _ = models.make_regression_inferences(
            model, test_loader, device, loss_fn
        )
        mean_pred = pred_std.cpu().numpy() * target_std + target_mean
        y_true = y_true_std.cpu().numpy() * target_std + target_mean
        std_pred = None
    else:
        # Any model whose `.predict` returns (mean, std): the GP heads and GNNEnsemble.
        mean_pred_std, std_pred_std, y_true_std = make_gp_predictions(
            model, test_loader, device
        )
        mean_pred = mean_pred_std.cpu().numpy() * target_std + target_mean
        std_pred = std_pred_std.cpu().numpy() * target_std
        y_true = y_true_std.cpu().numpy() * target_std + target_mean

    predictions = {"y_true": y_true, "mean": mean_pred}
    if std_pred is not None:
        predictions["std"] = std_pred

    rows = []
    for i, name in enumerate(target_names):
        row = {
            "target": name,
            "rmse": uq_metrics.rmse_score(y_true[:, i], mean_pred[:, i]),
            "mae": uq_metrics.mae_score(y_true[:, i], mean_pred[:, i]),
            "r2": uq_metrics.r2_score(y_true[:, i], mean_pred[:, i]),
        }
        if std_pred is not None:
            row.update(
                {
                    "nll": uq_metrics.gaussian_nll(
                        y_true[:, i], mean_pred[:, i], std_pred[:, i]
                    ),
                    "ece": uq_metrics.compute_ece(
                        y_true[:, i], mean_pred[:, i], std_pred[:, i]
                    ),
                    "ama": uq_metrics.compute_ama(
                        y_true[:, i], mean_pred[:, i], std_pred[:, i]
                    ),
                    "rusc": uq_metrics.compute_rusc(
                        y_true[:, i], mean_pred[:, i], std_pred[:, i]
                    ),
                    "cv": uq_metrics.compute_cv(std_pred[:, i]),
                    "sharpness": uq_metrics.compute_sharpness(std_pred[:, i]),
                }
            )
        rows.append(row)
    return rows, predictions


def make_cv_jobs(
    n: int,
    seeds: list[int],
    n_folds: int,
    cluster_labels: dict[str, np.ndarray] | None,
    max_baseline_folds: int | None = 5,
    *,
    min_cluster_size: int = 3,
    split_types: tuple[str, ...] | list[str] = ("iid", "loco"),
    only_cluster: str = "",
    collab_protocol: bool = False,
) -> list[dict]:
    """Builds the full list of CV jobs: IID folds x seeds, plus LOCO + matched-IID-baseline for datasets with clusters.

    Args:
        max_baseline_folds: Cap on size-matched IID baseline folds per (cluster, seed).
            Without it a small held-out cluster implies many folds and the baselines
            dominate the run (a 15-of-256 cluster yields 17 model fits per seed). 0 skips
            the baselines entirely.
        min_cluster_size: Smallest cluster eligible to be held out.
        split_types: Which job families to build, "iid" and/or "loco".
        only_cluster: If set, the one LOCO cluster key (e.g. "Donor=PBDB-T") to build.
        collab_protocol: Build "loco" the collaborator's way, from a single cluster
            column: every cluster is held out (no size floor, both halves of a two-cluster
            split), and each is matched by one random split of its size per
            `splits.COLLAB_IID_SEEDS` instead of K-fold. Those splits ignore the training
            seed and `max_baseline_folds`, except that 0 still skips them.
    """
    unknown = set(split_types) - {"iid", "loco"}
    if unknown:
        raise ValueError(f"Unknown split_types {sorted(unknown)}; expected 'iid' and/or 'loco'.")

    jobs = []
    if "iid" in split_types:
        for seed in seeds:
            for fold_id, split in enumerate(splits.get_regression_kfold_splits(n, n_folds, seed)):
                jobs.append(
                    {
                        "split_type": "iid",
                        "seed": seed,
                        "fold_id": str(fold_id),
                        "train_idx": split.train,
                        "test_idx": split.test,
                    }
                )

    if cluster_labels and "loco" in split_types:
        if collab_protocol:
            if len(cluster_labels) != 1:
                raise ValueError("The collaborator protocol holds out groups of one column.")
            ((column, labels),) = cluster_labels.items()
            loco_splits = {
                f"{column}={label}": split
                for label, split in splits.get_leave_one_group_out_splits(labels).items()
            }
            # A group's matched-IID seeds depend on its position among ALL groups, so the
            # positions are fixed here, before only_cluster narrows the set.
            group_index = {key: i for i, key in enumerate(loco_splits)}
        else:
            loco_splits = splits.get_loco_splits(cluster_labels, min_size=min_cluster_size)
        if only_cluster:
            if only_cluster not in loco_splits:
                raise ValueError(
                    f"only_cluster={only_cluster!r} is not an eligible cluster; "
                    f"choose from {sorted(loco_splits)}."
                )
            loco_splits = {only_cluster: loco_splits[only_cluster]}
        for seed in seeds:
            for cluster_name, split in loco_splits.items():
                jobs.append(
                    {
                        "split_type": "loco",
                        "seed": seed,
                        "fold_id": cluster_name,
                        "train_idx": split.train,
                        "test_idx": split.test,
                    }
                )
                if max_baseline_folds == 0:
                    baseline_folds = {}
                elif collab_protocol:
                    # Keyed by the split seed, so fold ids name their exact counterpart.
                    baseline_folds = {
                        str(base_seed): splits.get_size_matched_random_split(
                            n, len(split.test), group_index[cluster_name], base_seed
                        )
                        for base_seed in splits.COLLAB_IID_SEEDS
                    }
                else:
                    baseline_folds = dict(enumerate(splits.get_matched_iid_baseline_splits(
                        n, len(split.test), seed, max_folds=max_baseline_folds
                    )))
                for b_id, b_split in baseline_folds.items():
                    jobs.append(
                        {
                            "split_type": f"matched_iid[{cluster_name}]",
                            "seed": seed,
                            "fold_id": str(b_id),
                            "train_idx": b_split.train,
                            "test_idx": b_split.test,
                        }
                    )
    return jobs


def run_regression_cv(
    config: TrainConfig, device: torch.device, predictions_dir: Path | None = None
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Runs the full CV evaluation for one (problem, model_type) using tuned hyperparameters.

    Returns:
        The per-target metrics, and a long-format timing table with one row per
        (stage, job). Timings are kept out of the metrics frame because a single fold
        emits one metrics row per target, which would duplicate its duration and make any
        naive sum over-count.
    """
    spec = regression_data.REGRESSION_PROBLEMS[config.problem]
    smiles, x_features, y, target_names, cluster_labels = regression_data.load_regression_data(
        spec
    )
    output_dim = y.shape[1]
    if config.ood_clusters:
        cluster_labels = {
            config.ood_clusters: regression_data.load_collab_clusters(
                config.problem, config.ood_clusters, y
            )
        }

    hparams = resolve_hparams(config)

    graph_dataset = None
    featurize_seconds = 0.0
    if config.model_type in GRAPH_MODELS:
        logger.info(f"Featurizing {len(smiles)} SMILES into graphs using template '{config.graph}'...")
        featurize_start = time.perf_counter()
        graph_dataset = build_graph_dataset(
            smiles, x_features, y, graph_type=config.graph, use_3d=spec.use_3d
        )
        featurize_seconds = time.perf_counter() - featurize_start
        logger.info(
            f"Featurized in {featurize_seconds:.1f}s. This runs once and every CV job "
            f"reuses the same graphs, so it is not a per-seed cost."
        )


    jobs = make_cv_jobs(
        len(y),
        config.seeds,
        config.n_folds,
        cluster_labels,
        max_baseline_folds=config.matched_baseline_folds,
        min_cluster_size=regression_data.min_cluster_size(spec, len(y)),
        split_types=config.split_types,
        only_cluster=config.only_cluster,
        collab_protocol=bool(config.ood_clusters),
    )
    if not jobs:
        raise ValueError(
            f"No CV jobs for {config.problem} with split_types={config.split_types}: "
            f"'loco' needs a problem with cluster_cols."
        )
    logger.info(f"Running {len(jobs)} CV jobs for {config.problem}/{config.model_type}...")

    all_rows = []
    timings = [{
        "stage": "featurize", "problem": config.problem, "model_type": config.model_type,
        "seed": -1, "fold_id": "-", "split_type": "-", "seconds": featurize_seconds,
    }]
    for job_num, job in enumerate(jobs, start=1):
        set_seed(job["seed"])
        # Both fit on the train fold only, so the held-out fold never enters the statistics.
        y_std, target_mean, target_std = standardize_columns(y, job["train_idx"])
        x_std, x_mean, x_scale = standardize_columns(x_features, job["train_idx"])

        train_loader, test_loader = make_fold_loaders(
            config.model_type,
            graph_dataset,
            x_std,
            y_std,
            job["train_idx"],
            job["test_idx"],
            hparams["optimizer"]["batch_size"],
        )

        fit_start = time.perf_counter()
        model = build_model(config.model_type, hparams, output_dim, device)
        train_one_fold(
            config.model_type,
            model,
            hparams,
            train_loader,
            None,
            device,
            len(job["train_idx"]),
        )
        fit_seconds = time.perf_counter() - fit_start

        loss_fn = (
            models.get_loss_fn(hparams["optimizer"]["loss_type"])
            if config.model_type == "gnn_mlp"
            else None
        )
        eval_start = time.perf_counter()
        rows, predictions = evaluate_fold(
            config.model_type,
            model,
            test_loader,
            device,
            target_mean,
            target_std,
            target_names,
            loss_fn,
        )
        eval_seconds = time.perf_counter() - eval_start
        for stage, seconds in (("fit", fit_seconds), ("evaluate", eval_seconds)):
            timings.append({
                "stage": stage, "problem": config.problem,
                "model_type": config.model_type, "seed": job["seed"],
                "fold_id": job["fold_id"], "split_type": job["split_type"],
                "seconds": seconds,
            })
        for row in rows:
            row.update(
                {
                    "problem": config.problem,
                    "model_type": config.model_type,
                    "seed": job["seed"],
                    "fold_id": job["fold_id"],
                    "split_type": job["split_type"],
                }
            )
        all_rows.extend(rows)

        tag = f"{job['split_type']}_seed{job['seed']}_fold{job['fold_id']}"
        safe_tag = "".join(c if c.isalnum() or c in "-_" else "_" for c in tag)
        if config.save_models:
            save_fold_model(
                config.results_dir / "models" / safe_tag,
                model,
                {"seed": job["seed"], "fold_id": job["fold_id"], "split_type": job["split_type"]},
                hparams,
                output_dim,
                len(y),
                target_names,
                list(spec.x_cols),
                (target_mean, target_std, x_mean, x_scale),
                job["train_idx"],
                job["test_idx"],
            )
        if predictions_dir is not None:
            np.savez_compressed(
                predictions_dir / f"{safe_tag}.npz",
                test_idx=job["test_idx"],
                target_names=np.array(target_names),
                **predictions,
            )

        logger.info(
            f"[{job_num}/{len(jobs)}] {job['split_type']} seed={job['seed']} "
            f"fold={job['fold_id']} done in {fit_seconds:.1f}s fit + {eval_seconds:.1f}s eval"
        )

    return pd.DataFrame(all_rows), pd.DataFrame(timings)


def save_fold_model(
    out_dir: Path,
    model: nn.Module,
    job: dict,
    hparams: dict,
    output_dim: int,
    n_rows: int,
    target_names: list[str],
    x_cols: list[str],
    stats: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    train_idx: np.ndarray,
    test_idx: np.ndarray,
) -> None:
    """Write one CV job's model where the attribution scripts can load it.

    The attribution scripts read `hparams.yaml` for `hparams` and `output_dim`, and
    `standardization.npz` for the statistics the model was trained under and for
    `fit_idx` / `held_out_idx`, which here are the job's train and test folds.
    """
    target_mean, target_std, x_mean, x_std = stats
    out_dir.mkdir(parents=True, exist_ok=True)
    save_file(model.state_dict(), out_dir / _P.model_file)
    np.savez(
        out_dir / "standardization.npz",
        target_mean=target_mean, target_std=target_std,
        x_mean=x_mean, x_std=x_std,
        target_names=np.array(target_names, dtype=object),
        x_cols=np.array(x_cols, dtype=object),
        fit_idx=train_idx, held_out_idx=test_idx,
    )
    with open(out_dir / _P.hparams_file, "w") as f:
        yaml.safe_dump(
            {"cv_job": job, "hparams": hparams, "output_dim": output_dim,
             "n_rows": n_rows},
            f, sort_keys=False,
        )


def log_timing_summary(timings_df: pd.DataFrame) -> None:
    """Report fit time per seed, which is what a `seconds per seed` table column wants.

    Featurization is reported on its own line and never divided across seeds: it runs once
    before the CV loop and all jobs share the resulting graphs, so charging it per seed
    would bill one cost `len(seeds)` times. Quote the per-seed fit number against other
    models' `fit()` timings, and the featurization separately as a one-off setup cost.
    """
    fits = timings_df[timings_df.stage == "fit"]
    if fits.empty:
        return
    per_seed = fits.groupby("seed")["seconds"].agg(["sum", "count"])
    logger.info("Timing summary:")
    for seed, row in per_seed.iterrows():
        logger.info(
            f"  seed {seed}: {row['sum']:.1f}s of fitting over {int(row['count'])} folds"
        )
    logger.info(
        f"  mean per seed: {per_seed['sum'].mean():.1f}s "
        f"(sd {per_seed['sum'].std():.1f}s over {len(per_seed)} seeds)"
    )
    featurize = float(timings_df.loc[timings_df.stage == "featurize", "seconds"].sum())
    logger.info(f"  featurization, once and shared by every job: {featurize:.1f}s")
    logger.info(
        f"  evaluation: {timings_df.loc[timings_df.stage == 'evaluate', 'seconds'].sum():.1f}s total"
    )


@draccus.wrap()
def main(config: TrainConfig) -> None:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    config.results_dir.mkdir(parents=True, exist_ok=True)
    predictions_dir = config.results_dir / "predictions"
    predictions_dir.mkdir(exist_ok=True)

    run_config = {"cv": draccus.encode(config), "hparams": resolve_hparams(config)}
    with open(config.results_dir / _P.hparams_file, "w") as f:
        yaml.safe_dump(run_config, f, sort_keys=False)


    metrics_df, timings_df = run_regression_cv(
        config, device, predictions_dir=predictions_dir
    )
    metrics_df.to_csv(config.results_dir / _P.metrics_file, index=False)
    timings_df.to_csv(config.results_dir / "timings.csv", index=False)
    log_timing_summary(timings_df)

    # Aggregate across seeds/folds within each target, so differently-scaled targets are never pooled.
    metric_cols = [
        c
        for c in metrics_df.select_dtypes(include="number").columns
        if c not in ("seed", "fold_id")
    ]
    summary_df = (
        metrics_df.groupby(["target", "split_type"])[metric_cols]
        .agg(["mean", "std"])
        .reset_index()
    )
    summary_df.to_csv(config.results_dir / "metrics_summary.csv", index=False)

    logger.info(f"Results saved to {config.results_dir}")


if __name__ == "__main__":
    main()
