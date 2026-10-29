"""Models sub-package exposing neural network classes and training helpers.

Note:
    `polykernels.models.gp` is deliberately NOT re-exported here. It imports gpytorch,
    which is only needed by the regression pipeline -- re-exporting it would make every
    consumer of this package (including the classification scripts) require gpytorch.
    Import it directly instead: `from polykernels.models.gp import GNNWithGPHead`.
"""

from polykernels.models.ensemble import BetaNLLLoss, GNNEnsemble
from polykernels.models.gnn import (
    GNN,
    GlobalUpdate,
    GNNConfig,
    ModelSummaryWrapper,
    summarize_gnn,
)
from polykernels.models.losses import (
    ClassificationLossType,
    LossType,
    RegressionLossType,
    get_loss_fn,
)
from polykernels.models.mlp import (
    MLP,
    BlockType,
    GatedBlock,
    GatedMLP,
    MLPConfig,
    RMSNorm,
    build_mlp,
)
from polykernels.models.training import (
    Checkpoint,
    OptimizerConfig,
    make_inferences,
    make_regression_inferences,
    setup_optimizer_and_scheduler,
    train_model,
    train_regression_model,
    unwrap_model,
)

__all__ = [
    "BlockType",
    "MLPConfig",
    "GNNConfig",
    "GlobalUpdate",
    "OptimizerConfig",
    "LossType",
    "ClassificationLossType",
    "RegressionLossType",
    "MLP",
    "RMSNorm",
    "GatedBlock",
    "GatedMLP",
    "GNN",
    "GNNEnsemble",
    "BetaNLLLoss",
    "ModelSummaryWrapper",
    "build_mlp",
    "Checkpoint",
    "make_inferences",
    "make_regression_inferences",
    "setup_optimizer_and_scheduler",
    "train_model",
    "train_regression_model",
    "unwrap_model",
    "get_loss_fn",
    "summarize_gnn",
]
