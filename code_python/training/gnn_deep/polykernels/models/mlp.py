"""MLP variants (vanilla, gated, spectral-norm residual) and their configuration.

References:
- Gated MLP (SwiGLU): GLU Variants Improve Transformer (arXiv:2002.05202)
- Spectral-norm residual MLP: SNGP (arXiv:2006.10108), DUE (arXiv:2102.11409)

The hyperparameter search space is in `scripts/bayesopt_mlp.py`.
"""

from dataclasses import dataclass
from enum import StrEnum

import draccus
import torch
import torch.nn as nn
import torch.nn.functional as F  # noqa: N812

from polykernels.models import spectral


class BlockType(StrEnum):
    """Supported MLP block types."""

    GATED = "gated"
    VANILLA = "vanilla"


draccus.encode.register(BlockType, lambda x: x.value)


class RMSNorm(nn.Module):
    """Root Mean Square Layer Normalization."""

    def __init__(self, dim: int, eps: float = 1e-6) -> None:
        super().__init__()
        self.eps = eps
        self.weight = nn.Parameter(torch.ones(dim))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        variance = x.pow(2).mean(-1, keepdim=True)
        return (x * torch.rsqrt(variance + self.eps)) * self.weight


class SwiGLU(nn.Module):
    """SwiGLU activation helper chunking the input dimension."""

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x1, x2 = x.chunk(2, dim=-1)
        return F.silu(x1) * x2


class GatedBlock(nn.Module):
    """SwiGLU gated transformer-style feedforward block with RMSNorm."""

    def __init__(self, hidden_dim: int, dropout: float = 0.0) -> None:
        super().__init__()
        inner_dim = int(hidden_dim * (2.0 / 3.0))
        corrected_inner_dim = inner_dim + (inner_dim % 2)

        self.norm = RMSNorm(hidden_dim)
        self.in_proj = nn.Linear(hidden_dim, corrected_inner_dim * 2, bias=False)
        self.act = SwiGLU()
        self.drop = nn.Dropout(dropout)
        self.out_proj = nn.Linear(corrected_inner_dim, hidden_dim, bias=False)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x + self.out_proj(self.drop(self.act(self.in_proj(self.norm(x)))))


class GatedMLP(nn.Module):
    """Gated MLP architecture stacking GatedBlocks."""

    def __init__(
        self,
        hidden_dim: int,
        output_dim: int,
        num_layers: int,
        dropout: float = 0.0,
    ) -> None:
        super().__init__()
        blocks = [nn.LazyLinear(hidden_dim, bias=False)]
        blocks.extend([GatedBlock(hidden_dim, dropout) for _ in range(num_layers)])
        blocks.append(RMSNorm(hidden_dim))

        self.backbone = nn.Sequential(*blocks)
        self.pred_head = nn.Linear(hidden_dim, output_dim, bias=False)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.pred_head(self.backbone(x))

    def embed(self, x: torch.Tensor) -> torch.Tensor:
        return self.backbone(x)


class MLP(nn.Module):
    """Vanilla sequential MLP architecture with LayerNorm and SiLU.

    Note:
        The total number of linear projection layers is num_layers + 1. For example,
        setting num_layers=1 creates a 2-layer MLP (one hidden projection layer in the
        backbone followed by the final prediction head linear layer).
    """

    def __init__(
        self,
        hidden_dim: int,
        output_dim: int,
        num_layers: int,
        dropout: float = 0.0,
    ) -> None:
        super().__init__()
        blocks = []
        for _ in range(num_layers):
            blocks.extend(
                [
                    nn.LazyLinear(hidden_dim, bias=False),
                    nn.LayerNorm(hidden_dim),
                    nn.SiLU(),
                    nn.Dropout(dropout),
                ]
            )

        self.backbone = nn.Sequential(*blocks)
        self.pred_head = nn.Linear(hidden_dim, output_dim, bias=True)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.pred_head(self.backbone(x))

    def embed(self, x: torch.Tensor) -> torch.Tensor:
        return self.backbone(x)


class SNResidualMLP(nn.Module):
    """Distance-preserving MLP feature extractor for a GP head (SNGP / DUE).

    A fixed random input projection, then `num_layers` residual blocks
    h <- h + Dropout(ReLU(W h)) whose weights are spectral-norm bounded, with no
    LayerNorm: a norm layer rescales every embedding and would undo the bound. With
    ||W|| <= bound < 1 each block is bi-Lipschitz, so distances between inputs survive into
    the latent space instead of collapsing. The input projection stays frozen because a
    trained non-residual layer could still collapse input directions (SNGP's reference
    implementation freezes it for the same reason).

    The bound and the freeze are applied by `apply_spectral_norm()` once the lazy input
    projection is materialized -- the training warmup calls it before the GP's inducing
    points are seeded from these embeddings.
    """

    def __init__(
        self,
        hidden_dim: int,
        output_dim: int,
        num_layers: int,
        dropout: float = 0.0,
        spectral_norm_bound: float = 0.9,
    ) -> None:
        super().__init__()
        self.spectral_norm_bound = spectral_norm_bound
        self.in_proj = nn.LazyLinear(hidden_dim, bias=False)
        self.blocks = nn.ModuleList(nn.Linear(hidden_dim, hidden_dim) for _ in range(num_layers))
        self.dropout = nn.Dropout(dropout)
        self.pred_head = nn.Linear(hidden_dim, output_dim, bias=True)

    def apply_spectral_norm(self) -> None:
        """Freeze the input projection and bound every block's weight."""
        self.in_proj.requires_grad_(False)
        spectral.apply_spectral_norm(self, self.spectral_norm_bound)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.pred_head(self.embed(x))

    def embed(self, x: torch.Tensor) -> torch.Tensor:
        h = self.in_proj(x)
        for block in self.blocks:
            h = h + self.dropout(F.relu(block(h)))
        return h


@dataclass
class MLPConfig:
    """Configuration class for MLP architecture."""

    type: str = "mlp"
    block_type: BlockType = BlockType.VANILLA
    hidden_dim: int = 128
    output_dim: int = 256
    num_layers: int = 2
    dropout: float = 0.1
    # > 0 builds `SNResidualMLP` with this spectral-norm bound instead of `block_type`.
    spectral_norm_bound: float = 0.0


def build_mlp(config: MLPConfig) -> nn.Module:
    """Helper factory function to build the selected MLP typology."""
    if config.spectral_norm_bound > 0:
        return SNResidualMLP(
            hidden_dim=config.hidden_dim,
            output_dim=config.output_dim,
            num_layers=config.num_layers,
            dropout=config.dropout,
            spectral_norm_bound=config.spectral_norm_bound,
        )
    if config.block_type == BlockType.GATED:
        return GatedMLP(
            hidden_dim=config.hidden_dim,
            output_dim=config.output_dim,
            num_layers=config.num_layers,
            dropout=config.dropout,
        )
    elif config.block_type == BlockType.VANILLA:
        return MLP(
            hidden_dim=config.hidden_dim,
            output_dim=config.output_dim,
            num_layers=config.num_layers,
            dropout=config.dropout,
        )
    raise ValueError(f"Unsupported block_type: {config.block_type}")
