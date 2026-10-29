"""Graph Neural Network (GNN) architectures using MetaLayer and GAT components."""

from dataclasses import dataclass
from enum import StrEnum
import time

import draccus
from jaxtyping import Float, Int
from loguru import logger
import torch
from torch import Tensor
import torch.nn as nn
import torch_geometric as pyg
from torch_geometric.nn.dense.linear import Linear as PygLinear
import torchinfo

from polykernels.models import spectral
from polykernels.models.mlp import MLP


class GlobalUpdate(StrEnum):
    """Supported global update pooling mechanisms."""

    SUM_MAX = "sum_max"
    MEAN = "mean"
    PNA = "pna"


draccus.encode.register(
    GlobalUpdate, lambda x: str(x.value) if hasattr(x, "value") else str(x)
)


@dataclass
class GNNConfig:
    """Configuration class for GNN model architecture."""

    type: str = "gnn"
    node_dim: int = 128
    edge_dim: int = 32
    global_dim: int = 64
    emb_dim: int = 128
    output_dim: int = 138

    global_update: GlobalUpdate = GlobalUpdate.PNA
    depth: int = 4
    dropout: float = 0.05
    num_att_heads: int = 4
    # False feeds the graph globals `u` (the dataset's tabular x_cols) into the encoder.
    # True instead seeds the global state from pooled node embeddings and drops `u`
    # entirely, so any tabular features are silently ignored.
    seed_globals: bool = False
    weight_tied: bool = True
    # > 0 bounds every weight's spectral norm by this value and makes the encoder
    # distance preserving for a GP head (SNGP / DUE); see `_make_distance_preserving`.
    spectral_norm_bound: float = 0.0


def make_proj_block(in_dim: int, out_dim: int, dropout: float = 0.0) -> nn.Sequential:
    """Helper to build a Linear -> LayerNorm -> SiLU -> Dropout sequence."""
    return nn.Sequential(
        nn.Linear(in_dim, out_dim, bias=False),
        nn.LayerNorm(out_dim),
        nn.SiLU(),
        nn.Dropout(dropout),
    )


class GraphGlobalSeed(nn.Module):
    """Seed initial global latent state from node mean and log1p atom count."""

    def __init__(self, node_dim: int, global_dim: int) -> None:
        super().__init__()
        self.proj = nn.Sequential(
            nn.Linear(node_dim + 1, global_dim, bias=False),
            nn.LayerNorm(global_dim),
            nn.SiLU(),
        )

    def forward(
        self,
        x: Float[Tensor, "nodes node_dim"],
        batch: Int[Tensor, "nodes"],  # noqa: F821
    ) -> Float[Tensor, "graphs global_dim"]:
        x_mean = pyg.nn.global_mean_pool(x, batch)
        num_atoms = pyg.nn.global_add_pool(
            x.new_ones((x.size(0), 1)),
            batch,
        )
        log_num_atoms = torch.log1p(num_atoms)
        return self.proj(torch.cat([x_mean, log_num_atoms], dim=-1))


class GraphEncoder(nn.Module):
    """Encodes raw node, edge, and global features into hidden dimensions.

    Uses LazyBatchNorm1d to handle unnormalized inputs and LayerNorm + SiLU
    activations. If seed_globals is True, or if input `u` has 0 features, seeds
    the global representation from node embeddings via `GraphGlobalSeed`.
    """

    def __init__(
        self,
        node_dim: int,
        edge_dim: int,
        global_dim: int,
        seed_globals: bool = False,  # noqa: FBT001, FBT002
    ) -> None:
        super().__init__()
        self.global_dim = global_dim
        self.seed_globals = seed_globals
        self.node_proj = nn.Sequential(
            nn.LazyBatchNorm1d(),
            nn.LazyLinear(node_dim, bias=False),
            nn.LayerNorm(node_dim),
            nn.SiLU(),
        )
        self.edge_proj = nn.Sequential(
            nn.LazyBatchNorm1d(),
            nn.LazyLinear(edge_dim, bias=False),
            nn.LayerNorm(edge_dim),
            nn.SiLU(),
        )
        if seed_globals:
            self.global_proj = GraphGlobalSeed(node_dim, global_dim)
        else:
            self.global_proj = nn.Sequential(
                nn.LazyBatchNorm1d(),
                nn.LazyLinear(global_dim, bias=False),
                nn.LayerNorm(global_dim),
                nn.SiLU(),
            )
            self.fallback_seed = GraphGlobalSeed(node_dim, global_dim)

    def forward(
        self,
        x: Float[Tensor, "nodes node_in"],
        edge_attr: Float[Tensor, "edges edge_in"],
        u: Float[Tensor, "graphs global_in"] | None = None,
        batch: Int[Tensor, "nodes"] | None = None,  # noqa: F821
    ) -> tuple[
        Float[Tensor, "nodes node_dim"],
        Float[Tensor, "edges edge_dim"],
        Float[Tensor, "graphs global_dim"],
    ]:
        x = self.node_proj(x)
        edge_attr = self.edge_proj(edge_attr)
        if batch is None:
            batch = x.new_zeros(x.size(0), dtype=torch.long)

        if self.seed_globals:
            u_out = self.global_proj(x, batch)
        elif u is None or u.size(-1) == 0:
            n_graphs = int(batch.max().item()) + 1 if batch.numel() > 0 else 1
            # Unused on this path; called once so its lazy layers get weights.
            _ = self.global_proj(x.new_zeros((n_graphs, 1)))
            u_out = self.fallback_seed(x, batch)
        else:
            u_out = self._project_globals(u)

        return x, edge_attr, u_out

    def _project_globals(self, u: torch.Tensor) -> torch.Tensor:
        """Encode the tabular globals, tolerating a batch that holds a single graph.

        `u` has one row per graph, so in train mode a one-graph batch makes the leading
        BatchNorm raise "Expected more than 1 value per channel". That batch uses the
        running statistics instead, as eval mode does.
        """
        if not (self.training and u.shape[0] < 2):
            return self.global_proj(u)
        norm = self.global_proj[0]
        was_training = norm.training
        norm.eval()
        try:
            return self.global_proj(u)
        finally:
            norm.train(was_training)


class AdaLN(nn.Module):
    """Parameter-efficient global context injection (FiLM via LayerNorm)."""

    # True conditions additively, x + beta(u): the normalize-and-scale form rescales x by
    # an unbounded factor, which a spectral-norm-bounded encoder cannot allow.
    shift_only: bool = False

    def __init__(self, global_dim: int, feature_dim: int) -> None:
        super().__init__()
        self.ln = nn.LayerNorm(feature_dim, elementwise_affine=False)
        self.context_proj = nn.Linear(global_dim, feature_dim * 2)

        nn.init.zeros_(self.context_proj.weight)
        nn.init.zeros_(self.context_proj.bias)

    def forward(
        self, x: torch.Tensor, u: torch.Tensor, batch: torch.Tensor
    ) -> torch.Tensor:
        gamma, beta = self.context_proj(u[batch]).chunk(2, dim=-1)
        if self.shift_only:
            return x + beta
        return self.ln(x) * (1.0 + gamma) + beta


class EdgeUpdate(nn.Module):
    """Edge feature update network using local conditioning and AdaLN."""

    def __init__(
        self,
        node_dim: int,
        edge_dim: int,
        global_dim: int,
        dropout: float = 0.0,
    ) -> None:
        super().__init__()
        self.edge_mlp = MLP(edge_dim, edge_dim, num_layers=1, dropout=dropout)
        self.adaln = AdaLN(global_dim, edge_dim)

    def forward(
        self,
        src: torch.Tensor,
        dst: torch.Tensor,
        edge_attr: torch.Tensor,
        u: torch.Tensor,
        batch: torch.Tensor,
    ) -> torch.Tensor:
        local_cond = torch.cat([src, dst, edge_attr], dim=-1)
        out = edge_attr + self.edge_mlp(local_cond)
        return self.adaln(out, u, batch)


class NodeAttention(nn.Module):
    """Node feature update network using Graph Attention Networks (GATv2)."""

    def __init__(
        self,
        node_dim: int,
        edge_dim: int,
        global_dim: int,
        num_heads: int = 8,
        dropout: float = 0.0,
    ) -> None:
        super().__init__()
        self.self_attn = pyg.nn.GATv2Conv(
            in_channels=node_dim,
            out_channels=node_dim // num_heads,
            heads=num_heads,
            concat=True,
            edge_dim=edge_dim,
            dropout=dropout,
        )
        self.output_mlp = MLP(node_dim, node_dim, num_layers=1, dropout=dropout)
        self.dropout_layer = nn.Dropout(dropout)
        self.norm1 = AdaLN(global_dim, node_dim)
        self.norm2 = nn.LayerNorm(node_dim)

    def forward(
        self,
        x: torch.Tensor,
        edge_index: torch.Tensor,
        edge_attr: torch.Tensor,
        u: torch.Tensor,
        batch: torch.Tensor,
    ) -> torch.Tensor:
        attn = self.self_attn(x, edge_index, edge_attr=edge_attr)
        out = self.norm1(x + self.dropout_layer(attn), u, batch)
        out = self.norm2(out + self.dropout_layer(self.output_mlp(out)))
        return out


class PooledGlobal(nn.Module):
    """Global update: pool node and edge states per graph, then update `u` residually.

    Subclasses pick the pooling (`make_pool`) and how many statistics it returns per
    feature (`n_stats`).
    """

    n_stats: int

    def __init__(
        self,
        node_dim: int,
        edge_dim: int,
        global_dim: int,
        dropout: float = 0.0,
    ) -> None:
        super().__init__()
        self.node_pool = self.make_pool()
        self.edge_pool = self.make_pool()
        self.proj = make_proj_block(
            global_dim + self.n_stats * (node_dim + edge_dim), global_dim, dropout
        )
        self.out_norm = nn.LayerNorm(global_dim)

    @staticmethod
    def make_pool() -> nn.Module:
        raise NotImplementedError

    def forward(
        self,
        x: torch.Tensor,
        edge_index: torch.Tensor,
        edge_attr: torch.Tensor,
        u: torch.Tensor,
        batch: torch.Tensor,
    ) -> torch.Tensor:
        edge_batch = batch[edge_index[0]]
        node_aggr = self.node_pool(x, batch, dim_size=u.size(0))
        edge_aggr = self.edge_pool(edge_attr, edge_batch, dim_size=u.size(0))
        features = torch.cat([u, node_aggr, edge_aggr], dim=-1)
        delta = self.proj(features)
        return self.out_norm(u + delta)


class SumMaxGlobal(PooledGlobal):
    """Global update using Sum/Max aggregation over nodes and edges."""

    n_stats = 2

    @staticmethod
    def make_pool() -> nn.Module:
        return pyg.nn.aggr.MultiAggregation(["sum", "max"])


class MeanGlobal(PooledGlobal):
    """Global update using Mean aggregation over nodes and edges."""

    n_stats = 1

    @staticmethod
    def make_pool() -> nn.Module:
        return pyg.nn.aggr.MeanAggregation()


class PNAGlobal(PooledGlobal):
    """Global update using PNA aggregation over nodes and edges."""

    n_stats = 4

    @staticmethod
    def make_pool() -> nn.Module:
        return pyg.nn.aggr.MultiAggregation(["mean", "var", "min", "max"])


_GLOBAL_MODEL_REGISTRY = {
    GlobalUpdate.SUM_MAX: SumMaxGlobal,
    GlobalUpdate.MEAN: MeanGlobal,
    GlobalUpdate.PNA: PNAGlobal,
}


class Readout(nn.Module):
    """Readout layer combining pooled node states and inner global state."""

    def __init__(
        self,
        node_dim: int,
        global_inner_dim: int,
        output_dim: int,
        dropout: float = 0.0,
    ) -> None:
        super().__init__()
        self.pool = pyg.nn.aggr.MeanAggregation()
        self.readout = make_proj_block(node_dim + global_inner_dim, output_dim, dropout)

    def forward(
        self,
        x: torch.Tensor,
        u_inner: torch.Tensor,
        batch: torch.Tensor,
    ) -> torch.Tensor:
        x_pooled = self.pool(x, batch, dim_size=u_inner.size(0))
        fused = torch.cat([x_pooled, u_inner], dim=-1)
        return self.readout(fused)


class LearnedStepEmbedding(nn.Module):
    """Learned step context embedding with magnitude scaling for tied architectures."""

    def __init__(self, depth: int, dim: int) -> None:
        super().__init__()
        self.embedding = nn.Embedding(depth, dim)
        self.scaler = nn.Parameter(torch.ones(depth, 1))
        nn.init.zeros_(self.embedding.weight)

    def forward(self, step: int, u: torch.Tensor) -> torch.Tensor:
        step_ctx = self.embedding.weight[step] * self.scaler[step]
        return u + step_ctx.to(dtype=u.dtype)


class GraphNets(nn.Module):
    """Tied-weight or stacked weight Graph Network iterating MetaLayer blocks."""

    def __init__(
        self,
        node_dim: int,
        edge_dim: int,
        global_dim: int,
        global_update: GlobalUpdate | str = GlobalUpdate.SUM_MAX,
        depth: int = 4,
        dropout: float = 0.0,
        num_heads: int = 4,
        weight_tied: bool = True,
    ) -> None:
        super().__init__()
        self.depth = depth
        self.weight_tied = weight_tied
        global_update_enum = GlobalUpdate(global_update)

        if weight_tied:
            global_model = _GLOBAL_MODEL_REGISTRY[global_update_enum](
                node_dim, edge_dim, global_dim, dropout
            )
            self.shared_layer = pyg.nn.MetaLayer(
                edge_model=EdgeUpdate(node_dim, edge_dim, global_dim, dropout),
                node_model=NodeAttention(
                    node_dim,
                    edge_dim,
                    global_dim,
                    num_heads=num_heads,
                    dropout=dropout,
                ),
                global_model=global_model,
            )
            self.step_context = LearnedStepEmbedding(depth, global_dim)
        else:
            self.layers = nn.ModuleList(
                [
                    pyg.nn.MetaLayer(
                        edge_model=EdgeUpdate(
                            node_dim, edge_dim, global_dim, dropout
                        ),
                        node_model=NodeAttention(
                            node_dim,
                            edge_dim,
                            global_dim,
                            num_heads=num_heads,
                            dropout=dropout,
                        ),
                        global_model=_GLOBAL_MODEL_REGISTRY[global_update_enum](
                            node_dim, edge_dim, global_dim, dropout
                        ),
                    )
                    for _ in range(depth)
                ]
            )

    def forward(
        self,
        x: torch.Tensor,
        edge_index: torch.Tensor,
        edge_attr: torch.Tensor,
        u: torch.Tensor,
        batch: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        if self.weight_tied:
            for step in range(self.depth):
                u_conditioned = self.step_context(step, u)
                x, edge_attr, u = self.shared_layer(
                    x, edge_index, edge_attr, u_conditioned, batch
                )
        else:
            for layer in self.layers:
                x, edge_attr, u = layer(x, edge_index, edge_attr, u, batch)
        return x, edge_attr, u


# Backwards compatibility alias
TiedGraphNets = GraphNets


class GNN(nn.Module):
    """End-to-end Graph Neural Network model combining GraphNets and prediction head."""

    def __init__(
        self,
        config: GNNConfig,
    ) -> None:
        super().__init__()
        self.config = config
        self.global_update = config.global_update

        self.encoder = GraphEncoder(
            node_dim=config.node_dim,
            edge_dim=config.edge_dim,
            global_dim=config.global_dim,
            seed_globals=config.seed_globals,
        )
        self.gnn = GraphNets(
            node_dim=config.node_dim,
            edge_dim=config.edge_dim,
            global_dim=config.global_dim,
            global_update=config.global_update,
            depth=config.depth,
            dropout=config.dropout,
            num_heads=config.num_att_heads,
            weight_tied=config.weight_tied,
        )

        self.readout = Readout(
            node_dim=config.node_dim,
            global_inner_dim=config.global_dim,
            output_dim=config.emb_dim,
            dropout=config.dropout,
        )

        self.pred_head = MLP(
            hidden_dim=config.emb_dim,
            output_dim=config.output_dim,
            num_layers=1,
            dropout=config.dropout,
        )

        self.spectral_norm_bound = config.spectral_norm_bound
        if self.spectral_norm_bound > 0:
            self._make_distance_preserving()

    @property
    def latent_dim(self) -> int:
        """Width of `embed()`: `emb_dim`, or the unprojected readout under spectral norm."""
        if self.spectral_norm_bound > 0:
            return self.config.node_dim + self.config.global_dim
        return self.config.emb_dim

    def _make_distance_preserving(self) -> None:
        """Rewire the encoder so spectral normalization can preserve input distances.

        A LayerNorm rescales its input to unit spread, which both cancels any bound on the
        preceding weight and pins every embedding to the same shell, so all of them go.
        The input projections become a fixed standardization followed by a plain linear map
        (frozen later, in `apply_spectral_norm`), the readout keeps `[mean_pool(x), u]`
        whole rather than projecting it down, and AdaLN conditions additively.

        PyG's `Linear` (inside GATv2) is swapped for an equivalent `nn.Linear`: its custom
        `__deepcopy__` cannot copy a parametrized weight, which `AveragedModel` needs.
        """
        encoder = self.encoder
        for proj in (encoder.node_proj, encoder.edge_proj, encoder.global_proj):
            if isinstance(proj, nn.Sequential):  # BatchNorm, LazyLinear, LayerNorm, SiLU
                proj[0] = nn.LazyBatchNorm1d(affine=False)
                proj[3] = nn.Identity()
        self.readout.readout = nn.Identity()

        for parent in list(self.modules()):
            if isinstance(parent, AdaLN):
                parent.shift_only = True
            for name, child in list(parent.named_children()):
                if isinstance(child, nn.LayerNorm):
                    setattr(parent, name, nn.Identity())
                elif isinstance(child, PygLinear):
                    linear = nn.Linear(
                        child.in_channels, child.out_channels, bias=child.bias is not None
                    )
                    with torch.no_grad():
                        linear.weight.copy_(child.weight)
                        if child.bias is not None:
                            linear.bias.copy_(child.bias)
                    setattr(parent, name, linear)

    def apply_spectral_norm(self) -> None:
        """Freeze the input projections and bound every other weight's spectral norm.

        Called once the lazy layers are materialized, before training: the GP warmup does
        it before seeding inducing points, so those see the bounded encoder.
        """
        for proj in (self.encoder.node_proj, self.encoder.edge_proj, self.encoder.global_proj):
            if isinstance(proj, nn.Sequential):  # proj[1] is the (materialized) Linear
                proj[1].requires_grad_(False)
        spectral.apply_spectral_norm(self, self.spectral_norm_bound)

    def forward(self, data: pyg.data.Data) -> torch.Tensor:
        emb = self.embed(data)
        return self.pred_head(emb)

    def embed(self, data: pyg.data.Data) -> torch.Tensor:
        x, edge_index, edge_attr, u, batch = (
            data.x,
            data.edge_index,
            data.edge_attr,
            data.u,
            data.batch,
        )
        if batch is None:
            batch = torch.zeros(x.size(0), dtype=torch.long, device=x.device)

        x, edge_attr, u_pre = self.encoder(x, edge_attr, u, batch)
        x, edge_attr, u_inner = self.gnn(x, edge_index, edge_attr, u_pre, batch)
        return self.readout(x, u_inner, batch)


class ModelSummaryWrapper(nn.Module):
    """Wrapper to enable model summary for PyG models using torchinfo."""

    def __init__(self, model: nn.Module) -> None:
        super().__init__()
        self.model = model

    def forward(
        self,
        x: Float[Tensor, "num_nodes node_dim"],
        edge_index: Int[Tensor, "2 num_edges"],
        edge_attr: Float[Tensor, "num_edges edge_dim"],
        u: Float[Tensor, "num_graphs global_dim"],
        batch: Int[Tensor, "num_nodes"],  # noqa: F821
    ) -> Float[Tensor, "num_graphs output_dim"]:
        data = pyg.data.Data(
            x=x, edge_index=edge_index, edge_attr=edge_attr, u=u, batch=batch
        )
        return self.model(data)


def summarize_gnn(
    model: nn.Module,
    sample_graph: pyg.data.Data,
    device: torch.device | None = None,
) -> str:
    """Generate a detailed architectural summary of a GNN model using torchinfo."""
    if device is None:
        device = next(model.parameters()).device

    model.to(device)
    dummy_graph = sample_graph.clone().to(device)

    was_training = model.training
    model.eval()
    try:
        with torch.no_grad():
            _ = model(dummy_graph)

        t0 = time.perf_counter()
        with torch.no_grad():
            _ = model(dummy_graph)
        elapsed_ms = (time.perf_counter() - t0) * 1000
        logger.info(
            f"Model single-batch dummy inference time on {device}: {elapsed_ms:.3f} ms"
        )

        dummy_batch = torch.zeros(
            dummy_graph.x.size(0), dtype=torch.long, device=device
        )
        summary = torchinfo.summary(
            ModelSummaryWrapper(model),
            input_data=[
                dummy_graph.x,
                dummy_graph.edge_index,
                dummy_graph.edge_attr,
                dummy_graph.u,
                dummy_batch,
            ],
            depth=5,
            verbose=0,
        )
    finally:
        if was_training:
            model.train()
    return str(summary)
