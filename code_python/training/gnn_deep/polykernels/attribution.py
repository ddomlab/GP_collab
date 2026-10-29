"""Integrated Gradients attribution for the regression GNNs.

Splits a prediction into per-atom, per-node-feature, and per-tabular-feature credit, with
the completeness guarantee that the parts sum back to `F(input) - F(baseline)`.

The reference implementation is `scripts/Che1148_CodeLab8_IntegratedGradients.ipynb`. Two
adaptations were needed for this repo's models:

* The notebook's model takes `(x, edge_index, edge_attr, batch)`; `GNN.forward` takes a
  single `pyg.data.Data` that also carries a global feature vector `u` (the dataset's
  tabular `x_cols`, injected by `overwrite_dataset_globals`). Attribution therefore covers
  both `x` and `u`, and `Readout` concatenates their two pathways -- see
  `polykernels/models/gnn.py::Readout` -- so the split between them is meaningful.
* The notebook re-runs the model once per interpolation step. Here all `steps` points are
  packed into one `pyg.data.Batch`, so a molecule costs a single forward/backward.

Caveat worth knowing before reading any result: under `GNNConfig.seed_globals=True` (the
default in `get_default_hparams`), `GraphEncoder.forward` takes the
`self.global_proj(x, batch)` branch and never reads `u` -- the global state is synthesized
from mean-pooled node embeddings and log atom count instead. The tabular features are then
not an input to the model at all, `ig_u` is identically zero, and `IGResult.globals_connected`
is False. Verified empirically: scaling `u` by 0, 100, -5, or replacing it with noise moves
the output by exactly 0.0.
"""

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import torch
import torch_geometric as pyg

from polykernels.molgraph.templates import Template

# In-ring-size buckets are the one node feature whose column meanings are not recoverable
# from the feature object (no `values`, no `features`, empty `to_config`).
_IN_RING_SIZE_BUCKETS = ("3", "4", "5", "6", "7to11", "12plus")


def node_feature_names(template: Template) -> list[str]:
    """Human-readable label for every column of `data.x`, in hstack order.

    Derived from the template rather than hardcoded, so it stays correct if the feature
    set changes: block widths come from `GraphFeature.ndim`, and category names from
    `.values` / `.features` / `to_config()` depending on the feature kind.
    """
    names: list[str] = []
    for f in template.node_features:
        names.extend(_labels_for_feature(f))
    total = sum(f.ndim for f in template.node_features)
    if len(names) != total:
        raise ValueError(f"Built {len(names)} labels for {total} node feature columns.")
    return names


def _labels_for_feature(f) -> list[str]:   # f is a GraphFeature; untyped to avoid a cycle
    cls = f.__class__.__name__
    stem = f.name()

    if cls == "GraphwaveEmbeddings":
        cfg = f.to_config()
        per_part = cfg["ndim"] // 2
        times = [i * cfg["step_size"] for i in range(per_part)]
        # Layout is scale-major, then [real block, imaginary block] -- see
        # polykernels/molgraph/graphwave.py.
        return [
            f"graphwave_s{scale:g}_{part}_t{t:g}"
            for scale in cfg["scales"]
            for part in ("re", "im")
            for t in times
        ]
    if cls == "InRingSizes":
        return [f"{stem}={b}" for b in _IN_RING_SIZE_BUCKETS]

    values = getattr(f, "values", None)                      # CategoricalFeature
    if values:
        return [f"{stem}={v}" for v in values]
    class_features = getattr(f.__class__, "features", None)   # multi-hot class attribute
    if class_features:
        return [f"{stem}={v}" for v in class_features]

    cfg = f.to_config()
    if f.ndim > 1 and "min_value" in cfg and "max_value" in cfg:   # IntOneHotFeature
        return [f"{stem}={v}" for v in range(int(cfg["min_value"]), int(cfg["max_value"]) + 1)]
    if f.ndim == 1:
        return [stem]
    return [f"{stem}[{i}]" for i in range(f.ndim)]


def feature_blocks(template: Template) -> list[tuple[str, int, int]]:
    """(name, start, stop) column span of each node feature, for block-level aggregation."""
    blocks, offset = [], 0
    for f in template.node_features:
        blocks.append((f.name(), offset, offset + f.ndim))
        offset += f.ndim
    return blocks


def constant_columns(x: torch.Tensor) -> np.ndarray:
    """Indices of node feature columns that never vary across atoms.

    Some are structural (Graphwave's t=0 term is exp(0) = 1 + 0i on every atom); most are
    one-hot values a dataset never exhibits otherwise (`atom_radical_electrons=0`,
    `chiral_tag=CHI_UNSPECIFIED`). An always-zero column gets identically zero IG. An
    always-nonzero one is credited by a zero baseline even though it cannot separate any
    two molecules, and routinely takes most of the node attribution. Callers must zero
    these before reporting anything; `scripts/analyze_attributions.py` does.
    """
    return np.flatnonzero(x.std(dim=0).cpu().numpy() == 0.0)


@dataclass
class IGResult:
    """Integrated Gradients for one (molecule, target)."""

    ig_x: np.ndarray        # [n_atoms, n_node_features]
    ig_u: np.ndarray        # [n_x_cols]
    f_input: float
    f_baseline: float
    # False when the model's output does not depend on `u` at all, so `ig_u` is exactly
    # zero by construction rather than by measurement. True for `GNNConfig.seed_globals`,
    # which makes `GraphEncoder` synthesize the global state from node embeddings and
    # never read the tabular features -- see the module docstring.
    globals_connected: bool = True

    @property
    def total(self) -> float:
        """Attribution mass, which completeness says should equal f_input - f_baseline."""
        return float(self.ig_x.sum() + self.ig_u.sum())

    @property
    def completeness_abs_error(self) -> float:
        return abs((self.f_input - self.f_baseline) - self.total)

    @property
    def completeness_rel_error(self) -> float:
        """Absolute error relative to the attributed difference; NaN when that is ~0."""
        denom = abs(self.f_input - self.f_baseline)
        return float("nan") if denom < 1e-12 else self.completeness_abs_error / denom

    def signed_node_scores(self) -> np.ndarray:
        """Per-atom signed attribution -- notebook `signed_node_scores_from_ig`."""
        return self.ig_x.sum(axis=1)

    def magnitude_node_scores(self) -> np.ndarray:
        """Per-atom absolute attribution -- notebook `magnitude_node_scores_from_ig`."""
        return np.abs(self.ig_x).sum(axis=1)

    def chemistry_share(self) -> float:
        """Fraction of total |attribution| carried by node features rather than x_cols."""
        chem, proc = float(np.abs(self.ig_x).sum()), float(np.abs(self.ig_u).sum())
        return float("nan") if chem + proc == 0 else chem / (chem + proc)


ScalarFn = Callable[[pyg.data.Batch], torch.Tensor]


def scalar_output_fn(model: torch.nn.Module, model_type: str, target_idx: int) -> ScalarFn:
    """A differentiable [n_graphs] -> scalar-per-graph view of a model, for one target.

    `gnn_gp` and `gnn_gpsn` return a list of `MultivariateNormal`s from `forward`, so their
    predictive mean is taken through `predict()`. That must NOT be wrapped in
    `torch.no_grad` the way `polykernels.models.gp.make_gp_predictions` does, or no gradient
    reaches the input.
    """
    if model_type == "gnn_mlp":
        return lambda batch: model(batch)[:, target_idx]
    if model_type in ("gnn_gp", "gnn_gpsn"):
        import gpytorch

        def _gp_mean(batch: pyg.data.Batch) -> torch.Tensor:
            # gpytorch memoizes the variational strategy across calls in eval mode. The
            # first backward frees that cached graph, so a second attribution raises
            # "Trying to backward through the graph a second time". Clearing the cache
            # per call keeps each attribution independent.
            for gp in model.gp_head.gps:
                gp.variational_strategy._clear_cache()
            # Same jitter as training; inducing points sit on a thin manifold and the
            # default 1e-6 ceiling raises NotPSDError.
            with gpytorch.settings.cholesky_jitter(float_value=1e-3, double_value=1e-6):
                mean, _ = model.predict(batch)
            return mean[:, target_idx]

        return _gp_mean
    raise ValueError(
        f"Node-level attribution needs a graph model; got model_type={model_type!r}."
    )


def integrated_gradients(
    fn: ScalarFn,
    data: pyg.data.Data,
    device: torch.device,
    *,
    steps: int = 64,
) -> IGResult:
    """Integrated Gradients over node features and graph globals jointly.

    Baseline is the all-zero input for both `x` and `u`, so completeness covers the whole
    input and the node/global split accounts for the entire prediction.

    The model must be in eval mode: the encoder uses `LazyBatchNorm1d`, and in train mode
    the interpolation steps -- which share one batch here -- would normalize against each
    other, coupling points that must be independent.

    The path integral uses the trapezoid rule. The notebook's plain mean over a
    linspace that includes both endpoints gives each endpoint a full 1/steps weight
    instead of half -- an O(1/steps) bias. Sweeps at 64 and 512 steps cut completeness
    error almost exactly 8x, and the bias inflated the node-feature share in every run.
    """
    if steps < 2:
        raise ValueError(f"steps must be >= 2 to span both endpoints; got {steps}.")
    x0 = data.x.to(device)
    u0 = data.u.to(device)
    n_atoms, n_feat = x0.shape

    proto = pyg.data.Data(
        x=x0, edge_index=data.edge_index.to(device),
        edge_attr=data.edge_attr.to(device), u=u0,
    )
    batch = pyg.data.Batch.from_data_list([proto.clone() for _ in range(steps)]).to(device)

    # linspace includes both endpoints, so the same forward pass yields F(x) and F(x').
    alphas = torch.linspace(0.0, 1.0, steps, device=device)
    x_path = (alphas.view(-1, 1, 1) * x0.unsqueeze(0)).reshape(steps * n_atoms, n_feat)
    u_path = alphas.view(-1, 1) * u0
    x_path = x_path.detach().requires_grad_()
    u_path = u_path.detach().requires_grad_()
    batch.x, batch.u = x_path, u_path

    out = fn(batch)
    # Each graph's output depends only on its own slice, so one backward on the sum gives
    # every step's gradient at once.
    out.sum().backward()

    weights = torch.full((steps,), 1.0 / (steps - 1), device=device)
    weights[0] = weights[-1] = 0.5 / (steps - 1)
    grads_x = (weights @ x_path.grad.view(steps, -1)).view(n_atoms, n_feat)
    # `u_path.grad` is None when the globals are not in the computation graph at all --
    # the `seed_globals` path discards them. Attribution is then exactly zero, and
    # completeness still holds over the node features alone because the output genuinely
    # does not depend on `u`.
    globals_connected = u_path.grad is not None
    ig_u = (
        (u0.squeeze(0) * (weights @ u_path.grad)).detach().cpu().numpy()
        if globals_connected
        else np.zeros(u0.shape[-1], dtype=np.float32)
    )

    return IGResult(
        ig_x=(x0 * grads_x).detach().cpu().numpy(),
        ig_u=ig_u,
        f_input=float(out[-1].detach()),
        f_baseline=float(out[0].detach()),
        globals_connected=globals_connected,
    )
