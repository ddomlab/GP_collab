"""SHAP feature attribution over the tabular `x_cols`, for any of the regression models.

Players are the tabular `x_cols` for every model type, so this is the analysis that says
which process conditions drive each property. For the graph models (`gnn_*`) the molecule
is fixed *context*: each explained row keeps its own graph while its conditions are masked,
so a value reads as "contribution of this condition, for this molecule", and the base value
is that molecule's expected prediction over the background conditions. Structure itself
stays the job of Integrated Gradients (`scripts/analyze_attributions.py`). GP-SHAP
(`scripts/analyze_gpshap.py`) differs on exactly this point: it never evaluates the
explained row's own molecule.

Uses the model-agnostic permutation explainer, which needs nothing from the model beyond a
numpy-in/numpy-out callable -- convenient here, since the GP head is not differentiable in
the way a gradient explainer would want. A graph model is explained one row at a time,
through a callable that holds that row's molecule fixed and swaps only its globals `u`.

Explains a saved model in `<results_dir>/<problem>_<model_type>/`.
`hpc/killarney/attribution_cv.sh` copies each CV fold's model there (from a sweep run with
`--save_models true`) and runs this script on it.

Examples:
    python scripts/analyze_shap.py --problem dataset4_pervaporation

    python scripts/analyze_shap.py --problem dataset4_pervaporation --model_type gnn_mlp

    # Ground-truth cross-check: exact SHAP is tractable at d<=13
    python scripts/analyze_shap.py --problem dataset4_pervaporation --algorithm exact
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import draccus
import gpytorch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from safetensors.torch import load_file
from scripts.train_gnn_regression import (
    GRAPH_MODELS,
    ModelType,
    build_graph_dataset,
    build_model,
    make_fold_loaders,
    warmup_lazy_modules,
)
import shap
import torch
import torch_geometric as pyg
import yaml

from polykernels import paths_and_constants, regression_data
from polykernels.datasets import overwrite_dataset_globals
from polykernels.molgraph import GraphType
from polykernels.utils import get_logger, set_seed
from polykernels.vis import (
    plot_attribution_heatmap,
    save_figure,
    set_visualization_style,
)

_P = paths_and_constants.Paths()
logger = get_logger()


@dataclass
class ShapConfig:
    """Configuration for a SHAP attribution run on one fitted model."""

    problem: str = "dataset4_pervaporation"
    model_type: ModelType = "mlp_gp"
    graph: GraphType = GraphType.CUSTOM
    results_dir: Path = _P.results_dir / "shap"
    algorithm: Literal["permutation", "exact"] = "permutation"
    # Permutation rounds, as a total evaluation budget. Measured against exact SHAP on
    # dataset4 (d=6, 12 rows), max |error| as a percentage of the largest SHAP value:
    #     2d+1 (13):  13.2%  -- and the GLOBAL RANKING disagreed with exact
    #     5x   (65):   6.7%
    #     20x (260):   3.0%
    # The 2d+1 minimum is far too coarse to report a feature ranking from, so the default
    # is shap's own 500 (18x-38x the minimum here, since d is 6-13). Additivity does NOT
    # detect this: permutation SHAP satisfies it by construction at any budget.
    max_evals: int = 500
    # Rows drawn from the training split for shap.maskers.Independent.
    background_size: int = 100
    # 0 explains every held-out row; anything else truncates (smoke tests).
    max_explain_rows: int = 0
    # Rows per forward pass. Permutation SHAP hands `f` very wide batches.
    chunk: int = 4096
    # The same for a graph model, where every row is a whole copy of the molecule.
    graph_chunk: int = 1024
    seed: int = 6


def load_model(
    model_type: ModelType, run_dir: Path, hparams: dict, output_dim: int, graph_dataset,
    x_std, y_std, fit_idx, device: torch.device,
) -> torch.nn.Module:
    """Rebuild, materialize lazy layers, then load trained weights -- in that order.

    Every extractor starts with lazy layers, so `load_state_dict` cannot populate parameter
    shapes that a forward pass has not created yet. The same warmup re-registers a `*_gpsn`
    extractor's spectral-norm parametrization, which its saved weights are keyed under.
    """
    train_loader, _ = make_fold_loaders(
        model_type, graph_dataset, x_std, y_std, fit_idx, fit_idx,
        hparams["optimizer"]["batch_size"],
    )
    model = build_model(model_type, hparams, output_dim, device)
    warmup_lazy_modules(model, train_loader, device)
    model.load_state_dict(load_file(run_dir / _P.model_file, device=str(device)))
    model.eval()   # every extractor carries Dropout; also puts a GP in predictive mode
    return model


def predictive_mean(model, model_type: ModelType, inputs) -> torch.Tensor:
    """Standardized predictive mean `[n, n_targets]`: a GP head's mean, or `gnn_mlp`'s output."""
    if model_type == "gnn_mlp":
        return model(inputs)
    # The default 1e-6 jitter ceiling raises NotPSDError on these inducing points.
    with gpytorch.settings.cholesky_jitter(float_value=1e-3, double_value=1e-6):
        mean, _ = model.predict(inputs)
    return mean


def make_predict_fn(model, model_type: ModelType, device: torch.device, target_mean,
                    target_std, chunk: int):
    """numpy [n, d] -> numpy [n, n_targets], in each target's physical units.

    Un-standardizing the output is affine, so SHAP values come out directly in physical
    units rather than standard deviations.
    """
    def f(x_batch: np.ndarray) -> np.ndarray:
        preds = []
        for start in range(0, len(x_batch), chunk):
            xb = torch.from_numpy(
                np.ascontiguousarray(x_batch[start : start + chunk], dtype=np.float32)
            ).to(device)
            with torch.no_grad():
                preds.append(predictive_mean(model, model_type, xb).cpu().numpy())
        return np.concatenate(preds, axis=0) * target_std + target_mean

    return f


def make_graph_predict_fn(model, model_type: ModelType, graph: pyg.data.Data,
                          device: torch.device, target_mean, target_std, chunk: int):
    """numpy [n, d] of x_cols -> numpy [n, n_targets] for ONE molecule, held fixed.

    Each input row becomes a copy of `graph` whose globals `u` are that row, so only the
    tabular features vary. The copies are laid out back to back once, `chunk` of them, and
    every call takes a prefix: shap calls `f` thousands of times per explained row, and
    collating the copies anew on each call would repeat the same work.
    """
    n_nodes, n_edges = graph.num_nodes, graph.edge_index.shape[1]
    copies = torch.arange(chunk, device=device)
    x = graph.x.to(device).repeat(chunk, 1)
    edge_attr = graph.edge_attr.to(device).repeat(chunk, 1)
    edge_index = graph.edge_index.to(device).repeat(1, chunk) + (
        copies * n_nodes
    ).repeat_interleave(n_edges)
    batch = copies.repeat_interleave(n_nodes)

    def f(x_batch: np.ndarray) -> np.ndarray:
        preds = []
        for start in range(0, len(x_batch), chunk):
            u = torch.from_numpy(
                np.ascontiguousarray(x_batch[start : start + chunk], dtype=np.float32)
            ).to(device)
            k = len(u)
            data = pyg.data.Data(
                x=x[: k * n_nodes], edge_index=edge_index[:, : k * n_edges],
                edge_attr=edge_attr[: k * n_edges], u=u, batch=batch[: k * n_nodes],
            )
            with torch.no_grad():
                preds.append(predictive_mean(model, model_type, data).cpu().numpy())
        return np.concatenate(preds, axis=0) * target_std + target_mean

    return f


def explain_graph_rows(config: ShapConfig, model, graphs, x_std, explain_idx, background,
                       x_cols, stats, device: torch.device, kwargs: dict):
    """SHAP for a graph model, one row at a time, each with its own molecule held fixed.

    Returns values `[n_rows, d, n_targets]`, base values and model outputs
    `[n_rows, n_targets]`. The base value differs per row: it is that molecule's mean
    prediction over the background conditions.
    """
    masker = shap.maskers.Independent(x_std[background])
    values, base, f_x = [], [], []
    for n, row in enumerate(explain_idx, start=1):
        predict = make_graph_predict_fn(
            model, config.model_type, graphs[int(row)], device,
            stats["target_mean"], stats["target_std"], config.graph_chunk,
        )
        # One seed per row: each explainer reseeds numpy's global RNG (see `main`), so a
        # shared seed would replay the same permutations on every row.
        explainer = shap.Explainer(
            model=predict, algorithm=config.algorithm, masker=masker, feature_names=x_cols,
            seed=config.seed + n,
        )
        explanation = explainer(x_std[[row]], silent=True, **kwargs)
        values.append(as_3d(np.asarray(explanation.values))[0])
        base.append(np.asarray(explanation.base_values).reshape(-1))
        f_x.append(predict(x_std[[row]])[0])
        if n % 10 == 0 or n == len(explain_idx):
            logger.info(f"  {n}/{len(explain_idx)} rows explained")
    return np.stack(values), np.stack(base), np.stack(f_x)


def check_fixed_molecule(model, model_type: ModelType, graphs, explain_idx, f_x, stats,
                         device: torch.device) -> None:
    """The fixed-molecule callable must reproduce the model on each row's own graph.

    Additivity cannot catch a wrong callable -- the values sum to whatever `f` returns --
    so the replicated batch is checked against an ordinary collated one here. The two agree
    exactly on CPU; on GPU, reduction order still varies with the batch slightly.
    """
    batch = pyg.data.Batch.from_data_list([graphs[int(i)] for i in explain_idx]).to(device)
    with torch.no_grad():
        direct = predictive_mean(model, model_type, batch).cpu().numpy()
    direct = direct * stats["target_std"] + stats["target_mean"]
    gap = float(np.abs(direct - f_x).max())
    logger.info(f"Fixed-molecule callable vs direct prediction: max |gap| {gap:.3e}")
    if gap > 1e-3 * max(1.0, float(np.abs(direct).max())):
        logger.warning(
            "The fixed-molecule callable does not reproduce the model's own prediction, "
            "so these SHAP values explain a different function. Check the replicated batch."
        )


def as_3d(values: np.ndarray) -> np.ndarray:
    """Normalize SHAP output to [n_rows, n_features, n_targets] for uniform handling."""
    return values[:, :, None] if values.ndim == 2 else values


@draccus.wrap()
def main(config: ShapConfig) -> None:
    set_seed(config.seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    run_dir = config.results_dir / f"{config.problem}_{config.model_type}"
    if not (run_dir / _P.model_file).exists():
        raise FileNotFoundError(
            f"No weights at {run_dir / _P.model_file}. Copy a saved CV fold model there; "
            f"hpc/killarney/attribution_cv.sh does this for a sweep run with --save_models."
        )

    with open(run_dir / _P.hparams_file) as f_in:
        record = yaml.safe_load(f_in)
    hparams, output_dim = record["hparams"], record["output_dim"]
    stats = np.load(run_dir / "standardization.npz", allow_pickle=True)
    x_cols = [str(c) for c in stats["x_cols"]]
    target_names = [str(t) for t in stats["target_names"]]
    fit_idx, held_out = stats["fit_idx"], stats["held_out_idx"]

    spec = regression_data.REGRESSION_PROBLEMS[config.problem]
    smiles, x_features, y, _, _ = regression_data.load_regression_data(spec)
    # SHAP must operate in the standardized space the model was trained on.
    x_std = (x_features - stats["x_mean"]) / stats["x_std"]
    y_std = (y - stats["target_mean"]) / stats["target_std"]

    if len(held_out) == 0:
        logger.warning(
            "This model was fitted on all rows, so background and explained rows overlap "
            "and the SHAP values describe an in-sample fit. Explain a CV fold model instead."
        )
        held_out = np.arange(len(y))
    explain_idx = held_out[: config.max_explain_rows] if config.max_explain_rows else held_out

    graph_dataset = None
    if config.model_type in GRAPH_MODELS:
        if hparams["model"].get("seed_globals"):
            logger.warning(
                "This model was fitted with seed_globals=True, under which it never reads "
                "the tabular x_cols: every SHAP value will be exactly zero."
            )
        logger.info(f"Featurizing {len(smiles)} SMILES...")
        graph_dataset = build_graph_dataset(
            smiles, x_std, y, graph_type=config.graph, use_3d=spec.use_3d
        )
        if graph_dataset.graphs[0].u.shape[1] > 0:
            overwrite_dataset_globals(graph_dataset, np.ascontiguousarray(x_std))

    model = load_model(
        config.model_type, run_dir, hparams, output_dim, graph_dataset, x_std, y_std,
        fit_idx, device,
    )

    n_features = x_std.shape[1]
    max_evals = config.max_evals or (2 * n_features + 1)
    rng = np.random.default_rng(config.seed)
    background = fit_idx[rng.permutation(len(fit_idx))[: config.background_size]]
    logger.info(
        f"{config.algorithm} SHAP: {len(explain_idx)} rows x {n_features} features x "
        f"{output_dim} target(s), background {len(background)}, max_evals {max_evals}."
    )

    kwargs = {"max_evals": max_evals} if config.algorithm == "permutation" else {}
    if graph_dataset is None:
        predict = make_predict_fn(
            model, config.model_type, device, stats["target_mean"], stats["target_std"],
            config.chunk,
        )
        # The permutation explainer calls np.random.seed(seed) when it is built, so leaving
        # seed=None reseeds from the OS and undoes set_seed: the run is not reproducible.
        explainer = shap.Explainer(
            model=predict,
            algorithm=config.algorithm,
            masker=shap.maskers.Independent(x_std[background]),
            feature_names=x_cols,
            seed=config.seed,
        )
        explanation = explainer(x_std[explain_idx], **kwargs)
        values = as_3d(np.asarray(explanation.values))
        base = np.asarray(explanation.base_values).reshape(len(explain_idx), -1)
        f_x = predict(x_std[explain_idx])
    else:
        values, base, f_x = explain_graph_rows(
            config, model, graph_dataset.graphs, x_std, explain_idx, background, x_cols,
            stats, device, kwargs,
        )
        check_fixed_molecule(
            model, config.model_type, graph_dataset.graphs, explain_idx, f_x, stats, device
        )

    run_dir.mkdir(parents=True, exist_ok=True)
    _write_tables(run_dir, values, base, f_x, explain_idx, x_cols, target_names)
    np.savez_compressed(
        run_dir / "shap_values.npz",
        values=values, base_values=base, data=x_std[explain_idx],
        explain_idx=explain_idx, background_idx=background,
        feature_names=np.array(x_cols, dtype=object),
        target_names=np.array(target_names, dtype=object),
    )
    _write_figures(config, run_dir, values, base, x_std[explain_idx], x_cols, target_names)
    logger.info(f"SHAP artifacts written to {run_dir}")


def _write_tables(run_dir, values, base, f_x, explain_idx, x_cols, target_names) -> None:
    # The requested headline statistic: np.abs(shap_values.values).mean(axis=0).
    mean_abs = np.abs(values).mean(axis=0)          # [n_features, n_targets]
    mean_signed = values.mean(axis=0)
    rows = []
    for t_idx, t_name in enumerate(target_names):
        order = np.argsort(-mean_abs[:, t_idx])
        rank = {col: r + 1 for r, col in enumerate(order)}
        for col, name in enumerate(x_cols):
            rows.append({
                "target": t_name, "x_col": name,
                "mean_abs_shap": float(mean_abs[col, t_idx]),
                "mean_signed_shap": float(mean_signed[col, t_idx]),
                "rank": rank[col],
            })
    pd.DataFrame(rows).sort_values(["target", "rank"]).to_csv(
        run_dir / "shap_global_importance.csv", index=False
    )

    local = [
        {"row": int(explain_idx[i]), "target": t_name, "x_col": name,
         "shap_value": float(values[i, col, t_idx])}
        for i, _ in enumerate(explain_idx)
        for t_idx, t_name in enumerate(target_names)
        for col, name in enumerate(x_cols)
    ]
    pd.DataFrame(local).to_csv(run_dir / "shap_local.csv", index=False)

    # Local accuracy: base + sum(shap) must reproduce the model output. This is SHAP's
    # analogue of the completeness property, and the primary correctness test.
    total = base + values.sum(axis=1)
    f_x = np.asarray(f_x).reshape(total.shape)
    add = [
        {"row": int(explain_idx[i]), "target": target_names[t],
         "base_value": float(base[i, t]), "sum_shap": float(values[i, :, t].sum()),
         "f_x": float(f_x[i, t]), "abs_error": float(abs(total[i, t] - f_x[i, t])),
         "rel_error": float(abs(total[i, t] - f_x[i, t]) / max(abs(f_x[i, t]), 1e-12))}
        for i in range(values.shape[0])
        for t in range(values.shape[2])
    ]
    frame = pd.DataFrame(add)
    frame.to_csv(run_dir / "additivity.csv", index=False)
    logger.info(
        f"Additivity: median abs error {frame.abs_error.median():.3e}, "
        f"max {frame.abs_error.max():.3e}"
    )
    # Relative above 1, since targets come in their own units (LMH/bar, %, log units).
    if (frame.abs_error > 1e-3 * np.maximum(1.0, frame.f_x.abs())).any():
        logger.warning(
            "SHAP values do not sum to the model output. The prediction function and the "
            "explainer likely disagree about the feature space -- check standardization."
        )


def _write_figures(config, run_dir, values, base, data, x_cols, target_names) -> None:
    fig_dir = run_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tag = f"{config.problem} / {config.model_type}"

    set_visualization_style()
    save_figure(
        plot_attribution_heatmap(
            np.abs(values).mean(axis=0).T, target_names, x_cols,
            title=f"Mean |SHAP| per tabular feature -- {tag}",
            diverging=False, cbar_label="Mean |SHAP| (target units)",
        ),
        fig_dir / "shap_global_importance.png",
    )

    for t_idx, t_name in enumerate(target_names):
        sliced = shap.Explanation(
            values=values[:, :, t_idx], base_values=base[:, t_idx], data=data,
            feature_names=x_cols,
        )
        for kind, plotter in (("beeswarm", shap.plots.beeswarm), ("bar", shap.plots.bar)):
            try:
                plotter(sliced, show=False)
            except Exception as exc:  # noqa: BLE001 - a plotting failure must not lose data
                logger.warning(f"{kind} plot failed for {t_name}: {exc}")
                plt.close("all")
                continue
            fig = plt.gcf()
            fig.suptitle(f"{tag} -- {t_name}", fontsize=9)
            save_figure(fig, fig_dir / f"{kind}_{_slug(t_name)}.png")


def _slug(text: str) -> str:
    return "".join(c if c.isalnum() else "_" for c in text).strip("_")[:40]


if __name__ == "__main__":
    main()
