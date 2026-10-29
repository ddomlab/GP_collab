"""GP-SHAP: stochastic Shapley values with credible intervals, for the GP-headed models.

`scripts/analyze_shap.py` explains `mlp_gp` with permutation SHAP and gets point estimates;
`scripts/analyze_attributions.py` explains the graph models with Integrated Gradients and
also gets point estimates. Both discard the predictive posterior these models carry --
`analyze_shap.py` literally throws away the std it is handed.

This script keeps it. Following Chau et al. (arXiv:2305.15167), the explanations of a GP
are themselves multivariate Gaussian, so every attribution comes with a credible interval
and every pair of features with a covariance.

Scope, and its limits:

* Players are the tabular `x_cols` for every model type. For the graph models (`gnn_gp`,
  `gnn_gpsn`) the molecular structure is not an attributed player -- structure remains the
  job of Integrated Gradients. Nor is the explained row's own molecule held fixed: every
  payoff is built from the reference rows' posterior, so even the full coalition sees the
  reference rows' molecules, weighted by how close their conditions are. The conditional
  mean embedding marginalizes the removed conditions over the training distribution, which
  drags their associated molecules along; a graph-model value therefore reads as
  "contribution of this condition, averaged over the joint distribution of the other
  conditions and their molecules". `scripts/analyze_shap.py` keeps each row's molecule.
* The uncertainty is epistemic, from the GP posterior. A wide interval says the model is
  unsure what the explanation *is*, not that the feature is unimportant.
* `--ridge` stabilizes the conditional-mean-embedding inversion and has no principled
  value. If the intervals move materially with it, say so rather than reporting one run.

Requires a model that held rows out, such as a CV fold model, so the
conditional-mean-embedding reference set and the explained rows are disjoint.

Examples:
    python scripts/analyze_gpshap.py --problem dataset4_pervaporation --model_type mlp_gp
    python scripts/analyze_gpshap.py --problem dataset6_reaction --model_type gnn_gp \
        --bayes_coalitions 8
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import draccus
import numpy as np
import pandas as pd
from safetensors.torch import load_file
from scripts.train_gnn_regression import (
    TABULAR_MODELS,
    build_graph_dataset,
    build_model,
    make_fold_loaders,
    warmup_lazy_modules,
)
import torch
import yaml

from polykernels import gpshap as gps
from polykernels import paths_and_constants, regression_data
from polykernels.datasets import overwrite_dataset_globals
from polykernels.models.gp import MultiTargetGPHead
from polykernels.molgraph import GraphType
from polykernels.utils import get_logger, set_seed
from polykernels.vis import (
    plot_attribution_heatmap,
    plot_interval_bars,
    save_figure,
    set_visualization_style,
)

_P = paths_and_constants.Paths()
logger = get_logger()

# The *_gpsn variants carry the same GP head on a spectral-normalized extractor.
MODEL_TYPES = ("mlp_gp", "gnn_gp", "mlp_gpsn", "gnn_gpsn")
# 95% credible interval on a Gaussian explanation.
_Z95 = 1.959963985


@dataclass
class GPShapConfig:
    """Configuration for a GP-SHAP run."""

    problem: str = "dataset4_pervaporation"
    model_type: Literal["mlp_gp", "gnn_gp", "mlp_gpsn", "gnn_gpsn"] = "mlp_gp"
    results_dir: Path = _P.results_dir / "gpshap"
    graph: GraphType = GraphType.CUSTOM
    # Reference rows for the conditional mean embedding, chosen spread out over the
    # training split. This is the cost driver: each coalition solves an n_reference system.
    n_reference: int = 128
    # 0 enumerates every coalition, which drives the BayesGP-SHAP estimation term to zero.
    n_coalitions: int = 0
    # >0 additionally runs BayesGP-SHAP at that coalition budget, reproducing the paper's
    # Figure 1 ablation: subsampling should widen the intervals, exact enumeration not.
    bayes_coalitions: int = 0
    max_explain_rows: int = 64
    # Ridge on the CME inversion (lambda in Proposition 8).
    ridge: float = 1e-3
    n_exemplars: int = 2
    ell0: float = 1.0
    sigma0_sq: float = 1e-4
    seed: int = 6
    coalition_chunk: int = 256


def load_gp_model(
    config: GPShapConfig, run_dir: Path, hparams: dict, output_dim: int,
    x_std: np.ndarray, y_std: np.ndarray, fit_idx: np.ndarray,
    graph_dataset, device: torch.device,
):
    """Rebuild and reload a fitted GP-headed model.

    Order is mandatory: every extractor starts with LazyLinear, so the parameter shapes must
    be materialized by a forward pass before `load_state_dict` can populate them. The same
    warmup re-registers a `*_gpsn` extractor's spectral-norm parametrization, which its
    saved weights are keyed under.
    """
    train_loader, _ = make_fold_loaders(
        config.model_type, graph_dataset, x_std, y_std, fit_idx, fit_idx,
        hparams["optimizer"]["batch_size"],
    )
    model = build_model(config.model_type, hparams, output_dim, device)
    warmup_lazy_modules(model, train_loader, device)
    model.load_state_dict(load_file(run_dir / _P.model_file, device=str(device)))
    model.eval()
    return model


def embed_rows(config, model, x_std, graph_dataset, rows, device) -> torch.Tensor:
    """Feature-extractor embeddings for a set of dataset rows, `[len(rows), emb_dim]`.

    The GP kernel lives on these embeddings, which is why they cannot themselves be Shapley
    players -- but the posterior we need is indexed by row, so raw features and embeddings
    line up one to one.
    """
    with torch.no_grad():
        if config.model_type in TABULAR_MODELS:
            batch = torch.from_numpy(np.ascontiguousarray(x_std[rows])).float().to(device)
            return model.mlp.embed(batch)
        import torch_geometric as pyg

        graphs = [graph_dataset.graphs[int(i)] for i in rows]
        return model.gnn.embed(pyg.data.Batch.from_data_list(graphs).to(device))


def gp_posterior(model, features: torch.Tensor) -> list[tuple[torch.Tensor, torch.Tensor]]:
    """Per-target (mean, full covariance) at the given embeddings."""
    head: MultiTargetGPHead = model.gp_head
    with torch.no_grad():
        return head.posterior_at(features)


def run_gpshap(
    config: GPShapConfig, x_ref: torch.Tensor, x_query: torch.Tensor,
    posterior: list[tuple[torch.Tensor, torch.Tensor]], n_coalitions: int,
    rng: np.random.Generator,
) -> dict:
    """GP-SHAP for every target at one coalition budget.

    Returns mean `[n_query, d, n_targets]`, per-query covariance `[n_query, d, d, T]`, the
    payoff matrix and the coalition design, which BayesGP-SHAP needs.
    """
    d = x_ref.shape[1]
    z, w = gps.coalition_matrix(d, n_coalitions, rng=rng)
    bandwidth = MultiTargetGPHead._median_pairwise_distance(x_ref.float())
    logger.info(
        f"{z.shape[0]} coalitions over {d} features, CME bandwidth {bandwidth:.3f}, "
        f"{x_ref.shape[0]} reference rows, {x_query.shape[0]} explained rows."
    )
    b = gps.conditional_mean_projections(
        x_ref, x_query, z, bandwidth=bandwidth, ridge=config.ridge,
        coalition_chunk=config.coalition_chunk,
    )

    means, covs, payoff_list = [], [], []
    for m_ref, k_ref in posterior:
        phi, r, payoffs = gps.gpshap(
            m_ref.double().cpu(), k_ref.double().cpu(), b, z, w
        )
        means.append(phi.T)                          # [n_query, d]
        covs.append(gps.per_query_covariance(r))     # [n_query, d, d]
        payoff_list.append(payoffs)
    return {
        "mean": torch.stack(means, dim=-1).numpy(),
        "cov": torch.stack(covs, dim=-1).numpy(),
        "payoffs": payoff_list,
        "z": z,
        "w": w,
    }


def _slug(text: str) -> str:
    return "".join(c if c.isalnum() else "_" for c in text).strip("_")[:40]


def write_tables(config, run_dir, result, x_cols, target_names, explain_idx,
                 ref_idx) -> pd.DataFrame:
    """Write the value/importance/efficiency tables; return the efficiency frame.

    Everything arriving here is already in physical target units -- `main` rescales the
    means, covariances and payoffs together so unit conversion happens exactly once.
    """
    mean, cov = result["mean"], result["cov"]
    sd = np.sqrt(np.clip(np.einsum("nddt->ndt", cov), 0.0, None))

    np.savez_compressed(
        run_dir / "gpshap_values.npz",
        mean=mean, cov_per_query=cov, sd=sd,
        explain_idx=explain_idx, reference_idx=ref_idx,
        feature_names=np.array(x_cols, dtype=object),
        target_names=np.array(target_names, dtype=object),
    )

    rows = []
    for t_idx, t_name in enumerate(target_names):
        for f_idx, col in enumerate(x_cols):
            rows.append({
                "target": t_name, "x_col": col,
                # Mean of |SSV| accounts for explanation uncertainty; |mean SSV| does not.
                # The paper's Figure 2b contrasts them, and they can rank features
                # differently.
                "mean_abs_ssv": float(np.abs(mean[:, f_idx, t_idx]).mean()),
                "std_abs_ssv": float(sd[:, f_idx, t_idx].mean()),
                "abs_mean_ssv": float(abs(mean[:, f_idx, t_idx].mean())),
            })
    importance = pd.DataFrame(rows)
    importance["rank"] = (
        importance.groupby("target").mean_abs_ssv.rank(ascending=False, method="min")
        .astype(int)
    )
    importance.sort_values(["target", "rank"]).to_csv(
        run_dir / "gpshap_global_importance.csv", index=False
    )

    local = [
        {"row": int(explain_idx[i]), "target": t_name, "x_col": col,
         "ssv_mean": float(mean[i, f, t]), "ssv_std": float(sd[i, f, t])}
        for i in range(mean.shape[0])
        for t, t_name in enumerate(target_names)
        for f, col in enumerate(x_cols)
    ]
    pd.DataFrame(local).to_csv(run_dir / "gpshap_local.csv", index=False)

    # Efficiency: the values must sum to nu(full) - nu(empty). This is the axiom that
    # catches a mis-weighted coalition design, so it is checked per row rather than assumed.
    eff = []
    for t_idx, t_name in enumerate(target_names):
        payoffs = result["payoffs"][t_idx].numpy()
        span = payoffs[-1] - payoffs[0]
        total = mean[:, :, t_idx].sum(axis=1)
        for i in range(mean.shape[0]):
            eff.append({
                "row": int(explain_idx[i]), "target": t_name,
                "sum_ssv": float(total[i]), "payoff_span": float(span[i]),
                "abs_error": float(abs(total[i] - span[i])),
            })
    eff_df = pd.DataFrame(eff)
    eff_df.to_csv(run_dir / "efficiency.csv", index=False)
    return eff_df


def write_figures(config, run_dir, result, x_cols, target_names, explain_idx,
                  y_true, y_pred, bayes_sd=None) -> None:
    set_visualization_style()
    fig_dir = run_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tag = f"{config.problem} / {config.model_type}"
    mean, cov = result["mean"], result["cov"]
    sd = np.sqrt(np.clip(np.einsum("nddt->ndt", cov), 0.0, None))

    for t_idx, t_name in enumerate(target_names):
        slug = _slug(t_name)
        err = np.abs(y_pred[:, t_idx] - y_true[:, t_idx])
        picks = {"best": int(np.argmin(err))}
        if config.n_exemplars > 1:
            picks["worst"] = int(np.argmax(err))

        for label, i in picks.items():
            save_figure(
                plot_interval_bars(
                    mean[i, :, t_idx], _Z95 * sd[i, :, t_idx], x_cols,
                    title=(f"Stochastic Shapley values -- {tag}\n{t_name}, row "
                           f"{int(explain_idx[i])} ({label}), 95% credible interval"),
                    xlabel="Stochastic Shapley value (target units)",
                ),
                fig_dir / f"local_{slug}_{label}.png",
            )
            # Correlation between explanations, which the tractable covariance makes
            # available and ordinary SHAP cannot express at all.
            block = cov[i, :, :, t_idx]
            scale = np.sqrt(np.clip(np.diag(block), 1e-30, None))
            save_figure(
                plot_attribution_heatmap(
                    block / np.outer(scale, scale), x_cols, x_cols,
                    title=f"SSV correlation -- {tag}\n{t_name}, row {int(explain_idx[i])}",
                    cbar_label="Correlation", annot=len(x_cols) <= 8,
                ),
                fig_dir / f"ssv_correlation_{slug}_{label}.png",
            )

        save_figure(
            plot_interval_bars(
                np.abs(mean[:, :, t_idx]).mean(axis=0),
                sd[:, :, t_idx].mean(axis=0), x_cols,
                title=f"Global importance -- {tag}\n{t_name}",
                xlabel="Mean |stochastic Shapley value| (target units)",
                secondary=np.abs(mean[:, :, t_idx].mean(axis=0)),
                secondary_label="|mean SSV| (ignores explanation uncertainty)",
            ),
            fig_dir / f"global_importance_{slug}.png",
        )

        if bayes_sd is not None:
            save_figure(
                plot_interval_bars(
                    mean[picks["best"], :, t_idx],
                    _Z95 * bayes_sd[picks["best"], :, t_idx], x_cols,
                    title=(f"BayesGP-SHAP vs GP-SHAP -- {tag}\n{t_name}, row "
                           f"{int(explain_idx[picks['best']])}: bars are the subsampled "
                           f"interval, markers the exact one"),
                    xlabel="Stochastic Shapley value (target units)",
                    secondary=mean[picks["best"], :, t_idx]
                    + _Z95 * sd[picks["best"], :, t_idx],
                    secondary_label="GP-SHAP upper bound, exact coalitions",
                ),
                fig_dir / f"bayes_ablation_{slug}.png",
            )


@draccus.wrap()
def main(config: GPShapConfig) -> None:
    set_seed(config.seed)
    if config.model_type not in MODEL_TYPES:
        raise ValueError(f"model_type must be one of {MODEL_TYPES}.")
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
    x_std = (x_features - stats["x_mean"]) / stats["x_std"]
    y_std = (y - stats["target_mean"]) / stats["target_std"]

    if len(held_out) == 0:
        logger.warning(
            "Model was fitted on all rows, so the reference set and the explained rows "
            "overlap and these values describe an in-sample fit. Explain a CV fold model "
            "instead."
        )
        held_out = np.arange(len(y))
    explain_idx = (
        held_out[: config.max_explain_rows] if config.max_explain_rows else held_out
    )

    graph_dataset = None
    if config.model_type not in TABULAR_MODELS:
        logger.info(f"Featurizing {len(smiles)} SMILES...")
        graph_dataset = build_graph_dataset(
            smiles, x_std, y, graph_type=config.graph, use_3d=spec.use_3d
        )
        if graph_dataset.graphs[0].u.shape[1] > 0:
            overwrite_dataset_globals(graph_dataset, np.ascontiguousarray(x_std))

    model = load_gp_model(
        config, run_dir, hparams, output_dim, x_std, y_std, fit_idx, graph_dataset, device
    )

    # Reference rows spread over the training split, reusing the same farthest-point
    # traversal the GP head uses to place inducing points.
    pool = torch.from_numpy(np.ascontiguousarray(x_std[fit_idx])).float()
    n_ref = min(config.n_reference, len(fit_idx))
    if len(fit_idx) > n_ref:
        chosen = MultiTargetGPHead._spread_out_sample(pool, n_ref, scale=0.0)
        ref_idx = np.array([
            int(np.argmin(np.abs(x_std[fit_idx] - row.numpy()).sum(axis=1)))
            for row in chosen
        ])
        ref_idx = fit_idx[np.unique(ref_idx)]
    else:
        ref_idx = fit_idx

    posterior = gp_posterior(model, embed_rows(config, model, x_std, graph_dataset, ref_idx, device))
    x_ref = torch.from_numpy(np.ascontiguousarray(x_std[ref_idx])).double()
    x_query = torch.from_numpy(np.ascontiguousarray(x_std[explain_idx])).double()

    rng = np.random.default_rng(config.seed)
    result = run_gpshap(config, x_ref, x_query, posterior, config.n_coalitions, rng)

    # Scale explanations into physical target units; the GP works in standardized space.
    # Means, covariances and payoffs are converted together -- rescaling any of them a
    # second time downstream silently breaks the efficiency check.
    target_scale = np.asarray(stats["target_std"], dtype=float)
    result["mean"] = result["mean"] * target_scale
    result["cov"] = result["cov"] * (target_scale**2)
    result["payoffs"] = [p * target_scale[t] for t, p in enumerate(result["payoffs"])]

    bayes_sd = None
    if config.bayes_coalitions:
        sub = run_gpshap(config, x_ref, x_query, posterior, config.bayes_coalitions, rng)
        extra = []
        for t_idx in range(len(target_names)):
            phi = torch.from_numpy(sub["mean"][:, :, t_idx]).T
            extra.append(gps.bayes_extra_variance(
                sub["payoffs"][t_idx], phi, sub["z"], sub["w"],
                ell0=config.ell0, sigma0_sq=config.sigma0_sq, rng=rng,
            ).T.numpy())
        sub_sd = np.sqrt(np.clip(np.einsum("nddt->ndt", sub["cov"]), 0.0, None))
        bayes_sd = np.sqrt(sub_sd**2 + np.stack(extra, axis=-1)) * target_scale
        logger.info(
            f"BayesGP-SHAP at {sub['z'].shape[0]} coalitions: mean interval half-width "
            f"{(_Z95 * bayes_sd).mean():.4g} vs GP-SHAP exact "
            f"{(_Z95 * np.sqrt(np.einsum('nddt->ndt', result['cov']))).mean():.4g}."
        )

    with torch.no_grad():
        feats = embed_rows(config, model, x_std, graph_dataset, explain_idx, device)
        pred_mean, _ = model.gp_head.predict(feats)
    y_pred = pred_mean.cpu().numpy() * target_scale + np.asarray(stats["target_mean"])

    run_dir.mkdir(parents=True, exist_ok=True)
    eff = write_tables(
        config, run_dir, result, x_cols, target_names, explain_idx, ref_idx
    )
    logger.info(
        f"Efficiency: median |error| {eff.abs_error.median():.3e}, max "
        f"{eff.abs_error.max():.3e} (values must sum to the payoff span)."
    )
    if eff.abs_error.max() > 1e-6:
        logger.warning(
            f"Efficiency violated by up to {eff.abs_error.max():.3e}. The coalition design "
            f"or its weights are wrong -- the Shapley axioms should hold to float precision."
        )
    write_figures(
        config, run_dir, result, x_cols, target_names, explain_idx,
        y[explain_idx], y_pred, bayes_sd,
    )
    logger.info(f"GP-SHAP written to {run_dir}")


if __name__ == "__main__":
    main()
