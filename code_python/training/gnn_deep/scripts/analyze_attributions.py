"""Integrated Gradients attribution for a fitted regression GNN.

Runs IG for every row and every target, with a zero baseline over BOTH the node features
and the graph globals (the tabular `x_cols`). Because `Readout` concatenates the pooled
node states with the global state, the split between those two attribution masses is a
direct measure of how much of a prediction comes from chemistry versus process conditions.

Node-feature columns that are constant across the whole dataset are zeroed in every
reported number -- tables, `ig_raw.npz`, `chemistry_share` and figures. The zero baseline
credits them anyway, yet they cannot distinguish one molecule from another. Completeness is
still checked on the unmasked attributions, and `chemistry_share_all_columns` keeps the
unmasked share for comparison.

Explains a saved model in `<results_dir>/<problem>_<model_type>/`.
`hpc/killarney/attribution_cv.sh` copies each CV fold's model there (from a sweep run with
`--save_models true`) and runs this script on it.

Examples:
    # One configuration
    python scripts/analyze_attributions.py --problem dataset3_membrane --model_type gnn_mlp

    # Cross-dataset summary, after all configurations have been analyzed
    python scripts/analyze_attributions.py --summarize True

    # The same figure averaged over training seeds (after attribution_consistency.py)
    python scripts/analyze_attributions.py --summarize True \
        --consistency_csv results/attribution_seeds/consistency.csv
"""

import dataclasses
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import draccus
import numpy as np
import pandas as pd
from safetensors.torch import load_file
from scripts.train_gnn_regression import (
    ModelType,
    build_graph_dataset,
    build_model,
    make_fold_loaders,
    warmup_lazy_modules,
)
import torch
import yaml

from polykernels import paths_and_constants, regression_data
from polykernels.attribution import (
    constant_columns,
    feature_blocks,
    integrated_gradients,
    node_feature_names,
    scalar_output_fn,
)
from polykernels.datasets import overwrite_dataset_globals
from polykernels.molgraph import GraphType, get_template
from polykernels.utils import get_logger, set_seed
from polykernels.vis import (
    draw_molecule_with_scores,
    plot_attribution_heatmap,
    save_figure,
    set_visualization_style,
)

_P = paths_and_constants.Paths()
logger = get_logger()

PROBLEMS = list(regression_data.REGRESSION_PROBLEMS)
GRAPH_MODELS: tuple[ModelType, ...] = ("gnn_mlp", "gnn_gp", "gnn_gpsn")
# Models in the chemistry-vs-process summary. gnn_gpsn is left out: its zero-baseline IG
# gives the x_cols ~1% of the credit gnn_gp's does, yet swapping a row's x_cols for another
# row's moves both models' predictions by the same amount, so its chemistry share measures
# the baseline path through the unnormalized encoder rather than reliance on chemistry.
SUMMARY_MODELS: tuple[ModelType, ...] = ("gnn_mlp", "gnn_gp")

# What each problem's SMILES column actually enumerates, aligned with
# `DatasetSpec.smiles_cols`. Only the noun lives here -- the counts are read back from the
# source CSV in `_scale_label`, so a label can never drift from the data it describes.
ENTITY_NOUNS: dict[str, tuple[str, ...]] = {
    "dataset1_pce": ("donors", "acceptors"),
    "dataset2_flux": ("polymers",),
    "dataset2_separation": ("polymers",),
    "dataset3_membrane": ("polymers",),
    "dataset4_pervaporation": ("polymers",),
    "dataset5_rg": ("monomers",),
    "dataset6_reaction": ("substrates",),
}


@dataclass
class AttributionConfig:
    """Configuration for an Integrated Gradients attribution run."""

    problem: str = "dataset3_membrane"
    model_type: Literal["gnn_mlp", "gnn_gp", "gnn_gpsn"] = "gnn_mlp"
    graph: GraphType = GraphType.CUSTOM
    results_dir: Path = _P.results_dir / "attribution_globals"
    # Trapezoid-rule steps along the baseline->input path; completeness error falls
    # as this rises. 2048 is what the attribution_globals sweep used on every
    # problem, so a default run reproduces it instead of mixing in a coarser one.
    steps: int = 2048
    # Molecules drawn per target: the best- and worst-predicted rows.
    n_exemplars: int = 2
    seed: int = 6
    # Skip attribution; rebuild the cross-dataset figure from existing per-run CSVs.
    summarize: bool = False
    # With summarize: plot the chemistry share averaged over training seeds, read from
    # this table by scripts/attribution_consistency.py, instead of one sweep's
    # completeness.csv files. The figure is written beside the table.
    consistency_csv: str = ""
    # With consistency_csv: "folds" plots the mean and SD over the CV fold models instead
    # (attribution_consistency.py --folds), to chemistry_vs_process_folds.png.
    consistency_over: Literal["seeds", "folds"] = "seeds"
    # Attribute only the rows the model never trained on (`held_out_idx`), as
    # attribution_cv.sh does for each CV fold model, so IG explains the same network and the
    # same rows as SHAP / GP-SHAP.
    held_out_only: bool = False
    # Attribute only the first N of the selected rows; 0 means all of them. With
    # held_out_only, 64 picks the rows GP-SHAP explained (its max_explain_rows).
    # Otherwise it is a smoke test: IG cost scales linearly with rows.
    max_rows: int = 0
    # Atoms shown in the per-atom heatmap, ranked by total |IG|. A polymer repeat unit can
    # run to a couple of hundred atoms, and one row each makes a figure metres tall whose
    # labels are illegible once a viewer fits it to the screen. Full data is in
    # node_attribution.csv; 0 plots every atom.
    max_atoms_plotted: int = 40


def load_fitted_model(
    model_type: ModelType, hparams: dict, output_dim: int, run_dir: Path,
    train_loader, device: torch.device,
) -> torch.nn.Module:
    """Rebuild the architecture, materialize its lazy layers, then load trained weights.

    Order matters: the GNN uses `LazyLinear`, so `warmup_lazy_modules` must run a forward
    pass to create the parameter shapes before `load_state_dict` can populate them.
    """
    model = build_model(model_type, hparams, output_dim, device)
    warmup_lazy_modules(model, train_loader, device)
    model.load_state_dict(load_file(run_dir / _P.model_file, device=str(device)))
    model.eval()   # required: LazyBatchNorm1d must not pool statistics across IG steps
    return model


def attribute_one(config: AttributionConfig) -> None:
    """Run IG for every (row, target) of one fitted model and write all artifacts."""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    run_dir = config.results_dir / f"{config.problem}_{config.model_type}"
    if not (run_dir / _P.model_file).exists():
        raise FileNotFoundError(
            f"No weights at {run_dir / _P.model_file}. Copy a saved CV fold model there; "
            f"hpc/killarney/attribution_cv.sh does this for a sweep run with --save_models."
        )

    with open(run_dir / _P.hparams_file) as f:
        record = yaml.safe_load(f)
    hparams, output_dim = record["hparams"], record["output_dim"]
    stats = np.load(run_dir / "standardization.npz", allow_pickle=True)
    target_names = [str(t) for t in stats["target_names"]]
    x_cols = [str(c) for c in stats["x_cols"]]

    spec = regression_data.REGRESSION_PROBLEMS[config.problem]
    smiles, x_features, y, _, _ = regression_data.load_regression_data(spec)
    # Reuse the fitted standardization rather than recomputing it, so the model sees
    # exactly the inputs it was trained on.
    x_std = (x_features - stats["x_mean"]) / stats["x_std"]
    y_std = (y - stats["target_mean"]) / stats["target_std"]

    logger.info(f"Featurizing {len(smiles)} SMILES for {run_dir.name}...")
    graph_dataset = build_graph_dataset(
        smiles, x_std, y, graph_type=config.graph, use_3d=spec.use_3d
    )
    if graph_dataset.graphs[0].u.shape[1] > 0:
        overwrite_dataset_globals(graph_dataset, np.ascontiguousarray(x_std))
    # Found before any row selection: "constant" means constant across the whole dataset,
    # so a subset run masks the same columns as a full one.
    dead_cols = constant_columns(torch.cat([g.x for g in graph_dataset.graphs], dim=0))

    # Dataset row indices to attribute. Outputs report these, not positions in the
    # selection, so rows line up with SHAP's and with a full run's.
    row_idx = np.arange(len(y))
    if config.held_out_only:
        row_idx = np.asarray(stats["held_out_idx"], dtype=int)
        if len(row_idx) == 0:
            raise ValueError(
                f"{run_dir} was fit on every row, so it has no held-out rows. "
                "held_out_only needs a model that held rows out, such as a CV fold model."
            )
    if config.max_rows:
        row_idx = row_idx[: config.max_rows]
    if len(row_idx) < len(y):
        logger.info(
            f"Attributing {len(row_idx)} of {len(y)} rows "
            f"(held_out_only={config.held_out_only}, max_rows={config.max_rows})."
        )
        graph_dataset.graphs = [graph_dataset.graphs[i] for i in row_idx]
        smiles, y = smiles[row_idx], y[row_idx]
        x_std, y_std = x_std[row_idx], y_std[row_idx]
    all_idx = np.arange(len(y))
    train_loader, _ = make_fold_loaders(
        config.model_type, graph_dataset, x_std, y_std, all_idx, all_idx,
        hparams["optimizer"]["batch_size"],
    )
    model = load_fitted_model(
        config.model_type, hparams, output_dim, run_dir, train_loader, device
    )

    template = get_template(GraphType.CUSTOM, global_dim=x_features.shape[1])
    feat_names = node_feature_names(template)
    blocks = feature_blocks(template)
    live = np.ones(len(feat_names), dtype=bool)
    live[dead_cols] = False
    logger.info(
        f"{len(feat_names)} node features ({len(dead_cols)} constant, zeroed in every "
        f"reported attribution), {len(x_cols)} globals, {output_dim} target(s), "
        f"{len(y)} rows."
    )

    n_rows = len(y)
    feat_sum = np.zeros((n_rows, output_dim, len(feat_names)), dtype=np.float32)
    glob_ig = np.zeros((n_rows, output_dim, len(x_cols)), dtype=np.float32)
    f_in = np.zeros((n_rows, output_dim), dtype=np.float32)
    f_base = np.zeros((n_rows, output_dim), dtype=np.float32)
    completeness, node_rows, keep = [], [], {}

    logger.info(f"Running IG with {config.steps} steps...")
    for t_idx, t_name in enumerate(target_names):
        fn = scalar_output_fn(model, config.model_type, t_idx)
        # `pos` indexes the selected rows; `row` is the dataset row it came from.
        for pos, row in enumerate(row_idx.tolist()):
            data = graph_dataset.graphs[pos]
            result = integrated_gradients(fn, data, device, steps=config.steps)
            # A constant column holds the same value on every atom of every molecule, so
            # it cannot separate any two rows -- but the zero baseline still walks it from
            # 0 to that value and credits it, often with most of the node attribution.
            # Everything reported is computed on `shown`. Completeness stays on the raw
            # `result`: IG's sum rule only holds over every input.
            shown = dataclasses.replace(result, ig_x=result.ig_x * live)

            feat_sum[pos, t_idx] = shown.ig_x.sum(axis=0)
            glob_ig[pos, t_idx] = result.ig_u
            f_in[pos, t_idx], f_base[pos, t_idx] = result.f_input, result.f_baseline
            completeness.append({
                "row": row, "smiles": smiles[pos], "target": t_name,
                "y_true": float(y[pos, t_idx]),
                "y_pred": float(result.f_input * stats["target_std"][t_idx]
                                + stats["target_mean"][t_idx]),
                "f_input": result.f_input, "f_baseline": result.f_baseline,
                "sum_ig": result.total,
                "abs_error": result.completeness_abs_error,
                "rel_error": result.completeness_rel_error,
                "chemistry_share": shown.chemistry_share(),
                "chemistry_share_all_columns": result.chemistry_share(),
                "globals_connected": result.globals_connected,
            })
            signed, magnitude = shown.signed_node_scores(), shown.magnitude_node_scores()
            for atom in range(len(signed)):
                node_rows.append({
                    "row": row, "smiles": smiles[pos], "target": t_name,
                    "atom_idx": atom, "signed": float(signed[atom]),
                    "magnitude": float(magnitude[atom]),
                })
            keep[(row, t_idx)] = shown.ig_x
        logger.info(f"  {t_name}: done ({n_rows} rows)")

    if not all(rec["globals_connected"] for rec in completeness):
        logger.warning(
            "The model output does not depend on the tabular x_cols at all, so their "
            "attribution is identically zero. This is the `seed_globals=True` path in "
            "GraphEncoder.forward, which synthesizes the global state from node "
            "embeddings and discards `u`. Refit with model.seed_globals=false to give "
            "the graph model access to the process conditions."
        )

    _write_tables(run_dir, completeness, node_rows, feat_sum, glob_ig, feat_names,
                  x_cols, target_names, dead_cols)
    np.savez_compressed(
        run_dir / "ig_raw.npz",
        feature_sum=feat_sum, global_ig=glob_ig, f_input=f_in, f_baseline=f_base,
        feature_names=np.array(feat_names, dtype=object),
        # `feature_sum` already has these zeroed; recorded so readers can tell a masked
        # file from one written before masking existed.
        constant_columns=dead_cols,
        x_cols=np.array(x_cols, dtype=object),
        target_names=np.array(target_names, dtype=object),
        smiles=np.array(smiles, dtype=object),
        # Dataset row of each entry along the first axis of the arrays above.
        row_idx=row_idx,
    )
    _write_figures(config, run_dir, feat_sum, glob_ig, blocks, live, x_cols, target_names,
                   keep, pd.DataFrame(completeness))
    logger.info(f"Attribution artifacts written to {run_dir}")


def _write_tables(run_dir, completeness, node_rows, feat_sum, glob_ig, feat_names,
                  x_cols, target_names, dead_cols) -> None:
    comp = pd.DataFrame(completeness)
    comp.to_csv(run_dir / "completeness.csv", index=False)
    pd.DataFrame(node_rows).to_csv(run_dir / "node_attribution.csv", index=False)

    rows = []
    for t_idx, t_name in enumerate(target_names):
        for col, name in enumerate(feat_names):
            rows.append({
                "target": t_name, "feature": name, "column": col,
                "is_constant": col in set(dead_cols.tolist()),
                "mean_signed_ig": float(feat_sum[:, t_idx, col].mean()),
                "mean_abs_ig": float(np.abs(feat_sum[:, t_idx, col]).mean()),
            })
    pd.DataFrame(rows).to_csv(run_dir / "feature_attribution.csv", index=False)

    rows = []
    for t_idx, t_name in enumerate(target_names):
        for col, name in enumerate(x_cols):
            rows.append({
                "target": t_name, "x_col": name,
                "mean_signed_ig": float(glob_ig[:, t_idx, col].mean()),
                "mean_abs_ig": float(np.abs(glob_ig[:, t_idx, col]).mean()),
            })
    pd.DataFrame(rows).to_csv(run_dir / "global_attribution.csv", index=False)

    median_rel = comp["rel_error"].median()
    logger.info(
        f"Completeness: median relative error {median_rel:.2%}, "
        f"max {comp['rel_error'].max():.2%}"
    )
    if median_rel > 0.05:
        logger.warning(
            f"Median completeness error {median_rel:.2%} exceeds 5%. Attributions do not "
            f"account for the prediction; raise --steps."
        )


def _write_figures(config, run_dir, feat_sum, glob_ig, blocks, live, x_cols,
                   target_names, keep, comp) -> None:
    set_visualization_style()
    fig_dir = run_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tag = f"{config.problem} / {config.model_type}"

    # A block with no varying column carries no information in this dataset, so it is
    # dropped rather than drawn as zero, which would read as "tested and unimportant".
    n_dropped = sum(not live[a:b].any() for _, a, b in blocks)
    blocks = [(n, a, b) for n, a, b in blocks if live[a:b].any()]
    dropped_note = (f"\n{n_dropped} block(s) constant in this dataset omitted"
                    if n_dropped else "")
    # Mean over the block's varying columns only; averaging in the zeroed constant ones
    # would shrink a block in proportion to how many of its one-hot values never occur.
    block_matrix = np.stack([
        [np.abs(feat_sum[:, t, a:b][:, live[a:b]]).mean() for _, a, b in blocks]
        for t in range(len(target_names))
    ])
    save_figure(
        plot_attribution_heatmap(
            block_matrix, target_names, [n for n, _, _ in blocks],
            title=f"Mean |IG| per node-feature block -- {tag}{dropped_note}",
            diverging=False, cbar_label="Mean |IG| (standardized target units)",
        ),
        fig_dir / "feature_blocks.png",
    )

    if np.abs(glob_ig).max() > 0:
        save_figure(
            plot_attribution_heatmap(
                np.abs(glob_ig).mean(axis=0), target_names, x_cols,
                title=f"Mean |IG| per tabular feature (x_cols) -- {tag}",
                diverging=False, cbar_label="Mean |IG| per SD of feature",
            ),
            fig_dir / "global_features.png",
        )
    else:
        logger.info("Skipping global-feature figure: x_cols attribution is all zero.")

    for t_idx, t_name in enumerate(target_names):
        sub = comp[comp.target == t_name].copy()
        sub["err"] = (sub.y_pred - sub.y_true).abs()
        picks = pd.concat([sub.nsmallest(1, "err"), sub.nlargest(1, "err")])
        for label, (_, rec) in zip(("best", "worst"), picks.iterrows()):
            if len(picks) < config.n_exemplars and label == "worst":
                continue
            row = int(rec["row"])
            ig_x = keep[(row, t_idx)]
            safe = f"{_slug(t_name)}_{label}"
            title = (f"{tag}\n{t_name} | pred {rec.y_pred:.3g} vs true {rec.y_true:.3g}")
            for kind, scores, cmap, signed in (
                ("signed", ig_x.sum(axis=1), "coolwarm_r", True),
                ("magnitude", np.abs(ig_x).sum(axis=1), "Greens", False),
            ):
                mol_fig, mol_svg = draw_molecule_with_scores(
                    rec["smiles"], scores, cmap_name=cmap, signed=signed,
                    title=f"{kind.capitalize()} node attribution\n{title}",
                )
                save_figure(
                    mol_fig, fig_dir / f"mol_{safe}_{kind}.png", vector_overlay=mol_svg
                )
            atom_block = np.stack([ig_x[:, a:b].sum(axis=1) for _, a, b in blocks], axis=1)
            order = np.arange(atom_block.shape[0])
            subtitle = ""
            if 0 < config.max_atoms_plotted < len(order):
                # Rank by total |IG| so the atoms that carry the prediction stay visible.
                # Labels keep the original atom index, which matches the indices RDKit
                # writes on the molecule figures above.
                order = np.argsort(np.abs(atom_block).sum(axis=1))[::-1]
                order = order[: config.max_atoms_plotted]
                subtitle = (
                    f"\ntop {len(order)} of {atom_block.shape[0]} atoms by total |IG|"
                )
            save_figure(
                plot_attribution_heatmap(
                    atom_block[order], [str(i) for i in order],
                    [n for n, _, _ in blocks],
                    title=f"Per-atom signed IG by feature block -- {title}{subtitle}",
                    cbar_label="Signed IG",
                ),
                fig_dir / f"atoms_{safe}.png",
            )


def _slug(text: str) -> str:
    return "".join(c if c.isalnum() else "_" for c in text).strip("_")[:40]


def _scale_label(problem: str) -> str:
    """`n = 558 · 143 donors + 261 acceptors` -- the dataset's size and chemical diversity.

    A chemistry share is only readable against how much chemistry the model was shown:
    0.58 over 143 donors and 261 acceptors means something very different from 0.36 over
    22 monomers, where the graph branch can only ever distinguish 22 things.

    Counts are recomputed from the CSV under the same NaN-dropping rule as
    `load_regression_data`, so they describe the rows actually trained on.
    """
    spec = regression_data.REGRESSION_PROBLEMS[problem]
    df = pd.read_csv(spec.csv_path).dropna(
        subset=list(spec.smiles_cols) + list(spec.y_cols)
    )
    parts = [
        f"{df[col].nunique()} {noun}"
        for col, noun in zip(spec.smiles_cols, ENTITY_NOUNS[problem], strict=True)
    ]
    return f"n = {len(df)} · {' + '.join(parts)}"


def plot_chemistry_share(
    summary: pd.DataFrame, out_path: Path, subtitle: str, sd_col: str | None = None
) -> None:
    """Heatmap of `summary.chemistry_share` per (problem | target) row and model column.

    With `sd_col`, each cell reads "mean ± SD" instead of the bare value, to 3 decimals:
    seed SDs run down to ~0.003, which 2 decimals would print as 0.00.
    """
    index = ["problem", "target"]
    pivot = summary.pivot_table(index=index, columns="model", values="chemistry_share")
    annot: bool | np.ndarray = True
    if sd_col is not None:
        sd = summary.pivot_table(index=index, columns="model", values=sd_col)
        sd = sd.reindex(index=pivot.index, columns=pivot.columns)
        annot = np.array([
            [f"{m:.3f} ± {s:.3f}" for m, s in zip(mrow, srow, strict=True)]
            for mrow, srow in zip(pivot.to_numpy(), sd.to_numpy(), strict=True)
        ])
    # Second line per row: how many datapoints and how many distinct structures produced
    # that share. Repeated on every target of a multi-target problem so each row is
    # readable on its own rather than only in the context of its neighbours.
    scale = {p: _scale_label(p) for p in pivot.index.get_level_values("problem").unique()}
    row_labels = [f"{f'{p} | {t}'[:52]}\n{scale[p]}" for p, t in pivot.index]

    set_visualization_style()
    save_figure(
        plot_attribution_heatmap(
            pivot.to_numpy(),
            row_labels,
            list(pivot.columns),
            title=("Share of |attribution| carried by node features (vs tabular x_cols)"
                   f"\n{subtitle}"),
            diverging=False, cbar_label="Chemistry share of |IG|",
            annot=annot,
        ),
        out_path,
    )


def summarize_seeds(consistency_csv: Path, over: str = "seeds") -> None:
    """The chemistry-vs-process figure, averaged over the seeds in a consistency table.

    With `over="folds"`, averaged over its CV fold models instead, each the mean over the
    rows that fold model held out.
    """
    table = pd.read_csv(consistency_csv)
    if over == "folds":
        if "chem_share_fold_mean" not in table:
            raise ValueError(f"{consistency_csv} has no fold-level chemistry share; rerun "
                             "scripts/attribution_consistency.py with --folds.")
        n_folds = table.n_rankings.iloc[0]
        n_seeds = table.n_seeds.iloc[0]
        mean, sd, name = "chem_share_fold_mean", "chem_share_fold_sd", "chemistry_vs_process_folds"
        what = f"{n_folds} CV fold models ({n_folds // n_seeds} folds × {n_seeds} seeds)"
    else:
        lo, hi = table.n_seeds.min(), table.n_seeds.max()
        mean, sd, name = "chem_share_mean", "chem_share_sd", "chemistry_vs_process"
        what = f"{lo if lo == hi else f'{lo}-{hi}'} training seeds"
    plot_chemistry_share(
        table.rename(columns={mean: "chemistry_share"}),
        consistency_csv.parent / f"{name}.png",
        subtitle=f"mean ± SD over {what}; constant node-feature columns excluded",
        sd_col=sd,
    )
    logger.info(f"Seed-averaged summary written to {consistency_csv.parent}")


def summarize_all(config: AttributionConfig) -> None:
    """Cross-dataset figure: how much of each prediction is chemistry vs process."""
    if config.consistency_csv:
        summarize_seeds(Path(config.consistency_csv), config.consistency_over)
        return
    rows = []
    for problem in PROBLEMS:
        for model_type in SUMMARY_MODELS:
            path = config.results_dir / f"{problem}_{model_type}" / "completeness.csv"
            if not path.exists():
                logger.warning(f"Missing {path}; skipping.")
                continue
            df = pd.read_csv(path)
            for target, grp in df.groupby("target"):
                rows.append({
                    "problem": problem, "model": model_type, "target": target,
                    "chemistry_share": grp.chemistry_share.mean(),
                    "median_rel_error": grp.rel_error.median(),
                })
    if not rows:
        raise FileNotFoundError(f"No completeness.csv found under {config.results_dir}.")

    summary = pd.DataFrame(rows)
    summary.to_csv(config.results_dir / "chemistry_vs_process.csv", index=False)
    plot_chemistry_share(
        summary,
        config.results_dir / "chemistry_vs_process.png",
        subtitle="constant node-feature columns excluded",
    )
    logger.info(f"Summary written to {config.results_dir}")


@draccus.wrap()
def main(config: AttributionConfig) -> None:
    set_seed(config.seed)
    if config.summarize:
        summarize_all(config)
        return
    if config.model_type not in GRAPH_MODELS:
        raise ValueError(f"Need a graph model, one of {GRAPH_MODELS}.")
    attribute_one(config)


if __name__ == "__main__":
    main()
