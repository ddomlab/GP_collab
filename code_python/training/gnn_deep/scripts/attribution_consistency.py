"""Seed consistency of the Integrated Gradients attributions.

IG is deterministic given a model's weights, and the step count is already converged (see
completeness.csv), so the variation left worth measuring is between *models*: a conclusion
that holds for one network but not for the same architecture retrained under another seed
describes that network, not the data. This compares attribution sweeps fitted under
different training seeds -- each a directory laid out like results/attribution_globals --
and reports, per (problem, model, target), how far the conclusions move:

  dataset level   chemistry share (mean and SD across seeds); agreement of the node-feature
                  and x_col importance rankings (Kendall's W with a permutation p-value, and
                  mean pairwise Spearman of mean |IG|) and overlap of the top node features
  molecule level  per-molecule Kendall's W and Spearman of the signed per-atom
                  attribution, overlap of the top atoms by magnitude, and magnitude-weighted
                  sign agreement -- each summarized by the median over molecules

Kendall's W is the concordance of all seeds' rankings at once, 0 (none) to 1 (identical).
Without ties it equals ((m-1) * mean pairwise Spearman + 1) / m for m seeds, so it carries
no information the Spearman columns lack; its use is a single standard statistic with a
significance test. Its chance level is 1/m, not 0 -- about 0.33 for three seeds -- so read
a W next to its p-value, never alone.

Seeds isolate training randomness only if everything else matches, so sweeps whose saved
hyperparameters differ are refused rather than silently compared.

For CV attributions, `--folds` takes the per-fold sweeps (`iid_seed<s>_fold<f>`), and the
dataset-level rankings are then compared across all those fold models -- 15 for 5 folds x
3 seeds, so chance is 1/15 -- each ranking by the rows it held out. The molecule-level
metrics need sweeps that explain the same molecules, so they stay across `--runs`, the
seeds' merged sweeps (`atom_w_chance`). The chemistry share is then also summarized over
the fold models (`chem_share_fold_mean`, `chem_share_fold_sd`), each averaging its own rows,
and completeness by the worst fold model (`worst_completeness_fold_median`).

Examples:
    python scripts/attribution_consistency.py
    python scripts/attribution_consistency.py --runs results/attribution_globals \
        results/attribution_seeds/seed13 results/attribution_seeds/seed21
    python scripts/attribution_consistency.py --runs results/ig_cv/merged/seed{6,13,42} \
        --folds results/ig_cv/iid_seed{6,13,42}_fold{0,1,2,3,4} \
        --out results/attribution_cv_summary/ig_consistency.csv
"""

from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

TOP_K_FEATURES = 10
TOP_K_ATOMS = 5
# Null draws for Kendall's W. The chi-square approximation is poor for the 4-13 x_cols
# these datasets have, so the p-value comes from permutation; 10k puts its floor at 1e-4.
N_PERMUTATIONS = 10_000


def _tie_term(ranks: np.ndarray) -> float:
    """sum over raters of sum(t^3 - t) over each rater's groups of t tied ranks."""
    total = 0.0
    for r in ranks:
        _, t = np.unique(r, return_counts=True)
        total += float((t**3 - t).sum())
    return total


def kendall_w(scores: pd.DataFrame, rng: np.random.Generator) -> tuple[float, float]:
    """Tie-corrected Kendall's W across seeds (columns) of item scores (rows).

    The p-value is a permutation test of "the seeds rank independently": every seed but
    the first is shuffled, which keeps each seed's tie structure and so the correction.
    """
    scores = scores.dropna()
    n, m = scores.shape
    if m < 2 or n < 3:
        return float("nan"), float("nan")
    ranks = scores.rank().to_numpy().T                          # (m seeds, n items)
    denom = m**2 * (n**3 - n) - m * _tie_term(ranks)
    if denom <= 0:
        return float("nan"), float("nan")

    def w_of(rank_sums: np.ndarray) -> np.ndarray:
        return 12 * ((rank_sums - m * (n + 1) / 2) ** 2).sum(-1) / denom

    w = float(w_of(ranks.sum(0)))
    shuffled = rng.permuted(np.repeat(ranks[None, 1:], N_PERMUTATIONS, axis=0), axis=-1)
    null = w_of(ranks[0] + shuffled.sum(1))
    return w, float((1 + np.sum(null >= w - 1e-12)) / (1 + N_PERMUTATIONS))


def per_molecule_w(nodes: pd.DataFrame, labels: list[str]) -> pd.Series:
    """Kendall's W of the signed per-atom attribution across all seeds, per molecule."""
    m = len(labels)
    g = nodes.groupby("row", sort=False)
    rank_sum = sum(g[f"signed_{s}"].rank() for s in labels)
    n_atoms = g["row"].transform("size")
    s = pd.DataFrame({"row": nodes["row"], "n": n_atoms,
                      "dev2": (rank_sum - m * (n_atoms + 1) / 2) ** 2}
                     ).groupby("row", sort=False).agg(n=("n", "first"), S=("dev2", "sum"))
    ties = sum(
        nodes.groupby(["row", f"signed_{lab}"]).size().pipe(lambda t: t**3 - t)
        .groupby(level=0).sum().reindex(s.index, fill_value=0)
        for lab in labels
    )
    denom = m**2 * (s.n**3 - s.n) - m * ties
    return (12 * s.S / denom).where(denom > 0)


def spearman(a: pd.Series, b: pd.Series) -> float:
    """Rank correlation of two vectors, aligned by index; NaN when either is constant."""
    ra, rb = a.rank(), b.reindex(a.index).rank()
    if ra.std() == 0 or rb.std() == 0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def top_k_overlap(a: pd.Series, b: pd.Series, k: int) -> float:
    """Fraction of the k largest entries of `a` that are also among the k largest of `b`."""
    k = min(k, len(a))
    return len(set(a.nlargest(k).index) & set(b.nlargest(k).index)) / k


def per_molecule(nodes: pd.DataFrame, a: str, b: str) -> pd.DataFrame:
    """Per-molecule agreement between two seeds' per-atom attributions, vectorized.

    `nodes` holds one row per (molecule, atom) of one target, with `signed_<seed>` and
    `magnitude_<seed>` columns. Spearman is Pearson on within-molecule ranks, computed from
    grouped sums so a sweep of thousands of molecules takes seconds, not minutes.
    """
    g = nodes.groupby("row", sort=False)
    ra, rb = g[f"signed_{a}"].rank(), g[f"signed_{b}"].rank()
    # Top atoms by magnitude; `first` breaks ties so exactly k atoms are selected.
    top_a = g[f"magnitude_{a}"].rank(ascending=False, method="first") <= TOP_K_ATOMS
    top_b = g[f"magnitude_{b}"].rank(ascending=False, method="first") <= TOP_K_ATOMS
    sa, sb = nodes[f"signed_{a}"], nodes[f"signed_{b}"]
    # Near-zero atoms flip sign on noise alone; weighting by size counts the atoms a
    # reader would actually look at.
    w = (sa.abs() + sb.abs()) / 2
    t = pd.DataFrame({
        "row": nodes["row"], "ra": ra, "rb": rb, "ra2": ra * ra, "rb2": rb * rb,
        "rab": ra * rb, "top_a": top_a, "both": top_a & top_b,
        "w": w, "w_agree": w * (np.sign(sa) == np.sign(sb)),
    })
    s = t.groupby("row", sort=False).agg(
        n=("ra", "size"), sa=("ra", "sum"), sb=("rb", "sum"), saa=("ra2", "sum"),
        sbb=("rb2", "sum"), sab=("rab", "sum"), top=("top_a", "sum"),
        both=("both", "sum"), w=("w", "sum"), w_agree=("w_agree", "sum"),
    )
    var_a, var_b = s.saa - s.sa**2 / s.n, s.sbb - s.sb**2 / s.n
    cov = s.sab - s.sa * s.sb / s.n
    ok = (var_a > 0) & (var_b > 0)
    return pd.DataFrame({
        "atom_spearman": (cov / np.sqrt(var_a * var_b)).where(ok),
        "atom_top_overlap": s.both / s.top,
        "sign_agreement": (s.w_agree / s.w).where(s.w > 0),
    })


def training_config(run_dir: Path) -> dict | None:
    """Everything that defines the fitted model except the seed; None if unrecorded."""
    path = run_dir / "hparams.yaml"
    if not path.exists():
        return None
    rec = yaml.safe_load(path.read_text())
    fit = {k: v for k, v in rec.get("fit", {}).items()
           if k not in ("seed", "results_dir", "seed_globals_override")}
    return {"fit": fit, "hparams": rec["hparams"], "n_rows": rec.get("n_rows")}


def compare_run(name: str, dirs: dict[str, Path],
                rank_dirs: dict[str, Path] | None = None) -> list[dict]:
    """All consistency metrics for one (problem, model) across the sweeps that have it.

    `rank_dirs`, when given, replaces `dirs` for the dataset-level rankings (the CV fold
    models); the molecule-level metrics always compare `dirs`.
    """
    rank_dirs = rank_dirs or dirs
    configs = {label: training_config(d) for label, d in {**dirs, **rank_dirs}.items()}
    known = {label: c for label, c in configs.items() if c is not None}
    if len({yaml.safe_dump(c, sort_keys=True) for c in known.values()}) > 1:
        raise SystemExit(
            f"{name}: hyperparameters differ between {sorted(known)}, so seed-to-seed "
            f"differences would be confounded with the config change. Refit them alike."
        )
    missing = sorted(set(configs) - set(known))
    if missing:
        print(f"  warning: {name}: no hparams.yaml in {missing}; training config unverified")

    labels, rank_labels = list(dirs), list(rank_dirs)
    comp = {s: pd.read_csv(d / "completeness.csv") for s, d in dirs.items()}
    feats = {s: pd.read_csv(d / "feature_attribution.csv") for s, d in rank_dirs.items()}
    fold_comp = None if rank_dirs is dirs else {
        s: pd.read_csv(d / "completeness.csv", usecols=["target", "chemistry_share", "rel_error"])
        for s, d in rank_dirs.items()}
    globs = {s: pd.read_csv(d / "global_attribution.csv") for s, d in rank_dirs.items()}

    cols = ["row", "target", "atom_idx", "signed", "magnitude"]
    nodes = None
    for s, d in dirs.items():
        df = pd.read_csv(d / "node_attribution.csv", usecols=cols).rename(
            columns={"signed": f"signed_{s}", "magnitude": f"magnitude_{s}"})
        nodes = df if nodes is None else nodes.merge(df, on=cols[:3], how="inner")
    n_ref = len(pd.read_csv(dirs[labels[0]] / "node_attribution.csv", usecols=["row"]))
    if len(nodes) != n_ref:
        raise SystemExit(f"{name}: sweeps attribute different atoms; were the data changed?")

    problem, model = name.rsplit("_", 2)[0], "_".join(name.rsplit("_", 2)[1:])
    out = []
    for target in comp[labels[0]]["target"].unique():
        shares = [comp[s].loc[comp[s].target == target, "chemistry_share"].mean()
                  for s in labels]
        fold_stats = {}
        if fold_comp is not None:
            per_fold = [fold_comp[s][fold_comp[s].target == target] for s in rank_labels]
            shares_f = [f.chemistry_share.mean() for f in per_fold]
            fold_stats = {"chem_share_fold_mean": float(np.mean(shares_f)),
                          "chem_share_fold_sd": float(np.std(shares_f, ddof=1)),
                          "worst_completeness_fold_median": float(
                              max(f.rel_error.median() for f in per_fold))}
        worst_rel = max(comp[s].loc[comp[s].target == target, "rel_error"].median()
                        for s in labels)

        imp = {s: feats[s][(feats[s].target == target) & ~feats[s].is_constant]
               .set_index("feature")["mean_abs_ig"] for s in rank_labels}
        gimp = {s: globs[s][globs[s].target == target].set_index("x_col")["mean_abs_ig"]
                for s in rank_labels}

        node_t = nodes[nodes.target == target]
        pairs = list(combinations(labels, 2))
        rank_pairs = list(combinations(rank_labels, 2))
        mol = sum(per_molecule(node_t, a, b) for a, b in pairs) / len(pairs)
        # Seeded per target so a rerun reproduces the same p-values.
        rng = np.random.default_rng(0)
        feat_w, feat_p = kendall_w(pd.DataFrame(imp), rng)
        xcol_w, xcol_p = kendall_w(pd.DataFrame(gimp), rng)

        out.append({
            "problem": problem, "model": model, "target": target,
            "n_seeds": len(labels), "sweeps": " ".join(labels),
            "n_rankings": len(rank_labels),
            "chem_share_mean": float(np.mean(shares)),
            "chem_share_sd": float(np.std(shares, ddof=1)),
            **fold_stats,
            "w_chance": 1 / len(rank_labels),
            "feature_kendall_w": feat_w, "feature_w_p": feat_p,
            "xcol_kendall_w": xcol_w, "xcol_w_p": xcol_p,
            "atom_w_chance": 1 / len(labels),
            "atom_kendall_w_median": float(per_molecule_w(node_t, labels).median()),
            "feature_spearman":np.mean([spearman(imp[a], imp[b]) for a, b in rank_pairs]),
            f"feature_top{TOP_K_FEATURES}_overlap": np.mean(
                [top_k_overlap(imp[a], imp[b], TOP_K_FEATURES) for a, b in rank_pairs]),
            # Two or fewer x_cols cannot be meaningfully rank-correlated.
            "xcol_spearman": (np.mean([spearman(gimp[a], gimp[b]) for a, b in rank_pairs])
                              if len(gimp[rank_labels[0]]) > 2 else float("nan")),
            "atom_spearman_median": float(mol["atom_spearman"].median()),
            f"atom_top{TOP_K_ATOMS}_overlap_median": float(mol["atom_top_overlap"].median()),
            "sign_agreement_median": float(mol["sign_agreement"].median()),
            "n_molecules": len(mol),
            "worst_completeness_median": float(worst_rel),
        })
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", type=Path, nargs="+", default=None,
                    help="sweep directories, one per seed (default: attribution_globals "
                         "plus every results/attribution_seeds/seed*)")
    ap.add_argument("--folds", type=Path, nargs="+", default=None,
                    help="per-fold CV sweeps to compare the dataset-level rankings across "
                         "instead of --runs (molecule-level metrics stay on --runs)")
    ap.add_argument("--out", type=Path,
                    default=Path("results/attribution_seeds/consistency.csv"))
    args = ap.parse_args()

    runs = args.runs or [Path("results/attribution_globals"),
                         *sorted(Path("results/attribution_seeds").glob("seed*"))]
    labels = [r.name for r in runs]
    if len(set(labels)) != len(labels):
        raise SystemExit(f"Sweep directory names must be distinct labels; got {labels}.")
    folds = args.folds or []
    if len({f.name for f in folds}) != len(folds):
        raise SystemExit(f"Fold directory names must be distinct labels; got {folds}.")

    names = sorted({d.name for r in runs for d in r.iterdir()
                    if (d / "completeness.csv").exists()})
    rows = []
    for name in names:
        dirs = {r.name: r / name for r in runs if (r / name / "completeness.csv").exists()}
        if len(dirs) < 2:
            print(f"  skipping {name}: present in only {list(dirs)}")
            continue
        rank_dirs = None
        if folds:
            rank_dirs = {f.name: f / name for f in folds
                         if (f / name / "completeness.csv").exists()}
            if len(rank_dirs) < len(folds):
                print(f"  skipping {name}: {len(rank_dirs)} of {len(folds)} folds finished")
                continue
        rows.extend(compare_run(name, dirs, rank_dirs))
    if not rows:
        raise SystemExit(f"No run is present in two or more of {runs}.")

    table = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, index=False)
    show = ["problem", "model", "target", "n_seeds", "n_rankings", "chem_share_mean",
            "chem_share_sd", "feature_kendall_w", "feature_w_p", "xcol_kendall_w", "xcol_w_p",
            "atom_kendall_w_median", "sign_agreement_median"]
    with pd.option_context("display.width", 200, "display.max_columns", None,
                           "display.float_format", "{:.3f}".format):
        print(table[show].to_string(index=False))
    print(f"\nKendall's W at chance is 1/n_rankings for the feature and x_col rankings "
          f"(w_chance) and 1/n_seeds for atoms (atom_w_chance); p-values floor at "
          f"{1 / (1 + N_PERMUTATIONS):.0e}.")
    print(f"wrote {args.out}  ({len(table)} rows)")


if __name__ == "__main__":
    main()
