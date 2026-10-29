"""How far SHAP, GP-SHAP and IG agree with each other, and each with itself across models.

Reads the per-fold CV attribution runs of `hpc/killarney/attribution_cv.sh`,
`<method>_cv/iid_seed<s>_fold<f>/<problem>_<model>/`, where each fold model explains its own
held-out rows, and compares them on the tabular x_cols -- the inputs all three methods
attribute (IG also attributes node features; `attribution_consistency.py` covers those).
Two tables, one row per (problem, model, target):

  fold_consistency.csv   per method, across all the fold models of every training seed
                         (15 for 5 folds x 3 seeds)
  method_agreement.csv   per fold model, across the methods that cover the model
                         (+ method_agreement_summary.csv, averaged over the fold models)

Each reports Kendall's W of the global importance ranking (mean |attribution| per x_col over
the rows the fold model explains) with its permutation p-value, mean pairwise Spearman,
top-3 overlap, whether the top x_col is the same everywhere, and the median over rows of
the per-row W of the signed local values. W's chance level is 1/m (m = fold models or
methods): read it next to `w_chance` and `p`. A row is held out by one fold model per seed,
so across fold models the per-row W compares the seeds, at its own chance `row_w_chance`.
Rankings do not depend on units, so IG's standardized values compare directly.

    python scripts/attribution_agreement.py
"""

import argparse
from itertools import combinations
from pathlib import Path
import re

import numpy as np
import pandas as pd
from scripts.attribution_consistency import (
    kendall_w,
    per_molecule_w,
    spearman,
    top_k_overlap,
)

TOP_K = 3
METHODS = ("shap", "gpshap", "ig")
LABELS = {"shap": "SHAP", "gpshap": "GP-SHAP", "ig": "IG"}
NPZ = {"shap": "shap_values.npz", "gpshap": "gpshap_values.npz", "ig": "ig_raw.npz"}
FOLD_DIR = re.compile(r"iid_seed(\d+)_fold(\d+)")


def load_local(method: str, run_dir: Path) -> tuple[np.ndarray, list[str], list[str], np.ndarray]:
    """Rows, x_cols, targets and signed local values `[rows, x_cols, targets]` of one run."""
    z = np.load(run_dir / NPZ[method], allow_pickle=True)
    if method == "shap":
        rows, cols, values = z["explain_idx"], z["feature_names"], z["values"]
    elif method == "gpshap":
        rows, cols, values = z["explain_idx"], z["feature_names"], z["mean"]
    else:
        rows, cols, values = z["row_idx"], z["x_cols"], z["global_ig"].transpose(0, 2, 1)
    return rows, [str(c) for c in cols], [str(t) for t in z["target_names"]], values


def rank_agreement(imp: pd.DataFrame, *, pair_columns: bool, rng_seed: int = 0) -> dict:
    """Agreement of the raters' (columns) global importance of the x_cols (rows)."""
    labels = list(imp)
    w, p = kendall_w(imp, np.random.default_rng(rng_seed))
    pairs = list(combinations(labels, 2))
    out = {
        "n_raters": len(labels), "w_chance": 1 / len(labels), "kendall_w": w, "w_p": p,
        "spearman_mean": np.mean([spearman(imp[a], imp[b]) for a, b in pairs]),
        f"top{TOP_K}_overlap": np.mean([top_k_overlap(imp[a], imp[b], TOP_K) for a, b in pairs]),
        "same_top1": imp.idxmax().nunique() == 1,
        "n_x_cols": imp.shape[0],
    }
    if pair_columns:
        out.update({f"spearman_{a}_vs_{b}": spearman(imp[a], imp[b]) for a, b in pairs})
    return out


def aligned(local: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    """The raters' `[rows, x_cols]` frames in one order; they must cover the same of each."""
    first = next(iter(local.values()))
    if any(set(v.columns) != set(first.columns) or set(v.index) != set(first.index)
           for v in local.values()):
        raise ValueError("raters attribute different x_cols or rows")
    return {k: v.loc[first.index, first.columns] for k, v in local.items()}


def row_w_median(local: dict[str, pd.DataFrame]) -> float:
    """Median over rows of the W of each row's signed values, across raters of the same rows."""
    local = aligned(local)
    first = next(iter(local.values()))
    nodes = pd.DataFrame({
        "row": np.repeat(first.index.to_numpy(), first.shape[1]),
        **{f"signed_{k}": v.to_numpy().ravel() for k, v in local.items()},
    })
    return float(per_molecule_w(nodes, list(local)).median())


def agreement(local: dict[str, pd.DataFrame], rng_seed: int = 0) -> dict:
    """Agreement of several raters of the same rows (methods), each a `[rows, x_cols]` frame."""
    local = aligned(local)
    imp = pd.DataFrame({k: v.abs().mean(axis=0) for k, v in local.items()})
    return {**rank_agreement(imp, pair_columns=len(local) > 2, rng_seed=rng_seed),
            "row_kendall_w_median": row_w_median(local),
            "n_rows": len(next(iter(local.values())))}


def fold_agreement(folds: dict[tuple[int, int], pd.DataFrame], rng_seed: int = 0) -> dict:
    """Agreement of one method's fold models `{(seed, fold): [rows, x_cols]}`.

    Every fold model ranks the x_cols by its own held-out rows, so all of them are rated
    together. The folds of a seed partition the rows, so stacking them gives each seed one
    value per row, from the fold model that held it out, and the per-row W compares seeds.
    """
    imp = pd.DataFrame({f"seed{s}_fold{f}": v.abs().mean(axis=0) for (s, f), v in folds.items()})
    by_seed: dict[str, list[pd.DataFrame]] = {}
    for (s, _), v in folds.items():
        by_seed.setdefault(f"seed{s}", []).append(v)
    stacked = {k: pd.concat(v).sort_index() for k, v in by_seed.items()}
    return {**rank_agreement(imp, pair_columns=False, rng_seed=rng_seed),
            "row_raters": len(stacked), "row_w_chance": 1 / len(stacked),
            "row_kendall_w_median": row_w_median(stacked),
            "n_rows": len(next(iter(stacked.values())))}


def collect(roots: dict[str, Path]) -> dict:
    """{(method, seed, fold, run): (rows, x_cols, targets, values)} for every fold run."""
    runs = {}
    for method, root in roots.items():
        for fold_dir in sorted(root.glob("iid_seed*_fold*")):
            seed, fold = map(int, FOLD_DIR.fullmatch(fold_dir.name).groups())
            for run_dir in sorted(fold_dir.iterdir()):
                if (run_dir / NPZ[method]).exists():
                    runs[(method, seed, fold, run_dir.name)] = load_local(method, run_dir)
    return runs


def frame(entry, t_idx: int) -> pd.DataFrame:
    rows, cols, _, values = entry
    return pd.DataFrame(values[:, :, t_idx], index=rows, columns=cols).sort_index()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out_dir", type=Path, default=Path("results/attribution_cv_summary"))
    args = ap.parse_args()
    roots = {m: Path(f"results/{m}_cv") for m in METHODS}
    runs = collect(roots)
    seeds = sorted({s for _, s, _, _ in runs})
    folds = sorted({f for _, _, f, _ in runs})
    names = sorted({r for *_, r in runs})

    fold_rows, method_rows = [], []
    for run in names:
        problem, model = run.rsplit("_", 2)[0], "_".join(run.rsplit("_", 2)[1:])
        base = runs[next(k for k in runs if k[3] == run)]
        for t_idx, target in enumerate(base[2]):
            key = {"problem": problem, "model": model, "target": target}
            for method in METHODS:
                have = [(s, f) for s in seeds for f in folds if (method, s, f, run) in runs]
                if not have:
                    continue
                # The per-row W needs every seed to explain every row once.
                if len(have) < len(seeds) * len(folds):
                    print(f"  skip {LABELS[method]} {run}: has {len(have)} of "
                          f"{len(seeds) * len(folds)} fold models")
                    continue
                stats = fold_agreement({(s, f): frame(runs[(method, s, f, run)], t_idx)
                                        for s, f in have})
                fold_rows.append({"method": LABELS[method], **key,
                                  "seeds": " ".join(map(str, seeds)), **stats})
            for s in seeds:
                for f in folds:
                    have = [m for m in METHODS if (m, s, f, run) in runs]
                    if len(have) >= 2:
                        local = {LABELS[m]: frame(runs[(m, s, f, run)], t_idx) for m in have}
                        method_rows.append({"seed": s, "fold": f, **key,
                                            "methods": " ".join(local), **agreement(local)})

    args.out_dir.mkdir(parents=True, exist_ok=True)
    if fold_rows:
        pd.DataFrame(fold_rows).to_csv(args.out_dir / "fold_consistency.csv", index=False)
    by_fold = pd.DataFrame(method_rows)
    by_fold.to_csv(args.out_dir / "method_agreement.csv", index=False)
    stat_cols = [c for c in by_fold.columns if c not in
                 ("seed", "fold", "problem", "model", "target", "methods")]
    summary = (by_fold.groupby(["problem", "model", "target", "methods"], sort=False)
               .agg(n_folds=("fold", "size"),
                    **{f"{c}_mean": (c, "mean") for c in stat_cols
                       if c not in ("n_raters", "w_chance", "n_x_cols", "n_rows")},
                    kendall_w_sd=("kendall_w", "std"), w_chance=("w_chance", "first"))
               .reset_index())
    summary.to_csv(args.out_dir / "method_agreement_summary.csv", index=False)
    print(f"seeds {seeds}, folds {folds}; {len(fold_rows)} fold-consistency rows, "
          f"{len(by_fold)} method-agreement rows -> {args.out_dir}")


if __name__ == "__main__":
    main()
