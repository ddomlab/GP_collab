"""Dataset configuration and loading for the polymer/membrane regression problems."""

from dataclasses import dataclass
import math
from pathlib import Path

import numpy as np
import pandas as pd

from polykernels.paths_and_constants import Constants, Paths

_P = Paths()
_C = Constants()


@dataclass(frozen=True)
class DatasetSpec:
    """Configuration for one regression problem: SMILES columns, features, and targets.

    All column collections are tuples so a spec is genuinely immutable and hashable --
    several specs share the same underlying `Constants.x_cols` entry, and a mutable list
    would let one problem's edits leak into another's.
    """

    csv_path: Path
    smiles_cols: tuple[str, ...]
    x_cols: tuple[str, ...]
    y_cols: tuple[str, ...]
    impute_mean_cols: tuple[str, ...] = ()
    cluster_cols: tuple[str, ...] = ()
    # Smallest held-out cluster, as a fraction of the dataset's rows (the floor is always at
    # least 3 rows). dataset5 keeps 0.0, i.e. every cluster with 3+ members, as it was run;
    # the structural groupings added for the other datasets have long tails of tiny groups,
    # so they hold out only groups with at least 5% of rows.
    cluster_min_fraction: float = 0.0
    # Whether the GNN paths should embed real 3D conformers (default) or fall back to
    # fast, always-succeeding 2D coordinates. False for donor/acceptor-style problems
    # whose combined SMILES is two physically unrelated fragments -- 3D embedding has
    # no real joint geometry to solve for there and can be extremely slow or hang.
    use_3d: bool = True


# Per-spec details, kept here rather than in Constants: they describe one dataset's
# quirks (which feature columns need imputing, which columns define OOD clusters)
# rather than anything global.
_DATASET3_IMPUTE_MEAN_COLS = (
    "P_MW",
    "surface tension (mN/m)",
    "pore maker molecular weight (Da)",
    "solubility parameter (MPa1/2)",
    "organic compound size (Da)",
)
# Primary column contributes one LOCO split per eligible cluster; the finer Side Chain
# column contributes one split for its rarest eligible cluster.
_DATASET5_CLUSTER_COLS = ("substructure cluster", "Side Chain Cluster")
# The remaining datasets have no precomputed clusters, so each is grouped by the identity of
# its molecule: donor polymer (dataset1 acceptors are all non-fullerene, leaving no coarser
# acceptor split), membrane polymer (datasets 2-4), or substrate (dataset6).
_STRUCTURAL_MIN_FRACTION = 0.05


REGRESSION_PROBLEMS: dict[str, DatasetSpec] = {
    "dataset1_pce": DatasetSpec(
        csv_path=_P.dataset1_csv,
        smiles_cols=(_C.dataset1_donor_smiles, _C.dataset1_acceptor_smiles),
        x_cols=tuple(_C.x_cols[1]),
        y_cols=tuple(_C.y_cols[1]),
        use_3d=False,
        cluster_cols=("Donor",),
        cluster_min_fraction=_STRUCTURAL_MIN_FRACTION,
    ),
    "dataset2_flux": DatasetSpec(
        csv_path=_P.dataset2_1_csv,
        smiles_cols=(_C.dataset2_smiles,),
        x_cols=tuple(_C.x_cols[2]),
        y_cols=("log (Total flux)",),
        cluster_cols=("Abbrevation",),
        cluster_min_fraction=_STRUCTURAL_MIN_FRACTION,
    ),
    "dataset2_separation": DatasetSpec(
        csv_path=_P.dataset2_2_csv,
        smiles_cols=(_C.dataset2_smiles,),
        x_cols=tuple(_C.x_cols[2]),
        y_cols=("log (Separation factor)",),
        cluster_cols=("Abbrevation",),
        cluster_min_fraction=_STRUCTURAL_MIN_FRACTION,
    ),
    "dataset3_membrane": DatasetSpec(
        csv_path=_P.dataset3_csv,
        smiles_cols=(_C.dataset3_smiles,),
        x_cols=tuple(_C.x_cols[3]),
        y_cols=tuple(_C.y_cols[3]),
        impute_mean_cols=_DATASET3_IMPUTE_MEAN_COLS,
        cluster_cols=("Polymer",),
        cluster_min_fraction=_STRUCTURAL_MIN_FRACTION,
    ),
    "dataset4_pervaporation": DatasetSpec(
        csv_path=_P.dataset4_csv,
        smiles_cols=(_C.dataset4_smiles,),
        x_cols=tuple(_C.x_cols[4]),
        y_cols=tuple(_C.y_cols[4]),
        cluster_cols=("Polymer",),
        cluster_min_fraction=_STRUCTURAL_MIN_FRACTION,
    ),
    "dataset5_rg": DatasetSpec(
        csv_path=_P.dataset5_csv,
        smiles_cols=(_C.dataset5_smiles,),
        x_cols=tuple(_C.x_cols[5]),
        y_cols=tuple(_C.y_cols[5]),
        cluster_cols=_DATASET5_CLUSTER_COLS,
    ),
    "dataset6_reaction": DatasetSpec(
        csv_path=_P.dataset6_csv,
        smiles_cols=(_C.dataset6_smiles,),
        x_cols=tuple(_C.x_cols[6]),
        y_cols=tuple(_C.y_cols[6]),
        # Substrates are single small molecules with real 3D geometry, unlike dataset1's
        # two-fragment donor.acceptor strings.
        use_3d=True,
        cluster_cols=("Building Block",),
        cluster_min_fraction=_STRUCTURAL_MIN_FRACTION,
    ),
}


# The collaborator's two OOD groupings (github.com/ddomlab/GP_collab), imported by
# scripts/import_collab_clusters.py: `structure_cluster` clusters each unique molecule's
# ECFP6 counts, `continuous_cluster` clusters rows by their standardized numeric features.
COLLAB_CLUSTER_COLUMNS = ("structure_cluster", "continuous_cluster")
COLLAB_CLUSTERS_DIR = _P.data_dir / "processed" / "collab_clusters"


def load_collab_clusters(problem: str, column: str, y: np.ndarray) -> np.ndarray:
    """One of the collaborator's cluster columns, one label per loaded row.

    The labels are positional -- row i belongs to row i of `load_regression_data` -- so
    the file also carries the targets they were imported against, and `y` (as loaded
    now) must still match them. A dataset CSV that gains, loses or reorders rows fails
    here instead of silently pairing rows with another row's cluster.
    """
    if column not in COLLAB_CLUSTER_COLUMNS:
        raise ValueError(f"Unknown cluster column {column!r}; expected one of {COLLAB_CLUSTER_COLUMNS}.")
    path = COLLAB_CLUSTERS_DIR / f"{problem}.csv"
    if not path.exists():
        raise FileNotFoundError(f"No {path}. Run scripts/import_collab_clusters.py first.")
    table = pd.read_csv(path)
    imported_y = table[list(REGRESSION_PROBLEMS[problem].y_cols)].to_numpy(dtype=np.float32)
    if imported_y.shape != y.shape or not np.allclose(imported_y, y, equal_nan=True):
        raise ValueError(
            f"{path.name} was imported against different rows than {problem} loads now. "
            f"Re-run scripts/import_collab_clusters.py."
        )
    # Native label type (int, or str for dataset5's named structure clusters), so sorting
    # the labels gives the collaborator's group order.
    return table[column].to_numpy()


def min_cluster_size(spec: DatasetSpec, n_rows: int) -> int:
    """Smallest cluster eligible to be held out: 3 rows, or `cluster_min_fraction` of the data."""
    return max(3, math.ceil(spec.cluster_min_fraction * n_rows))


def load_regression_data(
    spec: DatasetSpec,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str], dict[str, np.ndarray] | None]:
    """Load one regression problem's SMILES, tabular features, targets, and cluster labels.

    Rows missing any SMILES or target value are dropped (a row must have every target
    present, since multi-target problems train one shared model across their columns).
    `impute_mean_cols` are mean-imputed from the (already NaN-dropped) column itself.

    Returns:
        smiles: [n] array of SMILES strings (donor and acceptor joined by "." for dataset1).
        x_features: [n, len(x_cols)] tabular feature matrix.
        y: [n, len(y_cols)] target matrix.
        target_names: `spec.y_cols`, aligned with y's columns.
        cluster_labels: {cluster_col: [n] array of labels}, or None if spec has no cluster_cols.
    """
    df = pd.read_csv(spec.csv_path)
    df = df.dropna(
        subset=list(spec.smiles_cols) + list(spec.y_cols)
    ).reset_index(drop=True)

    if len(spec.smiles_cols) == 2:
        donor, acceptor = spec.smiles_cols
        smiles = (df[donor].astype(str) + "." + df[acceptor].astype(str)).to_numpy()
    else:
        smiles = df[spec.smiles_cols[0]].to_numpy(dtype=str)

    x_df = df[list(spec.x_cols)].copy()
    for col in spec.impute_mean_cols:
        x_df[col] = x_df[col].fillna(x_df[col].mean())
    x_features = x_df.to_numpy(dtype=np.float32)

    # Fail loudly rather than let NaNs reach the model as a silent NaN loss: any feature
    # column that gains missing values must be added to `impute_mean_cols` deliberately.
    if np.isnan(x_features).any():
        bad = [c for c in spec.x_cols if x_df[c].isna().any()]
        raise ValueError(
            f"{spec.csv_path.name}: feature columns {bad} still contain missing values. "
            f"Add them to the spec's impute_mean_cols."
        )

    y = df[list(spec.y_cols)].to_numpy(dtype=np.float32)

    cluster_labels = (
        {col: df[col].astype(str).to_numpy() for col in spec.cluster_cols}
        if spec.cluster_cols
        else None
    )

    return smiles, x_features, y, list(spec.y_cols), cluster_labels
