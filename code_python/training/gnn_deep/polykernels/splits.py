"""Train, validation, and test split management and serialization."""

from collections.abc import Sequence
import dataclasses
from pathlib import Path

import jaxtyping as jt
import numpy as np
from sklearn.model_selection import KFold

from polykernels.io import load_npz, save_npz
from polykernels.stratification import (
    get_multilabel_train_test_split,
    iterative_stratification,
    map_fold_to_indices,
)
from polykernels.types import IntArray
from polykernels.utils import cast_random_state


@dataclasses.dataclass(frozen=True)
class SplitTuple:
    """Basic train, validation, test split, assumes non-overlapping indices."""

    train: IntArray[" n"]
    val: IntArray[" j"]
    test: IntArray[" k"]

    def __post_init__(self) -> None:
        """Raises error if split parts overlap."""
        for a, b in [("train", "val"), ("val", "test"), ("test", "train")]:
            overlap = np.intersect1d(getattr(self, a), getattr(self, b))
            if overlap.shape[0]:
                raise ValueError(
                    f"Found {len(overlap)} overlapping items between {a} and {b}."
                )

    def __len__(self) -> int:
        return len(self.train) + len(self.val) + len(self.test)

    def trainval(self) -> IntArray[" trainval"]:
        return np.concatenate([self.train, self.val], axis=0)

    def indices(self) -> IntArray[" indices"]:
        return np.concatenate([self.train, self.val, self.test], axis=0)

    def as_dict(self) -> dict[str, IntArray[" *"]]:
        return dataclasses.asdict(self)

    @classmethod
    def from_lists(
        cls,
        train: Sequence[int],
        val: Sequence[int],
        test: Sequence[int],
    ) -> "SplitTuple":
        return cls(
            train=np.asarray(train),
            val=np.asarray(val),
            test=np.asarray(test),
        )


def serialize_splits(
    path: str | Path,
    splits: dict[str, SplitTuple],
    index: jt.Shaped[np.ndarray, " n"],
) -> None:
    """Saves named splits and a master index as a compressed npz file."""
    data = {"index": index}
    for name, split in splits.items():
        data[f"{name}_train"] = split.train
        data[f"{name}_val"] = split.val
        data[f"{name}_test"] = split.test
    save_npz(Path(path), **data)


def deserialize_splits(
    path: str | Path,
) -> tuple[dict[str, SplitTuple], jt.Shaped[np.ndarray, " n"]]:
    """Loads named splits and a master index from a compressed npz file."""
    data = load_npz(Path(path))
    index = data["index"]

    split_names = set()
    for key in data.keys():
        if key == "index":
            continue
        if key.endswith(("_train", "_val", "_test")):
            name = key.rsplit("_", 1)[0]
            split_names.add(name)

    splits = {}
    for name in sorted(split_names):
        splits[name] = SplitTuple(
            train=data[f"{name}_train"],
            val=data[f"{name}_val"],
            test=data[f"{name}_test"],
        )
    return splits, index


def deserialize_split(
    path: str | Path,
    name: str = "tvt_split",
) -> tuple[SplitTuple, jt.Shaped[np.ndarray, " n"]]:
    """Loads a single named SplitTuple and the master index from an npz file."""
    splits, index = deserialize_splits(path)
    if name not in splits:
        raise ValueError(
            f"Split '{name}' not found. Available splits: {list(splits.keys())}"
        )
    return splits[name], index


def get_multilabel_train_val_test_split(
    labels: jt.Float[np.ndarray, " b l"],
    val_size: float,
    test_size: float,
    order: int = 1,
    random_state: int | np.random.RandomState | None = None,
) -> SplitTuple:
    """Gets indices for a multilabel train/val/test split.

    Args:
        labels: 2D target matrix of shape (n_samples, n_labels).
        val_size: Proportion of the dataset to include in the validation set.
        test_size: Proportion of the dataset to include in the test set.
        order: Stratification order (1 or 2).
        random_state: RandomState instance or integer seed.

    Returns:
        SplitTuple containing train, val, and test split indices.
    """
    rng = cast_random_state(random_state)
    trainval_index, test_index = get_multilabel_train_test_split(
        labels, test_size, order=order, random_state=rng
    )
    train_size = 1.0 - (val_size + test_size)
    train_val_ratio = val_size / (train_size + val_size)
    tmp_train_index, tmp_val_index = get_multilabel_train_test_split(
        labels[trainval_index], train_val_ratio, order=order, random_state=rng
    )
    train_index = trainval_index[tmp_train_index]
    val_index = trainval_index[tmp_val_index]
    return SplitTuple(train=train_index, val=val_index, test=test_index)


def get_multilabel_crossvalidation_splits(
    labels: jt.Float[np.ndarray, " b l"],
    n_splits: int,
    order: int = 1,
    random_state: int | np.random.RandomState | None = None,
) -> list[tuple[IntArray[" n"], IntArray[" k"]]]:
    """Gets indices for stratified train/test cross-validation folds.

    Args:
        labels: 2D target matrix of shape (n_samples, n_labels).
        n_splits: Number of splits/folds.
        order: Stratification order (1 or 2).
        random_state: RandomState instance or integer seed.

    Returns:
        splits: A list of (train_indices, test_indices) tuples.
    """
    rng = cast_random_state(random_state)
    ratios = np.ones(n_splits) / n_splits
    folds = iterative_stratification(labels, ratios, order=order, random_state=rng)
    splits = [map_fold_to_indices(folds, i) for i in range(n_splits)]
    return splits


def get_regression_kfold_splits(
    n_samples: int,
    n_splits: int,
    random_state: int | np.random.RandomState | None = None,
) -> list[SplitTuple]:
    """Gets plain (non-stratified) shuffled K-fold splits for scalar/continuous regression targets.

    Args:
        n_samples: Number of samples to split.
        n_splits: Number of folds.
        random_state: RandomState instance or integer seed.

    Returns:
        One SplitTuple per fold (empty `val`, since these are train/test CV folds).
    """
    kfold = KFold(n_splits=n_splits, shuffle=True, random_state=cast_random_state(random_state))
    empty = np.empty(0, dtype=int)
    return [
        SplitTuple(train=train_idx, val=empty, test=test_idx)
        for train_idx, test_idx in kfold.split(np.arange(n_samples))
    ]


def _eligible_cluster_labels(
    labels: jt.Shaped[np.ndarray, " n"],
    min_size: int,
    *,
    rarest_only: bool = False,
) -> np.ndarray:
    """Cluster labels eligible for a leave-one-cluster-out test fold.

    Clusters with fewer than `min_size` members are skipped. If there are exactly 2
    clusters total (or `rarest_only` is set), only the rarest eligible one is held out,
    avoiding two near-duplicate complementary splits.
    """
    unique, counts = np.unique(labels, return_counts=True)
    eligible_mask = counts >= min_size
    eligible, eligible_counts = unique[eligible_mask], counts[eligible_mask]
    if eligible.size == 0:
        return eligible
    if rarest_only or len(unique) == 2:
        return eligible[[np.argmin(eligible_counts)]]
    return eligible


def _leave_one_out_split(
    labels: jt.Shaped[np.ndarray, " n"], label: object
) -> SplitTuple:
    """Split holding out every sample whose cluster equals `label`."""
    return SplitTuple(
        train=np.flatnonzero(labels != label),
        val=np.empty(0, dtype=int),
        test=np.flatnonzero(labels == label),
    )


def get_loco_splits(
    cluster_labels: dict[str, jt.Shaped[np.ndarray, " n"]],
    min_size: int = 3,
) -> dict[str, SplitTuple]:
    """Leave-one-cluster-out (LOCO) splits for out-of-distribution evaluation.

    For the first (primary) cluster column, builds one split per eligible label
    (test = that cluster's members, train = everyone else). Every further column is
    treated as a finer sub-clustering -- the "hierarchical" case, e.g. a coarse
    "substructure cluster" plus a finer "Side Chain Cluster" -- and contributes one
    extra split for its rarest eligible label.

    Args:
        cluster_labels: {column_name: [n] array of cluster labels}, ordered
            primary-first (as returned by `polykernels.regression_data.load_regression_data`).
        min_size: Minimum cluster size to be eligible as a held-out test fold.

    Returns:
        {"{column_name}={label}": SplitTuple} for every eligible (column, label) split.
    """
    splits: dict[str, SplitTuple] = {}
    for position, (column, raw_labels) in enumerate(cluster_labels.items()):
        labels = np.asarray(raw_labels)
        # Primary column contributes every eligible cluster; finer columns only their rarest.
        for label in _eligible_cluster_labels(
            labels, min_size, rarest_only=position > 0
        ):
            splits[f"{column}={label}"] = _leave_one_out_split(labels, label)
    return splits


def get_matched_iid_baseline_splits(
    n_total: int,
    n_test: int,
    random_state: int | np.random.RandomState | None = None,
    max_folds: int | None = 5,
) -> list[SplitTuple]:
    """Size-matched IID K-fold baseline for fair comparison against one LOCO test cluster.

    Uses `n_splits = round(n_total / n_test)` so each baseline fold's test size roughly
    matches the LOCO cluster being compared against.

    Args:
        n_total: Total number of samples in the full dataset.
        n_test: Size of the LOCO test cluster being matched.
        random_state: RandomState instance or integer seed.
        max_folds: Keep at most this many of the size-matched folds. A small held-out
            cluster implies many folds (a 15-of-256 cluster gives 17), and training a
            model on every one of them costs far more than the single LOCO fold it is
            compared against, so the fold *size* is preserved while the *count* is
            capped. Pass None to keep all folds.

    Returns:
        One SplitTuple per retained size-matched fold.
    """
    n_splits = max(round(n_total / n_test), 2)
    folds = get_regression_kfold_splits(n_total, n_splits, random_state)
    return folds if max_folds is None else folds[:max_folds]


# The collaborator protocol (github.com/ddomlab/GP_collab, code_python/training/
# scoring_validation.py): every cluster is held out once, and each is paired with
# size-matched random IID splits drawn from these base seeds.
COLLAB_IID_SEEDS = (17, 29, 43, 71, 97)


def get_leave_one_group_out_splits(
    labels: jt.Shaped[np.ndarray, " n"],
) -> dict[object, SplitTuple]:
    """One split per group, holding out every group once -- no size floor, no exceptions.

    Matches sklearn's `LeaveOneGroupOut`, as the collaborator protocol uses: groups come in
    sorted label order, which matters because a group's position sets its matched-IID
    seed (`get_size_matched_random_split`). Unlike `get_loco_splits`, both groups of a
    two-group labelling are held out, so the larger one trains on the smaller.
    """
    labels = np.asarray(labels)
    return {label: _leave_one_out_split(labels, label) for label in np.unique(labels)}


def get_size_matched_random_split(
    n_total: int, n_test: int, group_index: int, base_seed: int
) -> SplitTuple:
    """Random IID split whose test set is the size of one held-out group.

    Reproduces the collaborator's `_fit_predict_score_iid` exactly, so for the same rows
    in the same order it selects the same test indices: the seed is
    `base_seed + group_index * 10_000`, with `group_index` the group's position in
    `get_leave_one_group_out_splits` order.
    """
    if not 0 < n_test < n_total:
        raise ValueError(f"Invalid IID test size {n_test} for {n_total} samples.")
    permutation = np.random.default_rng(base_seed + group_index * 10_000).permutation(n_total)
    return SplitTuple(
        train=np.sort(permutation[n_test:]),
        val=np.empty(0, dtype=int),
        test=np.sort(permutation[:n_test]),
    )
