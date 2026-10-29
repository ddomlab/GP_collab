"""Iterative stratification utilities for multi-label data.

This code has been migrated from polykernels.splits to encapsulate the stratification
logic.
"""

import jaxtyping as jt
import numpy as np

from polykernels.types import IntArray
from polykernels.utils import cast_random_state


def _split_count_by_ratios(
    count: int, ratios: jt.Float[np.ndarray, " f"]
) -> IntArray[" f"]:
    """Split an integer based on proportions (9, [0.5, 0.5] -> [5, 4])."""
    parts = np.floor(ratios * count).astype(int)
    remaining = count - parts.sum()
    sorted_indices = np.argsort(parts)
    for index in range(remaining):
        parts[sorted_indices[index]] += 1
    if np.sum(parts) != count:
        raise ValueError(f"Splitting {count} by ratios=({ratios}) failed")
    return parts


def _break_ties_randomly(
    indices: IntArray[" n"], random_state: np.random.RandomState
) -> int:
    """Choose randomly from indices to break ties."""
    if not indices.shape[0]:
        raise ValueError("Found empty indices!")
    return random_state.choice(indices) if len(indices) > 1 else indices[0]


def _verify_ratios(ratios: np.ndarray) -> None:
    """Raises ValueError if ratios not positive or do not sum to 1."""
    if not np.isclose(np.sum(ratios), 1.0):
        raise ValueError(f"Ratios should sum up to 1, got {ratios.sum()}")
    if not np.all(ratios > 0):
        raise ValueError(f"Ratios should be positive, got {ratios}")


def make_second_order_matrix(y: np.ndarray) -> np.ndarray:
    """Expand multilabel array to include second-order relationships (pairs of active labels).

    As described in Szymanski & Kajdanowicz (2017), "A Network Perspective on
    Stratification of Multi-Label Data" (arXiv:1704.08756).
    Empty pairs (all-zero columns) are filtered out.

    Args:
        y: 2D binary label target matrix of shape (n_samples, n_labels).

    Returns:
        Expanded 2D array of shape (n_samples, n_labels + n_nonempty_pairs).
    """
    n_samples, n_labels = y.shape
    triu_indices = np.triu_indices(n_labels, k=1)
    pairs = y[:, triu_indices[0]] * y[:, triu_indices[1]]
    non_empty = np.any(pairs, axis=0)
    return np.hstack([y, pairs[:, non_empty]])


def iterative_stratification(
    y: jt.Float[np.ndarray, " b l"],
    ratios: jt.Float[np.ndarray, " f"],
    order: int = 1,
    random_state: int | np.random.RandomState | None = None,
) -> IntArray[" b"]:
    """Stratify (balance) labels in an iterative fashion.

    Balanced splits contain approximately the same proportion of active labels
    on each split fold. Non-binary targets can be binarized and concatenated to
    a single multilabel array to stratify across multiple columns.

    Supports both first-order stratification (Sechidis et al., 2011) and
    second-order stratification (Szymanski & Kajdanowicz, 2017) by considering
    label pairs.

    Citations:
        - Sechidis, K., Tsoumakas, G., & Vlahavas, I. (2011). On the stratification
          of multi-label data. Machine Learning and Knowledge Discovery in Databases,
          Lecture Notes in Computer Science, vol 6913. Springer.
          DOI: https://doi.org/10.1007/978-3-642-23808-6_10
        - Szymański, P., & Kajdanowicz, T. (2017). A Network Perspective on
          Stratification of Multi-Label Data. First International Workshop on
          Learning with Imbalanced Domains: Theory and Applications, PMLR 74:22-35.
          arXiv: https://arxiv.org/abs/1704.08756

    Args:
        y: 2D binary label target matrix of shape (n_samples, n_labels).
        ratios: 1D array of proportions for each fold (must sum to 1.0).
        order: Stratification order. 1 for first-order, 2 for second-order.
        random_state: RandomState instance or integer seed for reproducibility.

    Returns:
        folds: 1D integer array mapping each sample index to its fold assignment.
    """
    _verify_ratios(ratios)
    if order == 2:
        y = make_second_order_matrix(y)
    elif order != 1:
        raise ValueError(f"Supported orders are 1 and 2, got {order}")

    labels = y.astype(bool)
    # 1. Setup counting variables for folds and labels.
    n_samples = labels.shape[0]
    rng = cast_random_state(random_state)
    folds = np.zeros(n_samples, dtype=int)
    fold_counts = _split_count_by_ratios(n_samples, ratios)
    label_counts = np.sum(labels, axis=0)
    label_mask = label_counts > 0
    nlabels_per_fold = np.outer(label_counts, ratios)
    not_processed_mask = np.ones(n_samples, dtype=bool)
    sample_indices = np.arange(n_samples)

    # 2. Balance samples for which there are positive labels.
    for _ in range(label_mask.sum()):
        # 2a. Identify sparsest label and unassigned samples with this label.
        lowest_count = np.min(label_counts[label_mask])
        label_indices = np.where(label_counts == lowest_count)[0]
        label_index = _break_ties_randomly(label_indices, rng)
        candidates = sample_indices[labels[:, label_index] & not_processed_mask]
        # 2b. For each sample, assign to fold with highest number of
        # unassigned examples.
        for sample_index in candidates:
            highest_unassigned = nlabels_per_fold[label_index].max()
            fold_indices = np.where(
                nlabels_per_fold[label_index] == highest_unassigned
            )[0]
            if len(fold_indices) > 1:
                temp_fold_indices = np.where(
                    fold_counts[fold_indices] == fold_counts[fold_indices].max()
                )[0]
                temp_fold_index = _break_ties_randomly(temp_fold_indices, rng)
                fold_index = fold_indices[temp_fold_index]
            else:
                fold_index = fold_indices[0]
            folds[sample_index] = fold_index

            # Update fold counts and label counts.
            not_processed_mask[sample_index] = False
            nlabels_per_fold[labels[sample_index], fold_index] -= 1
            fold_counts[fold_index] -= 1
            label_counts -= labels[sample_index]

        # Update label mask.
        label_mask = label_counts > 0

        # If all labels or samples have been processed, exit early.
        if not (not_processed_mask.any() and label_mask.any()):
            break

    # 3. Randomly assign leftover samples (no positive labels) to each fold.
    leftover_indices = sample_indices[not_processed_mask]
    rng.shuffle(leftover_indices)
    for sample_index in leftover_indices:
        fold_index = np.argmax(fold_counts)
        folds[sample_index] = fold_index
        fold_counts[fold_index] -= 1

    return folds


def map_fold_to_indices(
    folds: IntArray[" b"], fold_index: int
) -> tuple[IntArray[" n"], IntArray[" k"]]:
    """Map a single fold to train/test split indices."""
    indices = np.arange(len(folds))
    test_mask = folds == fold_index
    return indices[~test_mask], indices[test_mask]


def get_multilabel_train_test_split(
    labels: jt.Float[np.ndarray, " b l"],
    test_size: float,
    order: int = 1,
    random_state: int | np.random.RandomState | None = None,
) -> tuple[IntArray[" n"], IntArray[" k"]]:
    """Gets indices for a multilabel train/test split.

    Args:
        labels: 2D target matrix of shape (n_samples, n_labels).
        test_size: Proportion of the dataset to include in the test set.
        order: Stratification order (1 or 2).
        random_state: RandomState instance or integer seed.

    Returns:
        train_indices, test_indices
    """
    rng = cast_random_state(random_state)
    ratios = np.array([1.0 - test_size, test_size])
    folds = iterative_stratification(labels, ratios, order=order, random_state=rng)
    return map_fold_to_indices(folds, 1)
