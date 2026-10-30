"""Reusable grouped GP-SHAP calculations.

This module is intentionally independent of GPyTorch and Pyro.  A model adapter
only needs to supply the latent posterior mean and covariance on a reference
set.  The functions here then construct conditional-mean-embedding payoffs,
project them to grouped SHAP values, and propagate posterior uncertainty.
"""

from dataclasses import dataclass
from itertools import combinations
import math
from typing import Sequence

import numpy as np
import torch


def _as_double_cpu(value) -> torch.Tensor:
    """Detach tensor-like input and place it on the shared CPU compute path."""
    if isinstance(value, torch.Tensor):
        return value.detach().to(device="cpu", dtype=torch.float64)
    return torch.as_tensor(value, dtype=torch.float64)


@dataclass(frozen=True)
class GPShapResult:
    """Query-major grouped GP-SHAP results.

    Attributes:
        values: Posterior-mean SHAP values, ``[n_query, n_players]``.
        standard_deviations: Marginal posterior standard deviations with the
            same shape as ``values``.
        covariance: Per-query player covariance,
            ``[n_query, n_players, n_players]``.
        payoffs: Coalition payoff means, ``[n_coalitions, n_query]``.
        coalitions: Binary coalition design, ``[n_coalitions, n_players]``.
        weights: Kernel-SHAP weights for the coalition design.
        efficiency_error: Maximum absolute efficiency residual.
    """

    values: torch.Tensor
    standard_deviations: torch.Tensor
    covariance: torch.Tensor
    payoffs: torch.Tensor
    coalitions: torch.Tensor
    weights: torch.Tensor
    efficiency_error: float


def grouped_feature_players(
    feature_names: Sequence[object],
    feature_group_indices: dict[str, Sequence[int]],
    *,
    structural_prefix: str = "fp_",
) -> tuple[list[str], list[list[int]]]:
    """Build interpretable players from structural groups and scalar columns.

    Every non-empty group whose name starts with ``structural_prefix`` becomes
    one player.  Columns in the ``count`` group, and any otherwise ungrouped
    columns, become individual players.
    """
    names_by_column = [str(name) for name in feature_names]
    n_columns = len(names_by_column)
    players: list[tuple[str, list[int]]] = []
    covered: set[int] = set()

    def add_player(name: str, columns: Sequence[int]) -> None:
        normalized = [int(column) for column in columns]
        invalid = [column for column in normalized if not 0 <= column < n_columns]
        if invalid:
            raise IndexError(
                f"GP-SHAP group '{name}' contains invalid column indices {invalid}."
            )
        overlap = covered.intersection(normalized)
        if overlap:
            raise ValueError(
                "GP-SHAP feature groups overlap at column indices "
                f"{sorted(overlap)}."
            )
        players.append((str(name), normalized))
        covered.update(normalized)

    for name in sorted(feature_group_indices):
        columns = feature_group_indices[name]
        if name.startswith(structural_prefix) and columns:
            add_player(name, columns)

    for column in sorted(feature_group_indices.get("count", [])):
        column = int(column)
        name = names_by_column[column] if 0 <= column < n_columns else "count"
        add_player(name, [column])

    for column, name in enumerate(names_by_column):
        if column not in covered:
            add_player(name, [column])

    player_names = [name for name, _ in players]
    if len(player_names) != len(set(player_names)):
        raise ValueError(
            f"GP-SHAP player names must be unique; got {player_names}."
        )
    return player_names, [columns for _, columns in players]


def coalition_matrix(
    n_players: int,
    n_coalitions: int = 0,
    *,
    rng: np.random.Generator | None = None,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Return a Kernel-SHAP coalition design and its weights.

    ``n_coalitions`` counts interior coalitions.  A non-positive value, or a
    value at least as large as the complete interior power set, enumerates the
    design exactly.  Empty and full coalitions are always included as exact
    efficiency constraints.  Sampled designs retain every singleton so the
    weighted regression has full column rank.
    """
    if n_players < 1:
        raise ValueError("GP-SHAP needs at least one player.")

    if rng is None:
        rng = np.random.default_rng()
    total_interior = (1 << n_players) - 2
    exact = n_coalitions <= 0 or n_coalitions >= total_interior
    if exact:
        interior = [
            coalition
            for size in range(1, n_players)
            for coalition in combinations(range(n_players), size)
        ]
    else:
        if n_coalitions < n_players:
            raise ValueError(
                "n_coalitions must be at least the number of GP-SHAP players "
                "so the sampled regression is identifiable."
            )
        selected = {(player,) for player in range(n_players)}
        while len(selected) < n_coalitions:
            mask = rng.random(n_players) < 0.5
            size = int(mask.sum())
            if 0 < size < n_players:
                selected.add(tuple(np.flatnonzero(mask).tolist()))
        interior = sorted(selected, key=lambda members: (len(members), members))

    rows = [()] + interior + [tuple(range(n_players))]
    coalitions = torch.zeros((len(rows), n_players), dtype=torch.float64)
    for row, members in enumerate(rows):
        if members:
            coalitions[row, list(members)] = 1.0

    sizes = coalitions.sum(dim=1)
    weights = torch.full((len(rows),), float("inf"), dtype=torch.float64)
    interior_mask = (sizes > 0) & (sizes < n_players)
    if interior_mask.any():
        interior_sizes = sizes[interior_mask]
        binomial = torch.tensor(
            [math.comb(n_players, int(size)) for size in interior_sizes],
            dtype=torch.float64,
        )
        weights[interior_mask] = (n_players - 1) / (
            binomial * interior_sizes * (n_players - interior_sizes)
        )
    return coalitions, weights


def shapley_projection(
    coalitions: torch.Tensor,
    weights: torch.Tensor,
) -> torch.Tensor:
    """Return the exact-efficiency linear payoff-to-SHAP projection."""
    n_players, n_rows = coalitions.shape[1], coalitions.shape[0]
    if n_players == 1:
        return torch.tensor([[-1.0, 1.0]], dtype=coalitions.dtype)

    interior = torch.isfinite(weights)
    z_interior, w_interior = coalitions[interior], weights[interior]
    zt_w = z_interior.T * w_interior
    normal = zt_w @ z_interior
    unconstrained = torch.linalg.solve(normal, zt_w)

    ones = torch.ones(n_players, dtype=coalitions.dtype)
    normal_inv_one = torch.linalg.solve(normal, ones)
    constraint_direction = normal_inv_one / (ones @ normal_inv_one)
    correction = torch.eye(n_players, dtype=coalitions.dtype) - torch.outer(
        constraint_direction, ones
    )

    projection = torch.zeros((n_players, n_rows), dtype=coalitions.dtype)
    projection[:, interior] = correction @ unconstrained
    projection[:, 0] -= (
        correction
        @ unconstrained
        @ torch.ones(z_interior.shape[0], dtype=coalitions.dtype)
    ) + constraint_direction
    projection[:, -1] += constraint_direction
    return projection


def stable_cholesky(
    matrix: torch.Tensor,
    *,
    initial_jitter: float = 1e-8,
    max_tries: int = 6,
) -> torch.Tensor:
    """Compute a Cholesky factor using bounded adaptive diagonal jitter."""
    eye = torch.eye(matrix.shape[-1], dtype=matrix.dtype, device=matrix.device)
    jitter = initial_jitter
    for _ in range(max_tries):
        factor, info = torch.linalg.cholesky_ex(matrix + jitter * eye)
        if int(info.max()) == 0:
            return factor
        jitter *= 10.0
    raise RuntimeError(
        "GP-SHAP Cholesky failed after adding jitter up to "
        f"{jitter / 10.0:g}."
    )


def grouped_conditional_mean_projections(
    x_reference: torch.Tensor,
    x_query: torch.Tensor,
    coalitions: torch.Tensor,
    player_columns: Sequence[Sequence[int]],
    *,
    ridge: float,
) -> torch.Tensor:
    """Return coalition CME weights, ``[coalition, reference, query]``."""
    if ridge <= 0:
        raise ValueError("ridge must be positive.")
    if len(player_columns) != coalitions.shape[1]:
        raise ValueError(
            "player_columns must contain one column group per coalition player."
        )

    x_reference = _as_double_cpu(x_reference)
    x_query = _as_double_cpu(x_query)
    n_reference, n_query = x_reference.shape[0], x_query.shape[0]
    projections = torch.empty(
        (coalitions.shape[0], n_reference, n_query), dtype=torch.float64
    )

    center = x_reference.mean(dim=0)
    scale = x_reference.std(dim=0, unbiased=False)
    scale = torch.where(scale > 1e-12, scale, torch.ones_like(scale))
    x_reference = (x_reference - center) / scale
    x_query = (x_query - center) / scale
    eye = torch.eye(n_reference, dtype=torch.float64)

    for coalition_index, player_mask in enumerate(coalitions):
        active_players = torch.flatnonzero(player_mask).tolist()
        if not active_players:
            projections[coalition_index] = 1.0 / n_reference
            continue

        active_columns = sorted(
            {
                int(column)
                for player in active_players
                for column in player_columns[player]
            }
        )
        reference_subset = x_reference[:, active_columns]
        query_subset = x_query[:, active_columns]
        pairwise = torch.pdist(reference_subset)
        positive = pairwise[pairwise > 0]
        bandwidth = float(positive.median()) if len(positive) else 1.0
        gamma = 1.0 / (2.0 * bandwidth**2)
        reference_kernel = torch.exp(
            -gamma * torch.cdist(reference_subset, reference_subset).square()
        )
        query_kernel = torch.exp(
            -gamma * torch.cdist(reference_subset, query_subset).square()
        )
        factor = stable_cholesky(reference_kernel + ridge * eye)
        projections[coalition_index] = torch.cholesky_solve(query_kernel, factor)
    return projections


def mixture_posterior_moments(
    mean_samples: torch.Tensor,
    covariance_samples: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Combine conditional Gaussian draws using total expectation/covariance."""
    mean_samples = _as_double_cpu(mean_samples)
    covariance_samples = _as_double_cpu(covariance_samples)
    if mean_samples.ndim != 2:
        raise ValueError(
            "mean_samples must have shape [n_draws, n_reference]."
        )
    if mean_samples.shape[0] < 1:
        raise ValueError("At least one posterior draw is required.")
    expected_covariance_shape = (
        mean_samples.shape[0],
        mean_samples.shape[1],
        mean_samples.shape[1],
    )
    if covariance_samples.shape != expected_covariance_shape:
        raise ValueError(
            "covariance_samples must have shape "
            f"{expected_covariance_shape}; got {tuple(covariance_samples.shape)}."
        )

    posterior_mean = mean_samples.mean(dim=0)
    conditional_covariance = covariance_samples.mean(dim=0)
    centered_means = mean_samples - posterior_mean
    parameter_covariance = centered_means.T @ centered_means / mean_samples.shape[0]
    posterior_covariance = conditional_covariance + parameter_covariance
    posterior_covariance = 0.5 * (
        posterior_covariance + posterior_covariance.T
    )
    return posterior_mean, posterior_covariance


def explain_posterior(
    posterior_mean: torch.Tensor,
    posterior_covariance: torch.Tensor,
    x_reference: torch.Tensor,
    x_query: torch.Tensor,
    player_columns: Sequence[Sequence[int]],
    *,
    n_coalitions: int = 0,
    ridge: float = 1e-3,
    rng: np.random.Generator | None = None,
    efficiency_tolerance: float = 1e-7,
) -> GPShapResult:
    """Explain a latent Gaussian posterior using grouped GP-SHAP.

    Posterior moments must be evaluated on the same rows as ``x_reference``.
    This function is model-library agnostic and performs no fitting.
    """
    posterior_mean = _as_double_cpu(posterior_mean).reshape(-1)
    posterior_covariance = _as_double_cpu(posterior_covariance)
    x_reference = _as_double_cpu(x_reference)
    x_query = _as_double_cpu(x_query)

    n_reference = x_reference.shape[0]
    if posterior_mean.shape != (n_reference,):
        raise ValueError(
            "posterior_mean must have one value per reference row; got "
            f"{tuple(posterior_mean.shape)} for {n_reference} rows."
        )
    if posterior_covariance.shape != (n_reference, n_reference):
        raise ValueError(
            "posterior_covariance must be square over the reference rows; got "
            f"{tuple(posterior_covariance.shape)}."
        )
    if x_reference.ndim != 2 or x_query.ndim != 2:
        raise ValueError("x_reference and x_query must be two-dimensional.")
    if x_reference.shape[1] != x_query.shape[1]:
        raise ValueError("Reference and query rows must have the same columns.")

    coalitions, weights = coalition_matrix(
        len(player_columns), n_coalitions, rng=rng
    )
    projection = shapley_projection(coalitions, weights)
    cme = grouped_conditional_mean_projections(
        x_reference,
        x_query,
        coalitions,
        player_columns,
        ridge=ridge,
    )

    payoffs = torch.einsum("crq,r->cq", cme, posterior_mean)
    player_major_values = projection @ payoffs
    posterior_factor = stable_cholesky(posterior_covariance)
    payoff_factor = torch.einsum("crq,rs->cqs", cme, posterior_factor)
    shap_factor = torch.einsum("pc,cqs->pqs", projection, payoff_factor)
    covariance = torch.einsum("pqs,rqs->qpr", shap_factor, shap_factor)
    standard_deviations = torch.sqrt(
        torch.clamp(torch.diagonal(covariance, dim1=1, dim2=2), min=0.0)
    )

    efficiency_error = torch.max(
        torch.abs(
            player_major_values.sum(dim=0) - (payoffs[-1] - payoffs[0])
        )
    ).item()
    payoff_scale = max(1.0, float(torch.max(torch.abs(payoffs))))
    if efficiency_error > efficiency_tolerance * payoff_scale:
        raise RuntimeError(
            "GP-SHAP efficiency check failed: feature attributions do not "
            "sum to the payoff span "
            f"(max error {efficiency_error:.3e})."
        )

    return GPShapResult(
        values=player_major_values.T,
        standard_deviations=standard_deviations,
        covariance=covariance,
        payoffs=payoffs,
        coalitions=coalitions,
        weights=weights,
        efficiency_error=efficiency_error,
    )


def mean_absolute_importance(
    values: torch.Tensor | np.ndarray,
    player_names: Sequence[str],
) -> dict[str, float]:
    """Aggregate query-major local SHAP values into global importances."""
    values_array = (
        _as_double_cpu(values).numpy()
        if isinstance(values, torch.Tensor)
        else np.asarray(values)
    )
    if values_array.ndim != 2 or values_array.shape[1] != len(player_names):
        raise ValueError(
            "values must have shape [n_query, len(player_names)]; got "
            f"{values_array.shape}."
        )
    importance = np.abs(values_array).mean(axis=0)
    return {
        str(name): float(value)
        for name, value in zip(player_names, importance, strict=True)
    }


__all__ = [
    "GPShapResult",
    "coalition_matrix",
    "explain_posterior",
    "grouped_conditional_mean_projections",
    "grouped_feature_players",
    "mean_absolute_importance",
    "mixture_posterior_moments",
    "shapley_projection",
    "stable_cholesky",
]
