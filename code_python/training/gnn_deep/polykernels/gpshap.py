"""Stochastic Shapley values for Gaussian process models (GP-SHAP).

Implements Chau et al., "Explaining the Uncertain: Stochastic Shapley Values for Gaussian
Process Models" (arXiv:2305.15167), Propositions 8/9 (GP-SHAP) and 11 (BayesGP-SHAP).

The idea: Shapley values are a *linear* functional of the coalition payoffs, and a GP
posterior stays Gaussian under conditional expectation. So when the payoffs come from a
GP, the explanations are themselves multivariate Gaussian with a closed-form covariance --
credible intervals per feature, plus the covariance *between* features, neither of which
ordinary SHAP can express.

Two kernels appear, and conflating them is the easiest way to get this wrong:

* `k_S`, on the raw feature subspace, estimates the conditional expectation that "removes"
  features outside a coalition. It is ours to choose (RBF here, with the median heuristic)
  and has nothing to do with the model.
* the GP's own posterior, which enters only as the mean `m_ref` and covariance `K_ref`
  evaluated at the reference set.

That separation is what lets this work against a deep-kernel GP whose kernel lives on a
learned embedding rather than on interpretable features: the embedding is never indexed by
feature subset, only the raw reference rows are.

Deviation from the paper's Algorithm 1, and the reason for it: that algorithm uses the
model's inducing points as the conditional-mean-embedding reference set, which requires
indexing them by feature subset. Inducing points for a deep-kernel GP live in embedding
space and cannot be subset that way. We therefore follow Propositions 8/9 verbatim with a
reference set of real data rows, which are in raw feature space and subsettable.

Everything here takes arrays rather than models, so it is testable without a checkpoint.
"""

from itertools import combinations
import math

import numpy as np
import torch
from torch import Tensor

# Weight placeholder for the empty and full coalitions. The paper sets these to infinity to
# enforce efficiency; we instead impose that constraint analytically in
# `shapley_projection`, which is exact rather than approximate and avoids a badly
# conditioned normal-equations solve. The value is never used in the fit.
_ANCHOR_WEIGHT = float("inf")


def coalition_matrix(
    n_features: int, n_coalitions: int = 0, *, rng: np.random.Generator | None = None
) -> tuple[Tensor, Tensor]:
    """Binary coalition matrix `Z` and Shapley kernel weights `W`.

    Args:
        n_coalitions: 0 enumerates all `2**n_features` coalitions exactly, which is what
            makes the estimation-uncertainty term of BayesGP-SHAP vanish. Any other value
            samples that many *interior* coalitions uniformly.

    Returns:
        Z `[n_coalitions, n_features]` of 0/1, and W `[n_coalitions]`. The empty and full
        coalitions are always present and always first -- without them the efficiency
        axiom is unenforceable, so they survive subsampling too.
    """
    if n_features < 1:
        raise ValueError(f"Need at least one feature, got {n_features}.")
    interior = [
        c for size in range(1, n_features) for c in combinations(range(n_features), size)
    ]
    exact = n_coalitions <= 0 or n_coalitions >= len(interior)
    if not exact:
        rng = rng or np.random.default_rng()
        picked = rng.choice(len(interior), size=n_coalitions, replace=False)
        interior = [interior[i] for i in sorted(picked)]

    rows = [()] + interior + [tuple(range(n_features))]
    z = torch.zeros(len(rows), n_features, dtype=torch.float64)
    for i, members in enumerate(rows):
        z[i, list(members)] = 1.0

    sizes = z.sum(dim=1)
    weights = torch.full((len(rows),), _ANCHOR_WEIGHT, dtype=torch.float64)
    mask = (sizes > 0) & (sizes < n_features)
    s = sizes[mask]
    # w(S) = (d - 1) / (C(d, |S|) * |S| * (d - |S|)) -- the Shapley kernel.
    binom = torch.tensor(
        [math.comb(n_features, int(v)) for v in s], dtype=torch.float64
    )
    weights[mask] = (n_features - 1) / (binom * s * (n_features - s))
    return z, weights


def shapley_projection(z: Tensor, weights: Tensor) -> Tensor:
    """Linear map `A` from the payoff vector to Shapley values, `[n_features, n_coalitions]`.

    The paper writes `A = (Z^T W Z)^-1 Z^T W` with infinite weight on the empty and full
    coalitions. Two adjustments make that usable in finite precision:

    * **Centering.** `Z` has no intercept column, so the regression fits
      `nu(S) ~ sum_{i in S} phi_i` and the all-zero empty row can only drive `nu(empty)`
      towards zero rather than being fitted. Payoffs are therefore centered by
      `nu(empty)`, which is exactly SHAP's base value.
    * **Exact efficiency.** Rather than approximating an infinite weight with a large one
      -- which leaves a residual and wrecks the conditioning -- the constraint
      `sum_i phi_i = nu(full) - nu(empty)` is imposed in closed form via its Lagrange
      multiplier.

    Both stay linear in the payoffs, so `A` remains a fixed matrix. That matters: the
    covariance propagation of Theorem 7 is `A V[v] A^T`, which is only valid for a linear
    map.
    """
    interior = torch.isfinite(weights)
    z_int, w_int = z[interior], weights[interior]
    zt_w = z_int.T * w_int
    unconstrained = torch.linalg.solve(zt_w @ z_int, zt_w)

    n_features, n_coalitions = z.shape[1], z.shape[0]
    ones = torch.ones(n_features, dtype=z.dtype)
    m_inv_one = torch.linalg.solve(zt_w @ z_int, ones)
    # g spreads any efficiency shortfall across features in the metric of the WLS problem.
    g = m_inv_one / (ones @ m_inv_one)
    correction = torch.eye(n_features, dtype=z.dtype) - torch.outer(g, ones)

    a = torch.zeros(n_features, n_coalitions, dtype=z.dtype)
    a[:, interior] = correction @ unconstrained
    # Centering by nu(empty) and the constraint value nu(full) - nu(empty), both linear.
    a[:, 0] -= (correction @ unconstrained @ torch.ones(z_int.shape[0], dtype=z.dtype)) + g
    a[:, -1] += g
    return a


def conditional_mean_projections(
    x_ref: Tensor, x_query: Tensor, z: Tensor, *, bandwidth: float, ridge: float,
    coalition_chunk: int = 256,
) -> Tensor:
    """`B[j] = (K_{ref_S ref_S} + lambda I)^-1 k_S(ref_S, query_S)` for each coalition.

    This is the conditional mean embedding of Proposition 8: it estimates
    `E[f(X) | X_S = x_S]` as a weighted combination of the reference rows, using only the
    features inside the coalition.

    Returns `[n_coalitions, n_ref, n_query]`. The empty coalition is left at a uniform
    `1/n_ref`, which is the unconditional mean -- the correct payoff for `nu(empty)`.
    """
    n_ref, n_query = x_ref.shape[0], x_query.shape[0]
    out = torch.empty(z.shape[0], n_ref, n_query, dtype=torch.float64)
    eye = torch.eye(n_ref, dtype=torch.float64) * ridge
    gamma = 1.0 / (2.0 * bandwidth**2)

    for start in range(0, z.shape[0], coalition_chunk):
        block = z[start : start + coalition_chunk]
        for i, mask in enumerate(block):
            active = mask.bool()
            if not active.any():
                # Conditioning on nothing: every reference row contributes equally.
                out[start + i] = 1.0 / n_ref
                continue
            ref_s, query_s = x_ref[:, active], x_query[:, active]
            k_rr = torch.exp(-gamma * torch.cdist(ref_s, ref_s) ** 2)
            k_rq = torch.exp(-gamma * torch.cdist(ref_s, query_s) ** 2)
            # Cholesky rather than a general solve: K + lambda I is symmetric positive
            # definite by construction, and this is the algorithm's dominant cost.
            out[start + i] = torch.cholesky_solve(
                k_rq, torch.linalg.cholesky(k_rr + eye)
            )
    return out


def gpshap(
    m_ref: Tensor, k_ref: Tensor, b: Tensor, z: Tensor, weights: Tensor
) -> tuple[Tensor, Tensor, Tensor]:
    """Mean stochastic Shapley values and the low-rank factor of their covariance.

    Args:
        m_ref: GP posterior mean at the reference rows, `[n_ref]`.
        k_ref: GP posterior covariance at the reference rows, `[n_ref, n_ref]`.
        b: conditional mean projections, `[n_coalitions, n_ref, n_query]`.

    Returns:
        `phi [n_features, n_query]` -- the mean SSVs, which by Proposition 5 equal the
        deterministic Shapley values of the mean game;
        `r [n_features, n_query, n_ref]` -- the low-rank factor, with
        `Cov = R R^T` contracted over the reference axis. Never form the full
        `[d, d, n, n]` covariance; contract on demand instead.
        `payoffs [n_coalitions, n_query]` -- the mean game, needed by BayesGP-SHAP.
    """
    a = shapley_projection(z, weights)
    payoffs = torch.einsum("ijk,j->ik", b, m_ref)
    phi = a @ payoffs
    # Cholesky of the posterior covariance, so that R R^T reproduces A B^T K B A^T
    # without ever forming the n_ref x n_ref product per coalition.
    chol = torch.linalg.cholesky(
        k_ref + 1e-6 * torch.eye(k_ref.shape[0], dtype=k_ref.dtype)
    )
    q = torch.einsum("ijk,jl->ikl", b, chol)
    r = torch.einsum("ij,jkl->ikl", a, q)
    return phi, r, payoffs


def per_query_covariance(r: Tensor) -> Tensor:
    """Feature-by-feature SSV covariance for each query row: `[n_query, d, d]`."""
    return torch.einsum("ikl,jkl->kij", r, r)


def bayes_extra_variance(
    payoffs: Tensor, phi: Tensor, z: Tensor, weights: Tensor, *,
    ell0: float, sigma0_sq: float, rng: np.random.Generator,
) -> Tensor:
    """Additional per-feature variance from estimating Shapley values (Proposition 11).

    This is a different quantity from the GP-SHAP covariance: it is uncertainty from
    estimating the Shapley values, not from the model's posterior.

    It shrinks as coalitions are added but does not vanish at exact enumeration, and the
    reason is worth knowing before reading the numbers. Equation 13 makes
    `s^2 = [weighted residual + ||phi||^2] / l`, so both terms decay only as `1/l`. The
    paper reports ~0.0002 at `l = 2^16`; with a handful of features `l` is far smaller
    (62 interior coalitions at d=6, 8190 at d=13) and the term stays visible. Measured
    here: `s^2` falls 0.063 -> 0.00086 going from d=6 to d=13, i.e. tracking 1/l.

    So a non-zero value at full enumeration is expected, not a defect; what should hold is
    that it decreases monotonically with the coalition budget while the GP-SHAP covariance
    stays flat.

    Returns `[n_features, n_query]` variances to add to the GP-SHAP diagonal.
    """
    # Interior coalitions only. The empty and full rows carry infinite weight -- they are
    # constraints, not observations -- and their residual is zero by construction, so
    # including them would evaluate inf * 0 and poison the whole thing with NaN.
    interior = torch.isfinite(weights)
    z_int, w_int = z[interior], weights[interior]
    n_coalitions = int(interior.sum())

    residual = payoffs[interior] - z_int @ phi
    weighted = torch.einsum("ik,i,ik->k", residual, w_int, residual)
    s_sq = (weighted + (phi * phi).sum(dim=0)) / n_coalitions

    # sigma^2 ~ Scaled-Inv-chi^2(ell0 + n_Z, (ell0 * sigma0^2 + n_Z * s^2) / (ell0 + n_Z)),
    # sampled per query row, then scaled by the WLS covariance (Z^T W Z)^-1.
    dof = ell0 + n_coalitions
    scale = (ell0 * sigma0_sq + n_coalitions * s_sq.numpy()) / dof
    draw = scale * dof / rng.chisquare(dof, size=s_sq.shape[0])
    wls_cov = torch.linalg.inv((z_int.T * w_int) @ z_int)
    return torch.outer(torch.diagonal(wls_cov), torch.from_numpy(draw).to(z.dtype))
