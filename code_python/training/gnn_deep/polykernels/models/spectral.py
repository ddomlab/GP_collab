"""Spectral normalization with an upper bound, for distance-preserving feature extractors.

SNGP (Liu et al., arXiv:2006.10108) and DUE (van Amersfoort et al., arXiv:2102.11409)
make a GP head distance aware by bounding the spectral norm of every weight inside
residual blocks, h <- h + f(h) with Lip(f) < 1: the block can then neither collapse two
inputs together nor blow them apart, so a test point far from the training data stays far
from it in the latent space the GP sees.

`SpectralNormBound` is the TensorFlow `SpectralNormalization(norm_multiplier=c)` rule,
W <- W * c / sigma(W) only when sigma(W) > c, rather than the always-normalize-to-1 rule
of `torch.nn.utils.parametrizations.spectral_norm`.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F  # noqa: N812
from torch.nn.utils import parametrize


class SpectralNormBound(nn.Module):
    """Parametrization W -> W / max(1, sigma(W) / bound).

    sigma is estimated by power iteration, one step per training forward (Miyato et al.,
    arXiv:1802.05957); eval reuses the stored singular vectors. The max(1, .) leaves a
    weight already inside the bound untouched, and keeps a zero-initialized weight (AdaLN's
    context projection) finite where the built-in W / sigma would be 0 / 0.
    """

    def __init__(
        self,
        weight: torch.Tensor,
        bound: float,
        n_power_iterations: int = 1,
        eps: float = 1e-12,
    ) -> None:
        super().__init__()
        if bound <= 0:
            raise ValueError(f"Spectral norm bound must be positive; got {bound}.")
        self.bound = bound
        self.n_power_iterations = n_power_iterations
        self.eps = eps
        w = weight.detach().flatten(1)
        self.register_buffer("_u", F.normalize(w.new_empty(w.shape[0]).normal_(), dim=0, eps=eps))
        self.register_buffer("_v", F.normalize(w.new_empty(w.shape[1]).normal_(), dim=0, eps=eps))
        # Converge the estimate up front, so sigma is already accurate on the first step.
        self._power_iteration(w, 15)

    @torch.no_grad()
    def _power_iteration(self, w: torch.Tensor, n_iterations: int) -> None:
        # A zero weight has no singular direction to converge to and would zero u and v
        # for good; keep the random start until the weight moves.
        if not torch.any(w):
            return
        for _ in range(n_iterations):
            self._v.copy_(F.normalize(w.t() @ self._u, dim=0, eps=self.eps))
            self._u.copy_(F.normalize(w @ self._v, dim=0, eps=self.eps))

    def forward(self, weight: torch.Tensor) -> torch.Tensor:
        w = weight.flatten(1)
        u, v = self._u, self._v
        if self.training:
            self._power_iteration(w.detach(), self.n_power_iterations)
            # Cloned because a weight-tied layer runs this several times per step, and the
            # next call's in-place update would corrupt vectors saved for this backward.
            u, v = self._u.clone(), self._v.clone()
        sigma = torch.dot(u, w @ v)
        return weight / torch.clamp(sigma / self.bound, min=1.0)


def apply_spectral_norm(module: nn.Module, bound: float) -> int:
    """Bound the spectral norm of every trainable `nn.Linear` weight under `module`.

    Lazy layers must already be materialized (they become `nn.Linear` on their first
    forward). Frozen weights, such as a fixed random input projection, are left alone, as
    are weights already parametrized, so calling this twice is harmless. PyG's own `Linear`
    is deliberately not matched: its custom `__deepcopy__` cannot copy a parametrized weight,
    which breaks `AveragedModel`; swap such layers for `nn.Linear` first. Returns the number
    of weights newly bounded.
    """
    targets = [
        m
        for m in module.modules()
        if isinstance(m, nn.Linear)
        and m.weight.requires_grad
        and not parametrize.is_parametrized(m, "weight")
    ]
    for m in targets:
        parametrize.register_parametrization(m, "weight", SpectralNormBound(m.weight, bound))
    return len(targets)


def spectral_norms(module: nn.Module) -> list[float]:
    """Largest singular value of every effective (post-parametrization) Linear weight.

    Call with the module in eval mode: in train mode, reading a bounded weight also
    advances its power iteration.
    """
    return [
        float(torch.linalg.matrix_norm(m.weight.detach(), ord=2))
        for m in module.modules()
        if isinstance(m, nn.Linear) and not isinstance(m.weight, nn.parameter.UninitializedParameter)
    ]
