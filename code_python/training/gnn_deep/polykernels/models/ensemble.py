"""Deep ensemble of GNN-MLPs with a per-target variance head (GNN-MLP-Ens).

This gives GNN-MLP a predictive std, based on Lakshminarayanan et al. (2017), "Simple and Scalable
Predictive Uncertainty Estimation using Deep Ensembles":

* every member is an unchanged `GNN` whose head is twice as wide: per target, a mean and
  a raw variance mapped through softplus (`split_mean_var`);
* members are trained independently on the same rows, differing only in initialization
  and data order, with the beta-NLL loss (`BetaNLLLoss`);
* the prediction is the equal-weight Gaussian mixture of the members
  (`mixture_mean_std`). Its variance is the members' average noise variance plus their
  disagreement, so it covers both noise in the data and model uncertainty.

`GNNEnsemble.predict` returns `(mean, std)` like the GP heads, so the regression pipeline
evaluates it through the same path as them.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F  # noqa: N812
import torch_geometric as pyg

from polykernels.models.gnn import GNN

N_MEMBERS = 5
# 0.5 is the value Seitzer et al. recommend; see `BetaNLLLoss`.
BETA = 0.5
# Floor on predicted variance in standardized target units, so a member cannot drive its
# variance to zero and the NLL to minus infinity.
VAR_FLOOR = 1e-6


def split_mean_var(out: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    """Split a member's `[..., 2T]` output into a `[..., T]` mean and positive variance."""
    if out.shape[-1] % 2:
        raise ValueError(f"Expected an even last dimension (mean, variance), got {out.shape}.")
    mean, raw_var = out.chunk(2, dim=-1)
    return mean, F.softplus(raw_var) + VAR_FLOOR


class BetaNLLLoss(nn.Module):
    """Gaussian NLL with each element weighted by its detached variance to the power beta.

    Seitzer et al. (2022), "On the Pitfalls of Heteroscedastic Uncertainty Estimation with
    Probabilistic Neural Networks". Plain NLL divides each residual's gradient by its
    predicted variance, so a network can stall by inflating the variance on hard points
    instead of fitting their mean. The detached `var**beta` weight undoes part of that
    scaling while leaving the optimum of both the mean and the variance unchanged: beta=0
    is plain NLL, beta=1 gives the mean the gradient of MSE.

    Written out rather than via `F.gaussian_nll_loss`, whose variance clamp blocks the
    gradient below its eps.
    """

    def __init__(self, beta: float = BETA) -> None:
        super().__init__()
        self.beta = beta

    def forward(self, out: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
        mean, var = split_mean_var(out)
        nll = 0.5 * (torch.log(var) + (target - mean) ** 2 / var)
        if self.beta:
            nll = nll * var.detach() ** self.beta
        return nll.mean()


def mixture_mean_std(
    means: torch.Tensor, variances: torch.Tensor
) -> tuple[torch.Tensor, torch.Tensor]:
    """Mean and std of an equal-weight Gaussian mixture, over the leading member axis.

    The variance is written in its centred form, mean of the variances plus the variance
    of the means, rather than E[mu^2] - E[mu]^2, which loses precision to cancellation.
    """
    mean = means.mean(dim=0)
    var = variances.mean(dim=0) + means.var(dim=0, unbiased=False)
    return mean, var.sqrt()


class GNNEnsemble(nn.Module):
    """Independently trained `GNN` members, each predicting a mean and variance per target.

    Note:
        The members are held as `members`, never `module`: `unwrap_model` would mistake
        an attribute of that name for an `AveragedModel` wrapper.
    """

    def __init__(self, members: list[GNN]) -> None:
        super().__init__()
        self.members = nn.ModuleList(members)

    def member_outputs(self, data: pyg.data.Data) -> torch.Tensor:
        """Every member's raw `[n_graphs, 2T]` output, stacked to `[M, n_graphs, 2T]`."""
        return torch.stack([member(data) for member in self.members], dim=0)

    def predict(self, data: pyg.data.Data) -> tuple[torch.Tensor, torch.Tensor]:
        """Mixture `(mean, std)`, each `[n_graphs, T]`, in standardized target units."""
        means, variances = split_mean_var(self.member_outputs(data))
        return mixture_mean_std(means, variances)

    def forward(self, data: pyg.data.Data) -> torch.Tensor:
        """Point prediction (the mixture mean). Runs every member, so one pass
        materializes all of their lazy layers."""
        return self.predict(data)[0]
