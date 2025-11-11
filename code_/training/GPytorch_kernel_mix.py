import gpytorch
import numpy as np
import pandas as pd
import torch
# import torch.nn as nn
from sklearn.base import BaseEstimator
from sklearn.model_selection import train_test_split
from torch.utils import benchmark
from gpytorch.constraints import Positive
from gpytorch.kernels import Kernel
from gpytorch.kernels.rbf_kernel import postprocess_rbf
from gpytorch.priors import NormalPrior
from torch.distributions import InverseGamma
from gpytorch.kernels import Kernel, MaternKernel, ScaleKernel
# from torch_geometric.data import Batch
# from torch_geometric.loader import DataLoader

# from pytorch_mpnn import DMPNNPredictor, RevIndexedData, smiles2data


def binary_batch_tanimoto_sim(
        x1: torch.Tensor, x2: torch.Tensor, eps: float = 1e-6
) -> torch.Tensor:
    """
    Tanimoto between two batched tensors, across last 2 dimensions.
    eps argument ensures numerical stability if all zero tensors are added.
    """
    # Tanimoto distance is proportional to (<x, y>) / (||x||^2 + ||y||^2 - <x, y>) where x and y are bit vectors
    assert x1.ndim >= 2 and x2.ndim >= 2
    dot_prod = torch.matmul(x1, torch.transpose(x2, -1, -2))
    x1_sum = torch.sum(x1 ** 2, dim=-1, keepdims=True)
    x2_sum = torch.sum(x2 ** 2, dim=-1, keepdims=True)
    return (dot_prod + eps) / (
            eps + x1_sum + torch.transpose(x2_sum, -1, -2) - dot_prod
    )


def weighted_batch_tanimoto_similarity(x1: torch.Tensor, x2: torch.Tensor) -> torch.Tensor:
    assert x1.ndim >= 2 and x2.ndim >= 2, "Input tensors must have at least 2 dimensions."

    pairwise_min = torch.min(x1.unsqueeze(-2), x2.unsqueeze(-3))
    pairwise_max = torch.max(x1.unsqueeze(-2), x2.unsqueeze(-3))

    min_sum = torch.sum(pairwise_min, dim=-1)
    max_sum = torch.sum(pairwise_max, dim=-1)
    similarity = torch.where(max_sum == 0, torch.zeros_like(max_sum), min_sum / max_sum)

    return similarity
    


def _weighted_tanimoto_similarity(x1, x2):
    """
    Compute 1 - Tanimoto similarity between x1 and x2.
    Supports batched inputs.
    """
    x1 = x1.unsqueeze(-2)  # (..., N, 1, D)
    x2 = x2.unsqueeze(-3)  # (..., 1, M, D)

    numerator = torch.min(x1, x2).sum(dim=-1)
    denominator = torch.max(x1, x2).sum(dim=-1)
    tanimoto_sim = torch.where(denominator == 0, torch.zeros_like(denominator), numerator / denominator)
    tanimoto_dist = 1.0 - tanimoto_sim
    return tanimoto_dist



class TanimotoRBF(Kernel):
    has_lengthscale = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)


    def forward(self, x1, x2, diag=False, **params):
        """
        Standard RBF structure: exp(-0.5 * D^2 / lengthscale^2),
        but using Tanimoto distance instead of Euclidean.
        """
        if diag:
            return x1.new_ones(x1.shape[:-1])  # diagonal is 1

        # Normalize by lengthscale (same as RBFKernel)
        x1_ = x1.div(self.lengthscale)
        x2_ = x2.div(self.lengthscale)

        # Compute Tanimoto distance
        tanimoto_dist = _weighted_tanimoto_similarity(x1_, x2_)

        # Apply RBF exponentiation (same as postprocess_rbf)
        covar = postprocess_rbf(tanimoto_dist)

        return covar




class MixingKernel(Kernel):
    """
    Product of TanimotoRBF (for fingerprint features)
    and Matern(3/2) (for continuous features).
    All lengthscales are learnable with priors.
    """

    def __init__(self, feat_idx: dict, **kwargs):
        """
        feat_idx: dict with keys 'fp' and 'count', each a list of column indices.
                  Example: {'fp': [0,1,2,3], 'count': [4,5,6]}
        """
        super().__init__(**kwargs)

        self.feat_idx = feat_idx

        # --- Priors ---
        l_prior = InverseGamma(5.0, 5.0)
        sigma_prior = NormalPrior(0.0, 1.0)
        sigma_noise_prior = NormalPrior(0.0, 1.0)

        # --- Fingerprint kernel (TanimotoRBF) ---
        self.fp_kernel = ScaleKernel(
            TanimotoRBF(ard_num_dims=len(feat_idx.get('fp', []))),
            outputscale_prior=sigma_prior,
        )
        self.fp_kernel.base_kernel.register_prior(
            "lengthscale_prior", l_prior, "lengthscale"
        )

        # --- Continuous kernel (Matern 3/2) ---
        self.cont_kernel = ScaleKernel(
            MaternKernel(nu=1.5, ard_num_dims=len(feat_idx.get('count', []))),
            outputscale_prior=sigma_prior,
        )
        self.cont_kernel.base_kernel.register_prior(
            "lengthscale_prior", l_prior, "lengthscale"
        )

        # Store priors for convenience if used later in model
        self.sigma_noise_prior = sigma_noise_prior

    def forward(self, x1, x2, diag=False, **params):
        """
        Compute product kernel between subsets of x1/x2.
        """
        fp_idx = self.feat_idx.get('fp', [])
        cont_idx = self.feat_idx.get('count', [])

        x1_fp = x1[..., fp_idx] if fp_idx else None
        x2_fp = x2[..., fp_idx] if fp_idx else None
        x1_cont = x1[..., cont_idx] if cont_idx else None
        x2_cont = x2[..., cont_idx] if cont_idx else None

        # Compute product only on available parts
        if x1_fp is not None and x1_cont is not None:
            return self.fp_kernel(x1_fp, x2_fp, diag=diag) * self.cont_kernel(x1_cont, x2_cont, diag=diag)
        elif x1_fp is not None:
            return self.fp_kernel(x1_fp, x2_fp, diag=diag)
        elif x1_cont is not None:
            return self.cont_kernel(x1_cont, x2_cont, diag=diag)
        else:
            raise ValueError("No valid feature indices found in feat_idx.")


class GPMix(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood, feat_idx):
        super(GPMix, self).__init__(train_x, train_y, likelihood, feat_idx)
        self.mean_module = gpytorch.means.ZeroMean()
        self.covar_module = MixingKernel(feat_idx=feat_idx)

        def forward(self, x):
            mean_x = self.mean_module(x)
            covar_x = self.covar_module(x)
            return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)
        



class GPMixRegressor(BaseEstimator):
    









class SimpleSincKernel(gpytorch.kernels.Kernel):
    has_lengthscale = True

    # this is the kernel function
    def forward(self, x1, x2, **params):
        # apply lengthscale
        x1_ = x1.div(self.lengthscale)
        x2_ = x2.div(self.lengthscale)
        # calculate the distance between inputs
        diff = self.covar_dist(x1_, x2_, **params)
        # prevent divide by 0 errors
        diff.where(diff == 0, torch.as_tensor(1e-20))
        # return sinc(diff) = sin(diff) / diff
        return torch.sin(diff).div(diff)



# gpytorch.kernels.RBFKernel()

class SimpleSincGPModel(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood):
        super().__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ConstantMean()
        self.covar_module = SimpleSincKernel()

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)

# initialize the new model
model = SimpleSincGPModel(train_x, train_y, likelihood)

# set to training mode and train
model.train()
likelihood.train()
train(model, likelihood)

# Get into evaluation (predictive posterior) mode and predict
model.eval()
likelihood.eval()
observed_pred = predict(model, likelihood)




