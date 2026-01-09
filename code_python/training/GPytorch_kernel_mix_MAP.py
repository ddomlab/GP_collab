import gpytorch
import numpy as np
import pandas as pd
import torch
# import torch.nn as nn
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.model_selection import train_test_split
from torch.utils import benchmark
from gpytorch.constraints import Positive
from gpytorch.kernels import Kernel
from gpytorch.kernels.rbf_kernel import postprocess_rbf
from gpytorch.priors import NormalPrior, GammaPrior
from gpytorch.kernels import Kernel, MaternKernel, ScaleKernel
from pyro.infer.mcmc import NUTS, MCMC
import pyro
from gpytorch.models import PyroGP # <-- Import PyroGP
from pyro.distributions import inverse_gamma as InvGammaPrior
# from torch_geometric.data import Batch
# from torch_geometric.loader import DataLoader
import gc
# from pytorch_mpnn import DMPNNPredictor, RevIndexedData, smiles2data



    


def weighted_tanimoto_distance(x1, x2, eps=1e-6):
    x1e = x1.unsqueeze(-2)
    x2e = x2.unsqueeze(-3)

    numerator = torch.min(x1e, x2e).sum(dim=-1)
    denominator = torch.max(x1e, x2e).sum(dim=-1)

    sim = (numerator + eps) / (denominator + eps)
    dist = 1.0 - sim
    return torch.clamp(dist, min=0.)

# class TanimotoRBF(gpytorch.kernels.Kernel):
#     has_lengthscale = True

#     def __init__(self, eps=1e-6, **kwargs):
#         super().__init__(has_lengthscale=True, **kwargs)
#         self.eps = eps

#     def forward(self, x1, x2, diag=False, **params):
#         batch_shape = self.batch_shape  # e.g. torch.Size([1, 1]) in your error

#         if diag:
#             # We need shape batch_shape + (n,)
#             n = x1.size(-2) if x1.dim() > 1 else x1.size(0)
#             diag_covar = x1.new_ones(*batch_shape, n)
#             return diag_covar

#         # ---------------------------------------------------------
#         # 1. Compute pairwise Tanimoto distance without batch dims
#         # ---------------------------------------------------------
#         # x1: n x d, x2: m x d (or possibly have extra data batch dims)
#         # You can use your helper here if you like
#         x1e = x1.unsqueeze(-2)  # ... x n x 1 x d
#         x2e = x2.unsqueeze(-3)  # ... x 1 x m x d

#         numerator = torch.min(x1e, x2e).sum(dim=-1)
#         denominator = torch.max(x1e, x2e).sum(dim=-1)

#         sim = (numerator + self.eps) / (denominator + self.eps)
#         dist = torch.clamp(1.0 - sim, min=0.0)  # shape: ... x n x m

#         # ---------------------------------------------------------
#         # 2. Add *parameter* batch dims so shape matches batch_shape
#         # ---------------------------------------------------------
#         # Right now dist has shape:   (...data_batch..., n, m)
#         # We need:                    (*batch_shape, n, m)
#         # so we add leading singleton dims until there are enough.
#         while dist.dim() < len(batch_shape) + 2:
#             dist = dist.unsqueeze(0)

#         # Now dist should have shape batch_shape + (n, m),
#         # or something that broadcasts to it.

#         # ---------------------------------------------------------
#         # 3. Apply lengthscale with correct broadcasting
#         # ---------------------------------------------------------
#         ls = self.lengthscale  # shape: batch_shape + (1, 1)

#         covar = torch.exp(-0.5 * (dist / ls) ** 2)
#         return covar

# class CountDistance(torch.nn.Module):
#     r"""
#     Distance module for bit vector test_kernels.
#     """

#     def __init__(self, postprocess_script=lambda x: x):
#         super().__init__()
#         self._postprocess = postprocess_script

#     def _sim(self, x1, x2, postprocess, x1_eq_x2=False, metric="tanimoto"):
#         r"""
#         Computes the similarity between x1 and x2
#         Args:
#             :attr: `x1`: (Tensor `n x d` or `b x n x d`):
#                 First set of data where b is a batch dimension
#             :attr: `x2`: (Tensor `m x d` or `b x m x d`):
#                 Second set of data where b is a batch dimension
#             :attr: `postprocess` (bool):
#                 Whether to apply a postprocess script (default is none)
#             :attr: `x1_eq_x2` (bool):
#                 Is x1 equal to x2
#             :attr: `metric` (str):
#                 String specifying the similarity metric. One of ['tanimoto']
#         Returns:
#             (:class:`Tensor`, :class:`Tensor) corresponding to the similarity matrix between `x1` and `x2`
#         """

#         # Branch for Tanimoto metric
#         if metric == "tanimoto":
#             res = weighted_tanimoto_distance(x1, x2)
#             res.clamp_min_(0)  # zero out negative values
#             return self._postprocess(res) if postprocess else res
#         else:
#             raise RuntimeError(
#                 "Similarity metric not supported. Available options are 'tanimoto'"
#             )




# class TanimotoRBF(gpytorch.kernels.Kernel):
#     has_lengthscale = True
    
#     def __init__(self, eps=1e-6,metric="tanimoto", **kwargs):
#         super().__init__(**kwargs)
#         self.eps = eps
#         self.metric = metric

#     def covar_dist(
#             self,
#             x1,
#             x2,
#             last_dim_is_batch=False,
#             dist_postprocess_func=lambda x: x,
#             postprocess=True,
#             **params,
#     ):
#         r"""
#         This is a helper method for computing the bit vector similarity between
#         all pairs of points in x1 and x2.
#         Args:
#             :attr:`x1` (Tensor `n x d` or `b1 x ... x bk x n x d`):
#                 First set of data.
#             :attr:`x2` (Tensor `m x d` or `b1 x ... x bk x m x d`):
#                 Second set of data.
#             :attr:`last_dim_is_batch` (tuple, optional):
#                 Is the last dimension of the data a batch dimension or not?
#         Returns:
#             (:class:`Tensor`, :class:`Tensor) corresponding to the distance matrix between `x1` and `x2`.
#             The shape depends on the kernel's mode
#             * `diag=False`
#             * `diag=False` and `last_dim_is_batch=True`: (`b x d x n x n`)
#             * `diag=True`
#             * `diag=True` and `last_dim_is_batch=True`: (`b x d x n`)
#         """
#         if last_dim_is_batch:
#             x1 = x1.transpose(-1, -2).unsqueeze(-1)
#             x2 = x2.transpose(-1, -2).unsqueeze(-1)

#         x1_eq_x2 = torch.equal(x1, x2)

#         # torch scripts expect tensors
#         postprocess = torch.tensor(postprocess)

#         res = None

#         # Cache the Distance object or else JIT will recompile every time
#         if (
#                 not self.distance_module
#                 or self.distance_module._postprocess != dist_postprocess_func
#         ):
#             self.distance_module = CountDistance(dist_postprocess_func)

#         res = self.distance_module._sim(
#             x1, x2, postprocess, x1_eq_x2, self.metric
#         )
#         return res

    
#     def forward(self, x1, x2, diag=False, **params):
#         if diag:
#             assert x1.size() == x2.size() and torch.equal(x1, x2)
#             return torch.ones(
#                 *x1.shape[:-2], x1.shape[-2], dtype=x1.dtype, device=x1.device
#             )

#         else:
#             dist = self.covar_dist(x1, x2, **params)
#             exp_component = torch.exp(-0.5 * (dist / self.lengthscale)**2)
#             return exp_component


class TanimotoRBF(Kernel):
    has_lengthscale = True

    def _tanimoto_dist(self, x1, x2, eps=1e-6):
        # x1: (..., N, D), x2: (..., M, D)
        x1e = x1.unsqueeze(-2)
        x2e = x2.unsqueeze(-3)

        numerator = torch.min(x1e, x2e).sum(dim=-1)
        denominator = torch.max(x1e, x2e).sum(dim=-1)

        sim = (numerator + eps) / (denominator + eps)
        return (1.0 - sim).clamp(min=0.)

    def forward(self, x1, x2, diag=False, **params):
        if diag:
            # If x1 and x2 are the same, the distance is 0 and exp(0) = 1.
            # If they are different (e.g. in cross-validation), we compute 
            # only the element-wise distance.
            if torch.equal(x1, x2):
                # Returns a vector of 1s with the shape of the batch/n_samples
                return torch.ones(*x1.shape[:-1], device=x1.device, dtype=x1.dtype)
            
            # Element-wise Tanimoto distance
            num = torch.min(x1, x2).sum(dim=-1)
            den = torch.max(x1, x2).sum(dim=-1)
            dist = 1.0 - (num + 1e-6) / (den + 1e-6)
            return torch.exp(-0.5 * (dist.clamp(min=0.) / self.lengthscale).pow(2)).squeeze(-1)

        # Full covariance matrix path
        dist = self._tanimoto_dist(x1, x2)
        return torch.exp(-0.5 * (dist / self.lengthscale).pow(2))




# class MixingKernel(gpytorch.kernels.Kernel):
#     def __init__(self, feat_idx):
#         super().__init__()
#         self.k_fp = TanimotoRBF(
#                         active_dims=torch.tensor(feat_idx.get("fp", []), dtype=torch.long)
#                         )
        
#         self.k_fp.register_prior(
#                                 "fp_lengthscale_prior",
#                                 gpytorch.priors.LogNormalPrior(0, 1),
#                                 "lengthscale"
#                                 )
        
#         self.k_cont = gpytorch.kernels.MaternKernel(
#             nu=1.5,
#             ard_num_dims=len(feat_idx.get("count")),
#             active_dims=torch.tensor(feat_idx.get("count", []), dtype=torch.long)
#             )
#         self.k_cont.register_prior(
#                                 "cont_lengthscale_prior",
#                                 gpytorch.priors.LogNormalPrior(0, 1),
#                                 "lengthscale"
#                                 )
#         self.base_kernel = self.k_fp * self.k_cont
        
#     def forward(self, x1, x2, diag=False, **params):
#         return self.base_kernel(x1, x2, diag=diag, **params)

#     @property
#     def is_stationary(self):
#         return self.base_kernel.is_stationary



class GPMixGPyTorch(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood, feat_idx):
        super().__init__(train_x, train_y, likelihood)
    
        self.mean_module = gpytorch.means.ZeroMean()

        self.k_fp = TanimotoRBF(
                        active_dims=torch.tensor(feat_idx.get("fp", []), dtype=torch.long)
                        )
        self.k_fp.register_prior(
                                "fp_lengthscale_prior",
                                gpytorch.priors.GammaPrior(5,5),
                                "lengthscale"
                                )
        self.k_cont = gpytorch.kernels.MaternKernel(
            nu=1.5,
            ard_num_dims=len(feat_idx.get("count")),
            active_dims=torch.tensor(feat_idx.get("count", []), dtype=torch.long)
            )
        self.k_cont.register_prior(
                                "cont_lengthscale_prior",
                                gpytorch.priors.GammaPrior(5,5),
                                "lengthscale"
                                )
        self.base_kernel = self.k_fp * self.k_cont

        self.covar_module = gpytorch.kernels.ScaleKernel(
            self.base_kernel
        )

        self.covar_module.register_prior(
            "variance_prior",
            gpytorch.priors.LogNormalPrior(0.0, 1.0),
            "outputscale",
        )
    def forward(self, x):
            mean_x = self.mean_module(x)
            covar_x = self.covar_module(x)
            return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)

    

class GPytorchMixMCMCRegressor(BaseEstimator, RegressorMixin):
    def __init__(self, 
                feat_group=None, 
                num_samples=200,
                warmup_steps=200,
                num_drawn_samples=100,
                num_chains=1,
                use_cuda=False,
                random_state=42
        ):
        
        self.feat_group = feat_group
        self.num_samples = num_samples
        self.warmup_steps = warmup_steps
        self.num_drawn_samples = num_drawn_samples
        self.num_chains = num_chains
        self.use_cuda = use_cuda
        self.model = None
        self.likelihood = None
        self.mcmc_samples = None
        self.random_state = random_state


    def fit(self, X, y):
        if isinstance(X, pd.DataFrame):
            self.feat_idx = {
                        'fp': [X.columns.get_loc(c) for c in self.feat_group.get("fp")] if self.feat_group.get("fp") else None,
                        'count': [X.columns.get_loc(c) for c in self.feat_group.get("count") ] if self.feat_group.get("count") else None
                    }
            self.count_feat_name_idx = {c: X.columns.get_loc(c) for c in self.feat_group.get("count")} if self.feat_group.get("count") else {}
            X = X.to_numpy()
        if self.feat_idx is None:
             raise ValueError("feat_idx must be provided")

        if isinstance(y, (pd.DataFrame, pd.Series)):
            y = y.values

        X = torch.as_tensor(X, dtype=torch.float32)
        y = torch.as_tensor(y, dtype=torch.float32).view(-1)

        if self.use_cuda and torch.cuda.is_available():
            X = X.cuda()
            y = y.cuda()

        self.train_x = X
        self.train_y = y

        self.likelihood = gpytorch.likelihoods.GaussianLikelihood(noise_constraint=gpytorch.constraints.Positive())
        self.model = GPMixGPyTorch(X, y, self.likelihood, self.feat_idx)
        self.likelihood.register_prior(
            "noise_prior", 
            gpytorch.priors.LogNormalPrior(0.0, 1.0), 
            "noise"
        )
        
        if self.use_cuda and torch.cuda.is_available():
            self.model = self.model.cuda()
            self.likelihood = self.likelihood.cuda()


        # mll = gpytorch.mlls.ExactMarginalLogLikelihood(self.likelihood, self.model)

        def pyro_model(x, y):
            with gpytorch.settings.fast_computations(False, False, False):
                with gpytorch.settings.cholesky_jitter(1e-2):
                    sampled_model = self.model.pyro_sample_from_prior()
                    output = sampled_model.likelihood(sampled_model(x))
                    pyro.sample("obs", output, obs=y)
            return y
        pyro.clear_param_store()
        if self.random_state is not None:
            pyro.set_rng_seed(self.random_state)
        def hook(kernel, samples, stage, i):
            if stage == "Warmup" and i % 50 == 0:
                print(i, samples["likelihood.noise"], samples["covar_module.variance"])
        nuts_kernel = NUTS(pyro_model, ignore_jit_warnings=True, jit_compile=False)
        mcmc = MCMC(
            nuts_kernel, 
            num_samples=self.num_samples, 
            warmup_steps=self.warmup_steps, 
            num_chains=self.num_chains,
            # disable_progbar=True,
            hook_fn=hook,
        )

        mcmc.run(X, y)
        self.mcmc_samples = mcmc.get_samples(num_samples=self.num_drawn_samples)
        return self

    def predict(self, X_test, return_std=False):
        if isinstance(X_test, pd.DataFrame):
            X_test = X_test.to_numpy()

        X_test = torch.as_tensor(X_test, dtype=torch.float32)
        if self.use_cuda and torch.cuda.is_available():
            X_test = X_test.cuda()

        self.model.eval()
        # self.likelihood.eval()

        self.model.pyro_load_from_samples(self.mcmc_samples)
        self.model.eval()
        num_samples = list(self.mcmc_samples.values())[0].shape[0]
        expanded_test_x = X_test.unsqueeze(0).expand(num_samples, *X_test.shape)

        with torch.no_grad():
            output = self.model(expanded_test_x)
            observed_output = self.likelihood(output) 

            means = observed_output.mean
            mean_prediction = means.mean(0) 
        return mean_prediction.cpu().numpy()
    
    
