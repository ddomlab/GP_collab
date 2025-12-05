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

    
class TanimotoRBF(gpytorch.kernels.Kernel):
    is_stationary = False
    has_lengthscale = True

    def __init__(self, eps=1e-6, **kwargs):
        super().__init__(**kwargs)
        self.eps = eps

    def forward(self, x1, x2, diag=False, **params):
        if diag:
            zero_dist = torch.zeros(x1.shape[:-1], device=x1.device, dtype=x1.dtype)
            scaled_dist = zero_dist.unsqueeze(-1).div(self.lengthscale).squeeze(-1)
            return torch.exp(-scaled_dist)

        dist = weighted_tanimoto_distance(x1, x2, eps=self.eps)
        exp_component = dist.div(self.lengthscale).pow(2).mul(-0.5)
        
        return exp_component.exp()



class GPMixGPyTorch(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood, feat_idx):
        super().__init__(train_x, train_y, likelihood)
        
        self.mean_module = gpytorch.means.ConstantMean()
        
        fp_idx = feat_idx.get("fp") or []
        cont_idx = feat_idx.get("count") or []
        
        kernels = []
        safe_constraint = gpytorch.constraints.GreaterThan(1e-6)
        if len(fp_idx) > 0:
            k_fp = TanimotoRBF(
                active_dims=torch.tensor(fp_idx, dtype=torch.long),
                batch_shape=torch.Size([]),
                # lengthscale_constraint=safe_constraint,
            )
            k_fp.register_prior(
                "lengthscale_prior", 
                gpytorch.priors.GammaPrior(concentration=2.0, rate=0.5), 
                "lengthscale"
            )
            kernels.append(k_fp)
            
        if len(cont_idx) > 0:
            k_cont = gpytorch.kernels.MaternKernel(
                nu=1.5,
                ard_num_dims=len(cont_idx),
                active_dims=torch.tensor(cont_idx, dtype=torch.long),
                # lengthscale_constraint=safe_constraint,
            )
            k_cont.register_prior(
                "lengthscale_prior", 
                gpytorch.priors.GammaPrior(concentration=2.0, rate=0.5), 
                "lengthscale"
            )
            kernels.append(k_cont)

        if not kernels:
            raise ValueError("No features provided")
            
        combined_kernel = gpytorch.kernels.ProductKernel(*kernels)
            
        self.covar_module = gpytorch.kernels.ScaleKernel(combined_kernel)
        
        self.covar_module.register_prior(
            "outputscale_prior",
            gpytorch.priors.LogNormalPrior(0.0, 1.0),
            "outputscale"
        )
        
        self.mean_module.register_prior(
            "mean_prior",
            gpytorch.priors.NormalPrior(0.0, 1.0),
            "constant"
        )

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)




class GPytorchMixMCMCRegressor(BaseEstimator, RegressorMixin):
    def __init__(self, feat_idx=None, 
                num_samples=1000,
                warmup_steps=1000,
                num_chains=1,
                use_cuda=False
        ):
        
        self.feat_idx = feat_idx
        self.num_samples = num_samples
        self.warmup_steps = warmup_steps
        self.num_chains = num_chains
        self.use_cuda = use_cuda
        self.model = None
        self.likelihood = None
        self.mcmc_samples = None

    def fit(self, X, y):
        if self.feat_idx is None:
             raise ValueError("feat_idx must be provided")

        if isinstance(X, pd.DataFrame):
            X = X.values
        if isinstance(y, (pd.DataFrame, pd.Series)):
            y = y.values

        X = torch.as_tensor(X, dtype=torch.float32)
        y = torch.as_tensor(y, dtype=torch.float32).view(-1)

        if self.use_cuda and torch.cuda.is_available():
            X = X.cuda()
            y = y.cuda()

        self.train_x = X
        self.train_y = y

        self.likelihood = gpytorch.likelihoods.GaussianLikelihood(
            # noise_constraint=gpytorch.constraints.GreaterThan(1e-6)
        )
        self.likelihood.register_prior(
            "noise_prior", 
            gpytorch.priors.GammaPrior(concentration=2.0, rate=1.0), 
            "noise"
        )

        self.model = GPMixGPyTorch(X, y, self.likelihood, self.feat_idx)

        if self.use_cuda and torch.cuda.is_available():
            self.model = self.model.cuda()
            self.likelihood = self.likelihood.cuda()

        self.model.train()
        self.likelihood.train()

        mll = gpytorch.mlls.ExactMarginalLogLikelihood(self.likelihood, self.model)

        def pyro_model(x, y):
            self.model.pyro_sample_from_prior()
            output = self.model(x)
            loss = mll(output, y)
            pyro.factor("gp_mll", loss)

        nuts_kernel = NUTS(pyro_model)
        mcmc = MCMC(
            nuts_kernel, 
            num_samples=self.num_samples, 
            warmup_steps=self.warmup_steps, 
            num_chains=self.num_chains
        )

        mcmc.run(X, y)
        self.mcmc_samples = mcmc.get_samples()
        return self

    def predict(self, X_test, return_std=False):
        if isinstance(X_test, pd.DataFrame):
            X_test = X_test.values

        X_test = torch.as_tensor(X_test, dtype=torch.float32)
        if self.use_cuda and torch.cuda.is_available():
            X_test = X_test.cuda()

        self.model.eval()
        self.likelihood.eval()

        # 1. Load ALL samples at once. 
        # This converts the model into a "Batch Model" (e.g. 300 internal models)
        self.model.pyro_load_from_samples(self.mcmc_samples)

        # 2. Get the batch size (Num_Samples)
        # We check the first sample to find the number of MCMC draws
        num_samples = list(self.mcmc_samples.values())[0].shape[0]

        # 3. Expand the Test Data to match the batch size
        # Input shape: (N_test, D) -> Expanded: (Num_Samples, N_test, D)
        expanded_test_x = X_test.unsqueeze(0).expand(num_samples, *X_test.shape)

        with torch.no_grad():
            # 4. Forward pass
            # The model is now a batch model, so it returns a batch of MultivariateNormals
            output = self.model(expanded_test_x)
            observed_output = self.likelihood(output) # Adds noise variance

            # 5. Extract statistics
            # means shape: (Num_Samples, N_test)
            means = observed_output.mean
            # variances shape: (Num_Samples, N_test)
            variances = observed_output.variance

            # 6. Aggregate
            mean_prediction = means.mean(0) # Average over samples
            
            # Total Variance = Mean(Variances) + Variance(Means)
            total_variance = variances.mean(0) + means.var(0)

        if return_std:
            return mean_prediction.cpu().numpy(), total_variance.sqrt().cpu().numpy()
        
        return mean_prediction.cpu().numpy()
    

