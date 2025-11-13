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
from gpytorch.priors import NormalPrior, GammaPrior
from gpytorch.kernels import Kernel, MaternKernel, ScaleKernel
from pyro.infer.mcmc import NUTS, MCMC
import pyro
from pyro.distributions import inverse_gamma as InvGammaPrior
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
    


def _weighted_tanimoto_distance(x1, x2):
    # (..., N, 1, D)
    x1_expanded = x1.unsqueeze(-2) 
    # (..., 1, M, D)
    x2_expanded = x2.unsqueeze(-3) 

    # Intersection: min(x1, x2)
    numerator = torch.min(x1_expanded, x2_expanded).sum(dim=-1)
    # Union: max(x1, x2)
    denominator = torch.max(x1_expanded, x2_expanded).sum(dim=-1)

    # If denominator is 0 (both vectors are zero), set S_T to 0, matching the R code.
    tanimoto_sim = torch.where(
        denominator == 0, 
        torch.zeros_like(denominator, dtype=x1.dtype),  # S_T = 0
        numerator / denominator
    )
    
    tanimoto_dist = tanimoto_sim
    
    return torch.clamp(tanimoto_dist, min=0.0)



class TanimotoRBF(gpytorch.kernels.Kernel):

    is_stationary = False  # MUST be non-stationary

    def __init__(self, **kwargs):
        # makes lengthscale a learnable parameter
        super().__init__(has_lengthscale=True, **kwargs)

    def forward(self, x1, x2, **params):
        # compute Tanimoto distance (NOT scaled)
        dist = _weighted_tanimoto_distance(x1, x2)

        # RBF-style scaling with learnable lengthscale:
        # k = exp( -0.5 * d / ell^2 )
        scaled_dist = dist / (self.lengthscale ** 2)

        return torch.exp(-0.5 * scaled_dist)




class MixingKernel(Kernel):
    """
    Product of TanimotoRBF (for fingerprint features)
    and Matern(3/2) (for continuous features).
    Priors and hyperparameters can be passed from outside.
    """

    def __init__(
            self,
            feat_idx: dict,
            l_prior=None,
            sigma_prior=None,
            **kwargs,
        ):

        super().__init__(**kwargs)

        # Normalize None values → []
        fp_idx = feat_idx.get("fp") or []
        cont_idx = feat_idx.get("count") or []

        self.feat_idx = {"fp": fp_idx, "count": cont_idx}

        # --- Fingerprint kernel (only if fp features exist) ---
        if len(fp_idx) > 0:
            self.fp_kernel = ScaleKernel(
                TanimotoRBF(
                    ard_num_dims=len(fp_idx),
                    # lengthscale_constraint=Positive()
                ),
                outputscale_prior=sigma_prior,
            )
            if l_prior is not None:
                self.fp_kernel.base_kernel.register_prior(
                    "lengthscale_prior", l_prior, "lengthscale"
                )
        else:
            self.fp_kernel = None

        # --- Continuous kernel (only if count features exist) ---
        if len(cont_idx) > 0:
            self.cont_kernel = ScaleKernel(
                MaternKernel(
                    nu=1.5,
                    ard_num_dims=len(cont_idx),
                ),
                outputscale_prior=sigma_prior,
            )
            if l_prior is not None:
                self.cont_kernel.base_kernel.register_prior(
                    "lengthscale_prior", l_prior, "lengthscale"
                )
        else:
            self.cont_kernel = None

    def forward(self, x1, x2, diag=False, **params):

        fp_idx = self.feat_idx["fp"]
        cont_idx = self.feat_idx["count"]

        x1_fp = x1[..., fp_idx] if fp_idx else None
        x2_fp = x2[..., fp_idx] if fp_idx else None

        x1_cont = x1[..., cont_idx] if cont_idx else None
        x2_cont = x2[..., cont_idx] if cont_idx else None

        # Both exist → product kernel
        if self.fp_kernel is not None and self.cont_kernel is not None:
            return self.fp_kernel(x1_fp, x2_fp, diag=diag) * self.cont_kernel(x1_cont, x2_cont, diag=diag)

        # Only FP
        if self.fp_kernel is not None:
            return self.fp_kernel(x1_fp, x2_fp, diag=diag)

        # Only count
        if self.cont_kernel is not None:
            return self.cont_kernel(x1_cont, x2_cont, diag=diag)

        raise ValueError("Both fp and count feature lists are empty.")
        

class GPMix(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood, feat_idx, l_prior=None, sigma_prior=None):
        super(GPMix, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ZeroMean()
        
        # Pass priors to the kernel at initialization
        self.covar_module = MixingKernel(
            feat_idx=feat_idx, 
            l_prior=l_prior, 
            sigma_prior=sigma_prior
        )
        # self.covar_module = gpytorch.kernels.ScaleKernel(gpytorch.kernels.RBFKernel(),
        #                                                  outputscale_prior=sigma_prior)
        # self.covar_module.base_kernel.register_prior(
        #             "lengthscale_prior", l_prior, "lengthscale"
        #         )
    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)
        


class GPMixRegressor(BaseEstimator):
    def __init__(
        self,
        feat_idx=None,
        learning_rate=0.1,
        num_epochs=100,
        use_cuda=False,
    ):
        self.feat_idx = feat_idx
        self.learning_rate = learning_rate
        self.num_epochs = num_epochs
        self.use_cuda = use_cuda and torch.cuda.is_available()

    def fit(self, X_train, Y_train):
        if isinstance(X_train, pd.DataFrame):
            X_train = X_train.to_numpy()
        if isinstance(Y_train, pd.DataFrame):
            Y_train = Y_train.to_numpy()

        X_train = torch.tensor(X_train, dtype=torch.float)
        Y_train = torch.tensor(Y_train, dtype=torch.float)

        if self.use_cuda:
            X_train = X_train.cuda()
            Y_train = Y_train.cuda()

        l_prior = GammaPrior(5.0, 5.0)
        sigma_prior = NormalPrior(0.0, 1.0)
        sigma_noise_prior = NormalPrior(0.0, 1.0)

        self.likelihood = gpytorch.likelihoods.GaussianLikelihood(
            noise_prior=sigma_noise_prior
        )
        self.model = GPMix(
            X_train, 
            Y_train, 
            self.likelihood, 
            self.feat_idx,
            l_prior=l_prior,
            sigma_prior=sigma_prior
        )
        
        if self.use_cuda:
            self.model = self.model.cuda()
            self.likelihood = self.likelihood.cuda()

        self.model.train()
        self.likelihood.train()

        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=self.learning_rate)

        self.mll =  gpytorch.mlls.ExactMarginalLogLikelihood(self.likelihood, self.model)

        # print(f"Running MAP optimization for {self.num_epochs} epochs...")
        
        for i in range(self.num_epochs):
            self.optimizer.zero_grad()
            output = self.model(X_train)
            
            loss = -self.mll(output, Y_train)
            
            loss.backward()
            
            # if (i + 1) % 10 == 0:
            #     print(f"Epoch {i+1}/{self.num_epochs} - Loss: {loss.item():.3f} "
            #           f"Noise: {self.model.likelihood.noise.item():.3f}")
            
            self.optimizer.step()

        # print("Optimization complete. Model is trained (fitted).")
        
        self.train_x = X_train
        self.train_y = Y_train
        return self

    # @torch.no_grad()
    def predict(self, X_test):
        if isinstance(X_test, pd.DataFrame):
            X_test = X_test.to_numpy()
        X_test = torch.tensor(X_test, dtype=torch.float)

        if self.use_cuda:
            X_test = X_test.cuda()

        self.model.eval()
        self.likelihood.eval()

        with torch.no_grad():
            pred_distribution = self.likelihood(self.model(X_test))

        mean_pred = pred_distribution.mean.detach().cpu().numpy()
        
        return mean_pred













# class GPMixRegressor(BaseEstimator):
#     def __init__(
#         self,
#         feat_idx=None,
#         num_samples=100,
#         warmup_steps=100,
#         use_cuda=False,
#     ):
#         self.feat_idx = feat_idx
#         self.num_samples = num_samples
#         self.warmup_steps = warmup_steps
#         self.use_cuda = use_cuda and torch.cuda.is_available()

#     # --- NEW: Helper method for X data ---

#     def fit(self, X_train, Y_train):
#         # --- Convert input to torch tensors safely ---
#         if isinstance(X_train, pd.DataFrame):
#             X_train = X_train.to_numpy()

#         if isinstance(Y_train, pd.DataFrame):
#             Y_train = Y_train.to_numpy()

#         X_train = torch.tensor(X_train, dtype=torch.float)
#         Y_train = torch.tensor(Y_train.ravel(), dtype=torch.float)
#         # --- Priors (as you specified) ---
#         l_prior = GammaPrior(5.0, 5.0)
#         sigma_prior = NormalPrior(0.0, 1.0)
#         sigma_noise_prior = NormalPrior(0.0, 1.0)

#         # --- Likelihood and model ---
#         self.likelihood = gpytorch.likelihoods.GaussianLikelihood()
        
#         # Pass priors to the model (from previous refactor)
#         self.model = GPMix(
#             X_train, 
#             Y_train, 
#             self.likelihood, 
#             self.feat_idx,
#             l_prior=l_prior,
#             sigma_prior=sigma_prior
#         )
#         self.likelihood.register_prior("noise_prior", sigma_noise_prior, "noise")

#         def pyro_model(x, y):
#             with gpytorch.settings.fast_computations(False, False, False):
#                 sampled_model = self.model.pyro_sample_from_prior()
#                 output = sampled_model.likelihood(sampled_model(x))
#                 pyro.sample("obs", output, obs=y)
#             return y

#         nuts_kernel = NUTS(pyro_model)
#         mcmc = MCMC(
#             nuts_kernel,
#             num_samples=self.num_samples,
#             warmup_steps=self.warmup_steps,
#             disable_progbar=True
#         )

#         print(f"Running NUTS with {self.num_samples} samples, {self.warmup_steps} warmup...")
#         mcmc.run(X_train, Y_train)
#         self.mcmc_samples = mcmc.get_samples()
#         print("Sampling complete. Posterior hyperparameters stored.")

#         self.train_x = X_train
#         self.train_y = Y_train
#         return self

#     @torch.no_grad()
#     def predict(self, X_test):
#         # --- Convert test data to tensor ---
#         if isinstance(X_test, pd.DataFrame):
#             X_test = X_test.to_numpy()
#         X_test = torch.tensor(X_test, dtype=torch.float)
#         if self.use_cuda:
#             self.model = self.model.cuda()
#             self.likelihood = self.likelihood.cuda()

#         self.model.eval()
#         self.likelihood.eval()

#         preds = []
#         for i in range(self.num_samples):
#             sample_state = {k: v[i] for k, v in self.mcmc_samples.items()}
#             self.model.load_state_dict(sample_state, strict=False)
#             pred = self.likelihood(self.model(X_test)).mean.detach().cpu().numpy()
#             preds.append(pred)

#         mean_pred = np.mean(preds, axis=0)
#         return mean_pred






# if __name__ == "__main__":
    
#     # --- 1. Define Data Shape ---
#     N_train = 20
#     N_test = 5
#     D_fp = 10     # 10 fingerprint features (integer counts)
#     D_cont = 100  # 100 continuous features
    
#     # --- 2. Define Feature Indices ---
#     # This is critical for the kernel
#     feat_idx = {
#         'fp': list(range(D_fp)),                   # Cols 0-9
#         'count': list(range(D_fp, D_fp + D_cont))  # Cols 10-109
#     }
    
#     # --- 3. Generate Synthetic Data ---
#     # Fingerprint data (integer counts: 1 to 14) <-- MODIFIED HERE
#     X_fp_train = np.random.randint(1, 15, size=(N_train, D_fp))
#     X_fp_test = np.random.randint(1, 15, size=(N_test, D_fp))

#     # Continuous data (float)
#     X_cont_train = np.random.randn(N_train, D_cont)
#     X_cont_test = np.random.randn(N_test, D_cont)

#     # Combine data
#     X_train_np = np.hstack([X_fp_train, X_cont_train])
#     X_test_np = np.hstack([X_fp_test, X_cont_test])
    
#     # Create a synthetic target (Y)
#     # Y = 2 * (sum of 1st cont feature) - 0.5 * (sum of fp features) + noise
#     Y_train_np = (
#         2.0 * X_cont_train[:, 0] + 
#         -0.5 * X_fp_train.sum(axis=1) + 
#         np.random.randn(N_train) * 0.1
#     )

#     # --- 4. Convert to Pandas (to test data handlers) ---
#     columns = [f"fp_{i}" for i in range(D_fp)] + [f"cont_{i}" for i in range(D_cont)]
    
#     X_train_df = pd.DataFrame(X_train_np, columns=columns)
#     Y_train_series = pd.Series(Y_train_np, name="target")
#     X_test_df = pd.DataFrame(X_test_np, columns=columns)
    
#     print("--- Generated Training Data (Head) ---")
#     print(X_train_df.head())
    
#     # --- 5. Initialize and Run Model ---
#     print("\n--- Initializing Model ---")
#     gp_model = GPMixRegressor(
#         feat_idx=feat_idx,
#         num_samples=50,
#         warmup_steps=50,
#         use_cuda=False
#     )

#     print("\n--- Fitting Model ---")
#     # Assuming GPMixRegressor, GPMix, and MixingKernel classes are defined above this block
#     gp_model.fit(X_train_df, Y_train_series)
    
#     # --- 6. Make Predictions ---
#     print("\n--- Predicting ---")
#     predictions = gp_model.predict(X_test_df)
    
#     print("\n--- Test Predictions ---")
#     print(predictions)
#     print(f"\nPredicted {len(predictions)} values.")
