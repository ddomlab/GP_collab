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

    
class TanimotoRBF(Kernel):
    def __init__(self, eps: float = 1e-6,**kwargs):
        super().__init__(**kwargs)
        self.eps = eps


    has_lengthscale = True
    def forward(self, x1: torch.Tensor, x2: torch.Tensor, diag: bool = False, **params) -> torch.Tensor:
        
        if diag:
            # The distance from a point to itself is 0.
            # exp(-0.5 * (0 / l)^2) = exp(0) = 1.
            return x1.new_ones(x1.shape[:-1])

        # Compute Tanimoto distance D_T(x1, x2)
        tanimoto_dist = weighted_tanimoto_distance(x1, x2, eps=self.eps)

        # Apply the RBF Kernel formula using the single lengthscale
        # This is equivalent to: exp( -square(D_T) / (2 * square(l)) )
        covar = torch.exp(-0.5 * (tanimoto_dist / self.lengthscale)**2)
        return covar



class GPMixGPyTorch(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood, feat_idx):
        super().__init__(train_x, train_y, likelihood)
        
        self.mean_module = gpytorch.means.ConstantMean()
        
        fp_idx = feat_idx.get("fp") or []
        cont_idx = feat_idx.get("count") or []
        
        kernels = []
        
        if len(fp_idx) > 0:
            k_fp = TanimotoRBF(
                active_dims=torch.tensor(fp_idx, dtype=torch.long),
                batch_shape=torch.Size([]),
            )
            k_fp.register_prior(
                "lengthscale_prior", 
                gpytorch.priors.LogNormalPrior(0.0, 1.0), 
                "lengthscale"
            )
            kernels.append(k_fp)
            
        if len(cont_idx) > 0:
            k_cont = gpytorch.kernels.MaternKernel(
                nu=1.5,
                ard_num_dims=len(cont_idx),
                active_dims=torch.tensor(cont_idx, dtype=torch.long)
            )
            k_cont.register_prior(
                "lengthscale_prior", 
                gpytorch.priors.LogNormalPrior(0.0, 1.0), 
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

        self.likelihood = gpytorch.likelihoods.GaussianLikelihood()
        self.likelihood.register_prior(
            "noise_prior", 
            gpytorch.priors.LogNormalPrior(0.0, 1.0), 
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
    

# class MixingKernel(Kernel):
#     def __init__(
#             self,
#             feat_idx: dict,
#             l_prior=None,
#             **kwargs,
#         ):

#         super().__init__(**kwargs)

#         fp_idx = feat_idx.get("fp") or []
#         cont_idx = feat_idx.get("count") or []

#         self.feat_idx = {"fp": fp_idx, "count": cont_idx}

#         if len(fp_idx) > 0:
#             # self.fp_kernel = TanimotoMatern32()
#             self.fp_kernel = TanimotoRBF()
#             self.fp_kernel.initialize(lengthscale=1)
#             if l_prior is not None:
#                 self.fp_kernel.register_prior(
#                     "lengthscale_prior", l_prior, "lengthscale"
#                 )
#         else:
#             self.fp_kernel = None

#         if len(cont_idx) > 0:
#             self.cont_kernel = MaternKernel(
#                     nu=1.5,
#                     ard_num_dims=len(cont_idx),
#                 )
#             init_vec = torch.full((len(cont_idx),), float(2))
#             self.cont_kernel.initialize(lengthscale=init_vec)
#             if l_prior is not None:
#                 self.cont_kernel.register_prior(
#                     "lengthscale_prior", l_prior, "lengthscale",
                    
#                 )
#         else:
#             self.cont_kernel = None

#     def forward(self, x1, x2, diag=False, **params):
#         fp_idx = self.feat_idx["fp"]
#         cont_idx = self.feat_idx["count"]

#         x1_fp = x1[..., fp_idx] if fp_idx else None
#         x2_fp = x2[..., fp_idx] if fp_idx else None

#         x1_cont = x1[..., cont_idx] if cont_idx else None
#         x2_cont = x2[..., cont_idx] if cont_idx else None

#         k_fp = x1.new_ones(x1.shape[:-1]) if diag else x1.new_ones(x1.shape[:-2] + (x1.shape[-2], x2.shape[-2]))
#         k_cont = k_fp.clone() 

#         if self.fp_kernel is not None:
#             k_fp = self.fp_kernel(x1_fp, x2_fp, diag=diag)

#         if self.cont_kernel is not None:
#             k_cont = self.cont_kernel(x1_cont, x2_cont, diag=diag)
            
#         if self.fp_kernel is None and self.cont_kernel is None:
#              raise ValueError("Both fp and count feature lists are empty.")

#         return k_fp * k_cont
        

# class GPMix(gpytorch.models.ExactGP):
#     def __init__(self, train_x, train_y, likelihood, feat_idx, l_prior=None, sigma_prior=None):
#         super().__init__(train_x, train_y, likelihood)
#         self.mean_module = gpytorch.means.ZeroMean()
        
#         core_kernel = MixingKernel(
#             feat_idx=feat_idx, 
#             l_prior=l_prior
#         )
        
#         self.covar_module = gpytorch.kernels.ScaleKernel(
#             core_kernel,
#             outputscale_prior=sigma_prior
#         )

#     def forward(self, x):
#         mean_x = self.mean_module(x)
#         covar_x = self.covar_module(x)
#         return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)
        

# class GPMixRegressor(BaseEstimator):
#     def __init__(
#         self,
#         feat_idx=None,
#         learning_rate=0.01,
#         num_epochs=200,
#         use_cuda=False,
#     ):
#         self.feat_idx = feat_idx
#         self.learning_rate = learning_rate
#         self.num_epochs = num_epochs
#         self.use_cuda = use_cuda and torch.cuda.is_available()

#     def fit(self, X_train, Y_train):
#         if isinstance(X_train, pd.DataFrame):
#             X_train = X_train.to_numpy()
#         if isinstance(Y_train, pd.DataFrame):
#             Y_train = Y_train.to_numpy()

#         X_train = torch.tensor(X_train, dtype=torch.float)
#         Y_train = torch.tensor(Y_train, dtype=torch.float)

#         if self.use_cuda:
#             X_train = X_train.cuda()
#             Y_train = Y_train.cuda()

#         l_prior = GammaPrior(5.0, 5.0)
#         sigma_prior = NormalPrior(0.0, 1.0)
#         sigma_noise_prior = NormalPrior(0.0, 1.0)

#         self.likelihood = gpytorch.likelihoods.GaussianLikelihood(
#             noise_prior=sigma_noise_prior
#         )
#         self.model = GPMix(
#             X_train, 
#             Y_train, 
#             self.likelihood, 
#             self.feat_idx,
#             l_prior=l_prior,
#             sigma_prior=sigma_prior
#         )
        
#         if self.use_cuda:
#             self.model = self.model.cuda()
#             self.likelihood = self.likelihood.cuda()

#         self.model.train()
#         self.likelihood.train()

#         self.optimizer = torch.optim.Adam(self.model.parameters(), lr=self.learning_rate)

#         self.mll =  gpytorch.mlls.ExactMarginalLogLikelihood(self.likelihood, self.model)

#         # with gpytorch.settings.cholesky_jitter(1e-2):
#         for i in range(self.num_epochs):
#             # print(f"Epoch {i+1}/{self.num_epochs}")
#             self.optimizer.zero_grad()
#             output = self.model(X_train)
            
#             loss = -self.mll(output, Y_train)
            
#             loss.backward()
            
#             # if (i + 1) % 10 == 0:
#             #     print(f"Epoch {i+1}/{self.num_epochs} - Loss: {loss.item():.3f} "
#             #           f"Noise: {self.model.likelihood.noise.item():.3f}")
            
#             self.optimizer.step()

#         # print("Optimization complete. Model is trained (fitted).")
        
#         self.train_x = X_train
#         self.train_y = Y_train
#         return self

#     @torch.no_grad()
#     def predict(self, X_test):
#         if isinstance(X_test, pd.DataFrame):
#             X_test = X_test.to_numpy()
#         X_test = torch.tensor(X_test, dtype=torch.float)

#         if self.use_cuda:
#             X_test = X_test.cuda()

#         self.model.eval()
#         self.likelihood.eval()

#         # with torch.no_grad():
#         pred_distribution = self.likelihood(self.model(X_test))

#         mean_pred = pred_distribution.mean.detach().cpu().numpy()
        
#         return mean_pred


# # (Assume MixingKernel, TanimotoExponential, etc., are defined above)

# class GPMixMCMC(BaseEstimator): # <-- 1. MUST inherit from PyroGP, not ExactGP
#     def __init__(self, train_x, train_y, likelihood, feat_idx, l_prior=None, sigma_prior=None):
        
#         # 2. This super().__init__() call is CRITICAL.
#         # It initializes the PyroGP parent, which gives 
#         # your model the .pyro_observe() method.
#         super().__init__(train_x, train_y, likelihood) 
        
#         self.mean_module = gpytorch.means.ZeroMean()
        
#         core_kernel = MixingKernel(
#             feat_idx=feat_idx, 
#             l_prior=l_prior
#         )
        
#         self.covar_module = gpytorch.kernels.ScaleKernel(
#             core_kernel,
#             outputscale_prior=sigma_prior
#         )

#     def model(self, x, y):
#         """This method is required by PyroGP for MCMC."""
#         self.pyro_sample_from_prior() 
#         output = self.forward(x)
        
#         # This line will now work because super().__init__() was called
#         self.pyro_observe(output, y)

#     def forward(self, x):
#         mean_x = self.mean_module(x)
#         covar_x = self.covar_module(x)
#         return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)


# class GPMixMCMCRegressor(BaseEstimator):
#     def __init__(
#         self,
#         feat_idx=None,
#         num_samples=2000,
#         warmup_steps=1000,
#         use_cuda=False,
#     ):
#         self.feat_idx = feat_idx
#         self.num_samples = num_samples
#         self.warmup_steps = warmup_steps
#         self.use_cuda = use_cuda and torch.cuda.is_available()

#     def fit(self, X_train, Y_train):
#         if isinstance(X_train, pd.DataFrame):
#             X_train = X_train.to_numpy()
#         if isinstance(Y_train, pd.DataFrame):
#             Y_train = Y_train.to_numpy()

#         X_train = torch.tensor(X_train, dtype=torch.float)
#         Y_train = torch.tensor(Y_train, dtype=torch.float).squeeze(-1)

#         if self.use_cuda:
#             X_train = X_train.cuda()
#             Y_train = Y_train.cuda()

#         # --- Priors ---
#         l_prior = GammaPrior(5.0, 5.0)
#         sigma_prior = NormalPrior(0.0, 1.0)
#         sigma_noise_prior = NormalPrior(0.0, 1.0)

#         self.likelihood = gpytorch.likelihoods.GaussianLikelihood(
#             noise_prior=sigma_noise_prior
#         )
#         self.model = GPMix(
#             X_train, 
#             Y_train, 
#             self.likelihood, 
#             self.feat_idx,
#             l_prior=l_prior,
#             sigma_prior=sigma_prior
#         )
        
#         if self.use_cuda:
#             self.model = self.model.cuda()
#             self.likelihood = self.likelihood.cuda()

#         self.model.train()
#         self.likelihood.train()

#         # --- MCMC Setup (following your example) ---
        
#         # 1. Define the pyro_model function locally
#         def pyro_model(x, y):
#             with gpytorch.settings.fast_computations(False, False, False):
#                 # Sample parameters from their priors
#                 sampled_model = self.model.pyro_sample_from_prior()
#                 # Get the model output (and sample from likelihood)
#                 output = sampled_model.likelihood(sampled_model(x))
#                 # Observe the data
#                 pyro.sample("obs", output, obs=y)
#             return y

#         # 2. Set up NUTS and MCMC
#         nuts_kernel = NUTS(pyro_model, jit_compile=True)
#         mcmc_run = MCMC(
#             nuts_kernel,
#             num_samples=self.num_samples,
#             warmup_steps=self.warmup_steps,
#         )
        
#         # 3. Run the sampler
#         mcmc_run.run(X_train, Y_train)

#         # 4. Load the samples into the model
#         samples = mcmc_run.get_samples()
#         self.model.pyro_load_from_samples(samples)
        
#         self.train_x = X_train
#         self.train_y = Y_train
#         return self

#     @torch.no_grad()
#     def predict(self, X_test):
#         if isinstance(X_test, pd.DataFrame):
#             X_test = X_test.to_numpy()
#         X_test = torch.tensor(X_test, dtype=torch.float)

#         if self.use_cuda:
#             X_test = X_test.cuda()

#         self.model.eval()
#         self.likelihood.eval()

#         pred_distribution = self.likelihood(self.model(X_test))

#         mean_pred = pred_distribution.mean.detach().cpu().numpy()
        
#         return mean_pred







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
