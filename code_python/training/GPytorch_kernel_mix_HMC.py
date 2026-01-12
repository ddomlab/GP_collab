import gpytorch
import numpy as np
import pandas as pd
import torch
# import torch.nn as nn
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.model_selection import train_test_split
from torch.utils import benchmark
from gpytorch.constraints import Positive
from gpytorch.kernels import Kernel, ProductKernel, AdditiveKernel
from gpytorch.kernels.rbf_kernel import postprocess_rbf
from gpytorch.priors import NormalPrior, GammaPrior
from gpytorch.kernels import Kernel, MaternKernel, ScaleKernel, RBFKernel
import pyro
from gpytorch.models import PyroGP # <-- Import PyroGP
from pyro.distributions import inverse_gamma as InvGammaPrior
from pyro.infer import MCMC, NUTS, Predictive
import pyro.distributions as dist
from tqdm import trange
from tqdm import tqdm
import os, sys
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


class TanimotoRBF(gpytorch.kernels.Kernel):
    is_stationary = False
    has_lengthscale=True
    def __init__(self, **kwargs):
        super().__init__(has_lengthscale=True, **kwargs)


    def forward(self, x1, x2, diag=False, **params):
        if diag:
            if torch.equal(x1, x2):
                return torch.ones(
                    *x1.shape[:-1],
                    device=x1.device,
                    dtype=x1.dtype
                )

            num = torch.min(x1, x2).sum(dim=-1)
            den = torch.max(x1, x2).sum(dim=-1)
            dist = 1.0 - (num + 1e-6) / (den + 1e-6)

            return torch.exp(
                -0.5 * (dist.clamp(min=0.) / self.lengthscale).pow(2)
            )

        dist = weighted_tanimoto_distance(x1, x2)
        return torch.exp(
            -0.5 * (dist / self.lengthscale).pow(2)
        )




kernel_factory = {
    "TanimotoRBF": TanimotoRBF,
    # "TanimotoMatern32": lambda **kw: TanimotoMatern(nu=1.5, **kw),
    # "TanimotoMatern52": lambda **kw: TanimotoMatern(nu=2.5, **kw),
    # "Tanimoto": Tanimoto,
    "RBF": RBFKernel,
    "Matern32": lambda **kw: MaternKernel(nu=1.5, **kw),
    "Matern52": lambda **kw: MaternKernel(nu=2.5, **kw),
}


class MixingKernel:
    def __init__(
        self,
        feat_idx: dict,
        mixing_method: str,
        kernel_method: dict,
        variance: float | None = None,
    ):
        self.feat_idx = feat_idx
        self.mixing_method = mixing_method
        self.kernel_method = kernel_method
        self.variance = variance if variance else 1.0

    def _make_fp_kernels(self):
        fp_kernels = []
        fp_keys = sorted(k for k in self.feat_idx if k.startswith("fp_"))

        for key in fp_keys:
            idx = self.feat_idx[key]
            if idx:
                # print(f"{idx}")
                # print("ard dim for ", key, ":", len(idx))
                k = kernel_factory[self.kernel_method["fp"]](
                    active_dims=idx,
                    # ard_num_dims=len(idx)
                )
                fp_kernels.append(k)
        # print(len(fp_kernels), "fp kernels created")
        return fp_kernels

    def _make_count_kernels(self):
        count_kernels = []
        count_idx = sorted(self.feat_idx.get("count", []))

        for dim in count_idx:
            k = kernel_factory[self.kernel_method["count"]](
                active_dims=[dim],
                # ard_num_dims=1,
            )
            count_kernels.append(k)
        return count_kernels
    
    def _apply_variance(self, kernel):
        if self.variance is not None:
            kernel = ScaleKernel(kernel)
            kernel.outputscale = self.variance
        return kernel

    def build(self):
        fp_kernels = self._make_fp_kernels()
        count_kernels = self._make_count_kernels()

        if self.mixing_method in ("sum", "product"):
            all_kernels = fp_kernels + count_kernels
            if len(all_kernels) < 2:
                raise ValueError(
                    f"{self.mixing_method} mixing requires at least two kernels"
                )

            if self.mixing_method == "sum":
                kernel = AdditiveKernel(*all_kernels)
            else:
                kernel = ProductKernel(*all_kernels)

            return kernel

        if self.mixing_method == "averageProduct":
            if len(fp_kernels) < 1:
                raise ValueError("averageProduct requires ≥1 fp kernel")
            if len(count_kernels) < 1:
                raise ValueError("averageProduct requires ≥1 count kernel")

            # sum over count kernels (optionally averaged)
            prod_kernel = ProductKernel(*fp_kernels)
            sum_kernel = AdditiveKernel(*count_kernels)

            kernel = ProductKernel(prod_kernel, sum_kernel)
            return kernel

        raise ValueError(f"Unknown mixing_method: {self.mixing_method}")



# k_prod = ProductKernel(
#     gpytorch.kernels.RBFKernel(ard_num_dims=3),
#     gpytorch.kernels.LinearKernel()
# )

# k_mult = AdditiveKernel(
#     gpytorch.kernels.RBFKernel(ard_num_dims=3),
#     gpytorch.kernels.LinearKernel()
# )
# k_total = ProductKernel(k_prod, k_mult)
# k_var = ScaleKernel(k_total)

# for name, sk in k_var.base_kernel.named_sub_kernels():
#     print("name\n",name)
#     print("sk\n",sk)

class GPMix(gpytorch.models.ExactGP):
    def __init__(self, X, y, feat_idx, mixing_method:str, kernel_method:dict, likelihood):
        super().__init__(X, y, likelihood)
        self.feat_idx = feat_idx
        self.kernel_method = kernel_method
        self.mean_module = gpytorch.means.ZeroMean()

        self.kernel_builder = MixingKernel(
            feat_idx=feat_idx,
            mixing_method=mixing_method,
            kernel_method=kernel_method,
        )
        base_kernel  = self.kernel_builder.build()
        self.covar_module = gpytorch.kernels.ScaleKernel(base_kernel)
        self.covar_module.register_prior(
            "variance_prior",
            gpytorch.priors.LogNormalPrior(0.0, 1.0),
            "outputscale",
        )
        self.likelihood.register_prior(
            "noise_prior",
            gpytorch.priors.LogNormalPrior(0.0, 1.0),
            "noise",
        )
        fp_keys = sorted(k for k in feat_idx if k.startswith("fp_"))

        root_kernel = self.covar_module.base_kernel
        if mixing_method in ("sum", "product"):
            for i, (name, sk) in enumerate(root_kernel.named_sub_kernels()):
                if i < len(fp_keys):
                    if kernel_method["fp"].lower() == "tanimoto":
                        continue
                    if "tanimoto" in self.kernel_method["fp"].lower():
                        # ard_dim = len(feat_idx[fp_key])
                        sk.register_prior(
                            f"fp_lengthscale_prior_{i}",
                            gpytorch.priors.GammaPrior(5.0, 5.0),
                            "lengthscale",
                            )
                        # sk.ard_num_dims = len(feat_idx[fp_keys[i]]) 
                        
                    else:
                        sk.register_prior(
                            f"fp_lengthscale_prior_{i}",
                            gpytorch.priors.GammaPrior(5.0, 5.0),
                            "lengthscale",
                            )
                        sk.ard_num_dims = len(feat_idx[fp_keys[i]])   

                else:
                    sk.register_prior(
                        f"count_lengthscale_prior_{i}",
                        gpytorch.priors.GammaPrior(5.0, 5.0),
                        "lengthscale",
                    )

        elif mixing_method == "averageProduct":
            fp_product_kernel = root_kernel.kernels[0]
            count_sum_kernel = root_kernel.kernels[1]
            for i, (name, sk) in enumerate(count_sum_kernel.named_sub_kernels()):
                sk.register_prior(
                            f"count_lengthscale_prior_{i}",
                            gpytorch.priors.GammaPrior(5.0, 5.0),
                            "lengthscale",
                            )
            fp_keys = sorted(k for k in feat_idx if k.startswith("fp_"))
            for i, (name, sk) in enumerate(fp_product_kernel.named_sub_kernels()):
                if kernel_method["fp"].lower() == "tanimoto":
                    continue
                fp_key = fp_keys[i]
                ard_dim = len(feat_idx[fp_key])
                sk.register_prior(
                            f"fp_lengthscale_prior_{i}",
                            gpytorch.priors.GammaPrior(5.0, 5.0),
                            "lengthscale",
                            )

                if ard_dim > 1:
                    sk.ard_num_dims = ard_dim
        else:
            raise ValueError(f"Unknown mixing_method: {mixing_method}")

    
    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)

    


# def run_inference(gp_model,
#                 num_samples,
#                 warmup_steps,
#                 num_chains,
#                 num_drawn_samples,
#                 random_state=42
#                 ):
#     pyro.clear_param_store()
#     if random_state is not None:
#         pyro.set_rng_seed(random_state)
#     with gpytorch.settings.fast_computations(False, False, False):
#         # with gpytorch.settings.cholesky_jitter(1e-2):
#             sampled_model = self.model.pyro_sample_from_prior()
#             output = sampled_model.likelihood(sampled_model(x))
#             pyro.sample("obs", output, obs=y)


#     nuts_kernel = NUTS(gp_model, ignore_jit_warnings=True, jit_compile=True)

#     mcmc = MCMC(
#         nuts_kernel,
#         num_samples=num_samples,
#         warmup_steps=warmup_steps,
#         num_chains=num_chains,
#         disable_progbar=True
#     )

#     mcmc.run()
#     return mcmc.get_samples(num_samples=num_drawn_samples)




class CVProgressBar:
    """CV- and multi-bar-safe progress bar for GP MAP training."""

    def __init__(self, total_steps, disable=False, position=0, desc="MAP training"):
        # Auto-disable in CI or pytest parallel runs
        disable = disable or "CI" in os.environ or "PYTEST_XDIST_WORKER" in os.environ

        self.pbar = tqdm(
            total=total_steps,
            desc=desc,
            leave=False,
            position=position,
            file=sys.stderr,
            disable=disable,
            bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}, {rate_fmt}{postfix}]",
        )

    def update(self, n=1, postfix=None):
        if postfix is not None:
            self.pbar.set_postfix(postfix)
        self.pbar.update(n)

    def close(self):
        self.pbar.close()


class GPytorchMAPRegressor(BaseEstimator, RegressorMixin):
    def __init__(
        self,
        feat_group:dict,
        lr=1e-2,
        n_epoch=400,
        use_cuda=False,
        random_state=42,
        kernel_mixing_method:str="product",
        kernel_type:dict={"fp":"TanimotoRBF", "count":"Matern32"},
        progbar:bool=False,
    ):
        self.feat_group = feat_group
        self.lr = lr
        self.n_epoch = n_epoch
        self.use_cuda = use_cuda
        self.random_state = random_state
        self.kernel_mixing_method = kernel_mixing_method
        self.kernel_type = kernel_type
        self.progbar = progbar
        
        
    def fit(self, X_train, y_train):
        if isinstance(X_train, pd.DataFrame):
            # Map all fp_{unit} and count to integer indices
            self.feat_idx = {}
            for key, cols in self.feat_group.items():
                if cols:
                    self.feat_idx[key] = [X_train.columns.get_loc(c) for c in cols]
            
            self.count_feat_name_idx = {
                c: X_train.columns.get_loc(c) 
                for c in self.feat_group.get("count", [])
            }
            X_train = X_train.to_numpy()

        if isinstance(y_train, (pd.DataFrame, pd.Series)):
            y_train = y_train.to_numpy()

        X_t = torch.as_tensor(np.asarray(X_train), dtype=torch.float32)
        y_t = torch.as_tensor(np.asarray(y_train), dtype=torch.float32).view(-1)

        if self.use_cuda and torch.cuda.is_available():
            device = torch.device("cuda")
            X_t = X_t.to(device)
            y_t = y_t.to(device)
        else:
            device = torch.device("cpu")

        self._X_train = X_t
        self._y_train = y_t

        # Pyro GP model with priors
        self._likelihood = gpytorch.likelihoods.GaussianLikelihood(noise_constraint=gpytorch.constraints.Positive())
        self._gp_model = GPMix(X_t, y_t, self.feat_idx, self.kernel_mixing_method, self.kernel_type, self._likelihood)
        
        optimizer = torch.optim.Adam(self._gp_model.parameters(), lr=self.lr)
        mll = gpytorch.mlls.ExactMarginalLogLikelihood(self._likelihood, self._gp_model)

        if self.use_cuda and torch.cuda.is_available():
            mll = mll.cuda()

        self._gp_model.train()
        self._likelihood.train()
        pbar = None
        if self.progbar:
            position=0
            pbar = CVProgressBar(total_steps=self.n_epoch, position=position)
        for _ in range(self.n_epoch):
            optimizer.zero_grad()

            y_pred = self._gp_model(self._X_train)
            loss = -mll(y_pred, self._y_train)

            loss.backward()
            optimizer.step()
            
            if pbar is not None:
                # postfix = {"loss": f"{loss.item():.3f}", "log_prior": f"{self._gp_model.log_prior().item():.3f}"}
                pbar.update(1)
        if pbar is not None:
            pbar.close()


    def predict(self, X_test, return_std=False):
        if isinstance(X_test, pd.DataFrame):
            X_test = X_test.to_numpy()

        X_t = torch.as_tensor(np.asarray(X_test), dtype=torch.float32)

        if torch.cuda.is_available():
            X_t = X_t.cuda()
            self._gp_model = self._gp_model.cuda()
            self._likelihood = self._likelihood.cuda()

        self._gp_model.eval()
        self._likelihood.eval()
        with torch.no_grad():
            posterior = self._likelihood(self._gp_model(X_t))
            y_pred = posterior.mean.cpu().numpy()
            if return_std:
                y_std = posterior.variance.sqrt().cpu().numpy()
                return y_pred, y_std

        return y_pred

























class GPytorchMCMCRegressor(BaseEstimator, RegressorMixin):
    def __init__(
        self,
        feat_group:dict,
        num_samples=200,
        warmup_steps=200,
        num_chains=1,
        num_drawn_samples=100,
        use_cuda=False,
        random_state=42,
        kernel_mixing_method:str="product",
        kernel_type:dict={"fp":"TanimotoRBF", "count":"Matern32"},
    ):
        self.feat_group = feat_group
        self.num_samples = num_samples
        self.warmup_steps = warmup_steps
        self.num_chains = num_chains
        self.use_cuda = use_cuda
        self.random_state = random_state
        self.num_drawn_samples = num_drawn_samples
        self.kernel_mixing_method = kernel_mixing_method
        self.kernel_type = kernel_type

        self._is_fitted = False
        self._gp_model = None
        self._samples = None
        self._predictive = None

    def fit(self, X_train, y_train):
        if isinstance(X_train, pd.DataFrame):
            # Map all fp_{unit} and count to integer indices
            self.feat_idx = {}
            for key, cols in self.feat_group.items():
                if cols:
                    self.feat_idx[key] = [X_train.columns.get_loc(c) for c in cols]
            
            self.count_feat_name_idx = {
                c: X_train.columns.get_loc(c) 
                for c in self.feat_group.get("count", [])
            }
            X_train = X_train.to_numpy()

        if isinstance(y_train, (pd.DataFrame, pd.Series)):
            y_train = y_train.to_numpy()

        X_t = torch.as_tensor(np.asarray(X_train), dtype=torch.float32)
        y_t = torch.as_tensor(np.asarray(y_train), dtype=torch.float32).view(-1)

        if self.use_cuda and torch.cuda.is_available():
            device = torch.device("cuda")
            X_t = X_t.to(device)
            y_t = y_t.to(device)
        else:
            device = torch.device("cpu")

        self._X_train = X_t
        self._y_train = y_t

        # Pyro GP model with priors
        self.likelihood = gpytorch.likelihoods.GaussianLikelihood(noise_constraint=gpytorch.constraints.Positive())
        self._gp_model = GPMix(X_t, y_t, self.feat_idx, self.kernel_mixing_method, self.kernel_type, self.likelihood)

        def pyro_model(x, y):
                    sampled_model = self._gp_model.pyro_sample_from_prior()
                    output = sampled_model.likelihood(sampled_model(x))
                    pyro.sample("obs", output, obs=y)
                    return y
        pyro.clear_param_store()
        if self.random_state is not None:
            pyro.set_rng_seed(self.random_state)


        nuts_kernel = NUTS(pyro_model, ignore_jit_warnings=True, jit_compile=False)
        mcmc = MCMC(
            nuts_kernel, 
            num_samples=self.num_samples, 
            warmup_steps=self.warmup_steps, 
            num_chains=self.num_chains,
            # disable_progbar=True,
            # hook_fn=hook,
        )
        mcmc.run(self._X_train, self._y_train)
        samples =mcmc.get_samples(num_samples=self.num_drawn_samples)

        # make sure samples live on the same device as the model
        self._samples = {k: v.to(device) for k, v in samples.items()}
        # Build Predictive object that will plug in posterior parameter samples
        self._is_fitted = True
        return self
    
    def _predictive_strategy(self, X_new):
            print("ready to go")
            """
            Method used by Pyro's Predictive class. 
            Must be a method of the class to allow pickling.
            """
            # GP forward returns predictive mean and variance (including noise if noiseless=False)
            f_loc, f_var = self._gp_model(X_new, full_cov=False, noiseless=False)
            
            # Stability: clamp variance to avoid sqrt of negative numbers due to float errors
            f_scale = f_var.sqrt()

            # Sample observations
            pyro.sample(
                "y_pred",
                dist.Normal(f_loc, f_scale).to_event(1),
            )
            
    def predict(self, X_test):
        if not self._is_fitted:
            raise RuntimeError("Call fit before predict.")

        if isinstance(X_test, pd.DataFrame):
            X_test = X_test.to_numpy()

        X_t = torch.as_tensor(np.asarray(X_test), dtype=torch.float32)

        if X_t.ndim == 1:
            X_t = X_t.unsqueeze(0)

        device = next(self._gp_model.parameters()).device
        X_t = X_t.to(device)
        predictive = Predictive(
            self._predictive_strategy,
            posterior_samples=self._samples,
            return_sites=("y_pred",),
            parallel=False,
        )
        # draw predictive samples
        with torch.no_grad():
            samples_pred = predictive(X_t)
            y_samples = samples_pred["y_pred"]  # shape [S, N] or [S, N, 1]

        if y_samples.ndim == 3:
            y_samples = y_samples.squeeze(-1)

        mean_pred = y_samples.mean(dim=0)
        # if later you want std:
        # std_pred = y_samples.std(dim=0)

        return mean_pred.detach().cpu().numpy()


    def _get_latent_values(self, var_names):
        if isinstance(var_names, str):
            var_names = [var_names]

        return {k: self._samples.get(k) for k in var_names}
    

    def _get_lengthscale(self):
            summary = {}
            fp_keys = sorted([k for k in self.feat_idx.keys() if k.startswith("fp_")])
            count_names = sorted(list(self.count_feat_name_idx.keys()))

            fp_has_tanimoto = "tanimoto" in str(self.kernel_type["fp"]).lower()
            if fp_has_tanimoto:
                fp_name = fp_keys
            else:
                fp_name = [
                    f"{key}[{dim}]"
                    for key in fp_keys
                    for dim in sorted(list(self.feat_idx.get(key) or []))
                ]
            all_keys = fp_name + count_names
            # Extract FP lengthscales using their specific unit names
            for i, key in enumerate(all_keys):
                param_name = f"kernel.kern{i}.lengthscale"
                if param_name in self._samples:
                    summary[key] = self._samples[param_name].float().cpu().numpy()

            return summary

# class GPytorchMixMCMCRegressor(BaseEstimator, RegressorMixin):
#     def __init__(self, 
#                 feat_group=None, 
#                 num_samples=200,
#                 warmup_steps=200,
#                 num_drawn_samples=100,
#                 num_chains=1,
#                 use_cuda=False,
#                 random_state=42
#         ):
        
#         self.feat_group = feat_group
#         self.num_samples = num_samples
#         self.warmup_steps = warmup_steps
#         self.num_drawn_samples = num_drawn_samples
#         self.num_chains = num_chains
#         self.use_cuda = use_cuda
#         self.model = None
#         self.likelihood = None
#         self.mcmc_samples = None
#         self.random_state = random_state


#     def fit(self, X, y):
#         if isinstance(X, pd.DataFrame):
#             self.feat_idx = {
#                         'fp': [X.columns.get_loc(c) for c in self.feat_group.get("fp")] if self.feat_group.get("fp") else None,
#                         'count': [X.columns.get_loc(c) for c in self.feat_group.get("count") ] if self.feat_group.get("count") else None
#                     }
#             self.count_feat_name_idx = {c: X.columns.get_loc(c) for c in self.feat_group.get("count")} if self.feat_group.get("count") else {}
#             X = X.to_numpy()
#         if self.feat_idx is None:
#              raise ValueError("feat_idx must be provided")

#         if isinstance(y, (pd.DataFrame, pd.Series)):
#             y = y.values

#         X = torch.as_tensor(X, dtype=torch.float32)
#         y = torch.as_tensor(y, dtype=torch.float32).view(-1)

#         if self.use_cuda and torch.cuda.is_available():
#             X = X.cuda()
#             y = y.cuda()

#         self.train_x = X
#         self.train_y = y

#         self.likelihood = gpytorch.likelihoods.GaussianLikelihood(noise_constraint=gpytorch.constraints.Positive())
#         self.model = GPMix(X, y, self.likelihood, self.feat_idx)
#         self.likelihood.register_prior(
#             "noise_prior", 
#             gpytorch.priors.LogNormalPrior(0.0, 1.0), 
#             "noise"
#         )
        
#         if self.use_cuda and torch.cuda.is_available():
#             self.model = self.model.cuda()
#             self.likelihood = self.likelihood.cuda()


#         # mll = gpytorch.mlls.ExactMarginalLogLikelihood(self.likelihood, self.model)

#         def pyro_model(x, y):
#             with gpytorch.settings.fast_computations(False, False, False):
#                 with gpytorch.settings.cholesky_jitter(1e-2):
#                     sampled_model = self.model.pyro_sample_from_prior()
#                     output = sampled_model.likelihood(sampled_model(x))
#                     pyro.sample("obs", output, obs=y)
#             return y
#         pyro.clear_param_store()
#         if self.random_state is not None:
#             pyro.set_rng_seed(self.random_state)
#         def hook(kernel, samples, stage, i):
#             if stage == "Warmup" and i % 50 == 0:
#                 print(i, samples["likelihood.noise"], samples["covar_module.variance"])
#         nuts_kernel = NUTS(pyro_model, ignore_jit_warnings=True, jit_compile=False)
#         mcmc = MCMC(
#             nuts_kernel, 
#             num_samples=self.num_samples, 
#             warmup_steps=self.warmup_steps, 
#             num_chains=self.num_chains,
#             # disable_progbar=True,
#             hook_fn=hook,
#         )

#         mcmc.run(X, y)
#         self.mcmc_samples = mcmc.get_samples(num_samples=self.num_drawn_samples)
#         return self

#     def predict(self, X_test, return_std=False):
#         if isinstance(X_test, pd.DataFrame):
#             X_test = X_test.to_numpy()

#         X_test = torch.as_tensor(X_test, dtype=torch.float32)
#         if self.use_cuda and torch.cuda.is_available():
#             X_test = X_test.cuda()

#         self.model.eval()
#         # self.likelihood.eval()

#         self.model.pyro_load_from_samples(self.mcmc_samples)
#         self.model.eval()
#         num_samples = list(self.mcmc_samples.values())[0].shape[0]
#         expanded_test_x = X_test.unsqueeze(0).expand(num_samples, *X_test.shape)

#         with torch.no_grad():
#             output = self.model(expanded_test_x)
#             observed_output = self.likelihood(output) 

#             means = observed_output.mean
#             mean_prediction = means.mean(0) 
#         return mean_prediction.cpu().numpy()
    
    
