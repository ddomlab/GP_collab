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
# from pyro.distributions import inverse_gamma as InvGammaPrior
from pyro.infer import MCMC, NUTS, Predictive
import pyro.distributions as dist
from tqdm import trange
from tqdm import tqdm
import os, sys
import gc
# from torch_geometric.data import Batch
# from torch_geometric.loader import DataLoader
# from pytorch_mpnn import DMPNNPredictor, RevIndexedData, smiles2data

#Torch modules
import torch
# from torch.distributions import Gamma
# from gpytorch.priors import Prior
# from gpytorch.module import Module as TModule
# from gpytorch.priors.utils import _bufferize_attributes

from Gpytorch_sskkernel import SubsequenceStringKernel

# class InverseGammaPrior(Prior, Gamma):
#     r"""
#     Inverse-Gamma prior parameterized by concentration (alpha) and rate (beta).

#     If X ~ InverseGamma(alpha, beta), then 1/X ~ Gamma(alpha, beta).

#     log p_X(x) = log p_Gamma(1/x) + log|d(1/x)/dx|
#                = log p_Gamma(1/x) - 2 log x
#     """

#     def __init__(self, concentration, rate, validate_args=False, transform=None):
#         TModule.__init__(self)
#         Gamma.__init__(self, concentration=concentration, rate=rate, validate_args=validate_args)
#         _bufferize_attributes(self, ("concentration", "rate"))
#         self._transform = transform

#     def expand(self, batch_shape):
#         batch_shape = torch.Size(batch_shape)
#         return InverseGammaPrior(
#             self.concentration.expand(batch_shape),
#             self.rate.expand(batch_shape),
#             transform=self._transform,
#         )

#     def log_prob(self, x):
#         # Apply any gpytorch transform first (same convention as Prior.log_prob)
#         x = self.transform(x)

#         inv_x = x.reciprocal()
#         # Gamma.log_prob on inv_x (use super(Prior, self) to jump to Gamma in the MRO)
#         gamma_lp = super(Prior, self).log_prob(inv_x)

#         # Jacobian term for inv_x = 1/x is |d(1/x)/dx| = 1/x^2  -> log = -2 log x
#         return gamma_lp - 2.0 * torch.log(x)

#     def rsample(self, sample_shape=torch.Size()):
#         # Sample from Gamma then invert
#         s = super(Prior, self).rsample(sample_shape)
#         return s.reciprocal()

#     def __call__(self, *args, **kwargs):
#         # Match gpytorch's torch_priors pattern
#         return super(Gamma, self).__call__(*args, **kwargs)


    


def weighted_tanimoto_distance(x1, x2, eps=1e-6):
    x1e = x1.unsqueeze(-2)
    x2e = x2.unsqueeze(-3)

    numerator = torch.min(x1e, x2e).sum(dim=-1)
    denominator = torch.max(x1e, x2e).sum(dim=-1)

    sim = (numerator + eps) / (denominator + eps)
    dist = 1.0 - sim
    return torch.clamp(dist, min=0.)


class Tanimoto(Kernel):
    has_lengthscale = False
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
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

            return torch.clamp(dist, min=0.)
        
        else:
            # if self.last_dim_is_batch:
                # x1 = x1.transpose(-1, -2).unsqueeze(-1)
                # x2 = x2.transpose(-1, -2).unsqueeze(-1)
            return weighted_tanimoto_distance(x1, x2)

class TanimotoRBF(Kernel):
    is_stationary = False
    has_lengthscale=True
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


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


class TanimotoMatern(Kernel):
    is_stationary = False
    has_lengthscale=True
    def __init__(self, nu,**kwargs):
        super().__init__(has_lengthscale=True, **kwargs)
        self.nu = nu

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
            r = dist / self.lengthscale
            if self.nu == 1.5:
                sqrt3_r = 3**0.5 * r
                K = (1.0 + sqrt3_r) * torch.exp(-sqrt3_r)
            elif self.nu == 2.5:
                sqrt5_r = 5**0.5 * r
                K = (1.0 + sqrt5_r + 5.0/3.0 * r**2) * torch.exp(-sqrt5_r)
            else:
                raise RuntimeError("nu expected to be 1.5 or 2.5")
            return K

        dist = weighted_tanimoto_distance(x1, x2)
        r = dist / self.lengthscale
        if self.nu == 1.5:
            sqrt3_r = 3**0.5 * r
            K = (1.0 + sqrt3_r) * torch.exp(-sqrt3_r)
        elif self.nu == 2.5:
            sqrt5_r = 5**0.5 * r
            K = (1.0 + sqrt5_r + 5.0/3.0 * r**2) * torch.exp(-sqrt5_r)
        else:
            raise RuntimeError("nu expected to be 1.5 or 2.5")
        return K
    
    
kernel_factory = {
    "TanimotoRBF": TanimotoRBF,
    "TanimotoMatern32": lambda **kw: TanimotoMatern(nu=1.5, **kw),
    "TanimotoMatern52": lambda **kw: TanimotoMatern(nu=2.5, **kw),
    "Tanimoto": Tanimoto,
    "RBF": RBFKernel,
    "Matern32": lambda **kw: MaternKernel(nu=1.5, **kw),
    "Matern52": lambda **kw: MaternKernel(nu=2.5, **kw),
    "SSK": SubsequenceStringKernel,
}


class MixingKernel:
    def __init__(
        self,
        feat_idx: dict,
        mixing_method: str,
        kernel_method: dict,
        variance: float | None = None,
        ssk_parameters: dict | None = None,
        cuda_avail: dict = None
    ):
        self.feat_idx = feat_idx
        self.mixing_method = mixing_method
        self.kernel_method = kernel_method
        self.variance = variance if variance else 1.0
        self.ssk_parameters = ssk_parameters
        self.cuda_avail = cuda_avail
    def _make_fp_kernels(self):
        fp_kernels = []
        fp_keys = sorted(k for k in self.feat_idx if k.startswith("fp_"))

        for key in fp_keys:
            idx = self.feat_idx[key]
            if idx:
                # print("ard dim for ", key, ":", len(idx))
                
                if "tanimoto" in self.kernel_method["fp"].lower():
                    k = kernel_factory[self.kernel_method["fp"]](
                        active_dims=idx,
                        # ard_num_dims=len(idx)
                    )

                else:

                    if self.kernel_method["fp"].lower() == "ssk":
                        ssk_params = self.ssk_parameters.get(key, {})
                        k = kernel_factory[self.kernel_method["fp"]](
                            active_dims=idx,
                            ard_num_dims=len(idx),
                            **ssk_params,
                            **self.cuda_avail
                        )
                    else:
                        k = kernel_factory[self.kernel_method["fp"]](
                            active_dims=idx,
                            ard_num_dims=len(idx)
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
                **self.cuda_avail
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
            avg_kernel = ScaleKernel(sum_kernel)
            avg_kernel.outputscale = 1.0 / len(count_kernels)
            kernel = ProductKernel(prod_kernel, avg_kernel)
            return kernel

        raise ValueError(f"Unknown mixing_method: {self.mixing_method}")



class GPMix(gpytorch.models.ExactGP):
    def __init__(self, X, y, feat_idx,
                  mixing_method:str, kernel_method:dict, likelihood, prior, ssk_parameters, cuda_avail):
        super().__init__(X, y, likelihood)
        # self.feat_idx = feat_idx
        # self.kernel_method = kernel_method
        # self.ssk_parameters = ssk_parameters
        self.mean_module = gpytorch.means.ZeroMean()

        self.kernel_builder = MixingKernel(
            feat_idx=feat_idx,
            mixing_method=mixing_method,
            kernel_method=kernel_method,
            ssk_parameters=ssk_parameters,
            cuda_avail=cuda_avail
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
        if prior:
            if mixing_method in ("sum", "product"):
                for i, (name, sk) in enumerate(root_kernel.named_sub_kernels()):
                    if i < len(fp_keys):
                        if kernel_method["fp"].lower() == "tanimoto":
                            continue
                        else:
                            # "tanimoto" in self.kernel_method["fp"].lower():
                            # ard_dim = len(feat_idx[fp_key])
                            sk.register_prior(
                                f"fp_lengthscale_prior_{i}",
                                gpytorch.priors.GammaPrior(5.0, 5.0),
                                "lengthscale",
                                )
                        

                    else:
                        sk.register_prior(
                            f"count_lengthscale_prior_{i}",
                            gpytorch.priors.GammaPrior(3.0, 6.0),
                            "lengthscale",
                        )

            elif mixing_method == "averageProduct":
                fp_product_kernel = root_kernel.kernels[0]
                count_sum_kernel = root_kernel.kernels[1].base_kernel
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
                    # if "tanimoto" in self.kernel_method["fp"].lower():
                    else:
                        # fp_key = fp_keys[i]
                        sk.register_prior(
                                    f"fp_lengthscale_prior_{i}",
                                    gpytorch.priors.GammaPrior(5.0, 5.0),
                                    "lengthscale",
                                    )
            else:
                raise ValueError(f"Unknown mixing_method: {mixing_method}")
            

    
    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)

    

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
        # use_cuda=False,
        random_state=42,
        kernel_mixing_method:str="product",
        kernel_type:dict={"fp":"TanimotoRBF", "count":"Matern32"},
        ssk_parameters:dict=None,
        progbar:bool=True,
        prior=False,
        # return_mll=False,
        # return_std=False
    ):
        self.feat_group = feat_group
        self.lr = lr
        self.n_epoch = n_epoch
        # self.use_cuda = use_cuda
        self.random_state = random_state
        self.kernel_mixing_method = kernel_mixing_method
        self.kernel_type = kernel_type
        self.ssk_parameters = ssk_parameters
        self.progbar = progbar
        self.prior = prior
        # self.return_mll = return_mll
        # self.return_std = return_std
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.cuda_avail = {"dtype": torch.float, "device": device} 
        
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

        # X_t = torch.as_tensor(np.asarray(X_train), dtype=torch.float32)
        # y_t = torch.as_tensor(np.asarray(y_train), dtype=torch.float32).view(-1)

        # if self.use_cuda and torch.cuda.is_available():
            # device = torch.device("cuda")
        X_train = torch.as_tensor(
            np.asarray(X_train),
            **self.cuda_avail
        )

        y_train = torch.as_tensor(
            np.asarray(y_train),    
            **self.cuda_avail
        ).view(-1)
                
        # Pyro GP model with priors
        self._likelihood = gpytorch.likelihoods.GaussianLikelihood(noise_constraint=gpytorch.constraints.Positive()).to(**self.cuda_avail)
        self._gp_model = GPMix(X_train, y_train, self.feat_idx,
                                self.kernel_mixing_method,
                                self.kernel_type,
                                self._likelihood,
                                prior=self.prior,
                                ssk_parameters=self.ssk_parameters,
                                cuda_avail=self.cuda_avail
                                ).to(**self.cuda_avail)
        
        # if self.use_cuda and torch.cuda.is_available():
        optimizer = torch.optim.Adam(self._gp_model.parameters(), lr=self.lr)
        mll = gpytorch.mlls.ExactMarginalLogLikelihood(self._likelihood, self._gp_model)

        # if torch.cuda.is_available():

        self._gp_model.train()
        self._likelihood.train()
        pbar = None
        if self.progbar:
            position=0
            pbar = CVProgressBar(total_steps=self.n_epoch, position=position)
        for _ in range(self.n_epoch):
            optimizer.zero_grad()

            y_pred = self._gp_model(X_train)
            loss = -mll(y_pred, y_train)

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

            X_t = torch.as_tensor(
                    np.asarray(X_test),
                    **self.cuda_avail
                    )

        self._gp_model.eval()
        self._likelihood.eval()
        with torch.no_grad():
            posterior = self._likelihood(self._gp_model(X_t))
            y_pred = posterior.mean.cpu().numpy()
            # if return_std:
            y_std = posterior.stddev.cpu().numpy() if return_std else None
            # y_mll = posterior.log_prob(y_test).mean().item() if return_mll and y_test is not None else None
            # test_mll = posterior.log_prob(self._gp_model._y_train).mean().item
            # if self.return_mll:
            #     # test_mll = posterior.log_prob(y_test).mean().item()

        return {
            "y_pred": np.asarray(y_pred).ravel(),
            "y_std": np.asarray(y_std).ravel() if y_std is not None else None,
        }

    def _get_lengthscale(self):
        summary = {}
        root_kernel = self._gp_model.covar_module.base_kernel
        fp_keys = sorted([k for k in self.feat_idx.keys() if k.startswith("fp_")])
        count_names = sorted(list(self.count_feat_name_idx.keys()))

        def _extract_ls(kernel, key_prefix):
            if not hasattr(kernel, "lengthscale") or kernel.lengthscale is None:
                raise RuntimeError(f"No lengthscale found for kernel '{key_prefix}'")
            ls = kernel.lengthscale.detach().cpu().numpy()
            # If ARD, return one key per dimension
            if ls.shape[-1] > 1:

                return {f"{key_prefix}[{i}]": ls_i for i, ls_i in enumerate(ls.squeeze(0))}
            else:
                return {key_prefix: ls}

        if self.kernel_mixing_method in ("sum", "product"):
            for i, (_, sk) in enumerate(root_kernel.named_sub_kernels()):
                if i < len(fp_keys):
                    fp_key = fp_keys[i]
                    if self.kernel_type["fp"].lower() == "tanimoto":
                        summary[fp_key] = None
                    else:
                        summary.update(_extract_ls(sk, fp_key))
                else:
                    count_key = count_names[i - len(fp_keys)]
                    summary.update(_extract_ls(sk, count_key))

        elif self.kernel_mixing_method == "averageProduct":
            fp_product_kernel = root_kernel.kernels[0]
            count_sum_kernel = root_kernel.kernels[1].base_kernel

            for i, (_, sk) in enumerate(count_sum_kernel.named_sub_kernels()):
                count_key = count_names[i]
                summary.update(_extract_ls(sk, count_key))

            # FP kernels
            for i, (_, sk) in enumerate(fp_product_kernel.named_sub_kernels()):
                fp_key = fp_keys[i]
                if self.kernel_type["fp"].lower() == "tanimoto":
                    summary[fp_key] = None
                else:
                    summary.update(_extract_ls(sk, fp_key))
                
        else:
            raise ValueError(f"Unknown mixing_method: {self.kernel_mixing_method}")
        
        return summary
    


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

