import pyro
import pyro.distributions as dist
import pyro.contrib.gp.kernels as pk
import torch
import pyro
import pyro.distributions as dist
import pyro.contrib.gp as gp
# from numpyro.infer import MCMC, NUTS
from pyro.infer import MCMC, NUTS, Predictive
from sklearn.base import BaseEstimator, RegressorMixin
import pandas as pd
import numpy as np
# import numpyro
from pyro.nn import PyroSample
from pyro.nn.module import PyroParam
from torch.distributions import constraints

def weighted_tanimoto_distance(x1, x2, eps=1e-6):
    x1e = x1.unsqueeze(-2)
    x2e = x2.unsqueeze(-3)

    numerator = torch.min(x1e, x2e).sum(dim=-1)
    denominator = torch.max(x1e, x2e).sum(dim=-1)

    sim = (numerator + eps) / (denominator + eps)
    dist = 1.0 - sim
    return torch.clamp(dist, min=0.)



class TanimotoRBF(pk.Kernel):
    def __init__(self, variance=None, lengthscale=None, active_dims=None):
        super().__init__(input_dim=len(active_dims), active_dims=active_dims)
        
        variance = torch.tensor(1.0) if variance is None else variance
        self.variance = PyroParam(variance, constraints.positive)
        lengthscale = torch.tensor(1.0) if lengthscale is None else lengthscale
        self.lengthscale = PyroParam(lengthscale, constraints.positive)

    def _diag(self, X):
        """
        Calculates the diagonal part of covariance matrix on active features.
        """
        return self.variance.expand(X.size(0))

    def forward(self, X, Z=None, diag=False):
        if diag:
            return self._diag(X)
        if Z is None:
            Z = X
        X = self._slice_input(X)
        Z = self._slice_input(Z)
        if X.size(1) != Z.size(1):
            raise ValueError("Inputs must have the same number of features.")

        D = weighted_tanimoto_distance(X, Z, eps=1e-6)
        return  self.variance * torch.exp(-0.5 * (D / self.lengthscale)**2)


class ProductWithVariance(pk.Product):
    def __init__(self, kern0, kern1, variance=None):
        super().__init__(kern0=kern0, kern1=kern1)

        variance = torch.tensor(1.0) if variance is None else variance
        self.variance = PyroParam(variance, constraints.positive)

    def forward(self, X, Z=None, diag=False):
        base_val = super().forward(X, Z, diag)
        return self.variance * base_val


class SumWithVariance(pk.Sum):
    def __init__(self, kern0, kern1, variance=None):
        super().__init__(kern0=kern0, kern1=kern1)

        variance = torch.tensor(1.0) if variance is None else variance
        self.variance = PyroParam(variance, constraints.positive)

    def forward(self, X, Z=None, diag=False):
        base_val = super().forward(X, Z, diag)
        return self.variance * base_val


class SumMultipleWithVariance(pk.Kernel):
    """
    Sum of multiple kernels with a single shared variance,
    in the style of Pyro's Combination kernel (kern0, kern1, kern2, ...).
    """
    def __init__(self, *kernels, variance=None):
        if len(kernels) < 2:
            raise ValueError("At least two kernels are required")

        # type checks
        for k in kernels:
            if not isinstance(k, pk.Kernel):
                raise TypeError("All components must be Kernel instances")

        # combine active_dims
        active_dims = set()
        for k in kernels:
            active_dims |= set(k.active_dims)
        active_dims = sorted(active_dims)

        input_dim = len(active_dims)
        super().__init__(input_dim=input_dim, active_dims=active_dims)

        for i, k in enumerate(kernels):
            setattr(self, f"kern{i}", k)
        self._kernels = kernels  # store in a list too

        variance = torch.tensor(1.0) if variance is None else variance
        self.variance = PyroParam(variance, constraints.positive)

    def forward(self, X, Z=None, diag=False):
        val = 0.0
        for k in self._kernels:
            val = val + k(X, Z, diag=diag)
        return self.variance * val


class MixingKernelPyro:
    def __init__(self, feat_idx, variance=None):
        self.fp_idx = feat_idx.get("fp1") or []
        self.cont_idx = feat_idx.get("count") or []
        self.variance = variance

    def build(self):
        
        # FP kernel: one shared LS
        if len(self.fp_idx) > 0:
            # k_fp = TanimotoRBF(
            #     active_dims=self.fp_idx  # use full fp block but single LS
            # )
            k_fp = pk.RBF(
                input_dim=len(self.fp_idx),
                active_dims=self.fp_idx
            )

        # Continuous kernel still ARD
        if len(self.cont_idx) > 0:
            k_cont = pk.Matern32(
                input_dim=len(self.cont_idx),
                active_dims=self.cont_idx
            )
        return SumMultipleWithVariance(k_fp, k_cont, variance=self.variance)


class GPMixPyro(gp.models.GPRegression):
    def __init__(self, X, y, feat_idx):
        self.fp_dim = len(feat_idx.get("fp1") or [])
        self.cont_dim = len(feat_idx.get("count") or [])
        self.kernel_builder = MixingKernelPyro(feat_idx)
        kernel = self.kernel_builder.build()
        super().__init__(X, y, kernel, jitter=1e-6)

        self.noise = PyroSample(dist.LogNormal(0.0, 1.0))
        self.kernel.variance = PyroSample(dist.LogNormal(0.0, 1.0))
        if self.fp_dim > 0 and self.cont_dim > 0:
            # fp block (scalar or 1-d)
            self.kernel.kern0.lengthscale = PyroSample(
                dist.InverseGamma(5.0, 5.0)
            )
            # continuous block, ARD
            self.kernel.kern1.lengthscale = PyroSample(
                dist.InverseGamma(5.0, 5.0)
                .expand([self.cont_dim])
                .to_event(1)
            )

        elif self.fp_dim > 0:
            # only fp kernel, one LS
            self.kernel.lengthscale = PyroSample(
                dist.InverseGamma(5.0, 5.0)
            )

        elif self.cont_dim > 0:
            # only cont kernel, ARD over cont_dim
            self.kernel.lengthscale = PyroSample(
                dist.InverseGamma(5.0, 5.0)
                .expand([self.cont_dim])
                .to_event(1)
            )


# ----------------------------------------------------------------------
# Inference Routine
# ----------------------------------------------------------------------

def run_inference(gp_model,
                num_samples,
                warmup_steps,
                num_chains,
                num_drawn_samples,
                random_state=None
                ):
    pyro.clear_param_store()
    if random_state is not None:
        pyro.set_rng_seed(random_state)

    nuts_kernel = NUTS(gp_model.model, ignore_jit_warnings=True, jit_compile=True)

    mcmc = MCMC(
        nuts_kernel,
        num_samples=num_samples,
        warmup_steps=warmup_steps,
        num_chains=num_chains,
        disable_progbar=True
    )

    mcmc.run()
    return mcmc.get_samples(num_samples=num_drawn_samples)



class GPMixMCMCRegressor(BaseEstimator, RegressorMixin):
    def __init__(
        self,
        feat_group=None,
        num_samples=200,
        warmup_steps=200,
        num_chains=1,
        num_drawn_samples=100,
        use_cuda=False,
        random_state=42,
    ):
        self.feat_group = feat_group
        self.num_samples = num_samples
        self.warmup_steps = warmup_steps
        self.num_chains = num_chains
        self.use_cuda = use_cuda
        self.random_state = random_state
        self.num_drawn_samples = num_drawn_samples

        self._is_fitted = False
        self._gp_model = None
        self._samples = None
        self._predictive = None
    def fit(self, X_train, y_train):
        if isinstance(X_train, pd.DataFrame):
            self.feat_idx = {
                        'fp1': [X_train.columns.get_loc(c) for c in self.feat_group.get("fp1")] if self.feat_group.get("fp1") else None,
                        'count': [X_train.columns.get_loc(c) for c in self.feat_group.get("count") ] if self.feat_group.get("count") else None
                    }

            self.count_feat_name_idx = {c: X_train.columns.get_loc(c) for c in self.feat_group.get("count")} if self.feat_group.get("count") else {}
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

        self._gp_model = GPMixPyro(X_t, y_t, self.feat_idx)

        # MCMC over gp_model.model
        samples = run_inference(
            self._gp_model,
            num_samples=self.num_samples,
            warmup_steps=self.warmup_steps,
            num_chains=self.num_chains,
            random_state=self.random_state,
            num_drawn_samples=self.num_drawn_samples,
        )

        # make sure samples live on the same device as the model
        self._samples = {k: v.to(device) for k, v in samples.items()}
        # Build Predictive object that will plug in posterior parameter samples
        self._is_fitted = True
        
        return self
    

    def _predictive_strategy(self, X_new):
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

        if "kernel.kern1.lengthscale" in self._samples:
            ls_count = (
                self._samples["kernel.kern1.lengthscale"]
                .float()
                .cpu()
                .numpy()
            ) 

            assert self.count_feat_name_idx is not None, "Feature indices not captured during fit."
            assert len(self.count_feat_name_idx) == ls_count.shape[1], \
                f"Count feature size mismatch: {len(self.count_feat_name_idx)} names vs {ls_count.shape[1]} dims."

            for name, idx in self.count_feat_name_idx.items():
                summary[name] = ls_count[:, idx]

        if "kernel.kern0.lengthscale" in self._samples:
            summary["fp1"] = (
                self._samples["kernel.kern0.lengthscale"]
                .float()
                .cpu()
                .numpy()
            )
        else:
            summary["fp1"] = None

        return summary
