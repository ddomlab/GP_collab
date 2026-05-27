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


def _as_torch_tensor(value, *, dtype=None, device=None):
    return torch.as_tensor(value, dtype=dtype or torch.float32, device=device)


def weighted_tanimoto(x1, x2, eps=1e-6,dist=True):
    x1e = x1.unsqueeze(-2)
    x2e = x2.unsqueeze(-3)

    numerator = torch.min(x1e, x2e).sum(dim=-1)
    denominator = torch.max(x1e, x2e).sum(dim=-1)

    sim = torch.where(
        denominator > 0,
        numerator / (denominator + eps),
        torch.zeros_like(denominator),
    )
    if dist:
        dist = 1.0 - sim
        return torch.clamp(dist, min=0.)
    return torch.clamp(sim, min=0.)


class TanimotoRBF(pk.Kernel):
    def __init__(self, variance=None, lengthscale=None, input_dim=None, active_dims=None):
        super().__init__(input_dim=input_dim, active_dims=active_dims)
        
        variance = _as_torch_tensor(1.0 if variance is None else variance)
        self.variance = PyroParam(variance, constraints.positive)
        lengthscale = _as_torch_tensor(1.0 if lengthscale is None else lengthscale)
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

        D = weighted_tanimoto(X, Z, eps=1e-6, dist=True)
        return  self.variance * torch.exp(-0.5 * (D / self.lengthscale)**2)



class TanimotoMatern(pk.Kernel):
    def __init__(self, variance=None, lengthscale=None, input_dim=None, active_dims=None, nu=1.5):
        super().__init__(input_dim=input_dim, active_dims=active_dims)
        
        variance = _as_torch_tensor(1.0 if variance is None else variance)
        self.variance = PyroParam(variance, constraints.positive)
        lengthscale = _as_torch_tensor(1.0 if lengthscale is None else lengthscale)
        self.lengthscale = PyroParam(lengthscale, constraints.positive)
        self.nu = nu

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

        D = weighted_tanimoto(X, Z, eps=1e-6, dist=True)
        r = D / self.lengthscale
        if self.nu == 1.5:
            sqrt3_r = 3**0.5 * r
            K = self.variance * (1.0 + sqrt3_r) * torch.exp(-sqrt3_r)
        elif self.nu == 2.5:
            sqrt5_r = 5**0.5 * r
            K = self.variance * (1.0 + sqrt5_r + (5.0/3.0)*r**2) * torch.exp(-sqrt5_r)
        else:
            raise ValueError("Unsupported nu value. Only 1.5 and 2.5 are supported.")
        return K


class Tanimoto(pk.Kernel):
    def __init__(self, input_dim=None, active_dims=None):
        super().__init__(input_dim=input_dim, active_dims=active_dims)

    def _diag(self, X):
        """
        Calculates the diagonal part of covariance matrix on active features.
        """
        return X.new_ones(X.size(0))

    def forward(self, X, Z=None, diag=False):
        if diag:
            return self._diag(X)
        if Z is None:
            Z = X
        X = self._slice_input(X)
        Z = self._slice_input(Z)
        if X.size(1) != Z.size(1):
            raise ValueError("Inputs must have the same number of features.")

        D = weighted_tanimoto(X, Z, eps=1e-6, dist=False)
        D.clamp_min_(0)
        return  D



class ProductWithVariance(pk.Product):
    def __init__(self, kern0, kern1, variance=None):
        super().__init__(kern0=kern0, kern1=kern1)

        variance = _as_torch_tensor(1.0 if variance is None else variance)
        self.variance = PyroParam(variance, constraints.positive)

    def forward(self, X, Z=None, diag=False):
        base_val = super().forward(X, Z, diag)
        return self.variance * base_val


class ProductMultipleWithVariance(pk.Kernel):
    """
    product of multiple kernels with a single shared variance,
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

        variance = _as_torch_tensor(1.0 if variance is None else variance)
        self.variance = PyroParam(variance, constraints.positive)

    def forward(self, X, Z=None, diag=False):
        val = 1.0
        for k in self._kernels:
            val = val * k(X, Z, diag=diag)
        return self.variance * val


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

        variance = _as_torch_tensor(1.0 if variance is None else variance)
        self.variance = PyroParam(variance, constraints.positive)

    def forward(self, X, Z=None, diag=False):
        val = 0.0
        for k in self._kernels:
            val = val + k(X, Z, diag=diag)
        return self.variance * val


class AverageProductMultipleWithVariance(pk.Kernel):
    """
    (sum over sum_kernels) * (product over product_kernels), scaled by a shared variance.

    Backward compatible with your old code:
      - exposes self._kernels (flattened)
      - registers kernels as kern0, kern1, ...
    Also provides aliases:
      - sum_kern0, sum_kern1, ...
      - prod_kern0, prod_kern1, ...
    """

    def __init__(self, sum_kernels, product_kernels, variance=None, average_sum=True):
        sum_kernels = list(sum_kernels) if sum_kernels is not None else []
        product_kernels = list(product_kernels) if product_kernels is not None else []

        if len(sum_kernels) < 1:
            raise ValueError("sum_kernels must contain at least one kernel")
        if len(product_kernels) < 1:
            raise ValueError("product_kernels must contain at least one kernel")

        for k in sum_kernels + product_kernels:
            if not isinstance(k, pk.Kernel):
                raise TypeError("All components must be Kernel instances")

        self._sum_kernels = sum_kernels
        self._product_kernels = product_kernels
        self.average_sum = bool(average_sum)

        # Union active_dims
        active_dims = set()
        for k in  product_kernels + sum_kernels:
            active_dims |= set(k.active_dims)
        active_dims = sorted(active_dims)

        super().__init__(input_dim=len(active_dims), active_dims=active_dims)

        # Backward-compatible flattened ordering: sum kernels first, then product kernels
        self._kernels = self._product_kernels + self._sum_kernels

        # Register under old names: kern0, kern1, ...
        for i, k in enumerate(self._kernels):
            setattr(self, f"kern{i}", k)

        # Optional aliases (helpful for debugging or explicit grouping)
        # for i, k in enumerate(self._sum_kernels):
        #     setattr(self, f"sum_kern{i}", k)
        # for j, k in enumerate(self._product_kernels):
        #     setattr(self, f"prod_kern{j}", k)

        variance = _as_torch_tensor(1.0 if variance is None else variance)
        self.variance = PyroParam(variance, constraints.positive)

    def forward(self, X, Z=None, diag=False):
        # product group
        prod_val = None
        for k in self._product_kernels:
            v = k(X, Z, diag=diag)
            prod_val = v if prod_val is None else (prod_val * v)
        
        # sum group
        sum_val = None
        for k in self._sum_kernels:
            v = k(X, Z, diag=diag)
            sum_val = v if sum_val is None else (sum_val + v)
        if self.average_sum:
            sum_val = sum_val / float(len(self._sum_kernels))


        return self.variance * prod_val *sum_val 



mixing_factory: dict = {
    "product": ProductMultipleWithVariance,
    "sum": SumMultipleWithVariance,
    "averageProduct": AverageProductMultipleWithVariance,
} 

kernel_factory: dict = {
    "TanimotoRBF": TanimotoRBF,
    "TanimotoMatern32": lambda **kwargs: TanimotoMatern(nu=1.5, **kwargs),
    "TanimotoMatern52": lambda **kwargs: TanimotoMatern(nu=2.5, **kwargs),
    "Tanimoto": Tanimoto,
    "RBF": pk.RBF,
    "Matern32": pk.Matern32,
    "Matern52": pk.Matern52,
}


class MixingKernelPyro:
    def __init__(self, 
                feat_idx,
                mixing_method: str,
                kernel_method: dict,
                variance=None,
                tensor_kwargs: dict | None = None,
                ):
        self.feat_idx = feat_idx
        self.variance = variance
        self.mixing_method = mixing_method
        self.kernel_method = kernel_method

        self.tensor_kwargs = tensor_kwargs or {}

    def _tensor(self, value):
        return _as_torch_tensor(value, **self.tensor_kwargs)

    def _kernel_kwargs(self, kernel_name, input_dim, active_dims):
        kwargs = {
            "input_dim": input_dim,
            "active_dims": active_dims,
        }

        if kernel_name.lower() != "tanimoto":
            kwargs["variance"] = self._tensor(1.0)
            kwargs["lengthscale"] = self._tensor(1.0)

        return kwargs

    def _variance_tensor(self):
        return self._tensor(1.0 if self.variance is None else self.variance)

    def build(self):
        fp_kernels = []
        count_kernels = []

        # Fingerprint kernels
        fp_keys = sorted([k for k in self.feat_idx.keys() if k.startswith("fp_")])
        for key in fp_keys:
            idx = self.feat_idx[key]
            if idx:
                k_fp = kernel_factory[self.kernel_method["fp"]](
                    **self._kernel_kwargs(
                        self.kernel_method["fp"],
                        input_dim=len(idx),
                        active_dims=idx,
                    )
                )
                fp_kernels.append(k_fp)

        # Count kernels
        count_idx = sorted(self.feat_idx.get("count", []))
        for dim in count_idx:
            k_c = kernel_factory[self.kernel_method["count"]](
                **self._kernel_kwargs(
                    self.kernel_method["count"],
                    input_dim=1,
                    active_dims=[dim],
                )
            )
            count_kernels.append(k_c)

        # Compose based on mixing method
        if self.mixing_method in ("sum", "product"):
            all_kernels = fp_kernels + count_kernels
            if len(all_kernels) < 2:
                raise ValueError(f"{self.mixing_method} mixing requires at least two kernels total")
            return mixing_factory[self.mixing_method](
                *all_kernels,
                variance=self._variance_tensor(),
            )

        if self.mixing_method == "averageProduct":
            if len(fp_kernels) < 1:
                raise ValueError("average-product requires at least one fp kernel (sum group)")
            if len(count_kernels) < 1:
                raise ValueError("average-product requires at least one count kernel (product group)")
            return AverageProductMultipleWithVariance(
                sum_kernels=count_kernels,
                product_kernels=fp_kernels,
                variance=self._variance_tensor(),
                average_sum=True,
            )

        raise ValueError(f"Unknown mixing_method: {self.mixing_method}")


class GPMixPyro(gp.models.GPRegression):
    def __init__(self, X, y, feat_idx, mixing_method:str, kernel_method:dict):
        self.feat_idx = feat_idx
        self.kernel_method = kernel_method
        self.tensor_kwargs = {
            "dtype": X.dtype,
            "device": X.device,
        }
        self.kernel_builder = MixingKernelPyro(
            feat_idx,
            mixing_method,
            self.kernel_method,
            tensor_kwargs=self.tensor_kwargs,
        )
        kernel = self.kernel_builder.build()
        super().__init__(X, y, kernel, jitter=1e-4)

        self.noise = self._lognormal_prior("noise_prior", self)
        self.kernel.variance = self._lognormal_prior("variance_prior", self.kernel)

        # Dynamically assign lengthscales to all sub-kernels
        fp_keys = sorted([k for k in feat_idx.keys() if k.startswith("fp_")])
        
        for i, _ in enumerate(self.kernel._kernels):
            target_kern = getattr(self.kernel, f"kern{i}")
            
            # Check if this kernel corresponds to one of the FP groups
            if i < len(fp_keys):
                if self.kernel_method["fp"].lower() == "tanimoto":
                    continue
                if "tanimoto" in self.kernel_method["fp"].lower():
                    target_kern.lengthscale = self._gamma_prior(
                        "lengthscale_prior",
                        target_kern,
                    )
                else:
                    ard_length = len(self.feat_idx[fp_keys[i]])
                    target_kern.lengthscale = self._gamma_prior(
                        "lengthscale_prior",
                        target_kern,
                        event_shape=(ard_length,),
                    )
            else:
                target_kern.lengthscale = self._gamma_prior(
                    "lengthscale_prior",
                    target_kern,
                )

    def _register_prior_buffer(self, module, name, value):
        tensor = _as_torch_tensor(value, **self.tensor_kwargs)
        module.register_buffer(name, tensor)
        return name

    def _lognormal_prior(self, prefix, module):
        loc_name = self._register_prior_buffer(module, f"_{prefix}_loc", 0.0)
        scale_name = self._register_prior_buffer(module, f"_{prefix}_scale", 1.0)

        return PyroSample(
            lambda m, loc_name=loc_name, scale_name=scale_name: dist.LogNormal(
                getattr(m, loc_name),
                getattr(m, scale_name),
            )
        )

    def _gamma_prior(
        self,
        prefix,
        module,
        concentration=5.0,
        rate=5.0,
        event_shape=None,
    ):
        concentration_name = self._register_prior_buffer(
            module,
            f"_{prefix}_concentration",
            concentration,
        )
        rate_name = self._register_prior_buffer(module, f"_{prefix}_rate", rate)

        if event_shape is None:
            return PyroSample(
                lambda m, concentration_name=concentration_name, rate_name=rate_name: dist.Gamma(
                    getattr(m, concentration_name),
                    getattr(m, rate_name),
                )
            )

        return PyroSample(
            lambda m,
            concentration_name=concentration_name,
            rate_name=rate_name,
            event_shape=event_shape: dist.Gamma(
                getattr(m, concentration_name),
                getattr(m, rate_name),
            ).expand(event_shape).to_event(len(event_shape))
        )


# ----------------------------------------------------------------------
# Inference Routine
# ----------------------------------------------------------------------

def run_inference(gp_model,
                num_samples,
                warmup_steps,
                num_chains,
                num_drawn_samples,
                random_state=42
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



class GpyroHMCRegressor:
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
        kernel_type:dict|None=None,
        normalize_y: bool = False,
    ):
        self.feat_group = feat_group
        self.num_samples = num_samples
        self.warmup_steps = warmup_steps
        self.num_chains = num_chains
        self.random_state = random_state
        self.num_drawn_samples = num_drawn_samples
        self.kernel_mixing_method = kernel_mixing_method
        self.kernel_type = (
            kernel_type
            if kernel_type is not None
            else {"fp": "TanimotoRBF", "count": "Matern32"}
        )
        self.normalize_y = normalize_y

        self.is_fitted_ = False
        self.gp_model_ = None
        self.samples_ = None

        self.use_cuda= use_cuda
        device = torch.device(
            "cuda" if self.use_cuda and torch.cuda.is_available() else "cpu"
        )
        self.cuda_avail = {"dtype": torch.float, "device": device} 

    @property
    def y(self):
        """
        Return training targets on the original, unnormalized scale.
        """
        if not hasattr(self, "_y_train"):
            raise AttributeError("Training targets do not exist. Call fit() first.")

        y_original = self._y_train.detach().cpu().numpy() * self.y_std_ + self.y_mean_
        return np.asarray(y_original).ravel()

    @y.setter
    def y(self, y):
        """
        Store normalized y internally if normalize_y=True.
        """
        if isinstance(y, (pd.DataFrame, pd.Series)):
            y = y.to_numpy()

        y = np.asarray(y, dtype=float).reshape(-1)

        if self.normalize_y:
            self.y_mean_ = float(np.mean(y))
            self.y_std_ = float(np.std(y))

            if self.y_std_ == 0:
                self.y_std_ = 1.0

            y_scaled = (y - self.y_mean_) / self.y_std_
        else:
            self.y_mean_ = 0.0
            self.y_std_ = 1.0
            y_scaled = y

        self._y_train = torch.as_tensor(
            y_scaled,
            **self.cuda_avail,
        ).view(-1)

    def _prepare_X(self, X: pd.DataFrame, fit: bool = False):
        """
        Convert a DataFrame X to a torch tensor.

        The sklearn wrapper supplies DataFrames, and feature groups are mapped
        to integer indices during fit.
        """
        if not isinstance(X, pd.DataFrame):
            raise TypeError(
                "GpyroHMCRegressor expects X to be a pandas DataFrame. "
                f"Got {type(X).__name__} instead."
            )

        if fit:
            self.feat_idx_ = {}
            for key, cols in self.feat_group.items():
                if cols:
                    missing_cols = [c for c in cols if c not in X.columns]
                    if missing_cols:
                        raise ValueError(
                            f"Columns for feature group '{key}' are missing "
                            f"from X: {missing_cols}"
                        )

                    self.feat_idx_[key] = [
                        X.columns.get_loc(c) for c in cols
                    ]

            count_cols = self.feat_group.get("count", [])
            missing_count_cols = [
                c for c in count_cols if c not in X.columns
            ]
            if missing_count_cols:
                raise ValueError(
                    f"Count columns are missing from X: {missing_count_cols}"
                )

            self.count_feat_name_idx_ = {
                c: X.columns.get_loc(c)
                for c in count_cols
            }
            self.feature_names_in_ = np.asarray(X.columns, dtype=object)

            # Backward-compatible aliases for existing analysis helpers.
            self.feat_idx = self.feat_idx_
            self.count_feat_name_idx = self.count_feat_name_idx_

        else:
            if not hasattr(self, "feature_names_in_"):
                raise RuntimeError("Model is not fitted. Call fit() first.")

            missing_cols = [
                c for c in self.feature_names_in_
                if c not in X.columns
            ]
            if missing_cols:
                raise ValueError(
                    f"Prediction X is missing columns seen during fit: {missing_cols}"
                )

            X = X.loc[:, self.feature_names_in_]

        return torch.as_tensor(
            X.to_numpy(),
            **self.cuda_avail,
        )

    def fit(self, X_train: pd.DataFrame, y_train):
        torch.manual_seed(self.random_state)
        np.random.seed(self.random_state)

        X_train = self._prepare_X(X_train, fit=True)
        self._X_train = X_train

        # Calls the y setter.
        self.y = y_train
        y_train = self._y_train

        # Pyro GP model with priors
        self.gp_model_ = GPMixPyro(
            X_train,
            y_train,
            self.feat_idx_,
            self.kernel_mixing_method,
            self.kernel_type,
        ).to(**self.cuda_avail)

        # MCMC over gp_model.model
        samples = run_inference(
            self.gp_model_,
            num_samples=self.num_samples,
            warmup_steps=self.warmup_steps,
            num_chains=self.num_chains,
            random_state=self.random_state,
            num_drawn_samples=self.num_drawn_samples,
        )

        # make sure samples live on the same device as the model
        self.samples_ = {
            k: v.to(self.cuda_avail["device"]) for k, v in samples.items()
        }
        self._samples = self.samples_
        self._gp_model = self.gp_model_

        # Build Predictive object that will plug in posterior parameter samples
        self.is_fitted_ = True
        self._is_fitted = self.is_fitted_
        
        return self
    

    def _predictive_strategy(self, X_new):
        """
        Method used by Pyro's Predictive class.
        Must be a method of the class to allow pickling.
        """
        # GP forward returns predictive mean and variance (including noise if noiseless=False)
        f_loc, f_var = self.gp_model_(X_new, full_cov=False, noiseless=False)

        # Stability: clamp variance to avoid sqrt of negative numbers due to float errors
        f_scale = f_var.clamp_min(1e-12).sqrt()

        # Sample observations
        pyro.sample(
            "y_pred",
            dist.Normal(f_loc, f_scale).to_event(1),
        )
            
    def _predictive_samples(self, X_test: pd.DataFrame):
        if not self.is_fitted_:
            raise RuntimeError("Call fit before predict.")

        X_t = self._prepare_X(X_test, fit=False)

        if X_t.ndim == 1:
            X_t = X_t.unsqueeze(0)

        predictive = Predictive(
            self._predictive_strategy,
            posterior_samples=self.samples_,
            return_sites=("y_pred",),
            parallel=False,
        )

        # draw predictive samples
        with torch.no_grad():
            samples_pred = predictive(X_t)
            y_samples = samples_pred["y_pred"]  # shape [S, N] or [S, N, 1]

        if y_samples.ndim == 3:
            y_samples = y_samples.squeeze(-1)

        return y_samples

    def predict(self, X_test: pd.DataFrame, return_std=False):
        y_samples = self._predictive_samples(X_test)

        mean_pred = y_samples.mean(dim=0)
        y_pred = mean_pred.detach().cpu().numpy() * self.y_std_ + self.y_mean_

        if return_std:
            std_pred = y_samples.std(dim=0, unbiased=False)
            y_std = std_pred.detach().cpu().numpy() * self.y_std_
        else:
            y_std = None

        return {
            "y_pred": np.asarray(y_pred).ravel(),
            "y_std": np.asarray(y_std).ravel() if y_std is not None else None,
        }


    def _get_latent_values(self, var_names):
        if isinstance(var_names, str):
            var_names = [var_names]

        return {k: self.samples_.get(k) for k in var_names}
    

    def _get_lengthscale(self):
        if not self.is_fitted_:
            raise RuntimeError("Model is not fitted. Call fit() first.")

        summary = {}
        fp_keys = sorted([k for k in self.feat_idx_.keys() if k.startswith("fp_")])
        count_names = [
            name
            for name, _ in sorted(
                self.count_feat_name_idx_.items(),
                key=lambda item: item[1],
            )
        ]

        fp_has_tanimoto = "tanimoto" in str(self.kernel_type["fp"]).lower()
        if fp_has_tanimoto:
            fp_name = fp_keys
        else:
            fp_name = [
                f"{key}[{dim}]"
                for key in fp_keys
                for dim in sorted(list(self.feat_idx_.get(key) or []))
            ]
        all_keys = fp_name + count_names
        # Extract FP lengthscales using their specific unit names
        for i, key in enumerate(all_keys):
            param_name = f"kernel.kern{i}.lengthscale"
            if param_name in self.samples_:
                summary[key] = self.samples_[param_name].float().cpu().numpy()

        return summary


class GpyroHMCsklearnRegressor(BaseEstimator, RegressorMixin):
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
        kernel_type:dict|None=None,
        normalize_y: bool = False,
    ):
        self.feat_group = feat_group
        self.num_samples = num_samples
        self.warmup_steps = warmup_steps
        self.num_chains = num_chains
        self.num_drawn_samples = num_drawn_samples
        self.use_cuda = use_cuda
        self.random_state = random_state
        self.kernel_mixing_method = kernel_mixing_method
        self.kernel_type = kernel_type
        self.normalize_y = normalize_y

    def fit(self, X, y):
        if not isinstance(X, pd.DataFrame):
            raise TypeError(
                "GpyroHMCsklearnRegressor expects X to be a pandas DataFrame. "
                f"Got {type(X).__name__} instead."
            )

        self.regressor_ = GpyroHMCRegressor(
            feat_group=self.feat_group,
            num_samples=self.num_samples,
            warmup_steps=self.warmup_steps,
            num_chains=self.num_chains,
            num_drawn_samples=self.num_drawn_samples,
            use_cuda=self.use_cuda,
            random_state=self.random_state,
            kernel_mixing_method=self.kernel_mixing_method,
            kernel_type=self.kernel_type,
            normalize_y=self.normalize_y,
        )

        self.regressor_.fit(X, y)

        self.feature_names_in_ = np.asarray(X.columns, dtype=object)
        self.n_features_in_ = X.shape[1]
        self.is_fitted_ = True

        return self

    def predict(self, X, return_std=False):
        if not isinstance(X, pd.DataFrame):
            raise TypeError(
                "GpyroHMCsklearnRegressor expects X to be a pandas DataFrame. "
                f"Got {type(X).__name__} instead."
            )

        return self.regressor_.predict(X, return_std=return_std)

    def _get_latent_values(self, var_names):
        return self.regressor_._get_latent_values(var_names)

    def _get_lengthscale(self):
        return self.regressor_._get_lengthscale()
