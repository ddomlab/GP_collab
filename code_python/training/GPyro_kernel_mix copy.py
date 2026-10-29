import pyro
import pyro.distributions as dist
import pyro.contrib.gp.kernels as pk
import torch
import pyro
import pyro.distributions as dist
import pyro.contrib.gp as gp
# from numpyro.infer import MCMC, NUTS
from pyro.infer import MCMC, NUTS
from sklearn.base import BaseEstimator, RegressorMixin
import pandas as pd
import numpy as np
# import numpyro
from pyro.nn import PyroSample
from pyro.nn.module import PyroParam
from torch.distributions import constraints
from GPytorch_kernel_mix import weighted_tanimoto


def _as_torch_tensor(value, *, dtype=None, device=None):
    return torch.as_tensor(value, dtype=dtype or torch.float32, device=device)


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



# class ProductWithVariance(pk.Product):
#     def __init__(self, kern0, kern1, variance=None):
#         super().__init__(kern0=kern0, kern1=kern1)

#         variance = _as_torch_tensor(1.0 if variance is None else variance)
#         self.variance = PyroParam(variance, constraints.positive)

#     def forward(self, X, Z=None, diag=False):
#         base_val = super().forward(X, Z, diag)
#         return self.variance * base_val


class ProductMultipleWithVariance(pk.Kernel):
    """
    product of multiple kernels with a single shared variance,
    in the style of Pyro's Combination kernel (kern0, kern1, kern2, ...).
    """
    def __init__(self, *kernels, variance=None, learn_variance=True):
        if len(kernels) < 1:
            raise ValueError("At least one kernel is required")

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
        if learn_variance:
            self.variance = PyroParam(variance, constraints.positive)
        else:
            self.register_buffer("variance", variance)

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
    def __init__(self, *kernels, variance=None, learn_variance=True, average=False):
        if len(kernels) < 1:
            raise ValueError("At least one kernel is required")

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
        self.average = bool(average)

        variance = _as_torch_tensor(1.0 if variance is None else variance)
        if learn_variance:
            self.variance = PyroParam(variance, constraints.positive)
        else:
            self.register_buffer("variance", variance)

    def forward(self, X, Z=None, diag=False):
        val = 0.0
        for k in self._kernels:
            val = val + k(X, Z, diag=diag)
        if self.average:
            val = val / float(len(self._kernels))
        return self.variance * val


class SumProductMultipleWithVariance(pk.Kernel):
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

    @staticmethod
    def _combine(kernels, kernel_cls, *, variance=1.0, learn_variance=False, **kwargs):
        if len(kernels) == 1 and not learn_variance:
            return kernels[0]
        return kernel_cls(
            *kernels,
            variance=variance,
            learn_variance=learn_variance,
            **kwargs,
        )

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

        if not fp_kernels and not count_kernels:
            raise ValueError(
                "No kernels could be created: provide at least one fingerprint "
                "or continuous ('count') feature."
            )

        # Compose based on mixing method
        if self.mixing_method in ("sum", "product"):
            all_kernels = fp_kernels + count_kernels
            return mixing_factory[self.mixing_method](
                *all_kernels,
                variance=self._variance_tensor(),
                learn_variance=True,
            )

        if self.mixing_method == "averageProduct":
            grouped_kernels = []
            if fp_kernels:
                grouped_kernels.append(self._combine(
                    fp_kernels,
                    ProductMultipleWithVariance,
                    variance=self._tensor(1.0),
                    learn_variance=False,
                ))
            if count_kernels:
                grouped_kernels.append(self._combine(
                    count_kernels,
                    SumMultipleWithVariance,
                    variance=self._tensor(1.0),
                    learn_variance=False,
                    average=True,
                ))
            return ProductMultipleWithVariance(
                *grouped_kernels,
                variance=self._variance_tensor(),
                learn_variance=True,
            )

        if self.mixing_method =="(count:+)x(fp:x)":
            grouped_kernels = []
            if fp_kernels:
                grouped_kernels.append(self._combine(
                    fp_kernels,
                    ProductMultipleWithVariance,
                    variance=self._tensor(1.0),
                    learn_variance=False,
                ))
            if count_kernels:
                grouped_kernels.append(self._combine(
                    count_kernels,
                    SumMultipleWithVariance,
                    variance=self._tensor(1.0),
                    learn_variance=False,
                ))
            return ProductMultipleWithVariance(
                *grouped_kernels,
                variance=self._variance_tensor(),
                learn_variance=True,
            )

        if self.mixing_method == "(count:+)x(fp:+)":
            grouped_kernels = []
            if fp_kernels:
                grouped_kernels.append(self._combine(
                    fp_kernels,
                    SumMultipleWithVariance,
                    variance=self._tensor(1.0),
                    learn_variance=False,
                ))
            if count_kernels:
                grouped_kernels.append(self._combine(
                    count_kernels,
                    SumMultipleWithVariance,
                    variance=self._tensor(1.0),
                    learn_variance=False,
                ))
            return ProductMultipleWithVariance(
                *grouped_kernels,
                variance=self._variance_tensor(),
                learn_variance=True,
            )

        if self.mixing_method == "(count:x)+(fp:x)":
            grouped_kernels = []
            if fp_kernels:
                grouped_kernels.append(self._combine(
                    fp_kernels,
                    ProductMultipleWithVariance,
                    variance=self._tensor(1.0),
                    learn_variance=False,
                ))
            if count_kernels:
                grouped_kernels.append(self._combine(
                    count_kernels,
                    ProductMultipleWithVariance,
                    variance=self._tensor(1.0),
                    learn_variance=False,
                ))
            return SumMultipleWithVariance(
                *grouped_kernels,
                variance=self._variance_tensor(),
                learn_variance=True,
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

        zero = X.new_tensor(0.0)
        one = X.new_tensor(1.0)
        five = X.new_tensor(5.0)

        self.noise = PyroSample(dist.LogNormal(zero, one))
        self.kernel.variance = PyroSample(dist.LogNormal(zero, one))
        has_fp = any(
            key.startswith("fp_") and indices
            for key, indices in feat_idx.items()
        )
        has_count = bool(feat_idx.get("count", []))
        if mixing_method == "averageProduct" and has_fp and has_count:
            # The count kernel computes the arithmetic mean of its component
            # kernels. Infer a multiplier for that mean, matching the trainable
            # inner ScaleKernel used by the GPyTorch implementation.
            self.kernel.kern1.variance = PyroSample(
                dist.LogNormal(zero, one)
            )

        # Dynamically assign lengthscales to all sub-kernels
        fp_keys = sorted(
            k for k in feat_idx.keys()
            if k.startswith("fp_") and feat_idx[k]
        )

        def _leaf_kernels(kernel):
            if hasattr(kernel, "_kernels"):
                leaves = []
                for sub_kernel in kernel._kernels:
                    leaves.extend(_leaf_kernels(sub_kernel))
                return leaves
            return [kernel]

        for i, target_kern in enumerate(_leaf_kernels(self.kernel)):
            if i < len(fp_keys):
                if self.kernel_method["fp"].lower() == "tanimoto":
                    continue
                if "tanimoto" in self.kernel_method["fp"].lower():
                    target_kern.lengthscale = PyroSample(
                        dist.Gamma(five, five)
                    )
                else:
                    ard_length = len(self.feat_idx[fp_keys[i]])
                    target_kern.lengthscale = PyroSample(
                        dist.Gamma(five, five)
                        .expand([ard_length])
                        .to_event(1)
                    )
            else:
                target_kern.lengthscale = PyroSample(
                    dist.Gamma(five, five)
                )


# ----------------------------------------------------------------------
# Inference Routine
# ----------------------------------------------------------------------

def run_inference(gp_model,
                num_samples,
                warmup_steps,
                num_chains,
                random_state=42,
                jit_compile=False,
                ):
    pyro.clear_param_store()
    if random_state is not None:
        pyro.set_rng_seed(random_state)

    nuts_kernel = NUTS(
        gp_model.model,
        ignore_jit_warnings=True,
        jit_compile=jit_compile,
    )

    mcmc = MCMC(
        nuts_kernel,
        num_samples=num_samples,
        warmup_steps=warmup_steps,
        num_chains=num_chains,
        disable_progbar=True
    )

    mcmc.run()
    return mcmc.get_samples()



class GpyroHMCRegressor:
    def __init__(
        self,
        feat_group:dict,
        num_samples=200,
        warmup_steps=200,
        num_chains=1,
        use_cuda=False,
        random_state=42,
        kernel_mixing_method:str="product",
        kernel_type:dict|None=None,
        normalize_y: bool = True,
        jit_compile: bool = True,
    ):
        self.feat_group = feat_group
        self.num_samples = num_samples
        self.warmup_steps = warmup_steps
        self.num_chains = num_chains
        self.random_state = random_state
        self.kernel_mixing_method = kernel_mixing_method
        self.kernel_type = (
            kernel_type
            if kernel_type is not None
            else {"fp": "TanimotoRBF", "count": "Matern32"}
        )
        self.normalize_y = normalize_y
        self.jit_compile = jit_compile

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

            self.count_feat_name_idx_ = {
                c: X.columns.get_loc(c)
                for c in count_cols
            }
            self.feature_names_in_ = np.asarray(X.columns, dtype=object)

            # Backward-compatible aliases for existing analysis helpers.
            self.feat_idx = self.feat_idx_
            self.count_feat_name_idx = self.count_feat_name_idx_

        else:
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
            jit_compile=self.jit_compile,
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
        Record conditional predictive moments for each HMC parameter draw.

        ``noiseless=False`` includes observation noise in the conditional
        variance, matching the GPyTorch MAP prediction path. Observation noise
        does not change the conditional mean.
        """
        f_loc, f_var = self.gp_model_(X_new, full_cov=False, noiseless=False)
        f_var = f_var.clamp_min(0.0)

        pyro.deterministic("predictive_mean", f_loc)
        pyro.deterministic(
            "predictive_variance",
            f_var,
        )

        return f_loc, f_var
            
    def _predictive_moments(self, X_test: pd.DataFrame):
        if not self.is_fitted_:
            raise RuntimeError("Call fit before predict.")

        X_t = self._prepare_X(X_test, fit=False)

        if X_t.ndim == 1:
            X_t = X_t.unsqueeze(0)

        if not self.samples_:
            raise RuntimeError("No posterior samples are available.")

        sample_counts = {sample.size(0) for sample in self.samples_.values()}
        if len(sample_counts) != 1:
            raise RuntimeError(
                "Posterior sample sites have inconsistent sample counts: "
                f"{sorted(sample_counts)}"
            )

        num_draws = sample_counts.pop()
        mean_samples = []
        variance_samples = []

        with torch.no_grad():
            # Pyro's Predictive first executes the strategy without conditioning
            # on posterior samples in order to infer plate nesting and output
            # shapes. For a GP, that extra prior draw can yield a nearly singular
            # training covariance and fail Cholesky before posterior prediction
            # begins. Condition each evaluation explicitly instead.
            for draw_idx in range(num_draws):
                draw = {
                    name: samples[draw_idx]
                    for name, samples in self.samples_.items()
                }
                conditioned_strategy = pyro.poutine.condition(
                    self._predictive_strategy,
                    data=draw,
                )
                draw_mean, draw_variance = conditioned_strategy(X_t)
                mean_samples.append(draw_mean)
                variance_samples.append(draw_variance)

            mean_samples = torch.stack(mean_samples, dim=0)
            variance_samples = torch.stack(variance_samples, dim=0)

        # Average the conditional GP distributions over HMC parameter draws.
        mean_pred = mean_samples.mean(dim=0)
        variance_pred = (
            (variance_samples + mean_samples.square()).mean(dim=0)
            - mean_pred.square()
        ).clamp_min(0.0)

        return mean_pred, variance_pred

    def predict(self, X_test: pd.DataFrame, return_std=False):
        mean_pred, variance_pred = self._predictive_moments(X_test)

        y_pred = mean_pred.detach().cpu().numpy() * self.y_std_ + self.y_mean_

        if return_std:
            std_pred = variance_pred.sqrt()
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
        fp_keys = sorted(
            k for k in self.feat_idx_.keys()
            if k.startswith("fp_") and self.feat_idx_[k]
        )
        count_names = [
            name
            for name, _ in sorted(
                self.count_feat_name_idx_.items(),
                key=lambda item: item[1],
            )
        ]

        def _leaf_kernel_paths(kernel, path="kernel"):
            if hasattr(kernel, "_kernels"):
                paths = []
                for i, sub_kernel in enumerate(kernel._kernels):
                    paths.extend(
                        _leaf_kernel_paths(sub_kernel, f"{path}.kern{i}")
                    )
                return paths
            return [(path, kernel)]

        def _add_lengthscale(key_prefix, param_name):
            if param_name in self.samples_:
                ls = self.samples_[param_name].float().cpu().numpy()
                if ls.ndim > 1 and ls.shape[-1] > 1:
                    for i in range(ls.shape[-1]):
                        summary[f"{key_prefix}[{i}]"] = ls[..., i]
                else:
                    summary[key_prefix] = ls

        leaf_paths = _leaf_kernel_paths(self.gp_model_.kernel)
        fp_leaf_paths = leaf_paths[:len(fp_keys)]
        count_leaf_paths = leaf_paths[len(fp_keys):]

        if len(count_leaf_paths) > len(count_names):
            raise IndexError(
                "More count kernels were found than expected from "
                "`count_feat_name_idx_`."
            )

        for fp_key, (path, _) in zip(fp_keys, fp_leaf_paths):
            if self.kernel_type["fp"].lower() == "tanimoto":
                summary[fp_key] = None
                continue
            _add_lengthscale(fp_key, f"{path}.lengthscale")

        for count_key, (path, _) in zip(count_names, count_leaf_paths):
            _add_lengthscale(count_key, f"{path}.lengthscale")

        return summary


"""
The predictive moments integrate the conditional GP distributions over HMC
parameter draws. ``y_pred`` is the average conditional mean, while ``y_std``
uses the law of total variance and includes:
1. variation in conditional means across HMC parameter draws,
2. conditional GP predictive variance at X_test, and
3. observation noise because ``noiseless=False``.

"""


class GpyroHMCsklearnRegressor(BaseEstimator, RegressorMixin):
    def __init__(
        self,
        feat_group:dict,
        num_samples=200,
        warmup_steps=200,
        num_chains=1,
        use_cuda=False,
        random_state=42,
        kernel_mixing_method:str="product",
        kernel_type:dict|None=None,
        normalize_y: bool = True,
        jit_compile: bool = True,
    ):
        self.feat_group = feat_group
        self.num_samples = num_samples
        self.warmup_steps = warmup_steps
        self.num_chains = num_chains
        self.use_cuda = use_cuda
        self.random_state = random_state
        self.kernel_mixing_method = kernel_mixing_method
        self.kernel_type = kernel_type
        self.normalize_y = normalize_y
        self.jit_compile = jit_compile

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
            use_cuda=self.use_cuda,
            random_state=self.random_state,
            kernel_mixing_method=self.kernel_mixing_method,
            kernel_type=self.kernel_type,
            normalize_y=self.normalize_y,
            jit_compile=self.jit_compile,
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
