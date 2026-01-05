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

    sim = torch.where(
        denominator > 0,
        numerator / (denominator + eps),
        torch.zeros_like(denominator),
    )
    dist = 1.0 - sim
    return torch.clamp(dist, min=0.)



class TanimotoRBF(pk.Kernel):
    def __init__(self, variance=None, lengthscale=None, input_dim=None, active_dims=None):
        super().__init__(input_dim=input_dim, active_dims=active_dims)
        
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



class TanimotoMatern(pk.Kernel):
    def __init__(self, variance=None, lengthscale=None, input_dim=None, active_dims=None, nu=1.5):
        super().__init__(input_dim=input_dim, active_dims=active_dims)
        
        variance = torch.tensor(1.0) if variance is None else variance
        self.variance = PyroParam(variance, constraints.positive)
        lengthscale = torch.tensor(1.0) if lengthscale is None else lengthscale
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

        D = weighted_tanimoto_distance(X, Z, eps=1e-6)
        r = D / self.lengthscale
        if self.nu == 1.5:
            sqrt3_r = 3**0.5 * r
            K = self.variance * (1.0 + sqrt3_r) * torch.exp(-sqrt3_r)
        elif self.nu == 2.5:
            sqrt5_r = 5**0.5 * r
            K = self.variance * (1.0 + sqrt5_r + (5.0/3.0)*r**2) * torch.exp(-sqrt5_r)
        else:
            raise ValueError("Unsupported nu value. Only 1.5 and 2.5 are supported.")
        return self.variance * K


class Tanimoto(pk.Kernel):
    def __init__(self, input_dim=None, active_dims=None):
        super().__init__(input_dim=input_dim, active_dims=active_dims)

    def _diag(self, X):
        """
        Calculates the diagonal part of covariance matrix on active features.
        """
        return torch.ones(X.size(0))

    def forward(self, X, Z=None, diag=False):
        if diag:
            return self._diag(X)
        if Z is None:
            Z = X
        X = self._slice_input(X)
        Z = self._slice_input(Z)
        if X.size(1) != Z.size(1):
            raise ValueError("Inputs must have the same number of features.")

        D = 1- weighted_tanimoto_distance(X, Z, eps=1e-6)
        D.clamp_min_(0)
        return  D



class ProductWithVariance(pk.Product):
    def __init__(self, kern0, kern1, variance=None):
        super().__init__(kern0=kern0, kern1=kern1)

        variance = torch.tensor(1.0) if variance is None else variance
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

        variance = torch.tensor(1.0) if variance is None else variance
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

        variance = torch.tensor(1.0) if variance is None else variance
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

        variance = torch.tensor(1.0) if variance is None else variance
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
                variance=None
                ):
        self.feat_idx = feat_idx
        self.variance = variance
        self.mixing_method = mixing_method
        self.kernel_method = kernel_method
    def build(self):
        fp_kernels = []
        count_kernels = []

        # Fingerprint kernels
        fp_keys = sorted([k for k in self.feat_idx.keys() if k.startswith("fp_")])
        fp_has_tanimoto = "tanimoto" in str(self.kernel_method["fp"]).lower()
        for key in fp_keys:
            idx = self.feat_idx[key]
            if fp_has_tanimoto:
                k_fp = kernel_factory[self.kernel_method["fp"]](
                        input_dim=len(idx),
                        active_dims=idx,
                    )
                fp_kernels.append(k_fp)
            else:
                for dim in idx:
                    k_fp = kernel_factory[self.kernel_method["fp"]](
                        input_dim=1,
                        active_dims=[dim],
                    )
                    fp_kernels.append(k_fp)

        # else:
        #     for dim in idx:
        #         k_fp = kernel_factory[self.kernel_method["fp"]](
        #                 input_dim=1,
        #                 active_dims=[dim],
        #             )
        #         fp_kernels.append(k_fp)
                
        # Count kernels
        count_idx = sorted(self.feat_idx.get("count"))
        for dim in count_idx:
            k_c = kernel_factory[self.kernel_method["count"]](
                input_dim=1,
                active_dims=[dim],
            )
            count_kernels.append(k_c)

        # Compose based on mixing method
        if self.mixing_method in ("sum", "product"):
            all_kernels = fp_kernels + count_kernels
            if len(all_kernels) < 2:
                raise ValueError(f"{self.mixing_method} mixing requires at least two kernels total")
            return mixing_factory[self.mixing_method](*all_kernels, variance=self.variance)

        if self.mixing_method == "averageProduct":
            if len(fp_kernels) < 1:
                raise ValueError("average-product requires at least one fp kernel (sum group)")
            if len(count_kernels) < 1:
                raise ValueError("average-product requires at least one count kernel (product group)")
            return AverageProductMultipleWithVariance(
                sum_kernels=count_kernels,
                product_kernels=fp_kernels,
                variance=self.variance,
                average_sum=True,
            )

        raise ValueError(f"Unknown mixing_method: {self.mixing_method}")



# class GPMixPyro(gp.models.GPRegression):
#     def __init__(self, X, y, feat_idx, mixing_method:str, kernel_method:dict):
#         self.feat_idx = feat_idx
#         self.kernel_method = kernel_method
#         self.kernel_builder = MixingKernelPyro(feat_idx, mixing_method, self.kernel_method)
#         kernel = self.kernel_builder.build()
#         super().__init__(X, y, kernel, jitter=1e-6)

#         self.noise = PyroSample(dist.LogNormal(0.0, 1.0))
#         self.kernel.variance = PyroSample(dist.LogNormal(0.0, 1.0))

#         # Dynamically assign lengthscales to all sub-kernels
#         fp_keys = sorted([k for k in feat_idx.keys() if k.startswith("fp_")])
        
#         for i, _ in enumerate(self.kernel._kernels):
#             target_kern = getattr(self.kernel, f"kern{i}")
            
#             # Check if this kernel corresponds to one of the FP groups
#             if i < len(fp_keys):
#                 if self.kernel_method["fp"].lower() == "tanimoto":
#                     continue
                
#                 target_kern.lengthscale = PyroSample(dist.InverseGamma(5.0, 5.0))
#             else:
#                 target_kern.lengthscale = PyroSample(dist.InverseGamma(5.0, 5.0))


class GPMixPyro(gp.models.GPRegression):
    def __init__(
        self,
        X,
        y,
        feat_idx,
        mixing_method: str,
        kernel_method: dict,
        # fp_lengthscale_prior=None,
        # count_lengthscale_prior=None,
    ):
        self.feat_idx = feat_idx
        self.kernel_method = kernel_method

        self.kernel_builder = MixingKernelPyro(feat_idx, mixing_method, self.kernel_method)
        kernel = self.kernel_builder.build()
        super().__init__(X, y, kernel, jitter=1e-6)

        self.noise = PyroSample(dist.LogNormal(0.0, 1.0))
        self.kernel.variance = PyroSample(dist.LogNormal(0.0, 1.0))

        # Default priors if user does not pass them
        # if fp_lengthscale_prior is None:
        #     fp_lengthscale_prior = dist.InverseGamma(5.0, 5.0)
        # if count_lengthscale_prior is None:
        #     count_lengthscale_prior = dist.InverseGamma(5.0, 5.0)

        fp_name = str(self.kernel_method["fp"]).strip()
        fp_name_l = fp_name.lower()

        # Exact "Tanimoto" means no lengthscale. "TanimotoRBF" (and others containing tanimoto) do need one.
        fp_is_plain_tanimoto = (fp_name_l == "tanimoto")
        fp_is_block_kernel = ("tanimoto" in fp_name_l)  # includes TanimotoRBF

        fp_keys = sorted([k for k in feat_idx.keys() if k.startswith("fp_")])

        # Number of FP kernels created at the front of kernel._kernels
        if fp_is_block_kernel:
            # One kernel per fp group (non-empty)
            n_fp = sum(1 for k in fp_keys if (feat_idx.get(k) or []))
        else:
            # One kernel per fp bit (dimension)
            n_fp = sum(len(feat_idx.get(k) or []) for k in fp_keys)

        # Number of count kernels (one per count dim in your build)
        n_count = len(feat_idx.get("count") or [])

        # Assign priors based on kernel position:
        # [0, n_fp) are FP kernels, [n_fp, n_fp + n_count) are count kernels
        for i, _ in enumerate(self.kernel._kernels):
            is_fp_kernel = (i < n_fp)
            is_count_kernel = (n_fp <= i < n_fp + n_count)
            target_kern = getattr(self.kernel, f"kern{i}")
            if is_fp_kernel:
                if fp_is_plain_tanimoto:
                    continue
                target_kern.lengthscale = PyroSample(dist.InverseGamma(5.0, 5.0))

            elif is_count_kernel:
                target_kern.lengthscale = PyroSample(dist.InverseGamma(5.0, 5.0))

            else:
                # If you ever add a third kernel family, decide here what prior it should receive.
                raise RuntimeError("Kernel index exceeds expected number of kernels.")



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
        disable_progbar=False
    )

    mcmc.run()
    return mcmc.get_samples(num_samples=num_drawn_samples)



class GPMixMCMCRegressor(BaseEstimator, RegressorMixin):
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

        self._gp_model = GPMixPyro(X_t, y_t, self.feat_idx, self.kernel_mixing_method, self.kernel_type)

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
    

    # def _get_lengthscale(self):
    #         summary = {}
    #         fp_keys = sorted([k for k in self.feat_idx.keys() if k.startswith("fp_")])
    #         count_keys = sorted(list(self.count_feat_name_idx.keys()))

    #         fp_has_tanimoto = "tanimoto" in str(self.kernel_method["fp"]).lower()

    #         all_keys = fp_keys + count_keys
    #         # Extract FP lengthscales using their specific unit names
    #         for i, key in enumerate(all_keys):
    #             param_name = f"kernel.kern{i}.lengthscale"
    #             if param_name in self._samples:
    #                 summary[key] = self._samples[param_name].float().cpu().numpy()

    #         return summary

    def _get_lengthscale(self):
        if self._samples is None:
            raise RuntimeError("No MCMC samples found. Call fit() first.")

        summary = {}

        fp_name = str(self.kernel_type["fp"]).strip()
        fp_name_l = fp_name.lower()

        # Exact "Tanimoto" => no lengthscale for FP kernels
        fp_is_plain_tanimoto = (fp_name_l == "tanimoto")

        # Any tanimoto-* kernel (e.g. TanimotoRBF) is a block kernel in your construction
        fp_is_block_kernel = ("tanimoto" in fp_name_l)

        fp_keys = sorted([k for k in self.feat_idx.keys() if k.startswith("fp_")])

        # Count keys in build order (recommended) or sorted; choose one and keep consistent with build()
        count_keys = list(self.feat_group.get("count", []) or [])
        # If you prefer your previous behavior:
        # count_keys = sorted(list(self.count_feat_name_idx.keys()))

        labels = []

        # FP labels (must match how FP kernels were built)
        if not fp_is_plain_tanimoto:
            if fp_is_block_kernel:
                # one kernel per fp group (non-empty)
                for g in fp_keys:
                    if self.feat_idx.get(g):
                        labels.append(g)
            else:
                # one kernel per fp bit/dimension (non-empty), ordered deterministically
                for g in fp_keys:
                    dims = sorted(list(self.feat_idx.get(g) or []))
                    for dim in dims:
                        labels.append(f"{g}[{dim}]")

        # Count labels (one kernel per count feature)
        labels.extend(count_keys)

        # Extract lengthscales by kernel index
        for i, label in enumerate(labels):
            site = f"kernel.kern{i}.lengthscale"
            if site in self._samples:
                summary[label] = (
                    self._samples[site]
                    .detach()
                    .float()
                    .cpu()
                    .numpy()
                )

        return summary
