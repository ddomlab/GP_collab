import pyro
import pyro.distributions as dist
import pyro.contrib.gp.kernels as pk
import torch
import pyro
import pyro.distributions as dist
import pyro.contrib.gp as gp
# from numpyro.infer import MCMC, NUTS
from pyro.infer import MCMC
from sklearn.base import BaseEstimator, RegressorMixin
import pandas as pd
import numpy as np
import numpyro

def weighted_tanimoto_distance(x1, x2, eps=1e-6):
    x1e = x1.unsqueeze(-2)
    x2e = x2.unsqueeze(-3)

    numerator = torch.min(x1e, x2e).sum(dim=-1)
    denominator = torch.max(x1e, x2e).sum(dim=-1)

    sim = (numerator + eps) / (denominator + eps)
    dist = 1.0 - sim
    return torch.clamp(dist, min=0.)



class TanimotoRBF(pk.Kernel):
    def __init__(self, active_dims, eps=1e-6):
        super().__init__(input_dim=len(active_dims), active_dims=active_dims)
        self.eps = eps

        # plain PyTorch tensor placeholder (no pyro.param)
        # self.lengthscale = torch.tensor(1.0)

    def forward(self, X, Z=None, diag=False):
        X = self._slice_input(X)
        Z = self._slice_input(Z) if Z is not None else X

        if diag:
            return X.new_ones(X.shape[0])

        D = weighted_tanimoto_distance(X, Z, eps=self.eps)

        ls = torch.clamp(self.lengthscale, min=1e-6)
        return torch.exp(-0.5 * (D / ls)**2)

    
    
class MixingKernelPyro:
    def __init__(self, feat_idx):
        self.fp_idx = feat_idx.get("fp") or []
        self.cont_idx = feat_idx.get("count") or []

    def build(self):
        kernels = []

        # FP kernel: one shared LS
        if len(self.fp_idx) > 0:
            k_fp = TanimotoRBF(
                active_dims=self.fp_idx  # use full fp block but single LS
            )
            kernels.append(k_fp)

        # Continuous kernel still ARD
        if len(self.cont_idx) > 0:
            k_cont = pk.Matern32(
                input_dim=len(self.cont_idx),
                active_dims=self.cont_idx
            )
            kernels.append(k_cont)

        if not kernels:
            raise ValueError("Both feature groups empty")

        k = kernels[0]
        for other in kernels[1:]:
            k = pk.Product(k, other)

        return k




class GPMixPyro:
    def __init__(self, X, y, feat_idx):
        self.X = X
        self.y = y
        self.fp_dim = len(feat_idx.get("fp") or [])
        self.cont_dim = len(feat_idx.get("count") or [])
        self.kernel_builder = MixingKernelPyro(feat_idx)

    def model(self):
        # priors
        # dist.InverseGamma
        sigma_signal = pyro.sample("sigma_signal", dist.LogNormal(0,1.0))
        outputscale_variance = sigma_signal.pow(2)
        sigma_noise = pyro.sample("obs_noise_std", dist.LogNormal(0,1.0))
        obs_noise_variance = sigma_noise.pow(2)

        if self.fp_dim > 0:
            fp_ls = pyro.sample(
                "fp_lengthscale",
                dist.InverseGamma(5, 5)
            )

        if self.cont_dim > 0:
            cont_ls = pyro.sample(
                "cont_lengthscale",
                dist.InverseGamma(5, 5).expand([self.cont_dim]).to_event(1)
            )

        # build fresh kernel
        kernel = self.kernel_builder.build()

        # assign sampled params
        if self.fp_dim > 0:
            kernel.kern0.lengthscale = fp_ls   # FP kernel is kern0

        if self.cont_dim > 0:
            # automatically kern1 for Product(Tanimoto, Matern)
            kernel.kern1.lengthscale = cont_ls  

        # kernel.variance = outputscale
        kernel.variance = outputscale_variance
        # GPRegression
        gpmodel = gp.models.GPRegression(
            self.X,
            self.y,
            kernel,
            noise=obs_noise_variance,
            jitter=1e-4
        )

        gpmodel.model()


# ----------------------------------------------------------------------
# Inference Routine
# ----------------------------------------------------------------------

def run_inference(
                gp_model, 
                num_samples=300,
                warmup_steps=200,
                num_chains=1, 
                num_drawn_samples=100,
                random_state=None):
    pyro.clear_param_store()
    if random_state is not None:
        pyro.set_rng_seed(random_state)

    nuts_kernel = NUTS(gp_model.model,
                       ignore_jit_warnings=True,
                       jit_compile=True)

    mcmc = MCMC(
        nuts_kernel,
        num_samples=num_samples,
        warmup_steps=warmup_steps,
        num_chains=num_chains
    )

    mcmc.run()
    return mcmc.get_samples(num_samples=num_drawn_samples)


# ----------------------------------------------------------------------
# Posterior Prediction
# ----------------------------------------------------------------------

def predict_posterior(X_train, y_train, X_test, samples, kernel_builder):
    mean_list = []
    std_list = []

    if isinstance(X_test, torch.Tensor) and X_test.ndim == 1:
        X_test = X_test.unsqueeze(0)

    num_draws = samples["sigma_signal"].shape[0]

    fp_dim = len(kernel_builder.fp_idx)
    cont_dim = len(kernel_builder.cont_idx)

    for i in range(num_draws):
        sig_signal = samples["sigma_signal"][i]
        outputscale_variance = sig_signal.pow(2)

        sig_noise = samples["obs_noise_std"][i]
        obs_noise_variance = sig_noise.pow(2)

        kernel = kernel_builder.build()

        if fp_dim > 0:
            fp_ls = samples["fp_lengthscale"][i]
            kernel.kern0.lengthscale = fp_ls

        if cont_dim > 0:
            if fp_dim > 0:
                cont_ls = samples["cont_lengthscale"][i]
                kernel.kern1.lengthscale = cont_ls
            else:
                cont_ls = samples["cont_lengthscale"][i]
                kernel.kern0.lengthscale = cont_ls

        kernel.variance = outputscale_variance

        gpmodel = gp.models.GPRegression(
            X_train,
            y_train,
            kernel,
            noise=obs_noise_variance
        )

        f_dist = gpmodel(X_test, full_cov=False)

        if isinstance(f_dist, tuple):
            mean, var = f_dist
        else:
            mean = f_dist.mean
            var = f_dist.variance

        mean_list.append(mean)
        std_list.append(var.sqrt())

    mean_stack = torch.stack(mean_list)
    std_stack = torch.stack(std_list)

    mean_pred = mean_stack.mean(dim=0)
    std_pred = std_stack.mean(dim=0)

    return mean_pred, std_pred


# ----------------------------------------------------------------------
# Sklearn Wrapper
# ----------------------------------------------------------------------

class GPMixMCMCRegressor(BaseEstimator, RegressorMixin):
    def __init__(
        self,
        feat_idx=None,
        num_samples=1000,
        warmup_steps=1000,
        num_chains=4,
        num_drawn_samples=500,
        use_cuda=False,
        random_state=42,
    ):
        # store all parameters exactly as sklearn expects
        self.feat_idx = feat_idx
        self.num_samples = num_samples
        self.warmup_steps = warmup_steps
        self.num_chains = num_chains
        self.use_cuda = use_cuda
        self.random_state = random_state
        self.num_drawn_samples = num_drawn_samples


    def fit(self, X_train, y_train):
        # convert DataFrame to numpy
        if isinstance(X_train, pd.DataFrame):
            X_train = X_train.to_numpy()
        if isinstance(y_train, pd.DataFrame):
            y_train = y_train.to_numpy()

        # convert to tensors
        X_t = torch.as_tensor(np.asarray(X_train), dtype=torch.float32)
        y_t = torch.as_tensor(np.asarray(y_train), dtype=torch.float32).view(-1)

        if self.use_cuda and torch.cuda.is_available():
            X_t = X_t.cuda()
            y_t = y_t.cuda()

        self._X_train = X_t
        self._y_train = y_t

        # build Pyro GP model
        self._gp_model = GPMixPyro(X_t, y_t, self.feat_idx)

        # run MCMC
        self._samples = run_inference(
            self._gp_model,
            num_samples=self.num_samples,
            warmup_steps=self.warmup_steps,
            num_chains=self.num_chains,
            random_state=self.random_state,
            num_drawn_samples=self.num_drawn_samples
        )

        # summary = {}
        # if "fp_lengthscale" in self._samples:
        #     fp = self._samples["fp_lengthscale"]
        #     summary["fp_mean"] = fp.mean(dim=0).cpu().numpy()
        #     summary["fp_std"] = fp.std(dim=0).cpu().numpy()

        # if "cont_lengthscale" in self._samples:
        #     cont = self._samples["cont_lengthscale"]
        #     summary["cont_mean"] = cont.mean(dim=0).cpu().numpy()
        #     summary["cont_std"] = cont.std(dim=0).cpu().numpy()

        # self.lengthscale_summary_ = summary

        self._is_fitted = True
        return self
    

    def predict(self, X_test):
        if isinstance(X_test, pd.DataFrame):
            X_test = X_test.to_numpy()

        X_t = torch.as_tensor(np.asarray(X_test), dtype=torch.float32)

        if X_t.ndim == 1:
            X_t = X_t.unsqueeze(0)

        if self.use_cuda and torch.cuda.is_available():
            X_t = X_t.cuda()

        mean_pred, std_pred = predict_posterior(
            self._X_train,
            self._y_train,
            X_t,
            self._samples,
            self._gp_model.kernel_builder
        )

        return mean_pred.detach().cpu().numpy()


