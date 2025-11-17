import pyro
import pyro.distributions as dist
import pyro.contrib.gp.kernels as pk
import torch
import pyro
import pyro.distributions as dist
import pyro.contrib.gp as gp
from pyro.infer import MCMC, NUTS
from sklearn.base import BaseEstimator, RegressorMixin
import pandas as pd
import numpy as np

def weighted_tanimoto_distance(x1, x2, eps=1e-6):
    x1e = x1.unsqueeze(-2)
    x2e = x2.unsqueeze(-3)

    numerator = torch.min(x1e, x2e).sum(dim=-1)
    denominator = torch.max(x1e, x2e).sum(dim=-1)

    sim = (numerator + eps) / (denominator + eps)
    dist = 1.0 - sim
    return torch.clamp(dist, min=0.)



class TanimotoRBF(pk.Kernel):
    def __init__(self, input_dim, lengthscale=None, eps=1e-6, active_dims=None):
        super().__init__(input_dim=input_dim, active_dims=active_dims)
        self.eps = eps
        
        if lengthscale is None:
            # ARD: one lengthscale per dimension
            lengthscale = torch.ones(input_dim)
        self.lengthscale = torch.nn.Parameter(lengthscale)

    def forward(self, X, Z=None, diag=False):
        X = self._slice_input(X)
        Z = self._slice_input(Z) if Z is not None else None
        
        if Z is None:
            Z = X

        if diag is True:
            return X.new_ones(X.size(0))

        # Compute Tanimoto distance matrix
        D = weighted_tanimoto_distance(X, Z, eps=self.eps)

        # RBF on that distance
        ls = self.lengthscale.mean()  # Pyro kernels use scalar LS; match your code
        K = torch.exp(-0.5 * (D / ls)**2)
        return K
    

class MixingKernelPyro:
    """
    Builds a Pyro kernel that matches your GPyTorch MixingKernel
    """
    def __init__(self, feat_idx):
        self.fp_idx = feat_idx.get("fp") or []
        self.cont_idx = feat_idx.get("count") or []

    def build(self):
        kernels = []

        # FP block
        if len(self.fp_idx) > 0:
            k_fp = TanimotoRBF(input_dim=len(self.fp_idx))
            kernels.append(k_fp)

        # Continuous block
        if len(self.cont_idx) > 0:
            k_cont = pk.Matern32(
                input_dim=len(self.cont_idx),
                lengthscale=torch.ones(len(self.cont_idx))  # ARD for cont dims
            )
            kernels.append(k_cont)

        if not kernels:
            raise ValueError("Both FP and CONT feature lists are empty")

        # Product of kernels
        k = kernels[0]
        for other in kernels[1:]:
            k = pk.Product(k, other)

        return k


class GPMixPyro:
    def __init__(self, X, y, feat_idx):
        self.X = X
        self.y = y
        self.kernel_builder = MixingKernelPyro(feat_idx)

    def model(self):
        # sample hyperparameters
        outputscale = pyro.sample("outputscale", dist.LogNormal(0.0, 1.0))
        obs_noise = pyro.sample("obs_noise", dist.LogNormal(0.0, 1.0))

        # build a completely new kernel each time
        kernel = self.kernel_builder.build()
        kernel.variance = outputscale

        # instantiate GPRegression fresh each time
        gpmodel = gp.models.GPRegression(
            self.X,
            self.y,
            kernel,
            noise=obs_noise
        )

        gpmodel.model()




def run_inference(gpmodel):
    nuts_kernel = NUTS(gpmodel.model)
    mcmc = MCMC(
        nuts_kernel,
        num_samples=300,
        warmup_steps=200,
        num_chains=1
    )
    mcmc.run()
    return mcmc.get_samples()



def predict_posterior(X_train, y_train, X_test, samples, kernel_builder):
    K_list = []
    mean_list = []

    for i in range(samples["outputscale"].shape[0]):
        noise = samples["noise"][i]
        oscale = samples["outputscale"][i]

        # build kernel with sampled hyperparameters
        kernel = kernel_builder.build()
        kernel = pk.Scale(kernel, outputscale=oscale)

        # GP regression object
        gpmodel = gp.models.GPRegression(
            X_train, y_train, kernel, noise=noise
        )

        # manually compute predictive mean and covariance
        posterior = gpmodel(X_test, full_cov=False)
        mean_list.append(posterior.mean)
        K_list.append(posterior.variance)

    mean_stack = torch.stack(mean_list)
    var_stack = torch.stack(K_list)

    mean_pred = mean_stack.mean(dim=0)
    std_pred = var_stack.mean(dim=0).sqrt()

    return mean_pred, std_pred


class GPMixMCMCRegressor(BaseEstimator, RegressorMixin):
    def __init__(
        self,
        feat_idx=None,
        num_samples=300,
        warmup_steps=200,
        num_chains=1,
        use_cuda=False,
        random_state=None,
    ):
        # store all hyperparameters exactly as provided
        self.feat_idx = feat_idx
        self.num_samples = num_samples
        self.warmup_steps = warmup_steps
        self.num_chains = num_chains
        self.use_cuda = use_cuda
        self.random_state = random_state

    def fit(self, X_train, y_train):

        # convert DataFrame to numpy
        if isinstance(X_train, pd.DataFrame):
            X_train = X_train.to_numpy()
        if isinstance(y_train, pd.DataFrame):
            y_train = y_train.to_numpy()

        # convert to tensors
        X_t = torch.as_tensor(np.asarray(X_train), dtype=torch.float32)
        y_t = torch.as_tensor(np.asarray(y_train), dtype=torch.float32).view(-1)



        self._X_train = X_t
        self._y_train = y_t

        self._gp_model = GPMixPyro(X_t, y_t, self.feat_idx)

        self._samples = run_inference(
            self._gp_model,
            # num_samples=self.num_samples,
            # warmup_steps=self.warmup_steps,
            # num_chains=self.num_chains,
            # random_state=self.random_state,
        )

        self._is_fitted = True
        return self

    def predict(self, X_test):

        if isinstance(X_test, pd.DataFrame):
            X_test = X_test.to_numpy()

        X_t = torch.as_tensor(np.asarray(X_test), dtype=torch.float32)


        mean_pred, std_pred = predict_posterior(
            self._X_train,
            self._y_train,
            X_t,
            self._samples,
            self._gp_model.kernel_builder,
        )

        return mean_pred.detach().cpu().numpy()

