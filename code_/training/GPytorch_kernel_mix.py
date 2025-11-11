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
    


def _weighted_tanimoto_similarity(x1, x2):
    """
    Compute 1 - Tanimoto similarity between x1 and x2.
    Supports batched inputs.
    """
    x1 = x1.unsqueeze(-2)  # (..., N, 1, D)
    x2 = x2.unsqueeze(-3)  # (..., 1, M, D)

    numerator = torch.min(x1, x2).sum(dim=-1)
    denominator = torch.max(x1, x2).sum(dim=-1)
    tanimoto_sim = torch.where(denominator == 0, torch.zeros_like(denominator), numerator / denominator)
    tanimoto_dist = 1.0 - tanimoto_sim
    return tanimoto_dist



class TanimotoRBF(Kernel):
    has_lengthscale = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)


    def forward(self, x1, x2, diag=False, **params):
        """
        Standard RBF structure: exp(-0.5 * D^2 / lengthscale^2),
        but using Tanimoto distance instead of Euclidean.
        """
        if diag:
            return x1.new_ones(x1.shape[:-1])  # diagonal is 1

        # Normalize by lengthscale (same as RBFKernel)
        x1_ = x1.div(self.lengthscale)
        x2_ = x2.div(self.lengthscale)

        # Compute Tanimoto distance
        tanimoto_dist = _weighted_tanimoto_similarity(x1_, x2_)

        # Apply RBF exponentiation (same as postprocess_rbf)
        covar = postprocess_rbf(tanimoto_dist)

        return covar




class MixingKernel(Kernel):
    """
    Product of TanimotoRBF (for fingerprint features)
    and Matern(3/2) (for continuous features).
    Priors and hyperparameters can be passed from outside.
    """

    def __init__(
        self,
        feat_idx: dict,
        l_prior= None,
        sigma_prior= None,
        sigma_noise_prior= None,
        **kwargs,
    ):
        """
        feat_idx: dict with keys 'fp' and 'count', each a list of column indices.
                  Example: {'fp': [0,1,2,3], 'count': [4,5,6]}
        l_prior: optional prior for kernel lengthscales
        sigma_prior: optional prior for kernel outputscales
        sigma_noise_prior: optional prior for noise (stored for later use)
        """
        super().__init__(**kwargs)

        self.feat_idx = feat_idx
        self.sigma_noise_prior = sigma_noise_prior  # just stored for model use

        # --- Fingerprint kernel (TanimotoRBF) ---
        self.fp_kernel = ScaleKernel(
            TanimotoRBF(ard_num_dims=len(feat_idx.get('fp', []))),
            outputscale_prior=sigma_prior,
        )
        if l_prior is not None:
            self.fp_kernel.base_kernel.register_prior(
                "lengthscale_prior", l_prior, "lengthscale"
            )

        # --- Continuous kernel (Matern 3/2) ---
        self.cont_kernel = ScaleKernel(
            MaternKernel(nu=1.5, ard_num_dims=len(feat_idx.get('count', []))),
            outputscale_prior=sigma_prior,
        )
        if l_prior is not None:
            self.cont_kernel.base_kernel.register_prior(
                "lengthscale_prior", l_prior, "lengthscale"
            )

    def forward(self, x1, x2, diag=False, **params):
        """
        Compute product kernel between subsets of x1/x2.
        """
        fp_idx = self.feat_idx.get('fp', [])
        cont_idx = self.feat_idx.get('count', [])

        x1_fp = x1[..., fp_idx] if fp_idx else None
        x2_fp = x2[..., fp_idx] if fp_idx else None
        x1_cont = x1[..., cont_idx] if cont_idx else None
        x2_cont = x2[..., cont_idx] if cont_idx else None

        if x1_fp is not None and x1_cont is not None:
            return self.fp_kernel(x1_fp, x2_fp, diag=diag) * self.cont_kernel(x1_cont, x2_cont, diag=diag)
        elif x1_fp is not None:
            return self.fp_kernel(x1_fp, x2_fp, diag=diag)
        elif x1_cont is not None:
            return self.cont_kernel(x1_cont, x2_cont, diag=diag)
        else:
            raise ValueError("No valid feature indices found in feat_idx.")
        

class GPMix(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood, feat_idx):
        super(GPMix, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ZeroMean()
        self.covar_module = MixingKernel(feat_idx=feat_idx)

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)
        



class GPMixRegressor(BaseEstimator):
    def __init__(
        self,
        feat_idx,
        num_samples=100,
        warmup_steps=100,
        use_cuda=False,
        smoke_test=False,
    ):
        self.feat_idx = feat_idx
        self.num_samples = 2 if smoke_test else num_samples
        self.warmup_steps = 2 if smoke_test else warmup_steps
        self.use_cuda = use_cuda and torch.cuda.is_available()

    def fit(self, X_train, Y_train):
        # --- Convert input to torch tensors safely ---
        if isinstance(X_train, pd.DataFrame):
            X_train = X_train.to_numpy(dtype=np.float32)
        elif isinstance(X_train, np.ndarray):
            X_train = X_train.astype(np.float32)
        X_train = torch.tensor(X_train, dtype=torch.float)

        if isinstance(Y_train, (pd.Series, pd.DataFrame)):
            Y_train = Y_train.to_numpy(dtype=np.float32).ravel()
        elif isinstance(Y_train, np.ndarray):
            Y_train = Y_train.astype(np.float32).ravel()
        Y_train = torch.tensor(Y_train, dtype=torch.float)

        if self.use_cuda:
            X_train = X_train.cuda()
            Y_train = Y_train.cuda()
        # --- Priors (as you specified) ---
        l_prior = GammaPrior(5.0, 5.0)
        sigma_prior = NormalPrior(0.0, 1.0)
        sigma_noise_prior = NormalPrior(0.0, 1.0)

        # --- Likelihood and model ---
        self.likelihood = gpytorch.likelihoods.GaussianLikelihood()
        self.model = GPMix(X_train, Y_train, self.likelihood, self.feat_idx)

        # assign priors to the model components
        self.model.covar_module.fp_kernel.base_kernel.register_prior(
            "lengthscale_prior", l_prior, "lengthscale"
        )
        self.model.covar_module.cont_kernel.base_kernel.register_prior(
            "lengthscale_prior", l_prior, "lengthscale"
        )
        self.model.covar_module.fp_kernel.register_prior(
            "outputscale_prior", sigma_prior, "outputscale"
        )
        self.model.covar_module.cont_kernel.register_prior(
            "outputscale_prior", sigma_prior, "outputscale"
        )
        self.likelihood.register_prior("noise_prior", sigma_noise_prior, "noise")

        # --- Pyro MCMC model ---
        def pyro_model(x, y):
            with gpytorch.settings.fast_computations(False, False, False):
                sampled_model = self.model.pyro_sample_from_prior()
                output = sampled_model.likelihood(sampled_model(x))
                pyro.sample("obs", output, obs=y)
            return y

        nuts_kernel = NUTS(pyro_model)
        mcmc = MCMC(
            nuts_kernel,
            num_samples=self.num_samples,
            warmup_steps=self.warmup_steps,
            disable_progbar=True
        )

        print(f"Running NUTS with {self.num_samples} samples, {self.warmup_steps} warmup...")
        mcmc.run(X_train, Y_train)
        self.mcmc_samples = mcmc.get_samples()
        print("Sampling complete. Posterior hyperparameters stored.")

        self.train_x = X_train
        self.train_y = Y_train
        return self

    @torch.no_grad()
    def predict(self, X_test):
        if isinstance(X_test, (pd.DataFrame, np.ndarray)):
            X_test = torch.tensor(X_test, dtype=torch.float)

        if self.use_cuda:
            X_test = X_test.cuda()
            self.model = self.model.cuda()
            self.likelihood = self.likelihood.cuda()

        self.model.eval()
        self.likelihood.eval()

        preds = []
        for i in range(self.num_samples):
            sample_state = {k: v[i] for k, v in self.mcmc_samples.items()}
            self.model.load_state_dict(sample_state, strict=False)
            pred = self.likelihood(self.model(X_test)).mean.detach().cpu().numpy()
            preds.append(pred)

        mean_pred = np.mean(preds, axis=0)
        return mean_pred







n_samples = 100
n_fp = 8         # fingerprint features (binary)
n_cont = 3       # continuous features

# fingerprints (binary molecular-like descriptors)
fp_feats = np.random.randint(0, 2, size=(n_samples, n_fp))

# continuous features
cont_feats = np.random.randn(n_samples, n_cont)

# combine all features
X = np.concatenate([fp_feats, cont_feats], axis=1)

# true function (nonlinear in continuous features)
y_true = np.sin(cont_feats[:, 0]) + 0.3 * cont_feats[:, 1]**2 + (fp_feats[:, 0] * 0.5)
y = y_true + 0.1 * np.random.randn(n_samples)  # add noise

# pandas DataFrame for column-based indexing
columns = [f"fp{i}" for i in range(n_fp)] + [f"cont{i}" for i in range(n_cont)]
X_df = pd.DataFrame(X, columns=columns)

# feature index mapping for the kernel
feat_idx = {
    'fp': [X_df.columns.get_loc(c) for c in X_df.columns if "fp" in c],
    'count': [X_df.columns.get_loc(c) for c in X_df.columns if "cont" in c],
}

# split data
X_train, X_test, y_train, y_test = train_test_split(X_df, y, test_size=0.2, random_state=42)

# -------------------------------------------------------
# 2️⃣  Initialize and train the model
# -------------------------------------------------------

model = GPMixRegressor(
    feat_idx=feat_idx,
    num_samples=20,     # fewer samples for demo speed
    warmup_steps=20,
    use_cuda=False,     # set True if GPU available
)

print("Fitting model (this will run NUTS sampling)...")
model.fit(X_train, y_train)

# -------------------------------------------------------
# 3️⃣  Predict on test set
# -------------------------------------------------------

print("Predicting...")
y_pred = model.predict(X_test)

# -------------------------------------------------------
# 4️⃣  Evaluate and visualize
# -------------------------------------------------------
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

r2 = r2_score(y_test, y_pred)
print(f"R² on test data: {r2:.3f}")

plt.figure(figsize=(6, 6))
plt.scatter(y_test, y_pred, alpha=0.7)
plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'k--', lw=2)
plt.xlabel("True y")
plt.ylabel("Predicted y")
plt.title(f"GP-Mix Regressor (R² = {r2:.3f})")
plt.show()

