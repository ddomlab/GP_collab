import gpytorch
import numpy as np
import pandas as pd
import torch
# import torch.nn as nn
from sklearn.base import BaseEstimator
from sklearn.model_selection import train_test_split
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

def weighted_batch_tanimoto_similarity(
        x1: torch.Tensor, x2: torch.Tensor, eps: float = 1e-6
) -> torch.Tensor:
    """
    Weighted Tanimoto similarity (generalized Jaccard similarity) between two batched tensors across the last 2 dimensions.
    eps argument ensures numerical stability to avoid division by zero.
    """
    assert x1.ndim >= 2 and x2.ndim >= 2, "Input tensors must have at least 2 dimensions."

    pairwise_min = torch.min(x1.unsqueeze(-2), x2.unsqueeze(-3))  
    min_sum = torch.sum(pairwise_min, dim=-1, keepdims=True)  

    pairwise_max = torch.max(x1.unsqueeze(-2), x2.unsqueeze(-3))  
    max_sum = torch.sum(pairwise_max, dim=-1, keepdims=True)  

    # Weighted Tanimoto similarity
    similarity = (min_sum + eps) / (max_sum + eps)
    return similarity
    


class TanimotoDistance(torch.nn.Module):
    r"""
    Distance module for bit vector test_kernels.
    """

    def __init__(self, postprocess_script=lambda x: x):
        super().__init__()
        self._postprocess = postprocess_script

    def _sim(self, x1, x2, postprocess, x1_eq_x2=False, metric="binary tanimoto"):
        r"""
        Computes the similarity between x1 and x2
        Args:
            :attr: `x1`: (Tensor `n x d` or `b x n x d`):
                First set of data where b is a batch dimension
            :attr: `x2`: (Tensor `m x d` or `b x m x d`):
                Second set of data where b is a batch dimension
            :attr: `postprocess` (bool):
                Whether to apply a postprocess script (default is none)
            :attr: `x1_eq_x2` (bool):
                Is x1 equal to x2
            :attr: `metric` (str):
                String specifying the similarity metric. One of ['tanimoto']
        Returns:
            (:class:`Tensor`, :class:`Tensor) corresponding to the similarity matrix between `x1` and `x2`
        """

        # Branch for Tanimoto metric
        if metric == "binary tanimoto":
            res = binary_batch_tanimoto_sim(x1, x2)
            res.clamp_min_(0)  # zero out negative values
            return self._postprocess(res) if postprocess else res
        if metric == "weighted tanimoto":
            res = weighted_batch_tanimoto_similarity(x1, x2)
            res.clamp_min_(0)  # zero out negative values
            return self._postprocess(res) if postprocess else res
        else:
            raise RuntimeError(
                "Similarity metric not supported. Available options are 'tanimoto'"
            )


class TanimotoKernel(gpytorch.kernels.Kernel):
    ''' Tanimoto kernel from GAUCHE
    https://github.com/leojklarner/gauche/blob/main/gauche/kernels/fingerprint_kernels/tanimoto_kernel.py
    '''

    def __init__(self, metric="tanimoto", **kwargs):
        super(TanimotoKernel, self).__init__(**kwargs)
        self.metric = metric

    def covar_dist(
            self,
            x1,
            x2,
            last_dim_is_batch=False,
            dist_postprocess_func=lambda x: x,
            postprocess=True,
            **params,
    ):
        r"""
        This is a helper method for computing the bit vector similarity between
        all pairs of points in x1 and x2.
        Args:
            :attr:`x1` (Tensor `n x d` or `b1 x ... x bk x n x d`):
                First set of data.
            :attr:`x2` (Tensor `m x d` or `b1 x ... x bk x m x d`):
                Second set of data.
            :attr:`last_dim_is_batch` (tuple, optional):
                Is the last dimension of the data a batch dimension or not?
        Returns:
            (:class:`Tensor`, :class:`Tensor) corresponding to the distance matrix between `x1` and `x2`.
            The shape depends on the kernel's mode
            * `diag=False`
            * `diag=False` and `last_dim_is_batch=True`: (`b x d x n x n`)
            * `diag=True`
            * `diag=True` and `last_dim_is_batch=True`: (`b x d x n`)
        """
        if last_dim_is_batch:
            x1 = x1.transpose(-1, -2).unsqueeze(-1)
            x2 = x2.transpose(-1, -2).unsqueeze(-1)

        x1_eq_x2 = torch.equal(x1, x2)

        # torch scripts expect tensors
        postprocess = torch.tensor(postprocess)

        res = None

        # Cache the Distance object or else JIT will recompile every time
        if (
                not self.distance_module
                or self.distance_module._postprocess != dist_postprocess_func
        ):
            self.distance_module = TanimotoDistance(dist_postprocess_func)

        res = self.distance_module._sim(
            x1, x2, postprocess, x1_eq_x2, self.metric
        )

        return res

    def forward(self, x1, x2, diag=False, **params):
        if diag:
            assert x1.size() == x2.size() and torch.equal(x1, x2)
            return torch.ones(
                *x1.shape[:-2], x1.shape[-2], dtype=x1.dtype, device=x1.device
            )
        else:
            return self.covar_dist(x1, x2, **params)


class GP(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood, **kwargs):
        super(GP, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ZeroMean()
        # print(kwargs['kernel'])
        if kwargs['kernel'] == 'rbf':
            # for numerical
            self.covar_module = gpytorch.kernels.ScaleKernel(
                gpytorch.kernels.RBFKernel(ard_num_dims=train_x.shape[-1],**kwargs)
            )
            # self.covar_module.base_kernel.lengthscale = kwargs['lengthscale']
        elif kwargs['kernel'] == 'tanimoto':
            # for ECFP
            self.covar_module = gpytorch.kernels.ScaleKernel(TanimotoKernel())
        elif kwargs['kernel'] == 'RQ':
            # for numerical 
            self.covar_module = gpytorch.kernels.ScaleKernel(
                gpytorch.kernels.RQKernel(ard_num_dims=train_x.shape[-1])
                )
        elif kwargs['kernel'] == 'matern':
            # for numeical
            # nu = kwargs.get('nu', 2.5)
            self.covar_module = gpytorch.kernels.ScaleKernel(
                gpytorch.kernels.MaternKernel(ard_num_dims=train_x.shape[-1],**kwargs)
                )
        else:
            raise ValueError('Invalid kernel')

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)


class GPRegressor(BaseEstimator):

    def __init__(
            self,
            kernel,
            # length_scale=None,
            lr=1e-2,
            n_epoch=100,
            lengthscale=1,
            nu=2.5,

    ):
        self.ll = gpytorch.likelihoods.GaussianLikelihood()
        self.kernel = kernel
        self.lr = lr
        self.n_epoch = n_epoch
        self.lengthscale = lengthscale 
        self.nu = nu



    def fit(self, X_train, Y_train):

        if isinstance(X_train, pd.DataFrame):
            X_train = X_train.to_numpy()

        X_train = torch.tensor(X_train, dtype=torch.float)
        Y_train = torch.tensor(Y_train.ravel(), dtype=torch.float)
        self.model = GP(X_train, Y_train.ravel(),
                        self.ll, 
                        kernel=self.kernel,
                        lengthscale=self.lengthscale,
                        nu = self.nu)
        # gpytorch.settings.cholesky_jitter(1e-3)               
        # train return loss (minimize)
        optimizer = torch.optim.Adam(self.model.parameters(), lr=self.lr)
        mll = gpytorch.mlls.ExactMarginalLogLikelihood(self.ll, self.model)

        if torch.cuda.is_available():
            X_train = X_train.cuda()
            Y_train = Y_train.cuda()
            self.model = self.model.cuda()
            mll = mll.cuda()

        self.model.train()
        self.ll.train()
        for _ in range(self.n_epoch):
            optimizer.zero_grad()
            y_pred = self.model(X_train)
            loss = -mll(y_pred, Y_train.ravel())
            # print(f'LOSS: {loss.item()}', end='\r')
            loss.backward()
            optimizer.step()

    def predict(self, X_test):
        if isinstance(X_test, pd.DataFrame):
            X_test = X_test.to_numpy()

        X_test = torch.tensor(X_test, dtype=torch.float)

        if torch.cuda.is_available():
            X_test = X_test.cuda()
            self.model = self.model.cuda()
            self.ll = self.ll.cuda()

        self.model.eval()
        self.ll.eval()
        with torch.no_grad():
            y_pred = self.ll(self.model(X_test)).mean.cpu().numpy()

        return y_pred



