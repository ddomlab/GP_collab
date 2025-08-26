import math
from sklearn.gaussian_process.kernels import (
                                            Kernel, 
                                            StationaryKernelMixin, 
                                            NormalizedKernelMixin, 
                                            Hyperparameter,Matern,
                                            _approx_fprime,
                                            Matern,
                                            RBF,
                                            _check_length_scale)
import numpy as np
from scipy.spatial.distance import pdist, cdist, squareform
from sklearn.utils import check_array
from scipy.special import kv, gamma
from math import sqrt



def weighted_jaccard(u, v, eps: float = 1e-6):
    """Weighted Jaccard distance between two vectors u, v."""
    min_sum = np.sum(np.minimum(u, v))
    max_sum = np.sum(np.maximum(u, v))
    return 1.0 - ((min_sum + eps) / (max_sum + eps))


class JaccardKernel(NormalizedKernelMixin, Kernel):
    """Pure Jaccard similarity kernel (no RBF transform)."""

    def __init__(self):
        pass  # no hyperparameters here

    def __call__(self, X, Y=None, eval_gradient=False):
        X = np.atleast_2d(X)
        if Y is None:
            dists = pdist(X, metric=weighted_jaccard)
            sims = 1.0 - squareform(dists)  # similarity = 1 - distance
            np.fill_diagonal(sims, 1.0)
        else:
            if eval_gradient:
                raise ValueError("Gradient can only be evaluated when Y is None.")
            Y = np.atleast_2d(Y)
            dists = cdist(X, Y, metric=weighted_jaccard)
            sims = 1.0 - dists

        if eval_gradient:
            # no hyperparameters → return empty gradient
            return sims, np.empty((X.shape[0], X.shape[0], 0))
        else:
            return sims

    def diag(self, X):
        # all self-similarities are 1
        return np.ones(X.shape[0])

    def is_stationary(self):
        # Jaccard is not stationary in the GP sense
        return False


class Jaccard_RBF(StationaryKernelMixin, NormalizedKernelMixin, Kernel):

    def __init__(self, length_scale=1.0, length_scale_bounds=(1e-5, 1e5)):
        self.length_scale = length_scale
        self.length_scale_bounds = length_scale_bounds

    @property
    def anisotropic(self):
        return np.iterable(self.length_scale) and len(self.length_scale) > 1

    @property
    def hyperparameter_length_scale(self):
        if self.anisotropic:
            return Hyperparameter(
                "length_scale",
                "numeric",
                self.length_scale_bounds,
                len(self.length_scale),
            )
        return Hyperparameter("length_scale", "numeric", self.length_scale_bounds)

    def __call__(self, X, Y=None, eval_gradient=False):
        X = np.atleast_2d(X)
        length_scale = _check_length_scale(X, self.length_scale)
        if Y is None:
            dists = pdist(X / length_scale, metric=weighted_jaccard)
            K = np.exp(-0.5 * dists)
            # convert from upper-triangular matrix to square matrix
            K = squareform(K)
            np.fill_diagonal(K, 1)
        else:
            if eval_gradient:
                raise ValueError("Gradient can only be evaluated when Y is None.")
            dists = cdist(X / length_scale, Y / length_scale, metric=weighted_jaccard)
            K = np.exp(-0.5 * dists)

        if eval_gradient:
            if self.hyperparameter_length_scale.fixed:
                # Hyperparameter l kept fixed
                return K, np.empty((X.shape[0], X.shape[0], 0))
            elif not self.anisotropic or length_scale.shape[0] == 1:
                K_gradient = (K * squareform(dists))[:, :, np.newaxis]
                return K, K_gradient
            elif self.anisotropic:
                # We need to recompute the pairwise dimension-wise distances
                K_gradient = (X[:, np.newaxis, :] - X[np.newaxis, :, :]) ** 2 / (
                    length_scale**2
                )
                K_gradient *= K[..., np.newaxis]
                return K, K_gradient
        else:
            return K

    def __repr__(self):
        if self.anisotropic:
            return "{0}(length_scale=[{1}])".format(
                self.__class__.__name__,
                ", ".join(map("{0:.3g}".format, self.length_scale)),
            )
        else:  # isotropic
            return "{0}(length_scale={1:.3g})".format(
                self.__class__.__name__, np.ravel(self.length_scale)[0]
            )



class Jaccard_Matern(RBF):

    def __init__(self, length_scale=1.0, length_scale_bounds=(1e-5, 1e5), nu=1.5):
        super().__init__(length_scale, length_scale_bounds)
        self.nu = nu

    def __call__(self, X, Y=None, eval_gradient=False):
        X = np.atleast_2d(X)
        length_scale = _check_length_scale(X, self.length_scale)
        if Y is None:
            dists = pdist(X / length_scale, metric=weighted_jaccard)
        else:
            if eval_gradient:
                raise ValueError("Gradient can only be evaluated when Y is None.")
            dists = cdist(X / length_scale, Y / length_scale, metric=weighted_jaccard)

        if self.nu == 0.5:
            K = np.exp(-dists)
        elif self.nu == 1.5:
            K = dists * math.sqrt(3)
            K = (1.0 + K) * np.exp(-K)
        elif self.nu == 2.5:
            K = dists * math.sqrt(5)
            K = (1.0 + K + K**2 / 3.0) * np.exp(-K)
        elif self.nu == np.inf:
            K = np.exp(-(dists**2) / 2.0)
        else:  # general case; expensive to evaluate
            K = dists
            K[K == 0.0] += np.finfo(float).eps  # strict zeros result in nan
            tmp = math.sqrt(2 * self.nu) * K
            K.fill((2 ** (1.0 - self.nu)) / gamma(self.nu))
            K *= tmp**self.nu
            K *= kv(self.nu, tmp)

        if Y is None:
            # convert from upper-triangular matrix to square matrix
            K = squareform(K)
            np.fill_diagonal(K, 1)

        if eval_gradient:
            if self.hyperparameter_length_scale.fixed:
                # Hyperparameter l kept fixed
                K_gradient = np.empty((X.shape[0], X.shape[0], 0))
                return K, K_gradient

            # We need to recompute the pairwise dimension-wise distances
            if self.anisotropic:
                D = (X[:, np.newaxis, :] - X[np.newaxis, :, :]) ** 2 / (length_scale**2)
            else:
                D = squareform(dists**2)[:, :, np.newaxis]

            if self.nu == 0.5:
                denominator = np.sqrt(D.sum(axis=2))[:, :, np.newaxis]
                divide_result = np.zeros_like(D)
                np.divide(
                    D,
                    denominator,
                    out=divide_result,
                    where=denominator != 0,
                )
                K_gradient = K[..., np.newaxis] * divide_result
            elif self.nu == 1.5:
                K_gradient = 3 * D * np.exp(-np.sqrt(3 * D.sum(-1)))[..., np.newaxis]
            elif self.nu == 2.5:
                tmp = np.sqrt(5 * D.sum(-1))[..., np.newaxis]
                K_gradient = 5.0 / 3.0 * D * (tmp + 1) * np.exp(-tmp)
            elif self.nu == np.inf:
                K_gradient = D * K[..., np.newaxis]
            else:
                # approximate gradient numerically
                def f(theta):  # helper function
                    return self.clone_with_theta(theta)(X, Y)

                return K, _approx_fprime(self.theta, f, 1e-10)

            if not self.anisotropic:
                return K, K_gradient[:, :].sum(-1)[:, :, np.newaxis]
            else:
                return K, K_gradient
        else:
            return K

    def __repr__(self):
        if self.anisotropic:
            return "{0}(length_scale=[{1}], nu={2:.3g})".format(
                self.__class__.__name__,
                ", ".join(map("{0:.3g}".format, self.length_scale)),
                self.nu,
            )
        else:
            return "{0}(length_scale={1:.3g}, nu={2:.3g})".format(
                self.__class__.__name__, np.ravel(self.length_scale)[0], self.nu
            )


class AdditiveRBF(Kernel):
    """Additive or average RBF kernel across input dimensions."""

    def __init__(self, length_scale=1.0, length_scale_bounds=(1e-5, 1e5), mode="sum"):
        self.length_scale = length_scale
        self.length_scale_bounds = length_scale_bounds
        if mode not in ("sum", "avg"):
            raise ValueError("mode must be 'sum' or 'avg'")
        self.mode = mode

    @property
    def hyperparameter_length_scale(self):
        if np.iterable(self.length_scale) and np.ndim(self.length_scale) > 0:
            return Hyperparameter("length_scale", "numeric",
                                  self.length_scale_bounds,
                                  len(self.length_scale))
        return Hyperparameter("length_scale", "numeric",
                              self.length_scale_bounds)

    def __call__(self, X, Y=None, eval_gradient=False):
        X = check_array(X)
        length_scale = _check_length_scale(X, self.length_scale)

        if Y is None:
            Y = X
            symmetric = True
        else:
            Y = check_array(Y)
            symmetric = False

        n_features = X.shape[1]
        K = np.zeros((X.shape[0], Y.shape[0]))
        if eval_gradient:
            if not symmetric:
                raise ValueError("Gradient only supported for Y=None.")
            if self.hyperparameter_length_scale.fixed:
                return K, np.empty((X.shape[0], X.shape[0], 0))
            K_gradient = np.zeros((X.shape[0], Y.shape[0], n_features))

        for d in range(n_features):
            diff = (X[:, [d]] - Y[:, [d]].T) / length_scale[d]
            K_d = np.exp(-0.5 * diff**2)
            K += K_d

            if eval_gradient:
                grad = K_d * (diff**2)
                K_gradient[:, :, d] = grad

        if self.mode == "avg":
            K /= n_features
            if eval_gradient:
                K_gradient /= n_features

        if eval_gradient:
            return K, K_gradient
        return K

    def diag(self, X):
        n_features = X.shape[1]
        return np.full(X.shape[0], 1.0 if self.mode=="avg" else n_features)

    def is_stationary(self):
        return False

    def __repr__(self):
        return f"AdditiveRBF(length_scale={self.length_scale}, mode={self.mode})"


import numpy as np
from scipy.special import kv, gamma
from sklearn.gaussian_process.kernels import Kernel, Hyperparameter
from sklearn.utils.validation import check_array
from sklearn.gaussian_process.kernels import _check_length_scale


class AdditiveMatern(Matern):
    def __init__(self, length_scale=1.0, length_scale_bounds=(1e-5, 1e5), nu=1.5):
        super().__init__(length_scale=length_scale, length_scale_bounds=length_scale_bounds, nu=nu)

    def __call__(self, X, Y=None, eval_gradient=False):
        X = check_array(X)
        length_scale = _check_length_scale(X, self.length_scale)
        if Y is None:
            Y = X
            symmetric = True
        else:
            Y = check_array(Y)
            symmetric = False

        K = np.zeros((X.shape[0], Y.shape[0]))
        if eval_gradient:
            if not symmetric:
                raise ValueError("Gradient can only be evaluated when Y is None.")
            if self.hyperparameter_length_scale.fixed:
                return K, np.empty((X.shape[0], X.shape[0], 0))
            K_gradient = np.zeros((X.shape[0], Y.shape[0], X.shape[1]))

        for d in range(X.shape[1]):
            X_d = X[:, [d]] / length_scale[d]
            Y_d = Y[:, [d]] / length_scale[d]
            dists = np.abs(X_d - Y_d.T)

            if self.nu == 0.5:
                K_d = np.exp(-dists)
            elif self.nu == 1.5:
                sqrt3 = np.sqrt(3.0)
                K_d = (1.0 + sqrt3 * dists) * np.exp(-sqrt3 * dists)
            elif self.nu == 2.5:
                sqrt5 = np.sqrt(5.0)
                K_d = (1.0 + sqrt5 * dists + 5.0 / 3.0 * dists**2) * np.exp(-sqrt5 * dists)
            elif self.nu == np.inf:
                K_d = np.exp(-0.5 * dists**2)
            else:
                eps = 1e-6
                dists[dists == 0.0] += eps  # avoid nan in kv
                tmp = np.sqrt(2 * self.nu) * dists
                coef = (2 ** (1.0 - self.nu)) / gamma(self.nu)
                K_d = coef * (tmp ** self.nu) * kv(self.nu, tmp)
                K_d = np.where(np.isnan(K_d), 1.0, K_d)  # replace nan (from kv(ν, 0)) with 1.0

            K += K_d

            if eval_gradient:
                if self.nu == 1.5:
                    sqrt3 = np.sqrt(3.0)
                    grad = 3.0 * (X[:, [d]] - Y[:, [d]].T) ** 2 / length_scale[d]**3
                    K_gradient[:, :, d] = grad * np.exp(-sqrt3 * dists)
                elif self.nu == 2.5:
                    sqrt5 = np.sqrt(5.0)
                    grad = 5.0 * (X[:, [d]] - Y[:, [d]].T) ** 2 / (3 * length_scale[d]**3)
                    tmp = sqrt5 * dists
                    K_gradient[:, :, d] = grad * (tmp + 1) * np.exp(-tmp)
                elif self.nu == np.inf:
                    grad = (X[:, [d]] - Y[:, [d]].T) ** 2 / (length_scale[d]**3)
                    K_gradient[:, :, d] = K_d * grad
                else:
                    # general case not supported analytically
                    return K, _approx_fprime(self.theta, lambda theta: self.clone_with_theta(theta)(X), 1e-10)

        # -------- only change: take average instead of sum --------
        K /= X.shape[1]
        if eval_gradient:
            K_gradient /= X.shape[1]

        if Y is None:
            np.fill_diagonal(K, 1.0)  # average → diag entries should be 1

        if eval_gradient:
            return K, K_gradient
        return K

    def diag(self, X):
        X = check_array(X)
        return np.full(X.shape[0], 1.0)  # average mode → diagonal = 1

    def is_stationary(self):
        return False  # Sum/avg of 1D kernels is not stationary

    def __repr__(self):
        if self.anisotropic:
            return "{0}(length_scale=[{1}], nu={2:.3g})".format(
                self.__class__.__name__,
                ", ".join(map("{0:.3g}".format, self.length_scale)),
                self.nu,
            )
        else:
            return "{0}(length_scale={1:.3g}, nu={2:.3g})".format(
                self.__class__.__name__, np.ravel(self.length_scale)[0], self.nu
            )





import numpy as np
from sklearn.gaussian_process.kernels import Kernel

class ProductKernel(Kernel):
    """Product kernel = K_cont * K_fp on disjoint feature groups.

    Parameters
    ----------
    kernel_cont : Kernel
        Kernel to apply on continuous features.
    kernel_fp : Kernel
        Kernel to apply on fingerprint features.
    cont_idx : array-like of int
        Column indices of continuous features.
    fp_idx : array-like of int
        Column indices of fingerprint features.
    """

    def __init__(self, kernel_cont, kernel_fp, cont_idx, fp_idx):
        self.kernel_cont = kernel_cont
        self.kernel_fp = kernel_fp
        self.cont_idx = cont_idx
        self.fp_idx = fp_idx

    # --- sklearn plumbing: params, hyperparams, theta, bounds ---
    def get_params(self, deep=True):
        return dict(
            kernel_cont=self.kernel_cont,
            kernel_fp=self.kernel_fp,
            cont_idx=self.cont_idx,
            fp_idx=self.fp_idx,
        )

    @property
    def hyperparameters(self):
        return [*self.kernel_cont.hyperparameters,
                *self.kernel_fp.hyperparameters]

    @property
    def theta(self):
        return np.hstack([self.kernel_cont.theta, self.kernel_fp.theta])

    @theta.setter
    def theta(self, theta):
        n1 = self.kernel_cont.theta.shape[0]
        self.kernel_cont.theta = theta[:n1]
        self.kernel_fp.theta = theta[n1:]

    @property
    def bounds(self):
        b1, b2 = self.kernel_cont.bounds, self.kernel_fp.bounds
        if b1.shape[0] == 0 and b2.shape[0] == 0:
            return np.empty((0, 2))
        if b1.shape[0] == 0:
            return b2
        if b2.shape[0] == 0:
            return b1
        return np.vstack([b1, b2])

    # --- main kernel logic ---
    def __call__(self, X, Y=None, eval_gradient=False):
        X_cont, X_fp = X[:, self.cont_idx], X[:, self.fp_idx]
        Y_cont = None if Y is None else Y[:, self.cont_idx]
        Y_fp   = None if Y is None else Y[:, self.fp_idx]

        if eval_gradient:
            Kc, grad_c = self.kernel_cont(X_cont, Y_cont, eval_gradient=True)
            Kf, grad_f = self.kernel_fp(X_fp, Y_fp, eval_gradient=True)
            K = Kc * Kf

            grads = []
            if grad_c.shape[2] > 0:
                grads.append(grad_c * Kf[:, :, np.newaxis])
            if grad_f.shape[2] > 0:
                grads.append(grad_f * Kc[:, :, np.newaxis])

            grad = (np.concatenate(grads, axis=2)
                    if grads else np.empty((X.shape[0], X.shape[0], 0)))
            return K, grad
        else:
            return self.kernel_cont(X_cont, Y_cont) * self.kernel_fp(X_fp, Y_fp)

    def diag(self, X):
        return (self.kernel_cont.diag(X[:, self.cont_idx])
                * self.kernel_fp.diag(X[:, self.fp_idx]))

    def is_stationary(self):
        return self.kernel_cont.is_stationary() and self.kernel_fp.is_stationary()

    # convenience forwards
    @property
    def length_scale_cont(self):
        return getattr(self.kernel_cont, "length_scale", None)

    @property
    def length_scale_fp(self):
        return getattr(self.kernel_fp, "length_scale", None)

    def __repr__(self):
        return (f"FeatureGroupProductKernel("
                f"{self.kernel_cont}, {self.kernel_fp}, "
                f"cont_idx={self.cont_idx}, fp_idx={self.fp_idx})")




class AdditiveKernel(Kernel):
    def __init__(self, kernel_cont, kernel_fp, cont_idx, fp_idx):
        self.kernel_cont = kernel_cont
        self.kernel_fp = kernel_fp
        self.cont_idx = cont_idx
        self.fp_idx = fp_idx

    # --- sklearn plumbing: params, hyperparams, theta, bounds ---
    def get_params(self, deep=True):
        return dict(
            kernel_cont=self.kernel_cont,
            kernel_fp=self.kernel_fp,
            cont_idx=self.cont_idx,
            fp_idx=self.fp_idx,
        )

    @property
    def hyperparameters(self):
        return [*self.kernel_cont.hyperparameters,
                *self.kernel_fp.hyperparameters]

    @property
    def theta(self):
        return np.hstack([self.kernel_cont.theta, self.kernel_fp.theta])

    @theta.setter
    def theta(self, theta):
        n1 = self.kernel_cont.theta.shape[0]
        self.kernel_cont.theta = theta[:n1]
        self.kernel_fp.theta = theta[n1:]

    @property
    def bounds(self):
        b1, b2 = self.kernel_cont.bounds, self.kernel_fp.bounds
        if b1.shape[0] == 0 and b2.shape[0] == 0:
            return np.empty((0, 2))
        if b1.shape[0] == 0:
            return b2
        if b2.shape[0] == 0:
            return b1
        return np.vstack([b1, b2])

    # --- main kernel logic ---
    def __call__(self, X, Y=None, eval_gradient=False):
        X_cont, X_fp = X[:, self.cont_idx], X[:, self.fp_idx]
        Y_cont = None if Y is None else Y[:, self.cont_idx]
        Y_fp   = None if Y is None else Y[:, self.fp_idx]

        if eval_gradient:
            Kc, grad_c = self.kernel_cont(X_cont, Y_cont, eval_gradient=True)
            Kf, grad_f = self.kernel_fp(X_fp, Y_fp, eval_gradient=True)
            K = Kc + Kf

            grads = []
            if grad_c.shape[2] > 0:
                grads.append(grad_c)
            if grad_f.shape[2] > 0:
                grads.append(grad_f)

            grad = (np.concatenate(grads, axis=2)
                    if grads else np.empty((X.shape[0], X.shape[0], 0)))
            return K, grad
        else:
            return self.kernel_cont(X_cont, Y_cont) + self.kernel_fp(X_fp, Y_fp)

    def diag(self, X):
        return (self.kernel_cont.diag(X[:, self.cont_idx])
                + self.kernel_fp.diag(X[:, self.fp_idx]))

    def is_stationary(self):
        return self.kernel_cont.is_stationary() and self.kernel_fp.is_stationary()

    # convenience forwards
    @property
    def length_scale_cont(self):
        return getattr(self.kernel_cont, "length_scale", None)

    @property
    def length_scale_fp(self):
        return getattr(self.kernel_fp, "length_scale", None)

    def __repr__(self):
        return (f"FeatureGroupAdditiveKernel("
                f"{self.kernel_cont}, {self.kernel_fp}, "
                f"cont_idx={self.cont_idx}, fp_idx={self.fp_idx})")



from scipy.optimize import minimize


def map_optimizer(obj_func, initial_theta, bounds):
    def objective_with_prior(theta):
        nll, grad = obj_func(theta)

        # Example: Gaussian prior on log-length-scales
        log_prior = 0.5 * np.sum(theta**2)   # -log p(theta)
        grad_prior = theta                   # ∂/∂theta
        return nll + log_prior, grad + grad_prior

    res = minimize(
        fun=lambda th: objective_with_prior(th),
        x0=initial_theta,
        jac=True,
        bounds=bounds,
        method="L-BFGS-B"
    )
    # sklearn expects (theta_opt, func_min)
    return res.x, res.fun



import numpy as np
from sklearn.gaussian_process.kernels import Kernel, Hyperparameter
from scipy.special import kv, gamma

class AdditiveContinuousKernel(Kernel):
    """Additive kernel over continuous features with choice of base kernel.

    Parameters
    ----------
    kernel_type : str, {"rbf", "matern"}
        Type of base kernel to use per dimension.
    length_scale : array-like of shape (n_dims,), default=1.0
        Length scale per dimension.
    nu : float, default=1.5
        Smoothness parameter for Matern kernel.
    """
    def __init__(self, kernel_type="rbf", length_scale=1.0, nu=1.5):
        self.kernel_type = kernel_type
        self.length_scale = np.atleast_1d(length_scale).astype(float)
        self.nu = nu

    # ---- sklearn plumbing ----
    @property
    def hyperparameters(self):
        return [Hyperparameter("length_scale", "numeric", (1e-5, 1e5), len(self.length_scale))]

    @property
    def theta(self):
        return np.log(self.length_scale)

    @theta.setter
    def theta(self, theta):
        self.length_scale = np.exp(theta)

    @property
    def bounds(self):
        return np.log(np.array([[1e-5, 1e5]] * len(self.length_scale)))

    # ---- main kernel ----
    def __call__(self, X, Y=None, eval_gradient=False):
        if Y is None:
            Y = X
        n_dims = X.shape[1]
        K = np.zeros((X.shape[0], Y.shape[0]))
        grads = []

        for d in range(n_dims):
            ls = self.length_scale[d]
            X_d = X[:, [d]] / ls
            Y_d = Y[:, [d]] / ls
            dists = np.abs(X_d - Y_d.T)

            if self.kernel_type == "rbf":
                K_d = np.exp(-0.5 * dists**2)

            elif self.kernel_type == "matern":
                if self.nu == 0.5:
                    K_d = np.exp(-dists)
                elif self.nu == 1.5:
                    sqrt3 = np.sqrt(3.0)
                    K_d = (1.0 + sqrt3 * dists) * np.exp(-sqrt3 * dists)
                elif self.nu == 2.5:
                    sqrt5 = np.sqrt(5.0)
                    K_d = (1.0 + sqrt5 * dists + 5.0/3.0*dists**2) * np.exp(-sqrt5 * dists)
                elif self.nu == np.inf:  # limit -> RBF
                    K_d = np.exp(-0.5 * dists**2)
                else:
                    eps = 1e-6
                    dists[dists == 0.0] += eps
                    tmp = np.sqrt(2 * self.nu) * dists
                    coef = (2**(1.0 - self.nu)) / gamma(self.nu)
                    K_d = coef * (tmp**self.nu) * kv(self.nu, tmp)
                    K_d = np.where(np.isnan(K_d), 1.0, K_d)

            else:
                raise ValueError(f"Unknown kernel_type {self.kernel_type}")

            K += K_d

            if eval_gradient and self.kernel_type == "rbf":
                # dK/d log(length_scale)
                grad_d = K_d * (dists**2)   # shape (n_samples_X, n_samples_Y)
                grad_d = grad_d[:, :, np.newaxis]
                grads.append(grad_d)

        K /= n_dims  # average
        if eval_gradient and self.kernel_type == "rbf":
            grad = np.concatenate(grads, axis=2) / n_dims
            return K, grad
        elif eval_gradient:
            # For Matern gradients: left as extension
            return K, np.empty((X.shape[0], Y.shape[0], len(self.length_scale)))
        return K

    def diag(self, X):
        return np.ones(X.shape[0])

    def is_stationary(self):
        return True

    def __repr__(self):
        return f"AdditiveContinuousKernel({self.kernel_type}, length_scale={self.length_scale}, nu={self.nu})"


class HybridKernel(Kernel):
    """Hybrid kernel = (average over continuous features) * (fingerprint kernel)."""

    def __init__(self, cont_kernel, fp_kernel, cont_idx, fp_idx):
        self.cont_kernel = cont_kernel
        self.fp_kernel = fp_kernel
        self.cont_idx = cont_idx
        self.fp_idx = fp_idx

    @property
    def hyperparameters(self):
        return [*self.cont_kernel.hyperparameters, *self.fp_kernel.hyperparameters]

    @property
    def theta(self):
        return np.hstack([self.cont_kernel.theta, self.fp_kernel.theta])

    @theta.setter
    def theta(self, theta):
        n1 = self.cont_kernel.theta.shape[0]
        self.cont_kernel.theta = theta[:n1]
        self.fp_kernel.theta = theta[n1:]

    @property
    def bounds(self):
        return np.vstack([self.cont_kernel.bounds, self.fp_kernel.bounds])

    def __call__(self, X, Y=None, eval_gradient=False):
        X_cont, X_fp = X[:, self.cont_idx], X[:, self.fp_idx]
        Y_cont = None if Y is None else Y[:, self.cont_idx]
        Y_fp   = None if Y is None else Y[:, self.fp_idx]

        if eval_gradient:
            Kc, grad_c = self.cont_kernel(X_cont, Y_cont, eval_gradient=True)
            Kf, grad_f = self.fp_kernel(X_fp, Y_fp, eval_gradient=True)
            K = Kc * Kf
            grads = []
            if grad_c.shape[2] > 0:
                grads.append(grad_c * Kf[:, :, np.newaxis])
            if grad_f.shape[2] > 0:
                grads.append(grad_f * Kc[:, :, np.newaxis])
            grad = np.concatenate(grads, axis=2) if grads else np.empty((X.shape[0], X.shape[0], 0))
            return K, grad
        else:
            return self.cont_kernel(X_cont, Y_cont) * self.fp_kernel(X_fp, Y_fp)

    def diag(self, X):
        return self.cont_kernel.diag(X[:, self.cont_idx]) * self.fp_kernel.diag(X[:, self.fp_idx])

    def is_stationary(self):
        return self.cont_kernel.is_stationary() and self.fp_kernel.is_stationary()

    def __repr__(self):
        return f"HybridKernel({self.cont_kernel}, {self.fp_kernel})"
