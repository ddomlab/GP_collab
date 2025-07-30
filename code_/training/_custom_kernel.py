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




class JaccardRBF(StationaryKernelMixin, NormalizedKernelMixin, Kernel):
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
        def weighted_jaccard(u, v, eps: float = 1e-6):
            min_sum = np.sum(np.minimum(u, v))
            max_sum = np.sum(np.maximum(u, v))
            return 1.0 - ((min_sum + eps) / (max_sum + eps))

        X = np.atleast_2d(X)
        if Y is None:
            dists = pdist(X, metric=weighted_jaccard)
            K = np.exp(-squareform(dists) / (2 * self.length_scale**2))
            np.fill_diagonal(K, 1.0)
        else:
            if eval_gradient:
                raise ValueError("Gradient can only be evaluated when Y is None.")
            Y = np.atleast_2d(Y)
            dists = cdist(X, Y, metric=weighted_jaccard)
            K = np.exp(-dists / (2 * self.length_scale**2))

        if eval_gradient:
            if self.hyperparameter_length_scale.fixed:
                return K, np.empty((X.shape[0], X.shape[0], 0))
            else:
                raise NotImplementedError("Gradient for weighted Jaccard is not implemented.")
        else:
            return K

    def __repr__(self):
        if self.anisotropic:
            return "{0}(length_scale=[{1}])".format(
                self.__class__.__name__,
                ", ".join(map("{0:.3g}".format, self.length_scale)),
            )
        else:
            return "{0}(length_scale={1:.3g})".format(
                self.__class__.__name__, np.ravel(self.length_scale)[0]
            )



class JaccardMatern(StationaryKernelMixin, NormalizedKernelMixin, Kernel):
    def __init__(self, length_scale=1.0, length_scale_bounds=(1e-5, 1e5), nu=1.5, jitter=1e-6):
        self.length_scale = length_scale
        self.length_scale_bounds = length_scale_bounds
        self.nu = nu
        self.jitter = jitter

    @property
    def anisotropic(self):
        return False

    @property
    def hyperparameter_length_scale(self):
        return Hyperparameter("length_scale", "numeric", self.length_scale_bounds)

    def __call__(self, X, Y=None, eval_gradient=False):
        def weighted_jaccard(u, v, eps: float = 1e-6):
            min_sum = np.sum(np.minimum(u, v))
            max_sum = np.sum(np.maximum(u, v))
            return 1.0 - ((min_sum + eps) / (max_sum + eps))

        X = np.atleast_2d(X)
        length_scale = float(self.length_scale)

        if Y is None:
            dists = pdist(X, metric=weighted_jaccard)
            dists = squareform(dists)
        else:
            if eval_gradient:
                raise ValueError("Gradient can only be evaluated when Y is None.")
            Y = np.atleast_2d(Y)
            dists = cdist(X, Y, metric=weighted_jaccard)

        dists /= length_scale

        if self.nu == 0.5:
            K = np.exp(-dists)
        elif self.nu == 1.5:
            K = sqrt(3) * dists
            K = (1.0 + K) * np.exp(-K)
        elif self.nu == 2.5:
            K = sqrt(5) * dists
            K = (1.0 + K + K**2 / 3.0) * np.exp(-K)
        elif self.nu == np.inf:
            K = np.exp(-(dists**2) / 2.0)
        else:
            dists[dists == 0.0] += np.finfo(float).eps
            tmp = np.sqrt(2 * self.nu) * dists
            K = (2**(1.0 - self.nu)) / gamma(self.nu) * tmp**self.nu * kv(self.nu, tmp)

        if Y is None:
            np.fill_diagonal(K, 1.0 + self.jitter)
        else:
            # No diagonal to adjust in K(X, Y)
            pass

        if eval_gradient:
            if self.hyperparameter_length_scale.fixed:
                return K, np.empty((X.shape[0], X.shape[0], 0))
            else:
                raise NotImplementedError("Gradient for weighted Jaccard is not implemented.")

        return K

    def __repr__(self):
        return f"{self.__class__.__name__}(length_scale={self.length_scale:.3g}, nu={self.nu:.3g})"



class AdditiveRBF(RBF):
    def __init__(self, length_scale=1.0, length_scale_bounds=(1e-5, 1e5)):
        super().__init__(length_scale=length_scale, length_scale_bounds=length_scale_bounds)

    def __call__(self, X, Y=None, eval_gradient=False):
        """Compute the additive RBF kernel and optionally its gradient."""
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
            X_d = X[:, d:d+1] / length_scale[d]
            Y_d = Y[:, d:d+1] / length_scale[d]
            diff = X_d - Y_d.T
            K_d = np.exp(-0.5 * diff ** 2)
            K += K_d

            if eval_gradient:
                orig_diff = (X[:, d:d+1] - Y[:, d:d+1].T)
                K_gradient[:, :, d] = K_d * (orig_diff ** 2) / (length_scale[d] ** 3)

        if eval_gradient:
            return K, K_gradient
        return K

    def diag(self, X):
        X = check_array(X)
        return np.full(X.shape[0], X.shape[1])  # Sum of exp(0) = 1 per dimension

    def is_stationary(self):
        return False  # Additive kernels are not stationary

    def __repr__(self):
        if self.anisotropic:
            return "{0}(length_scale=[{1}])".format(
                self.__class__.__name__,
                ", ".join(map("{0:.3g}".format, np.ravel(self.length_scale))),
            )
        else:
            return "{0}(length_scale={1:.3g})".format(
                self.__class__.__name__, np.ravel(self.length_scale)[0]
            )


def _approx_fprime(xk, f, epsilon, args=()):
    f0 = f(*((xk,) + args))
    grad = np.zeros((f0.shape[0], f0.shape[1], len(xk)), float)
    ei = np.zeros((len(xk),), float)
    for k in range(len(xk)):
        ei[k] = 1.0
        d = epsilon * ei
        grad[:, :, k] = (f(*((xk + d,) + args)) - f0) / d[k]
        ei[k] = 0.0
    return grad


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

        if Y is None:
            np.fill_diagonal(K, X.shape[1])  # sum of 1s from each dim when dist == 0

        if eval_gradient:
            return K, K_gradient
        return K

    def diag(self, X):
        X = check_array(X)
        return np.full(X.shape[0], X.shape[1])

    def is_stationary(self):
        return False  # Sum of 1D kernels is not stationary

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
