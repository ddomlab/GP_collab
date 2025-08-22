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





class JaccardKernel(NormalizedKernelMixin, Kernel):
    """Pure Jaccard similarity kernel (no RBF transform)."""

    def __init__(self):
        pass  # no hyperparameters here

    def __call__(self, X, Y=None, eval_gradient=False):
        def weighted_jaccard(u, v, eps: float = 1e-6):
            min_sum = np.sum(np.minimum(u, v))
            max_sum = np.sum(np.maximum(u, v))
            # return distance
            return 1.0 - ((min_sum + eps) / (max_sum + eps))

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



# import numpy as np
# from sklearn.gaussian_process.kernels import Kernel, StationaryKernelMixin, NormalizedKernelMixin, Hyperparameter
# from scipy.spatial.distance import pdist, cdist, squareform


def weighted_jaccard(u, v, eps: float = 1e-6):
    """Weighted Jaccard distance between two vectors u, v."""
    min_sum = np.sum(np.minimum(u, v))
    max_sum = np.sum(np.maximum(u, v))
    return 1.0 - ((min_sum + eps) / (max_sum + eps))


class Jaccard_gradient_RBF(StationaryKernelMixin, NormalizedKernelMixin, Kernel):
    """RBF kernel but using Jaccard distance instead of Euclidean distance."""

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
                "length_scale", "numeric", self.length_scale_bounds, len(self.length_scale)
            )
        return Hyperparameter("length_scale", "numeric", self.length_scale_bounds)

    def __call__(self, X, Y=None, eval_gradient=False):
        X = np.atleast_2d(X)
        length_scale = np.array(self.length_scale, ndmin=1)

        # Rescale features for anisotropy
        X_scaled = X / length_scale
        if Y is None:
            Y_scaled = None
        else:
            Y_scaled = np.atleast_2d(Y) / length_scale

        # ----- Distance computation -----
        if Y is None:
            dists = pdist(X_scaled, metric=lambda u, v: weighted_jaccard(u, v))
            K = np.exp(-0.5 * squareform(dists)**2)
            np.fill_diagonal(K, 1.0)
        else:
            if eval_gradient:
                raise ValueError("Gradient can only be evaluated when Y is None.")
            dists = cdist(X_scaled, Y_scaled, metric=lambda u, v: weighted_jaccard(u, v))
            K = np.exp(-0.5 * dists**2)

        # ----- Gradient computation -----
        if eval_gradient:
            if self.hyperparameter_length_scale.fixed:
                return K, np.empty((X.shape[0], X.shape[0], 0))

            # Gradient wrt length_scales
            if not self.anisotropic or length_scale.shape[0] == 1:
                # isotropic: derivative wrt log(l)
                # dK/dlog(l) = K * d(dists^2)/dlog(l)
                # but since distance is Jaccard, analytic gradient is tricky
                # fallback: numerical approx
                n = X.shape[0]
                K_gradient = np.zeros((n, n, 1))
                eps = 1e-5
                for i in range(n):
                    for j in range(n):
                        if i == j:
                            continue
                        d0 = weighted_jaccard(X_scaled[i], X_scaled[j])
                        # perturb l
                        l_eps = length_scale * np.exp(eps)
                        Xi = X[i] / l_eps
                        Xj = X[j] / l_eps
                        d1 = weighted_jaccard(Xi, Xj)
                        dd = (d1**2 - d0**2) / eps
                        K_gradient[i, j, 0] = K[i, j] * (-0.5 * dd)
                return K, K_gradient
            else:
                # anisotropic: one gradient per dimension
                n, d = X.shape
                K_gradient = np.zeros((n, n, d))
                eps = 1e-5
                for dim in range(d):
                    l_eps = length_scale.copy()
                    l_eps[dim] *= np.exp(eps)
                    X_scaled_eps = X / l_eps
                    dists_eps = pdist(X_scaled_eps, metric=lambda u, v: weighted_jaccard(u, v))
                    dists_eps = squareform(dists_eps)
                    dists0 = pdist(X_scaled, metric=lambda u, v: weighted_jaccard(u, v))
                    dists0 = squareform(dists0)
                    dd = (dists_eps**2 - dists0**2) / eps
                    K_gradient[:, :, dim] = K * (-0.5 * dd)
                return K, K_gradient

        return K

    def __repr__(self):
        if self.anisotropic:
            return "{0}(length_scale=[{1}])".format(
                self.__class__.__name__, ", ".join(map("{0:.3g}".format, self.length_scale))
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





from sklearn.gaussian_process.kernels import Kernel
import numpy as np
import pandas as pd

class CombinedKernel(Kernel):
    """
    A kernel that applies different kernels to named feature subsets
    and sums the results. Accepts DataFrame at init, but will work with
    numpy arrays when sklearn calls __call__.
    """

    def __init__(self, feature_groups, kernels, all_columns):
        """
        Parameters
        ----------
        feature_groups : dict
            Mapping { "group_name": [list_of_feature_names] }
        kernels : dict
            Mapping { "group_name": Kernel instance }
        all_columns : list
            List of all column names in the dataset (X DataFrame.columns)
        """
        assert feature_groups.keys() == kernels.keys(), \
            "feature_groups and kernels must have the same keys"
        self.feature_groups = feature_groups
        self.kernels = kernels
        self.all_columns = list(all_columns)

        # Map feature names → column indices
        self.feature_indices = {
            group: [self.all_columns.index(c) for c in cols]
            for group, cols in feature_groups.items()
        }

    def __call__(self, X, Y=None, eval_gradient=False):
        # sklearn will pass numpy arrays here
        K_total = None
        grads = []

        for group, indices in self.feature_indices.items():
            X_sub = X[:, indices]
            Y_sub = None if Y is None else Y[:, indices]
            kernel = self.kernels[group]

            if eval_gradient:
                K, grad = kernel(X_sub, Y_sub, eval_gradient=True)
                grads.append(grad)
            else:
                K = kernel(X_sub, Y_sub, eval_gradient=False)

            K_total = K if K_total is None else K_total + K

        if eval_gradient:
            return K_total, np.concatenate(grads, axis=2)
        else:
            return K_total

    def diag(self, X):
        return sum(
            self.kernels[group].diag(X[:, indices])
            for group, indices in self.feature_indices.items()
        )

    def is_stationary(self):
        return all(kernel.is_stationary() for kernel in self.kernels.values())

    def __repr__(self):
        parts = [
            f"{group}: {self.kernels[group]} on {self.feature_groups[group]}"
            for group in self.feature_groups
        ]
        return " + ".join(parts)
