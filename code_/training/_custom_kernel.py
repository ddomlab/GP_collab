from sklearn.gaussian_process.kernels import Kernel, StationaryKernelMixin, NormalizedKernelMixin, Hyperparameter
import numpy as np
from scipy.spatial.distance import pdist, cdist, squareform
from sklearn.utils import check_array


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
# from sklearn.gaussian_process.kernels import Kernel, Hyperparameter
# from sklearn.utils import check_array
from sklearn.gaussian_process.kernels import _check_length_scale

# import numpy as np
from sklearn.gaussian_process.kernels import RBF
# from sklearn.utils import check_array

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


