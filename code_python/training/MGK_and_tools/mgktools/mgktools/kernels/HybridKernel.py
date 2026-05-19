#!/usr/bin/env python
# -*- coding: utf-8 -*-
from typing import Dict, List, Literal, Tuple, Any, Optional
import copy
import numpy as np
from sklearn.gaussian_process.kernels import RBF, DotProduct, Matern
from mgktools.kernels.base import BaseKernelConfig


class HybridKernel:
    def __init__(
        self,
        kernel_list: List,
        composition: List[Tuple[int]],
        hybrid_rule: Literal["product", "sum"] = "product",
        kernel_names: List[str] = None,
        composition_names: List[Tuple[str, ...]] = None,
        lengthscale_names: List[Optional[Tuple[str, ...]]] = None,
    ):
        self.kernel_list = kernel_list
        self.composition = composition
        self.hybrid_rule = hybrid_rule
        self.kernel_names = list(kernel_names) if kernel_names is not None else [
            f"kernel_{i}" for i in range(len(kernel_list))
        ]
        self.composition_names = (
            [tuple(names) for names in composition_names]
            if composition_names is not None
            else None
        )
        self.lengthscale_names = (
            [None if names is None else tuple(names) for names in lengthscale_names]
            if lengthscale_names is not None
            else None
        )

    @property
    def nkernel(self) -> int:
        return len(self.kernel_list)

    def get_X_list(self, X: np.ndarray) -> List[np.ndarray]:
        # for i in range(X.shape[1]):
        #     print(f"col {i}: {type(X[0, i]).__name__}")
        
        def f(c):
            return X[:, c]

        X = self._format_X(X)
        X_list = list(map(f, self.composition))
        for i, kernel in enumerate(self.kernel_list):
            if kernel.__class__ in [RBF, DotProduct, Matern]:
                X_list[i] = X_list[i].astype("float64")
        return X_list

    @staticmethod
    def _format_X(X: np.ndarray) -> np.ndarray:
        if X.ndim == 1:
            return X.reshape(1, X.size)  # .tolist()
        else:
            return X

    def __call__(
        self, X: np.ndarray, Y: np.ndarray = None, eval_gradient: bool = False
    ):
        X_list = self.get_X_list(X)
        Y_list = self.get_X_list(Y) if Y is not None else None
        if self.hybrid_rule == "product":
            if eval_gradient:
                covariance_matrix = 1.
                gradient_matrix_list = list(map(int, np.ones(self.nkernel).tolist()))
                for i, kernel in enumerate(self.kernel_list):
                    Xi = X_list[i]
                    Yi = Y_list[i] if Y is not None else None
                    output = kernel(Xi, Y=Yi, eval_gradient=True)
                    covariance_matrix *= np.asarray(output[0], dtype=np.float64)
                    for j in range(self.nkernel):
                        if j == i:
                            gradient_matrix_list[j] = (
                                gradient_matrix_list[j] * output[1]
                            )
                        else:
                            shape = output[0].shape + (1,)
                            gradient_matrix_list[j] = gradient_matrix_list[j] * output[
                                0
                            ].reshape(shape)
                gradient_matrix = gradient_matrix_list[0]
                for i, gm in enumerate(gradient_matrix_list):
                    if i != 0:
                        gradient_matrix = np.c_[gradient_matrix, gradient_matrix_list[i]]
                return covariance_matrix, np.asarray(gradient_matrix, dtype=np.float64)
            else:
                covariance_matrix = 1.
                for i, kernel in enumerate(self.kernel_list):
                    Xi = X_list[i]
                    Yi = Y_list[i] if Y is not None else None
                    output = kernel(Xi, Y=Yi, eval_gradient=False)
                    if output.dtype == object:
                        output = np.asarray(output, dtype=np.float64)
                    covariance_matrix *= output
                return covariance_matrix
        elif self.hybrid_rule == "sum":
            if eval_gradient:
                covariance_matrix = 0.
                gradient_matrix = 0.
                for i, kernel in enumerate(self.kernel_list):
                    Xi = X_list[i]
                    Yi = Y_list[i] if Y is not None else None
                    output = kernel(Xi, Y=Yi, eval_gradient=False)
                    covariance_matrix += np.asarray(output[0], dtype=np.float64)
                    gradient_matrix += np.asarray(output[1], dtype=np.float64)
                return covariance_matrix, gradient_matrix
            else:
                covariance_matrix = 0.
                for i, kernel in enumerate(self.kernel_list):
                    Xi = X_list[i]
                    Yi = Y_list[i] if Y is not None else None
                    output = kernel(Xi, Y=Yi, eval_gradient=False)
                    if output.dtype == object:
                        output = np.asarray(output, dtype=np.float64)
                    covariance_matrix += output
                return covariance_matrix
        else:
            raise ValueError

    def diag(self, X) -> List[float]:
        X_list = self.get_X_list(X)
        diag_list = [
            self.kernel_list[i].diag(X_list[i]) for i in range(len(self.kernel_list))
        ]
        if self.hybrid_rule == "product":
            return np.prod(diag_list, axis=0)
        else:
            raise Exception("Unknown hybrid rule %s" % self.hybrid_rule)

    def is_stationary(self):
        return False

    @property
    def requires_vector_input(self):
        return False

    @property
    def n_dims_list(self) -> List[int]:
        """Numbers of hyperparameters."""
        return [len(kernel.theta) for kernel in self.kernel_list]

    @property
    def n_dims(self) -> int:
        return sum(self.n_dims_list)

    @property
    def hyperparameters(self):
        return np.exp(self.theta)

    @property
    def theta(self):
        return np.concatenate([kernel.theta for kernel in self.kernel_list])

    @theta.setter
    def theta(self, value):
        if len(value) != len(self.theta):
            raise Exception("The length of n_dims and theta must the same")
        s = 0
        e = 0
        for i, kernel in enumerate(self.kernel_list):
            e += self.n_dims_list[i]
            kernel.theta = value[s:e]
            s += self.n_dims_list[i]

    @property
    def bounds(self):
        for i, kernel in enumerate(self.kernel_list):
            if i == 0:
                bounds = self.kernel_list[0].bounds
            elif kernel.bounds.shape != (0,):
                bounds = np.r_[bounds, kernel.bounds]
        return bounds

    def clone_with_theta(self, theta):
        clone = copy.deepcopy(self)
        clone.theta = theta
        return clone

    def get_params(self, deep=False):
        return dict(
            kernel_list=self.kernel_list,
            composition=self.composition,
            hybrid_rule=self.hybrid_rule,
            kernel_names=self.kernel_names,
            composition_names=self.composition_names,
            lengthscale_names=self.lengthscale_names,
        )

    def load(self, result_dir):
        for i, kernel in enumerate(self.kernel_list):
            if hasattr(kernel, "PreCalculate"):
                kernel.load(result_dir)


    def get_lengthscale(self):
        out = {}

        for i, kernel in enumerate(self.kernel_list):
            if not hasattr(kernel, "get_params"):
                continue

            params = kernel.get_params()

            if "length_scale" not in params:
                continue

            kernel_name = self.kernel_names[i] if i < len(self.kernel_names) else f"kernel_{i}"
            length_scale = np.asarray(params["length_scale"], dtype=float).ravel()
            names = None

            if self.lengthscale_names is not None:
                names = self.lengthscale_names[i] if i < len(self.lengthscale_names) else None
                if names is None:
                    continue

            if length_scale.size == 1:
                name = names[0] if names is not None and len(names) == 1 else kernel_name
                out[f"{name}"] = float(length_scale[0])
                continue

            if names is None and self.composition_names is not None and i < len(self.composition_names):
                names = self.composition_names[i]
            elif names is None:
                names = tuple(f"{kernel_name}_{j}" for j in range(length_scale.size))

            for j, value in enumerate(length_scale):
                name = names[j] if j < len(names) else f"{kernel_name}_{j}"
                out[f"{name}"] = float(value)

        return out

class HybridKernelConfig(BaseKernelConfig):
    def __init__(
        self,
        kernel_configs: List[BaseKernelConfig],
        composition: List[Tuple[int]],
        hybrid_rule: Literal["product", "sum"] = "product",
        kernel_names: List[str] = None,
        composition_names: List[Tuple[str, ...]] = None,
        lengthscale_names: List[Optional[Tuple[str, ...]]] = None,
    ):
        assert len(kernel_configs) == len(composition) >= 2
        self.kernel_configs = kernel_configs
        self.composition = composition
        self.hybrid_rule = hybrid_rule
        self.kernel_names = list(kernel_names) if kernel_names is not None else [
            f"kernel_{i}" for i in range(len(kernel_configs))
        ]
        self.composition_names = (
            [tuple(names) for names in composition_names]
            if composition_names is not None
            else None
        )
        self.lengthscale_names = (
            [None if names is None else tuple(names) for names in lengthscale_names]
            if lengthscale_names is not None
            else None
        )
        self.update_kernel()

    def update_kernel(self):
        kernels = []
        for kernel_config in self.kernel_configs:
            kernel_config.update_kernel()
            kernels.append(kernel_config.kernel)
        self.kernel = HybridKernel(
            kernels,
            self.composition,
            self.hybrid_rule,
            kernel_names=self.kernel_names,
            composition_names=self.composition_names,
            lengthscale_names=self.lengthscale_names,
        )

    def get_trial(self, trial) -> Dict:
        return self.combine_dicts(
            [kernel_config.get_trial(trial) for kernel_config in self.kernel_configs]
        )

    def update_from_trial(self, trial: Dict[str, Any]):
        for kernel_config in self.kernel_configs:
            kernel_config.update_from_trial(trial)

    def update_from_theta(self):
        for kernel_config in self.kernel_configs:
            kernel_config.update_from_theta()

    def save(self, path: str):
        for i, kernel_config in enumerate(self.kernel_configs):
            kernel_config.save(path=path, name=f'kernel_{i}.json')
    

    # def get_feature_length_scale(self):
    #     out = {}

    #     for i, kernel_config in enumerate(self.kernel_configs):
    #         # FeatureKernelConfig has microkernels_feature.
    #         # GraphKernelConfig does not.
    #         if not hasattr(kernel_config, "microkernels_feature"):
    #             continue

    #         kernel_config.update_from_theta()

    #         for microkernel in kernel_config.microkernels_feature:
    #             kernel_type = list(microkernel.microdict.keys())[0]

    #             out[f"kernel_{i}.{kernel_type}"] = microkernel.value

    #     return out


