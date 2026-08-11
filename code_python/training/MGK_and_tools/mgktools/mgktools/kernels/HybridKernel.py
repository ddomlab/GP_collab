#!/usr/bin/env python
# -*- coding: utf-8 -*-
from typing import Dict, List, Literal, Tuple, Any, Optional
import copy
import numpy as np
from sklearn.gaussian_process.kernels import RBF, DotProduct, Matern
from mgktools.kernels.base import BaseKernelConfig


HybridRule = Literal[
    "product",
    "sum",
    "(count:+)x(graph:x)",
    "(count:+)x(graph:+)",
    "(count:x)+(graph:x)",
]


class HybridKernel:
    def __init__(
        self,
        kernel_list: List,
        composition: List[Tuple[int]],
        hybrid_rule: HybridRule = "product",
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
        self.graph_kernel_indices, self.count_kernel_indices = (
            self._get_kernel_groups()
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

    def _get_kernel_groups(self) -> Tuple[List[int], List[int]]:
        if self.lengthscale_names is None:
            return list(range(self.nkernel)), []

        graph_indices = []
        count_indices = []
        for i in range(self.nkernel):
            names = (
                self.lengthscale_names[i]
                if i < len(self.lengthscale_names)
                else None
            )
            if names is None:
                graph_indices.append(i)
            else:
                count_indices.append(i)
        return graph_indices, count_indices

    @staticmethod
    def _as_float_array(value):
        if getattr(value, "dtype", None) == object:
            return np.asarray(value, dtype=np.float64)
        return np.asarray(value, dtype=np.float64)

    @staticmethod
    def _align_matrix_dims(matrix, target_ndim):
        while np.ndim(matrix) < target_ndim:
            matrix = np.expand_dims(matrix, axis=-1)
        return matrix

    @classmethod
    def _multiply_matrices(cls, matrices: List[np.ndarray]):
        target_ndim = max(np.ndim(matrix) for matrix in matrices)
        out = 1.
        for matrix in matrices:
            matrix = cls._align_matrix_dims(matrix, target_ndim)
            out = out * matrix
        return out

    @classmethod
    def _sum_matrices(
        cls, matrices: List[np.ndarray], average: bool = False
    ):
        target_ndim = max(np.ndim(matrix) for matrix in matrices)
        out = 0.
        for matrix in matrices:
            matrix = cls._align_matrix_dims(matrix, target_ndim)
            out = out + matrix
        if average:
            out = out / float(len(matrices))
        return out

    def _kernel_outputs(self, X_list, Y_list, Y, eval_gradient):
        covariances = []
        gradients = []

        for i, kernel in enumerate(self.kernel_list):
            Xi = X_list[i]
            Yi = Y_list[i] if Y is not None else None
            output = kernel(Xi, Y=Yi, eval_gradient=eval_gradient)

            if eval_gradient:
                covariance, gradient = output
                covariances.append(self._as_float_array(covariance))
                gradients.append(self._as_float_array(gradient))
            else:
                covariances.append(self._as_float_array(output))

        return covariances, gradients

    def _validate_grouped_rule(self):
        if not self.graph_kernel_indices:
            raise ValueError(f"{self.hybrid_rule} requires at least one graph kernel")
        if not self.count_kernel_indices:
            raise ValueError(f"{self.hybrid_rule} requires at least one count kernel")

    def _combine_grouped_covariances(self, covariances):
        self._validate_grouped_rule()

        graph_covariances = [covariances[i] for i in self.graph_kernel_indices]
        count_covariances = [covariances[i] for i in self.count_kernel_indices]

        if self.hybrid_rule == "(count:+)x(graph:x)":
            return (
                self._sum_matrices(count_covariances)
                * self._multiply_matrices(graph_covariances)
            )

        if self.hybrid_rule == "(count:+)x(graph:+)":
            return (
                self._sum_matrices(count_covariances)
                * self._sum_matrices(graph_covariances)
            )

        if self.hybrid_rule == "(count:x)+(graph:x)":
            return (
                self._multiply_matrices(count_covariances)
                + self._multiply_matrices(graph_covariances)
            )

        raise ValueError(f"Unknown hybrid rule {self.hybrid_rule}")

    def _combine_covariances(self, covariances):
        if self.hybrid_rule == "product":
            return self._multiply_matrices(covariances)
        if self.hybrid_rule == "sum":
            return self._sum_matrices(covariances)
        return self._combine_grouped_covariances(covariances)

    def _combine_gradients(self, covariances, gradients):
        gradient_matrix_list = []

        for i, gradient in enumerate(gradients):
            gradient_matrix_list.append(
                self._kernel_gradient_covariance(i, covariances, gradient)
            )

        gradient_matrix = gradient_matrix_list[0]
        for i, gm in enumerate(gradient_matrix_list):
            if i != 0:
                gradient_matrix = np.c_[gradient_matrix, gm]

        return np.asarray(gradient_matrix, dtype=np.float64)

    def _kernel_gradient_covariance(self, kernel_index, covariances, gradient):
        other_indices = [i for i in range(self.nkernel) if i != kernel_index]

        if self.hybrid_rule == "product":
            return self._multiply_matrices(
                [gradient] + [covariances[i] for i in other_indices]
            )

        if self.hybrid_rule == "sum":
            return gradient

        self._validate_grouped_rule()
        graph_indices = self.graph_kernel_indices
        count_indices = self.count_kernel_indices
        graph_covariances = [covariances[i] for i in graph_indices]
        count_covariances = [covariances[i] for i in count_indices]

        if self.hybrid_rule == "(count:+)x(graph:x)":
            if kernel_index in count_indices:
                return self._multiply_matrices(
                    [
                        gradient,
                        self._multiply_matrices(graph_covariances),
                    ]
                )

            graph_others = [
                covariances[i] for i in graph_indices if i != kernel_index
            ]
            return self._multiply_matrices(
                [gradient, self._sum_matrices(count_covariances)]
                + graph_others
            )

        if self.hybrid_rule == "(count:+)x(graph:+)":
            if kernel_index in count_indices:
                return self._multiply_matrices(
                    [gradient, self._sum_matrices(graph_covariances)]
                )

            return self._multiply_matrices(
                [gradient, self._sum_matrices(count_covariances)]
            )

        if self.hybrid_rule == "(count:x)+(graph:x)":
            if kernel_index in count_indices:
                count_others = [
                    covariances[i] for i in count_indices if i != kernel_index
                ]
                return self._multiply_matrices([gradient] + count_others)

            graph_others = [
                covariances[i] for i in graph_indices if i != kernel_index
            ]
            return self._multiply_matrices([gradient] + graph_others)

        raise ValueError(f"Unknown hybrid rule {self.hybrid_rule}")

    def __call__(
        self, X: np.ndarray, Y: np.ndarray = None, eval_gradient: bool = False
    ):
        X_list = self.get_X_list(X)
        Y_list = self.get_X_list(Y) if Y is not None else None

        covariances, gradients = self._kernel_outputs(
            X_list,
            Y_list,
            Y,
            eval_gradient=eval_gradient,
        )
        covariance_matrix = self._combine_covariances(covariances)

        if not eval_gradient:
            return covariance_matrix

        gradient_matrix = self._combine_gradients(covariances, gradients)
        return covariance_matrix, gradient_matrix

    def diag(self, X) -> List[float]:
        X_list = self.get_X_list(X)
        diag_list = [
            self._as_float_array(self.kernel_list[i].diag(X_list[i]))
            for i in range(len(self.kernel_list))
        ]
        return self._combine_covariances(diag_list)

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
        hybrid_rule: HybridRule = "product",
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

