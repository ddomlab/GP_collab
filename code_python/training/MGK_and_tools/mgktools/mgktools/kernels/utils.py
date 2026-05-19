# !/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import pickle
from typing import Dict, List, Literal, Tuple
from mgktools.data.data import Dataset
from mgktools.kernels.FeatureKernel import FeatureKernelConfig
from mgktools.kernels.PreComputed import PreComputedKernelConfig
from mgktools.kernels.HybridKernel import HybridKernelConfig


def _fallback_names(names: List[str], prefix: str, count: int) -> List[str]:
    if names is not None and len(names) == count:
        return list(names)
    return [f"{prefix}_{i}" for i in range(count)]


def get_kernel_config(
    dataset: Dataset,
    graph_kernel_type: Literal["graph", "pre-computed", 'no'],
    # arguments for vectorized features.
    features_kernel_type: Literal["dot_product", "RBF", "Matern32", "Matern52"] = None,
    features_hyperparameters: List[float] = None,
    features_hyperparameters_bounds: Tuple[float, float] = None,
    features_hyperparameters_file: str = None,
    # arguments for marginalized graph kernel
    mgk_hyperparameters_files: List[str] = None,
    # arguments for pre-computed kernel
    kernel_dict: Dict = None,
    kernel_pkl: str = None,
    hybrid_rule: Literal["product", "sum"] = "product",
    feature_mode: Literal["joint", "per_feature"] = "joint"
):
    if kernel_pkl is not None and os.path.exists(kernel_pkl) and graph_kernel_type == "pre-computed":
        return pickle.load(open(kernel_pkl, "rb"))

    if graph_kernel_type == "pre-computed":
        n_features = dataset.N_features_add
        feature_names = dataset.X_features_add_names
    else:
        n_features = dataset.N_features_mol + dataset.N_features_add
        feature_names = dataset.X_feature_names

    if features_hyperparameters_file is not None:
        if not os.path.exists(features_hyperparameters_file):
            saved_features_hyperparameters_file = os.path.join(
                os.path.dirname(__file__), "../hyperparameters/configs", features_hyperparameters_file
            )
            if os.path.exists(saved_features_hyperparameters_file):
                features_hyperparameters_file = saved_features_hyperparameters_file
            else:
                raise FileNotFoundError(f"{features_hyperparameters_file} not found.")
        features_kernel_config = FeatureKernelConfig.load(
            path=".", name=features_hyperparameters_file, idx=0
        )
    else:
        if features_kernel_type is None:
            features_kernel_config = None
        else:
            assert n_features != 0

            if len(features_hyperparameters) == 1:
                features_hyperparameters = features_hyperparameters[0]
            else:
                assert len(features_hyperparameters) == n_features
                assert features_kernel_type == 'RBF'
            hyperdict = {
                "features_kernel": {
                    features_kernel_type: [
                        features_hyperparameters,
                        features_hyperparameters_bounds,
                        None,
                        None,
                    ]
                }
            }
            features_kernel_config = FeatureKernelConfig(hyperdict, idx=0)

    if graph_kernel_type == "graph":
        from mgktools.kernels.GraphKernel import GraphKernelConfig
        graph_kernel_configs = []
        for i, mgk_file in enumerate(mgk_hyperparameters_files):
            if not os.path.exists(mgk_file):
                saved_mgk_file = os.path.join(
                    os.path.dirname(__file__), "../hyperparameters/configs", mgk_file
                )
                if os.path.exists(saved_mgk_file):
                    mgk_file = saved_mgk_file
                else:
                    raise FileNotFoundError(f"{mgk_file} not found.")

            graph_kernel_config = GraphKernelConfig.load(path=".", name=mgk_file, idx=i)
            graph_kernel_configs.append(graph_kernel_config)
        n_configs = len(graph_kernel_configs)
        graph_names = _fallback_names(dataset.X_graph_names, "graph", n_configs)
        if features_kernel_config is not None:
            if feature_mode == "joint":
                feature_kernel_name = feature_names[0] if len(feature_names) == 1 else "features"
                hybrid_kernel_config = HybridKernelConfig(
                    kernel_configs=[*graph_kernel_configs, features_kernel_config],
                    composition=[(i,) for i in range(n_configs)] + [tuple(range(n_configs, n_configs + n_features))],
                    hybrid_rule=hybrid_rule,
                    kernel_names=[*graph_names, feature_kernel_name],
                    composition_names=[(name,) for name in graph_names] + [tuple(feature_names)],
                    lengthscale_names=[None] * n_configs + [tuple(feature_names)],
                )
            elif feature_mode == "per_feature":
                import copy
                feature_kernel_configs = []
                for j in range(n_features):
                    fk = FeatureKernelConfig(
                        copy.deepcopy(features_kernel_config.hyperdict),
                        idx=n_configs + j
                    )
                    feature_kernel_configs.append(fk)

                hybrid_kernel_config = HybridKernelConfig(
                    kernel_configs=[*graph_kernel_configs, *feature_kernel_configs],
                    composition=[(i,) for i in range(n_configs)] +
                                [(n_configs + j,) for j in range(n_features)],
                    hybrid_rule=hybrid_rule,
                    kernel_names=[*graph_names, *feature_names],
                    composition_names=[(name,) for name in graph_names] +
                                      [(name,) for name in feature_names],
                    lengthscale_names=[None] * n_configs +
                                      [(name,) for name in feature_names],
                )
            return hybrid_kernel_config
        elif n_configs >= 2:
            hybrid_kernel_config = HybridKernelConfig(
                kernel_configs=graph_kernel_configs,
                composition=[(i,) for i in range(n_configs)],
                hybrid_rule=hybrid_rule,
                kernel_names=graph_names,
                composition_names=[(name,) for name in graph_names],
                lengthscale_names=[None] * n_configs,
            )
            return hybrid_kernel_config
        else:
            return graph_kernel_configs[0]
    elif graph_kernel_type == "pre-computed":
        if kernel_dict is None:
            kernel_dict = pickle.load(open(kernel_pkl, "rb"))
        precomputed_kernel_config = PreComputedKernelConfig(kernel_dict=kernel_dict)

        if features_kernel_config is not None:
            precomputed_name = "precomputed"
            if feature_mode == "joint":
                feature_kernel_name = feature_names[0] if len(feature_names) == 1 else "features"
                kernel_configs = [precomputed_kernel_config, features_kernel_config]
                composition = [(0,), ] + [tuple(range(1, 1 + n_features))]
                kernel_names = [precomputed_name, feature_kernel_name]
                composition_names = [(precomputed_name,), tuple(feature_names)]
                lengthscale_names = [None, tuple(feature_names)]
            elif feature_mode == "per_feature":
                import copy
                feature_kernel_configs = [
                    FeatureKernelConfig(
                        copy.deepcopy(features_kernel_config.hyperdict),
                        idx=1 + j
                    )
                    for j in range(n_features)
                ]
                kernel_configs = [precomputed_kernel_config, *feature_kernel_configs]
                composition = [(0,), ] + [(1 + j,) for j in range(n_features)]
                kernel_names = [precomputed_name, *feature_names]
                composition_names = [(precomputed_name,), ] + [
                    (name,) for name in feature_names
                ]
                lengthscale_names = [None] + [(name,) for name in feature_names]
            else:
                raise ValueError(f"Unknown feature_mode: {feature_mode}")


            hybrid_kernel_config = HybridKernelConfig(
                kernel_configs=kernel_configs,
                composition=composition,
                hybrid_rule=hybrid_rule,
                kernel_names=kernel_names,
                composition_names=composition_names,
                lengthscale_names=lengthscale_names,
            )

            return hybrid_kernel_config
        else:
            return precomputed_kernel_config
    else:
        assert features_kernel_config is not None
        return features_kernel_config
