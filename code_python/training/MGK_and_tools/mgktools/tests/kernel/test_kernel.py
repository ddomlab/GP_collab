#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pytest
import pandas as pd
import numpy as np
from mgktools.data.data import Dataset
from mgktools.features_mol.features_generators import FeaturesGenerator
from mgktools.kernels.utils import get_kernel_config
from mgktools.hyperparameters import *


smiles_1 = ['CCCC', 'CCCCCO', 'c1ccccc1', 'CCNCCO']
smiles_2 = ['CCCCC', 'CCCCCCO', 'Cc1ccccc1', 'CCCNCCO']
targets_1 = [3.4, 4.5, 5.6, 6.7]
targets_2 = [3.1, 14.5, 25.6, 56.7]
features_add = [0.6, 0.1, 0.3, 0.9]
df = pd.DataFrame({'smiles_1': smiles_1, 'smiles_2': smiles_2, 'targets_1': targets_1, 'targets_2': targets_2, 'features_add': features_add})


@pytest.mark.parametrize('smiles_columns', [['smiles_1'], ['smiles_1', 'smiles_2']])
@pytest.mark.parametrize('features_columns', [None, ['features_add']])
@pytest.mark.parametrize('mgk_file', [additive, additive_norm, additive_pnorm, additive_msnorm,
                                      product, product_norm, product_pnorm, product_msnorm])
@pytest.mark.parametrize('fingerprints', [
    (None, None, None),
    (['morgan'], 'concat', 2048),
    (['rdkit_2d_normalized'], 'concat', 200),
    (['morgan', 'rdkit_2d_normalized'], 'mean', 2248),
    (['morgan', 'rdkit_2d_normalized'], 'concat', 2248),
])
@pytest.mark.parametrize('features_kernel_type', ['dot_product', 'rbf'])
@pytest.mark.parametrize('normalize_feature', [True, False])
def test_kernel(smiles_columns, features_columns, mgk_file, fingerprints, features_kernel_type, normalize_feature):
    features_generators, features_combination, n_features = fingerprints
    if features_generators is not None:
        fgs = [FeaturesGenerator(features_generator_name=fg) for fg in features_generators]
    else:
        fgs = None
    dataset = Dataset.from_df(df=df,
                              smiles_columns=smiles_columns,
                              features_columns=features_columns,
                              targets_columns=['targets_1'])
    dataset.set_status(graph_kernel_type='graph',
                       features_generators=fgs,
                       features_combination=features_combination)
    dataset.create_graphs(n_jobs=2)
    dataset.unify_datatype()
    if features_generators is not None:
        dataset.create_features_mol(n_jobs=2)
    if features_generators is None and features_columns is None:
        kernel_config = get_kernel_config(dataset=dataset,
                                          graph_kernel_type='graph',
                                          mgk_hyperparameters_files=[mgk_file] * len(smiles_columns))
    else:
        kernel_config = get_kernel_config(dataset=dataset,
                                          graph_kernel_type='graph',
                                          features_kernel_type=features_kernel_type,
                                          features_hyperparameters=[0.5],
                                          features_hyperparameters_bounds='fixed',
                                          mgk_hyperparameters_files=[mgk_file] * len(smiles_columns))
    N = len(dataset)
    K = kernel_config.kernel(dataset.X)
    assert K.shape == (N, N)
    # invertable
    assert np.linalg.det(K) > 10 ** -3
    for i in range(N):
        for j in range(i + 1, N):
            # symmetric
            assert K[i, j] == pytest.approx(K[j, i], 1e-5)
            # diagonal largest
            assert np.sqrt(K[i, i] * K[j, j]) > K[i, j]
