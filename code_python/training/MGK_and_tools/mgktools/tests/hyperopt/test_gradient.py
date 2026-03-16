#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pytest
import pandas as pd
import numpy as np
from mgktools.data.data import Dataset
from mgktools.kernels.utils import get_kernel_config
from mgktools.hyperparameters import *
from mgktools.models import GaussianProcessRegressor


smiles = ['CCCC', 'CCCCCO', 'c1ccccc1', 'CCNCCO', 'OCCCO']
targets = [3.1, 14.5, 25.6, 56.7, 12.3]
df = pd.DataFrame({'smiles': smiles, 'targets': targets})


@pytest.mark.parametrize('mgk_file', [additive, additive_norm, additive_pnorm, additive_msnorm,
                                      product, product_norm, product_pnorm, product_msnorm])
@pytest.mark.parametrize('loss_function', ['loocv', 'likelihood'])
@pytest.mark.parametrize('optimizer', ['L-BFGS-B', 'SLSQP'])
def test_gradient_Graph(mgk_file, loss_function, optimizer):
    dataset = Dataset.from_df(df=df,
                              smiles_columns=['smiles'],
                              targets_columns=['targets'])
    dataset.set_status(graph_kernel_type='graph', features_generators=None, features_combination=None)
    dataset.create_graphs(n_jobs=4)
    dataset.unify_datatype()
    kernel_config = get_kernel_config(dataset=dataset,
                                      graph_kernel_type='graph',
                                      mgk_hyperparameters_files=[mgk_file])
    gpr = GaussianProcessRegressor(kernel=kernel_config.kernel,
                                   optimizer=optimizer,
                                   alpha=0.01,
                                   normalize_y=True)
    gpr.fit(dataset.X, dataset.y, loss=loss_function, verbose=True)
