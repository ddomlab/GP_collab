#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pathlib import Path
# import pytest
import pandas as pd
import numpy as np
# from mgktools.data.data import Dataset
# from mgktools.kernels.utils import get_kernel_config
# from mgktools.hyperparameters import *
# from mgktools.models import GaussianProcessRegressor




from training.filter_data import _get_dataset_features, get_structural_info

HERE = Path(__file__).resolve().parent
DATASETS = HERE.parent.parent / "datasets" / "Validation datasets"
RESULTS = HERE.parent.parent / "results"


PAPER = "Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation"
# "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation"
dataset_name = "Rg data with clusters aging imputed"
# "cleaned_dataset_pervaporation_membranes_wang"
#non_imputed_dropped_nan_Rg_data
w_data, feats, all_targets, polymer_unit = _get_dataset_features(DATASETS, PAPER, dataset_name)



smiles = w_data[f"{polymer_unit} SMILES"].tolist()
targets = w_data[all_targets].tolist()

df = pd.DataFrame({'smiles': smiles, 'targets': targets})

print(df.head())
# @pytest.mark.parametrize('mgk_file', [
#                                       product])
# @pytest.mark.parametrize('loss_function', ['likelihood'])
# @pytest.mark.parametrize('optimizer', ['L-BFGS-B'])
# def test_gradient_Graph(mgk_file, loss_function, optimizer):
#     dataset = Dataset.from_df(df=df,
#                               smiles_columns=['smiles'],
#                               targets_columns=['targets'])
#     dataset.set_status(graph_kernel_type='graph', features_generators=None, features_combination=None)
#     dataset.create_graphs(n_jobs=4)
#     dataset.unify_datatype()
#     kernel_config = get_kernel_config(dataset=dataset,
#                                       graph_kernel_type='graph',
#                                       mgk_hyperparameters_files=[mgk_file])
#     gpr = GaussianProcessRegressor(kernel=kernel_config.kernel,
#                                    optimizer=optimizer,
#                                    alpha=0.01,
#                                    normalize_y=True)
#     gpr.fit(dataset.X, dataset.y, loss=loss_function, verbose=True)
