#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pathlib import Path
import pytest
import pandas as pd
import numpy as np
from mgktools.data.data import Dataset
from mgktools.kernels.utils import get_kernel_config
from mgktools.hyperparameters import *
from mgktools.models import GaussianProcessRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, root_mean_squared_error

from rdkit import Chem
from filter_data import _get_dataset_features, get_structural_info

HERE = Path(__file__).resolve().parent
DATASETS = HERE.parent.parent / "datasets" / "Validation datasets"
RESULTS = HERE.parent.parent / "results"


PAPER = "Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation"
# "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation"
dataset_name = "Rg data with clusters aging imputed"
# "cleaned_dataset_pervaporation_membranes_wang"
#non_imputed_dropped_nan_Rg_data
w_data, feats, all_targets, polymer_unit = _get_dataset_features(DATASETS, PAPER, dataset_name)



smiles = w_data[f"{polymer_unit[0]} SMILES"].apply(lambda x: Chem.CanonSmiles(x))
targets = w_data[all_targets]
df = pd.concat([smiles, w_data[feats],w_data[all_targets]], axis=1)


# print(df.head())
@pytest.mark.parametrize('mgk_file', [
                                      product
                                      ])
@pytest.mark.parametrize('loss_function', ['likelihood'])
@pytest.mark.parametrize('optimizer', ['L-BFGS-B'])
def test_gradient_Graph(mgk_file, loss_function, optimizer):
    dataset = Dataset.from_df(df=df,
                              smiles_columns=[f"{polymer_unit[0]} SMILES"],
                              features_columns=feats,
                              targets_columns=all_targets)
    dataset.set_status(graph_kernel_type='graph', features_generators=None, features_combination=None)
    dataset.create_graphs(n_jobs=4)
    dataset.unify_datatype()
    kernel_config = get_kernel_config(dataset=dataset,
                                      graph_kernel_type='graph',
                                      mgk_hyperparameters_files=[mgk_file],
                                      features_kernel_type='rbf',
                                    #   features_hyperparameters=[rbf],
                                    #   features_hyperparameters_bounds=,
                                      features_hyperparameters_file="rbf.json",
                                      hybrid_rule="product",
                                      feature_mode="per_feature",
                                      )
    gpr = GaussianProcessRegressor(kernel=kernel_config.kernel,
                                   optimizer=optimizer,
                                   alpha=0.01,
                                   normalize_y=True)
    
    X_train, X_test, y_train, y_test = train_test_split(dataset.X, dataset.y, test_size=0.2, random_state=42, shuffle=True)
    
    gpr.fit(X_train, y_train, loss=loss_function, verbose=True)
    # gpr.predict(X_test, return_std=True)

    r2_score_calc = r2_score(y_test, gpr.predict(X_test))
    rmse_calc = root_mean_squared_error(y_test, gpr.predict(X_test))
    print(f"R^2 Score: {r2_score_calc}")
    print(f"RMSE: {rmse_calc}")
