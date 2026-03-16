#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd
import pytest
from mgktools.data.data import Dataset
from mgktools.interpret.interpret import interpret_training_mols, interpret_atoms, get_interpreted_mols
from mgktools.hyperparameters import additive_pnorm, product_pnorm

smiles = ['CN(C)C(=O)c1ccc(cc1)OC', 'CS(=O)(=O)Cl', 'CC(C)C=C', 'CCc1cnccn1']
targets = [-11.01, -4.87, 1.83, -5.45]
df = pd.DataFrame({'smiles': smiles, 'targets_1': targets})


@pytest.mark.parametrize('testset', [
    (additive_pnorm),
    (product_pnorm),
])
def test_interpret_training_mols(testset):
    mgk_hyperparameters_file = testset
    y_pred, y_std, df_interpret = interpret_training_mols(
        smiles_to_be_interpret=['Cc1cc(cc(c1)O)C', 'CCC'],
        smiles_train=smiles,
        targets_train=targets,
        alpha=0.01,
        n_mol=10,
        output_order='sort_by_value',
        mgk_hyperparameters_file=mgk_hyperparameters_file,
        n_jobs=6)
    for i, df in enumerate(df_interpret):
        assert df['contribution_value'].sum() == pytest.approx(y_pred[i], 1e-4)


@pytest.mark.parametrize('testset', [
    (additive_pnorm),
    (product_pnorm),
])
def test_interpret_atoms(testset):
    mgk_hyperparameters_file = testset
    y_pred, y_std, mol = interpret_atoms(
        smiles_to_be_interpret='Cc1cc(cc(c1)O)C',
        smiles_train=smiles,
        targets_train=targets,
        alpha=0.01,
        mgk_hyperparameters_file=mgk_hyperparameters_file)
    y_sum = 0.
    for atom in mol.GetAtoms():
        y_sum += float(atom.GetProp('atomNote'))
    assert y_sum == pytest.approx(y_pred, 1e-5)


@pytest.mark.parametrize('testset', [
    (additive_pnorm),
    (product_pnorm),
])
@pytest.mark.parametrize('batch_size', [(1), (2), (3)])
def test_get_interpreted_mols(testset, batch_size):
    mgk_hyperparameters_file = testset
    batch_size = batch_size
    smiles_to_be_interpret = ['Cc1cc(cc(c1)O)C', 'CC(C)C(C)C', 'CCCO', 'C1CCCC1CCO']
    y_pred, y_std, mols = get_interpreted_mols(smiles_train=smiles,
                                               targets_train=targets,
                                               smiles_to_be_interpret=smiles_to_be_interpret,
                                               mgk_hyperparameters_file=mgk_hyperparameters_file,
                                               alpha=0.01,
                                               return_mols_only=False,
                                               batch_size=batch_size)
    for i, mol in enumerate(mols):
        y_sum = 0.
        for atom in mol.GetAtoms():
            y_sum += float(atom.GetProp('atomNote'))
        assert y_sum == pytest.approx(y_pred[i], 1e-5)
