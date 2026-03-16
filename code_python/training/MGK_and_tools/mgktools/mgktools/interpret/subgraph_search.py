#!/usr/bin/env python
# -*- coding: utf-8 -*-
from typing import List
from tqdm import tqdm
import numpy as np
import pandas as pd
import rdkit.Chem.AllChem as Chem
from mgktools.interpret.subgraph import get_fragments_smarts
from mgktools.interpret.interpret import interpret_substructures


def substructure_search(mols_sorted: List[Chem.Mol],
                        n_mols: int = 5,
                        nb_nodes_min: int = 1,
                        nb_nodes_max: int = 5,
                        n_occurrence: int = 50,
                        maximize: bool = True,
                        contribution_per_atom: bool = True):
    """Find all possible substructures, and rank them according to

    Parameters
    ----------
    mols_sorted: list of RDKit Mol object.
        interpreted RDKit Mol objects, make sure it is ordered by the target value, because only the first-n_mols
        molecules' substructures will be considered.
    n_mols: int
        the number of molecules' substructures to be considered.
    nb_nodes_min: int
        the minimum number of nodes of a substructure.
    nb_nodes_max: int
        the maximum number of nodes of a substructure.
    n_occurrence: int
        only substructures that occur more than this value will be considered.
    maximize: bool
        if True, try to find substructures with maximum interpreted value.
        If False, try to find substructures with minimum interpreted value.
    contribution_per_atom: bool
        If True, the output will be ordered by interpreted contribution per atom.
        If False, the output will be ordered by interpreted contribution.

    Returns
    -------
    A dataframe contains substructures and its interpretation.
    """
    df = pd.DataFrame({'smarts': [],
                       'mean_contribution': [],
                       'std_contribution': [],
                       'mean_contribution_per_atom': [],
                       'contributions': []})
    for i in tqdm(range(n_mols), total=n_mols):
        mol = mols_sorted[i]
        fragments_smarts = get_fragments_smarts(mol, nb_nodes_min=nb_nodes_min, nb_nodes_max=nb_nodes_max)
        for smarts in fragments_smarts:
            if smarts not in df['smarts'].tolist():
                contributions = interpret_substructures(smarts, mols_sorted)
                n_heavy = Chem.MolFromSmarts(smarts).GetNumAtoms()
                if len(contributions) >= n_occurrence:
                    df.loc[len(df)] = smarts, \
                                      np.mean(contributions), \
                                      np.std(contributions), \
                                      np.mean(contributions) / n_heavy, \
                                      contributions
    if contribution_per_atom:
        return df.sort_values('mean_contribution_per_atom', ascending=not maximize)
    else:
        return df.sort_values('mean_contribution', ascending=not maximize)
