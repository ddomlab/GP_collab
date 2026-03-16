#!/usr/bin/env python
# -*- coding: utf-8 -*-
from typing import List, Tuple
import pickle
import os
import rdkit.Chem.AllChem as Chem
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
from IPython.display import display


def save_mols_pkl(mols: List[Chem.Mol], path, filename='mols.pkl'):
    f_mols = os.path.join(path, filename)
    atomic_attribution = []
    for mol in mols:
        atomNote = []
        for atom in mol.GetAtoms():
            atomNote.append(atom.GetProp('atomNote'))
        atomic_attribution.append(atomNote)
    with open(f_mols, "wb") as f:
        pickle.dump([mols, atomic_attribution], f)


def load_mols_pkl(path, filename='mols.pkl'):
    f_mols = os.path.join(path, filename)
    with open(f_mols, "rb") as f:
        mols, atomic_attribution = pickle.load(f)
    for i, mol in enumerate(mols):
        for j, atom in enumerate(mol.GetAtoms()):
            atom.SetProp('atomNote', atomic_attribution[i][j])
    return mols


def display_mol(mol: Chem.Mol, highlight_threshold_upper: float = None, highlight_threshold_lower: float = None,
                highlight_color_upper: Tuple[float, float, float] = (1, 0, 0), 
                highlight_color_lower: Tuple[float, float, float] = (0, 1, 0), size: Tuple[float, float] = (500, 500),
                remove_number: bool = False
                ):
    mol_copy = Chem.Mol(mol)
    d = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
    highlight_atoms_upper = []
    highlight_atoms_lower = []
    highlight_bonds = []
    atom_color = {}
    for atom in mol_copy.GetAtoms():
        if highlight_color_upper is not None and float(atom.GetProp('atomNote')) > highlight_threshold_upper:
            highlight_atoms_upper.append(atom.GetIdx())
            atom_color[atom.GetIdx()] = highlight_color_upper
        elif highlight_color_lower is not None and float(atom.GetProp('atomNote')) < highlight_threshold_lower:
            highlight_atoms_lower.append(atom.GetIdx())
            atom_color[atom.GetIdx()] = highlight_color_lower
    for bond in mol_copy.GetBonds():
        if bond.GetBeginAtomIdx() in highlight_atoms_upper and bond.GetEndAtomIdx() in highlight_atoms_upper:
            highlight_bonds.append(bond.GetIdx())
        elif bond.GetBeginAtomIdx() in highlight_atoms_lower and bond.GetEndAtomIdx() in highlight_atoms_lower:
            highlight_bonds.append(bond.GetIdx())
    if remove_number:
        for atom in mol_copy.GetAtoms():
            atom.SetProp('atomNote', '')
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol_copy, highlightAtoms=highlight_atoms_upper + highlight_atoms_lower,
                                       highlightBonds=highlight_bonds, highlightAtomColors=atom_color)
    d.FinishDrawing()
    display(SVG(d.GetDrawingText()))
