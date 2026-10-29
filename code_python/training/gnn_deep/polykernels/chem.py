"""Chemistry utilities for GSLF dataset loading, SMILES featurization, and volatility calculations."""

from collections.abc import Sequence
from typing import Optional, Tuple

from descriptastorus.descriptors import MakeGenerator
from loguru import logger
import numpy as np
import pandas as pd
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import rdkit.Chem.Crippen as Crippen
import rdkit.Chem.Descriptors as Descriptors

from polykernels.io import TensorMap
from polykernels.paths_and_constants import Constants, Paths

_P = Paths()
_C = Constants()

BURNOP_VALIDITY_ELEMENTS = {1, 6, 7, 8, 9, 17, 35, 53}  # H, C, N, O, F, Cl, Br, I
IDEAL_GAS_CONSTANT = 0.08314  # in L * bar / K / mol


def load_gslf_data() -> tuple[np.ndarray, np.ndarray]:
    """Load GSLF dataset and return smiles, labels array, and label names."""
    df = pd.read_csv(_P.data_csv)
    smiles = df[_C.smiles_col].to_numpy(dtype=str)
    labels = df[_C.labels].to_numpy(dtype=np.float32)
    return smiles, labels


def load_xy() -> tuple[np.ndarray, np.ndarray]:
    """Load precomputed features and labels for all molecules in GSLF dataset."""
    smiles, labels = load_gslf_data()
    tmap = TensorMap.load(_P.features_npz)
    features = tmap(smiles)
    return features, labels


def featurize_rdkit2d(smiles: Sequence[str] | np.ndarray) -> np.ndarray:
    """Featurize a list of SMILES using normalized 2D RDKit descriptors."""
    generator = MakeGenerator(("rdkit2dhistogramnormalized",))
    features = []
    for smi in smiles:
        res = generator.process(smi)
        if res[0] is None:
            logger.warning(f"RDKit 2D featurization failed for SMILES: {smi}")
        features.append(res[1:])
    return np.array(features, dtype=np.float32)


# =====================================================================
# Boiling Point & Volatility calculations based on Burnop/Banks methods
# =====================================================================

def is_within_atom_set(valid_elements: set[int], mol: Chem.Mol) -> bool:
    """Check if all atoms in the molecule belong to the allowed atomic numbers set."""
    return all(atom.GetAtomicNum() in valid_elements for atom in mol.GetAtoms())


def get_ring_sizes(mol: Chem.Mol) -> list[int]:
    """Get the sizes of all rings in the molecule."""
    return [len(ring) for ring in mol.GetRingInfo().AtomRings()]


def predict_boiling_point_banks(mol: Chem.Mol) -> float:
    """Estimate boiling point using Banks (1939) method based on molecular weight.

    Domains of applicability: All chemicals.

    Section 2.4.1 from "Handbook of Property Estimation Method for Chemicals".
    Link:
    https://books.google.com/books?id=48YgCAAAQBAJ&pg=PA29&lpg=PA29&dq=banks+1939+boiling+point

    Args:
        mol: RDKit molecule

    Returns:
        bp: Estimated boiling point in Kelvin.
    """
    bp = 10**(2.98 - 4. / np.sqrt(Descriptors.ExactMolWt(mol)))
    return bp


def predict_boiling_point_burnop(mol: Chem.Mol) -> float:
    """Estimate boiling point using Burnop (1938) method of weighted atom/bond features.

    Domains of applicability: All chemicals with C, H, O, N, and halides.

    Section 2.4.2 from "Handbook of Property Estimation Method for Chemicals".
    Link:
    https://books.google.com/books?id=48YgCAAAQBAJ&pg=PA29&lpg=PA29&dq=banks+1939+boiling+point

    Args:
        mol: RDKit molecule

    Returns:
        bp: Estimated boiling point in Kelvin.
    """
    if not is_within_atom_set(BURNOP_VALIDITY_ELEMENTS, mol):
        raise ValueError('Burnop method is not applicable. '
                         'SMILES: {}'.format(Chem.MolToSmiles(mol)))
    
    # Work on a copy of the molecule to avoid in-place side effects
    processed_mol = Chem.Mol(mol)
    processed_mol = Chem.AddHs(processed_mol)
    Chem.Kekulize(processed_mol, clearAromaticFlags=True)
    
    atoms = processed_mol.GetAtoms()
    bonds = processed_mol.GetBonds()
    ring_sizes = get_ring_sizes(processed_mol)

    num_c = sum(23.2 for atom in atoms if atom.GetAtomicNum() == 6)
    num_h = sum(10.9 for atom in atoms if atom.GetAtomicNum() == 1)
    num_o = sum(51.0 for atom in atoms if atom.GetAtomicNum() == 8)
    num_n = sum(39.7 for atom in atoms if atom.GetAtomicNum() == 7)
    num_f = sum(68.0 for atom in atoms if atom.GetAtomicNum() == 9)
    num_cl = sum(121.0 for atom in atoms if atom.GetAtomicNum() == 17)
    num_br = sum(255.0 for atom in atoms if atom.GetAtomicNum() == 35)
    num_i = sum(398.0 for atom in atoms if atom.GetAtomicNum() == 53)
    
    num_6_member_ring = sum(17.6 for ring in ring_sizes if ring == 6)
    num_double_bond = sum(
        16.1 for bond in bonds
        if bond.GetBondType() == Chem.BondType.DOUBLE
    )
    num_triple_bond = sum(
        33.0 for bond in bonds
        if bond.GetBondType() == Chem.BondType.TRIPLE
    )

    weighted_features = (
        num_c + num_h + num_o + num_n + num_f + num_cl + num_br + num_i +
        num_6_member_ring + num_double_bond + num_triple_bond
    )
    molecular_weight = Descriptors.MolWt(mol)
    bp = 10**((weighted_features / molecular_weight - 8) / np.sqrt(molecular_weight))
    return bp


def predict_boiling_point(mol: Chem.Mol) -> float:
    """Estimate boiling point using the most specific applicable estimator.

    Falls back to predict_boiling_point_banks if the molecule contains
    unsupported elements not covered by Burnop's domain of applicability.

    Args:
        mol: RDKit molecule

    Returns:
        Estimated boiling point in Kelvin.
    """
    if is_within_atom_set(BURNOP_VALIDITY_ELEMENTS, mol):
        return predict_boiling_point_burnop(mol)
    return predict_boiling_point_banks(mol)


def predict_has_odor_debug(
    mol: Optional[Chem.Mol],
    exp_logp: Optional[float] = None,
    exp_bp: Optional[float] = None,
) -> Tuple[bool, float, float]:
    """Variant of predict_has_odor that also returns values used in decision."""
    if mol is None and (exp_bp is None or exp_logp is None):
        raise ValueError('`mol` is needed if experimental values are not provided.')
    bp = exp_bp or predict_boiling_point(mol)
    logp = exp_logp or Crippen.MolLogP(mol)
    has_odor = (bp - 273.15) < 300 and logp > 0
    return has_odor, logp, bp


def predict_has_odor(
    mol: Optional[Chem.Mol] = None,
    exp_logp: Optional[float] = None,
    exp_bp: Optional[float] = None,
) -> bool:
    """Apply a simple threshold model to predict if a molecule has an odor.

    Odor detection threshold is based on predicted boiling point in Celsius
    and octanol-water partition coefficient (LogP).

    Args:
        mol: Molecule to predict odor. Required if experimental logp/bp are not provided.
        exp_logp: Experimental partition coefficient (LogP), if available.
        exp_bp: Experimental boiling point in Kelvin, if available.

    Returns:
        bool: True if the molecule is predicted to have an odor, False otherwise.
    """
    has_odor, _, _ = predict_has_odor_debug(mol, exp_logp=exp_logp, exp_bp=exp_bp)
    return has_odor


def get_vapor_mass_g(
    mol: Chem.Mol,
    pressure_bar: float,
    volume_mL: float,
    temperature_K: float = 298.15,
) -> float:
    """Compute saturated vapor mass using the Ideal Gas Law.

    Uses PV = nRT to determine the number of moles and multiplies by molecular weight.

    Args:
        mol: RDKit molecule
        pressure_bar: Vapor pressure in bar
        volume_mL: Gas phase volume in mL
        temperature_K: Temperature in Kelvin (defaults to 298.15 K)

    Returns:
        Vapor mass in grams.
    """
    volume_L = volume_mL / 1000.0
    moles = (pressure_bar * volume_L) / (IDEAL_GAS_CONSTANT * temperature_K)
    mass = moles * Descriptors.MolWt(mol)
    return mass


def get_robust_geometry(
    mol: Chem.Mol,
    use_3d: bool = True,
    max_tries: int = 50,
) -> Chem.Mol:
    """Add 3D/2D coordinates to a molecule in a robust way.

    Molecules constructed from SMILES do not include all 3D geometry information.
    To get 3D coordinates we utilize RDKit's EmbedMolecule method which is based
    on doi:10.1021/acs.jcim.0c00025 and doi:10.1021/acs.jcim.5b00654.
    If standard embedding fails, several empirical strategies are tried.

    See: https://www.rdkit.org/docs/RDKit_Book.html
    See: https://www.rdkit.org/docs/GettingStartedInPython.html

    Args:
        mol: Molecule to add a geometric conformer.
        use_3d: If False, will compute 2D coordinates instead.
        max_tries: Number of attempts for embedding.

    Returns:
        RDKit molecule with coordinate conformer.
    """
    mol_copy = Chem.Mol(mol)
    if not use_3d:
        AllChem.Compute2DCoords(mol_copy)
        return mol_copy

    # Track if the input had explicit Hs to match original hydrogen representation
    has_explicit_hs = any(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms())
    
    mol_copy = Chem.AddHs(mol_copy)
    params = AllChem.ETKDGv3()
    params.useSmallRingTorsions = True
    params.maxIterations = max_tries
    params.randomSeed = 42

    res = -1
    try:
        res = AllChem.EmbedMolecule(mol_copy, params=params)
    except Exception:
        pass

    if res < 0:
        # Fallback 1: Try ETKDGv3 but turn off basic knowledge and experimental torsion preferences
        # (highly helpful for strained bridged bicyclic ring systems)
        try:
            params.useBasicKnowledge = False
            params.useExpTorsionAnglePrefs = False
            res = AllChem.EmbedMolecule(mol_copy, params=params)
        except Exception:  # noqa: BLE001
            pass

    if res < 0:
        # Fallback 2: Try embedding with random coordinates
        try:
            res = AllChem.EmbedMolecule(mol_copy, useRandomCoords=True, randomSeed=42)
        except Exception:
            pass

    if res >= 0:
        try:
            if AllChem.MMFFHasAllMoleculeParams(mol_copy):
                AllChem.MMFFOptimizeMolecule(mol_copy)
        except Exception:
            pass
    else:
        AllChem.Compute2DCoords(mol_copy)

    if not has_explicit_hs:
        mol_copy = Chem.RemoveHs(mol_copy)

    return mol_copy
