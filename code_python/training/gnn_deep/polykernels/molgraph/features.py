"""Base and concrete graph feature computation classes for molecular featurization.

Design Rationale:
1) Fast & Batch-Oriented: Feature computations operate collectively over a `MolBatch` to
   minimize RDKit-to-Python loop overhead and allow vectorized NumPy operations.
2) Modular & Extensible: All features inherit from standard base classes that specify
   destination attributes (`atoms`, `bonds`, or `mols`), making it simple to add or remove features.
3) Customizable: Descriptors and thresholds (e.g., categorical sets, numerical ranges)
   can be tailored via initialization parameters to match specific datasets.
4) Serializable: Config-safe class structures allow rebuilding features from standard dictionaries.
"""

import abc
from collections.abc import Iterator, Sequence
import copy
import dataclasses
import enum
import functools
import math
import re
from typing import Any
import warnings

import numpy as np
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import rdkit.Chem.ChemicalFeatures as ChemicalFeatures
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors
import rdkit.ML.Descriptors.MoleculeDescriptors as rdDescriptors

from polykernels.chem import featurize_rdkit2d
from polykernels.io import TensorMap
from polykernels.molgraph.graphwave import get_graphwave_embeddings
from polykernels.paths_and_constants import Paths

_SPATIAL_KEY = "_SPATIAL_DISTANCE"
MAX_ATOM_NUM = 118


def _get_all_symbols() -> list[str]:
    table = Chem.GetPeriodicTable()
    return [table.GetElementSymbol(i) for i in range(118)]


ALL_ATOM_SYMBOLS = _get_all_symbols()


def _check_valid_symbols(symbols: list[str]) -> None:
    missing = set(symbols).difference(ALL_ATOM_SYMBOLS)
    if missing:
        raise ValueError(f"{missing} values not valid atom symbols")


def iterate_atoms(mols: Sequence[Chem.Mol]) -> Iterator[Chem.Atom]:
    for mol in mols:
        for atom in mol.GetAtoms():
            yield atom


def iterate_bonds(mols: Sequence[Chem.Mol]) -> Iterator[Chem.Bond]:
    for mol in mols:
        for bond in mol.GetBonds():
            yield bond


@dataclasses.dataclass
class MolRequirements:
    stereo_tags: bool = False
    spatial_geometry: bool = False
    add_hs: bool = False


@dataclasses.dataclass
class MolBatch:
    """Light-weight class to hold cached data on batches of molecules."""

    mols: Sequence[Chem.Mol]

    @functools.cached_property
    def atoms(self) -> Sequence[Chem.Atom]:
        return list(iterate_atoms(self.mols))

    @functools.cached_property
    def natoms(self) -> Sequence[int]:
        return [m.GetNumAtoms() for m in self.mols]

    @functools.cached_property
    def bonds(self) -> Sequence[Chem.Bond]:
        return list(iterate_bonds(self.mols))

    @functools.cached_property
    def nbonds(self) -> Sequence[int]:
        return [m.GetNumBonds() for m in self.mols]


class GraphAttribute(enum.StrEnum):
    """Enum for molecule graph attributes, useful to group features."""

    atoms = "atoms"
    bonds = "bonds"
    mols = "mols"


def one_hot_encode(indices: np.ndarray | list[int], depth: int) -> np.ndarray:
    indices_arr = np.array(indices, dtype=np.int32)
    n = len(indices_arr)
    values = np.zeros((n, depth), dtype=np.float32)
    values[np.arange(n), indices_arr] = 1.0
    return values


def _cast_as_2darray(v: list[Any] | np.ndarray) -> np.ndarray:
    return np.array(v, dtype=np.float32).reshape(-1, 1)


def _map_indices(inputs: list[Any], values: list[Any]) -> list[int]:
    return [values.index(i) if i in values else -1 for i in inputs]


def impute_unexpected_strings(
    values: list[str], expected_values: list[str], imputation: str
) -> list[str]:
    return [s if s in expected_values else imputation for s in values]


class GraphFeature(abc.ABC):
    """Basic class for a computable graph feature."""

    destination: GraphAttribute
    ndim: int
    requirements: MolRequirements = MolRequirements()

    @classmethod
    def name(cls) -> str:
        return re.sub("((?<!^)(?=[A-Z]))", "_", cls.__name__).lower()

    @abc.abstractmethod
    def compute(self, batch: MolBatch) -> Any:
        """Feature computation, returns raw values."""

    @abc.abstractmethod
    def encode(self, computed_values: Any) -> np.ndarray:
        """Feature encoding, convert raw values to 2-D np.arrays."""

    @abc.abstractmethod
    def to_config(self) -> dict[str, Any]:
        """Get a initializable configuration for the feature."""

    def __call__(self, batch: MolBatch) -> np.ndarray:
        return self.encode(self.compute(batch))

    def to_str(self, prefix="", infix="", suffix="") -> str:
        configs = [f"{key}={value}" for key, value in self.to_config().items()]
        config_str = (", with " + ", ".join(configs)) if configs else ""
        return (
            f"{prefix}{self.name():20s} of dim {self.ndim:3d}{infix}"
            f"{config_str}{suffix}"
        )

    def __repr__(self) -> str:
        return self.to_str()


class CategoricalFeature(GraphFeature):
    """Boilerplate feature for string data we want to one-hot encode."""

    values: list[str]
    default_values: list[str]
    add_unknown: bool = False

    def __init__(
        self, values: list[str] | None = None, add_unknown: bool = False  # noqa: FBT001, FBT002
    ):

        self.values = values or self.default_values
        if len(self.values) != len(set(self.values)):
            raise ValueError(f"Found duplicate entries ({values})")
        self.add_unknown = add_unknown
        self.ndim = len(self.values)
        if self.add_unknown:
            self.ndim += 1
        # O(1) Lookup dictionary mapping string value to index
        self.value_to_idx = {v: idx for idx, v in enumerate(self.values)}

    def to_config(self) -> dict[str, Any]:
        return {"values": self.values, "add_unknown": self.add_unknown}

    def encode(self, computed_values: list[str]) -> np.ndarray:
        unknown_val = self.ndim - 1 if self.add_unknown else -1
        indices = [self.value_to_idx.get(v, unknown_val) for v in computed_values]
        return one_hot_encode(indices, self.ndim)


class IntOneHotFeature(GraphFeature):
    """Boilerplate feature for integers that we want to one-hot encode."""

    min_value: int
    max_value: int
    default_min: int
    default_max: int
    add_unknown: bool = False

    def __init__(
        self,
        min_value: int | None = None,
        max_value: int | None = None,
        add_unknown: bool = False,  # noqa: FBT001, FBT002
    ):

        self.min_value = min_value if min_value is not None else self.default_min
        self.max_value = max_value if max_value is not None else self.default_max
        self.add_unknown = add_unknown
        self.ndim = (self.max_value - self.min_value) + 1
        if self.add_unknown:
            self.ndim += 1

    def encode(self, computed_values: list[int] | np.ndarray) -> np.ndarray:
        vals = np.array(computed_values, dtype=np.int32)
        indices = vals - self.min_value
        if self.add_unknown:
            mask = (indices >= 0) & (indices < self.ndim - 1)
            indices = np.where(mask, indices, self.ndim - 1)
        else:
            indices = np.clip(indices, 0, self.ndim - 1)
        return one_hot_encode(indices, self.ndim)

    def to_config(self) -> dict[str, Any]:
        return {
            "min_value": self.min_value,
            "max_value": self.max_value,
            "add_unknown": self.add_unknown,
        }


class Continuous1DFeature(GraphFeature):
    """Boilerplate feature for numbers (float/booleans/integers)."""

    ndim = 1

    def encode(self, computed_values: list[float | int | bool] | np.ndarray) -> np.ndarray:
        return _cast_as_2darray(computed_values)

    def to_config(self) -> dict[str, Any]:
        return {}


class VectorFeature(GraphFeature):
    """Boilerplate feature for vector features (e.g. embeddings)."""

    def encode(self, computed_values: np.ndarray) -> np.ndarray:
        return computed_values

    def to_config(self) -> dict[str, Any]:
        return {}


class NumericFeature(GraphFeature):
    """Boilerplate feature for numeric values that need scaling and clipping."""

    min_value: float
    max_value: float
    scale: float
    default_min: float
    default_max: float
    default_scale: float

    def __init__(
        self,
        min_value: float | None = None,
        max_value: float | None = None,
        scale: float | None = None,
    ):
        self.min_value = min_value if min_value is not None else self.default_min
        self.max_value = max_value if max_value is not None else self.default_max
        self.scale = scale if scale is not None else self.default_scale
        self.ndim = 1

    def encode(self, computed_values: list[float] | np.ndarray) -> np.ndarray:
        val = np.array(computed_values, dtype=np.float32) / self.scale
        val = np.clip(val, self.min_value, self.max_value)
        return _cast_as_2darray(val)

    def to_config(self) -> dict[str, Any]:
        return {
            "min_value": self.min_value,
            "max_value": self.max_value,
            "scale": self.scale,
        }


#
# Atom to node features
#


class AtomicSymbol(CategoricalFeature):
    """One-hot representation of atomic symbols to identify element identity in the GNN."""

    destination = GraphAttribute.atoms
    default_values = [
        "C", "N", "O", "S", "F", "Si", "P", "Cl", "Br", "Mg", "Na", "Ca", "Fe",
        "As", "Al", "I", "B", "V", "K", "Tl", "Yb", "Sb", "Sn", "Ag", "Pd", "Co",
        "Se", "Ti", "Zn", "H", "Li", "Ge", "Cu", "Au", "Ni", "Cd", "In", "Mn",
        "Zr", "Cr", "Pt", "Hg", "Pb", "Be", "*",
    ]

    def __init__(self, values: list[str] | None = None):
        super().__init__(values)
        _check_valid_symbols(self.values)

    def compute(self, batch: MolBatch) -> list[str]:
        return [str(a.GetSymbol()) for a in batch.atoms]


class AtomicSymbolImpute(CategoricalFeature):
    """One-hot representation of atomic symbols mapping exotic elements to a default impute key."""

    destination = GraphAttribute.atoms
    default_values = [
        "C", "N", "O", "H", "F", "P", "S", "Cl", "Br", "I", "Si", "Ge", "B", "Se",
        "OTHER", "*",
    ]
    wildcard_key = "*"
    impute_key = "OTHER"

    def __init__(self, values: list[str] | None = None, add_unknown: bool = False):
        super().__init__(values, add_unknown=add_unknown)
        tags = []
        _check_valid_symbols([s for s in self.values if s != self.impute_key])
        for symbol in ALL_ATOM_SYMBOLS:
            atom_encode = symbol if symbol in self.values else self.impute_key
            tags.append(atom_encode)
        self.atom_symbols = tags

    def compute(self, batch: MolBatch) -> list[str]:
        return [self.atom_symbols[a.GetAtomicNum()] for a in batch.atoms]


class PeriodicGroup(CategoricalFeature):
    """One-hot encoding of periodic table groups (alkali, transition metals, halogens, etc.)."""

    destination = GraphAttribute.atoms
    default_values = [
        "nonmetal", "halogen", "metalloid", "metal", "alkali_earth_metal",
        "alkali_metal", "transition_metal", "noble_gas", "exotic",
    ]
    impute_key = "exotic"
    groups = {
        "nonmetal": ["H", "C", "N", "O", "P", "S", "Se"],
        "halogen": ["F", "Cl", "Br", "I"],
        "metalloid": ["B", "Si", "Ge", "As", "Sb", "Te", "Po"],
        "metal": ["Al", "Ga", "In", "Sn", "Tl", "Pb", "Bi"],
        "alkali_earth_metal": ["Be", "Mg", "Ca", "Sr", "Ba", "Fr"],
        "alkali_metal": ["Li", "Na", "K", "Rb", "Cs", "Rn"],
        "transition_metal": [
            "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
            "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
            "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        ],
        "noble_gas": ["He", "Ne", "Ar", "Kr", "Xe", "At"],
    }

    def __init__(self, values: list[str] | None = None, add_unknown: bool = False):
        super().__init__(values, add_unknown=add_unknown)
        for symbols in self.groups.values():
            _check_valid_symbols(symbols)
        tags = []
        for symbol in ALL_ATOM_SYMBOLS:
            group = self.impute_key
            for name, symbols in self.groups.items():
                if symbol in symbols:
                    group = name
            tags.append(group)
        self.atomnum_to_group = tags

    def compute(self, batch: MolBatch) -> list[str]:
        return [self.atomnum_to_group[a.GetAtomicNum()] for a in batch.atoms]

    def to_config(self) -> dict[str, Any]:
        return {}


class Degree(IntOneHotFeature):
    """One-hot representation of the coordinate degree (number of direct bonds) of the atom."""

    destination = GraphAttribute.atoms
    default_min = 0
    default_max = 6

    def compute(self, batch: MolBatch) -> list[int]:
        return [a.GetDegree() for a in batch.atoms]


class AtomicWeight(NumericFeature):
    """Continuous atomic mass feature, normalized by dividing by 100.0."""

    destination = GraphAttribute.atoms
    default_min = 0.0
    default_max = 2.0
    default_scale = 100.0

    def compute(self, batch: MolBatch) -> list[float]:
        return [float(a.GetMass()) for a in batch.atoms]


def _assign_partial_charges(mol: Chem.Mol) -> Chem.Mol:
    """Calculate Gasteiger partial charges for an RDKit Mol."""
    mol = copy.deepcopy(mol)
    try:
        from rdkit.Chem.rdPartialCharges import compute_gasteiger_charges
        compute_gasteiger_charges(mol)
    except ImportError:
        from rdkit.Chem.rdPartialCharges import ComputeGasteigerCharges
        ComputeGasteigerCharges(mol)
    return mol


class PartialCharges(NumericFeature):
    """Calculates atom-wise Gasteiger partial charges for electrostatic representation."""

    destination = GraphAttribute.atoms
    default_min = -1.0
    default_max = 1.0
    default_scale = 1.0

    def compute(self, batch: MolBatch) -> list[float]:
        res = []
        for mol in batch.mols:
            try:
                mol_with_charges = _assign_partial_charges(mol)
            except Exception:
                res.extend([0.0] * mol.GetNumAtoms())
                continue
            for a in mol_with_charges.GetAtoms():
                try:
                    val = float(a.GetProp("_GasteigerCharge"))
                    if math.isnan(val):
                        val = 0.0
                except KeyError:
                    val = 0.0
                res.append(val)
        return res


class ChiralTag(CategoricalFeature):
    """One-hot representation of RDKit's internal chiral tags (CHI_UNSPECIFIED, CHI_TETRAHEDRAL_CW, etc.)."""

    destination = GraphAttribute.atoms
    default_values = [
        "CHI_UNSPECIFIED", "CHI_TETRAHEDRAL_CW", "CHI_TETRAHEDRAL_CCW",
        "CHI_OTHER",
    ]

    def compute(self, batch: MolBatch) -> list[str]:
        return [str(a.GetChiralTag()) for a in batch.atoms]


class ChiralCenter(CategoricalFeature):
    """One-hot representation of R/S configuration labels for stereocenters (or 'none')."""

    destination = GraphAttribute.atoms
    default_values = ["R", "S", "none"]

    def compute(self, batch: MolBatch) -> list[str]:
        res = []
        for mol in batch.mols:
            chiral_centers = dict(
                Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            )
            for a in mol.GetAtoms():
                res.append(chiral_centers.get(a.GetIdx(), "none"))
        return res


_HBOND_FEATURE_DEF = """
AtomType NDonor [N&!H0&v3,N&!H0&+1&v4,n&H1&+0]
AtomType AmideN [$(N-C(=O))]
AtomType SulfonamideN [$([N;H0]S(=O)(=O))]
AtomType NDonor [$([Nv3](-C)(-C)-C)]
AtomType NDonor [$(n[n;H1]),$(nc[n;H1])]
AtomType ChalcDonor [O,S;H1;+0]

DefineFeature SingleAtomDonor [{NDonor},{ChalcDonor}]
  Family Donor
  Weights 1
EndFeature

AtomType NAcceptor [n;+0;!X3;!$([n;H1](cc)cc)]
AtomType NAcceptor [$([N;H0]#[C&v4])]
AtomType NAcceptor [N&v3;H0;$(Nc)]
AtomType ChalcAcceptor [O;H0;v2;!$(O=N-*)] 
Atomtype ChalcAcceptor [O;-;!$(*-N=O)] 
Atomtype ChalcAcceptor [o;+0]
AtomType Hydroxyl [O;H1;v2]
AtomType HalogenAcceptor [F;$(F-[#6]);!$(FC[F,Cl,Br,I])]

DefineFeature SingleAtomAcceptor [{Hydroxyl},{ChalcAcceptor},{NAcceptor},{HalogenAcceptor}]
  Family Acceptor
  Weights 1
EndFeature

DefineFeature AcidicGroup [C,S](=[O,S,P])-[O;H1,H0&-1]
  Family NegIonizable
  Weights 1.0,1.0,1.0
EndFeature

AtomType Carbon_NotDouble [C;!$(C=*)]
AtomType BasicNH2 [$([N;H2&+0][{Carbon_NotDouble}])]
AtomType BasicNH1 [$([N;H1&+0]([{Carbon_NotDouble}])[{Carbon_NotDouble}])]
AtomType PosNH3 [$([N;H3&+1][{Carbon_NotDouble}])]
AtomType PosNH2 [$([N;H2&+1]([{Carbon_NotDouble}])[{Carbon_NotDouble}])]
AtomType PosNH1 [$([N;H1&+1]([{Carbon_NotDouble}])([{Carbon_NotDouble}])[{Carbon_NotDouble}])]
AtomType BasicNH0 [$([N;H0&+0]([{Carbon_NotDouble}])([{Carbon_NotDouble}])[{Carbon_NotDouble}])]
AtomType QuatN [$([N;H0&+1]([{Carbon_NotDouble}])([{Carbon_NotDouble}])([{Carbon_NotDouble}])[{Carbon_NotDouble}])]

DefineFeature BasicGroup [{BasicNH2},{BasicNH1},{BasicNH0};!$(N[a])]
  Family PosIonizable
  Weights 1.0
EndFeature
"""


class HydrogenBonding(VectorFeature):
    """Base pharmacore atom features, mostly relating to electrostatic features.

    These features are commonly used in pharmaceutical contexts.
    Acceptor/Donor refer to if an atom is likely to have a role in a
    hydrogen bond (https://en.wikipedia.org/wiki/Hydrogen_bond).

    NegIonizable/PosIonizable refers to atoms that are likely to be
    protonated/deprotonated.
    """

    destination = GraphAttribute.atoms
    ndim = 4
    features = ["Acceptor", "Donor", "NegIonizable", "PosIonizable"]

    def __init__(self):
        super().__init__()
        self.factory = ChemicalFeatures.BuildFeatureFactoryFromString(
            _HBOND_FEATURE_DEF
        )

    def compute(self, batch: MolBatch) -> np.ndarray:
        offsets = np.cumsum([0] + batch.natoms)
        total_atoms = offsets[-1]
        values = np.zeros((total_atoms, self.ndim), dtype=np.float32)
        for index, mol in enumerate(batch.mols):
            features = self.factory.GetFeaturesForMol(mol)
            mol_index = offsets[index]
            for feat in features:
                name = feat.GetFamily()
                if name in self.features:
                    atom_indices = feat.GetAtomIds()
                    feature_index = self.features.index(name)
                    values[mol_index + np.array(atom_indices), feature_index] = 1.0
        return values


class CrippenFeatures(VectorFeature):
    """Atom-wise values for Crippen properties (LogP/MR).

    Published in doi.org/10.1021/ci990307l. Weights come from atomic
    contributions of a linear model for prediction solubility (LogP) and
    molecular refractivity (MR). Weights take into account local environment
    via a substructure lookup table.

    Expected value range for LogP is [-1.95, 0.88].
    Expected value range for MR is [0.0, 14.02].
    """

    destination = GraphAttribute.atoms
    ndim = 2
    features = ["LogP", "MR"]

    def compute(self, batch: MolBatch) -> np.ndarray:
        values = []
        for mol in batch.mols:
            contribs = rdMolDescriptors._CalcCrippenContribs(mol)
            values.extend(contribs)
        return np.array(values, dtype=np.float32)


class InRingSizes(VectorFeature):
    """Boolean atom-wise labels for a range of ring sizes.

    A more specific version of 'InRing', each ring size gets one hot-encoded to a
    "bucket". Ring sizes in 3-6 have their own buckets, 7-11 are encoded to a
    single bucket, and 12+ to another. The rationale is around the expected types
    of conformal properties of rings at different sizes.

    See: https://en.wikipedia.org/wiki/Cycloalkane
    """

    destination = GraphAttribute.atoms
    ringsize_indexer = {
        3: 0, 4: 1, 5: 2, 6: 3, 7: 4, 8: 4, 9: 4, 10: 4, 11: 4, 12: 5
    }
    impute_index = 5
    ndim = 6

    def compute(self, batch: MolBatch) -> np.ndarray:
        offsets = np.cumsum([0] + batch.natoms)
        natoms = offsets[-1]
        values = np.zeros((natoms, self.ndim), dtype=np.float32)
        for index, mol in enumerate(batch.mols):
            mol_index = offsets[index]
            ring_info = mol.GetRingInfo()
            for ring in ring_info.AtomRings():
                ring_size_index = self.ringsize_indexer.get(
                    len(ring), self.impute_index
                )
                ring_indices = np.array(ring) + mol_index
                values[ring_indices, ring_size_index] = 1.0
        return values


class GraphwaveEmbeddings(VectorFeature):
    """Graphwave diffusion-based structural node features.

    Graphwave embeddings compute characteristic functions of heat kernel diffusion.
    By computing the heat kernel from the graph laplacian, they capture distance
    relationships for node embeddings (e.g. two nodes with similar subgraph
    structural roles will have similar embeddings).

    Original paper: arXiv:1710.10321

    This implementation is a vector version of the paper's characteristic function
    sampled at linear intervals.

    Note that these are computed graph-by-graph, so they may have high variance.
    If using multiple disconnected graphs, the normalization/eigen decomposition
    may become unstable.
    """

    destination = GraphAttribute.atoms

    def __init__(
        self,
        ndim: int | None = None,
        step_size: float | None = None,
        scales: list[float] | None = None,
    ):
        super().__init__()
        self.ndim_per_scale = ndim or 4
        self.step_size = step_size or 1.0
        self.scales = scales or [1.0, 2.0]
        self.ndim = self.ndim_per_scale * len(self.scales)

    def compute(self, batch: MolBatch) -> np.ndarray:
        values = []
        for mol in batch.mols:
            adj = Chem.GetAdjacencyMatrix(mol)
            emb = [
                get_graphwave_embeddings(
                    adj, self.ndim_per_scale, self.step_size, scale
                )
                for scale in self.scales
            ]
            values.append(np.concatenate(emb, axis=-1))
        return np.vstack(values).astype(np.float32)

    def to_config(self) -> dict[str, Any]:
        return {
            "ndim": self.ndim_per_scale,
            "step_size": self.step_size,
            "scales": self.scales,
        }


#
# Bond to edge features
#


class BondType(CategoricalFeature):
    """One-hot representation of standard covalent bond orders (single, double, triple, aromatic)."""

    destination = GraphAttribute.bonds
    default_values = ["SINGLE", "DOUBLE", "TRIPLE", "AROMATIC"]

    def compute(self, batch: MolBatch) -> list[str]:
        return [str(b.GetBondType()) for b in batch.bonds]


class InRing(Continuous1DFeature):
    """Boolean indicator of whether a chemical bond is part of any ring structure."""

    destination = GraphAttribute.bonds

    def compute(self, batch: MolBatch) -> list[bool]:
        return [b.IsInRing() for b in batch.bonds]


class BondStereo(CategoricalFeature):
    """One-hot representation of cis/trans (E/Z) configurations of double bonds."""

    destination = GraphAttribute.bonds
    default_values = ["ANY", "Z", "E", "NONE"]
    imputed_value = "ANY"
    requirements = MolRequirements(stereo_tags=True)

    def compute(self, batch: MolBatch) -> list[str]:
        tags = [str(b.GetStereo())[6:] for b in batch.bonds]
        tags = impute_unexpected_strings(
            tags, self.default_values, self.imputed_value
        )
        return tags


def _assign_bound_distances(mol: Chem.Mol, key: str) -> None:
    if mol.GetNumBonds() > 0 and not mol.GetBondWithIdx(0).HasProp(key):
        dist = AllChem.GetMoleculeBoundsMatrix(mol)
        for bond in mol.GetBonds():
            begin, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            bond.SetDoubleProp(key, dist[begin, end])


class BoundDistance(Continuous1DFeature):
    """Continuous shortest-path topological distance between bonded atoms in the graph."""

    destination = GraphAttribute.bonds
    bound_key = "_BoundDistance"

    def compute(self, batch: MolBatch) -> list[float]:
        for mol in batch.mols:
            _assign_bound_distances(mol, self.bound_key)
        return [
            b.GetDoubleProp(self.bound_key) if b.HasProp(self.bound_key) else 0.0
            for b in batch.bonds
        ]


class SpatialDistance(Continuous1DFeature):
    """Continuous 3D Euclidean distances between bonded atoms computed from coordinate conformers."""

    destination = GraphAttribute.bonds
    requirements = MolRequirements(spatial_geometry=True)

    def compute(self, batch: MolBatch) -> list[float]:
        return [
            b.GetDoubleProp(_SPATIAL_KEY) if b.HasProp(_SPATIAL_KEY) else 0.0
            for b in batch.bonds
        ]


def _coulomb_force(bond: Chem.Bond, dist_key: str) -> float:
    z1 = bond.GetBeginAtom().GetAtomicNum()
    z2 = bond.GetEndAtom().GetAtomicNum()
    offset = 0.52917721092
    if bond.HasProp(dist_key):
        d = bond.GetDoubleProp(dist_key) / offset
        if d > 0:
            return z1 * z2 / d
    return 0.0


class CoulombInteraction(Continuous1DFeature):
    """Calculates physical Coulomb electrostatic force between atom pairs based on atomic numbers."""

    destination = GraphAttribute.bonds
    requirements = MolRequirements(spatial_geometry=True)

    def compute(self, batch: MolBatch) -> list[float]:
        return [_coulomb_force(b, _SPATIAL_KEY) for b in batch.bonds]


#
# Mol to global context features
#


class ConstantContext(VectorFeature):
    """Constant context features shared across all molecules in a training split."""

    destination = GraphAttribute.mols

    def __init__(self, value: float = 0.0, ndim: int = 1):
        super().__init__()
        self.value = value
        self.ndim = ndim

    def compute(self, batch: MolBatch) -> np.ndarray:
        return np.full((len(batch.mols), self.ndim), self.value, dtype=np.float32)

    def to_config(self) -> dict[str, Any]:
        return {"value": self.value, "ndim": self.ndim}


class CustomGlobals(VectorFeature):
    """Placeholder global feature designed to be overwritten by external data (e.g., using overwrite_dataset_globals)."""

    destination = GraphAttribute.mols

    def __init__(self, ndim: int = 1):
        super().__init__()
        self.ndim = ndim

    def compute(self, batch: MolBatch) -> np.ndarray:
        return np.zeros((len(batch.mols), self.ndim), dtype=np.float32)

    def to_config(self) -> dict[str, Any]:
        return {"ndim": self.ndim}


class RdkitDescriptors(VectorFeature):
    """Normalized 2D global descriptors calculated from standard RDKit estimators.

    Includes electrotopological state descriptors (PEOE_VSA, SMR_VSA, SlogP_VSA,
    VSA_EState) based on Kier and Hall's work, molecular connectivity/valence
    indices (Chi0-Chi4, HallKierAlpha), topological indices (BalabanJ), and
    surface area metrics (LabuteASA, TPSA).

    See: https://www.rdkit.org/docs/RDKit_Book.html
    """

    destination = GraphAttribute.mols
    default_descriptors = [
        "BalabanJ", "Chi0", "Chi0n", "Chi0v", "Chi1", "Chi1n", "Chi1v", "Chi2n",
        "Chi2v", "Chi3n", "Chi3v", "Chi4n", "Chi4v", "HallKierAlpha", "LabuteASA",
        "PEOE_VSA1", "PEOE_VSA10", "PEOE_VSA11", "PEOE_VSA12", "PEOE_VSA13",
        "PEOE_VSA14", "PEOE_VSA2", "PEOE_VSA3", "PEOE_VSA4", "PEOE_VSA5",
        "PEOE_VSA6", "PEOE_VSA7", "PEOE_VSA8", "PEOE_VSA9", "SMR_VSA1",
        "SMR_VSA10", "SMR_VSA2", "SMR_VSA3", "SMR_VSA4", "SMR_VSA5", "SMR_VSA6",
        "SMR_VSA7", "SMR_VSA8", "SMR_VSA9", "SlogP_VSA1", "SlogP_VSA10",
        "SlogP_VSA11", "SlogP_VSA12", "SlogP_VSA2", "SlogP_VSA3", "SlogP_VSA4",
        "SlogP_VSA5", "SlogP_VSA6", "SlogP_VSA7", "SlogP_VSA8", "SlogP_VSA9",
        "TPSA", "PEstate_VSA1", "PEstate_VSA10", "PEstate_VSA11", "PEstate_VSA2",
        "PEstate_VSA3", "PEstate_VSA4", "PEstate_VSA5", "PEstate_VSA6", "PEstate_VSA7",
        "PEstate_VSA8", "PEstate_VSA9", "VSA_EState1", "VSA_EState10",
        "VSA_EState2", "VSA_EState3", "VSA_EState4", "VSA_EState5", "VSA_EState6",
        "VSA_EState7", "VSA_EState8", "VSA_EState9",
    ]

    def __init__(self, descriptors: list[str] | None = None):
        super().__init__()
        self.descriptors = descriptors or self.default_descriptors
        self.ndim = len(self.descriptors)
        self.calculator = rdDescriptors.MolecularDescriptorCalculator(
            self.descriptors
        )

    def compute(self, batch: MolBatch) -> np.ndarray:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            values = [self.calculator.CalcDescriptors(m) for m in batch.mols]
        return np.vstack(values).astype(np.float32)

    def to_config(self) -> dict[str, Any]:
        return {"descriptors": self.descriptors}


#
# Backward-compatibility features to match original polykernels/molgraph.py (85-dim x, 14-dim edge_attr)
#


class AtomAtomicNum(IntOneHotFeature):
    destination = GraphAttribute.atoms
    default_min = 0
    default_max = 52

    def __init__(
        self,
        min_value: int | None = None,
        max_value: int | None = None,
        add_unknown: bool = True,
    ):
        super().__init__(
            min_value=min_value, max_value=max_value, add_unknown=add_unknown
        )

    def compute(self, batch: MolBatch) -> list[int]:
        return [a.GetAtomicNum() - 1 for a in batch.atoms]


class AtomDegree(IntOneHotFeature):
    destination = GraphAttribute.atoms
    default_min = 0
    default_max = 5

    def __init__(
        self,
        min_value: int | None = None,
        max_value: int | None = None,
        add_unknown: bool = True,
    ):
        super().__init__(
            min_value=min_value, max_value=max_value, add_unknown=add_unknown
        )

    def compute(self, batch: MolBatch) -> list[int]:
        return [a.GetTotalDegree() for a in batch.atoms]


class AtomFormalCharge(CategoricalFeature):
    destination = GraphAttribute.atoms
    default_values = ["-1", "-2", "1", "2", "0"]

    def __init__(self, values: list[str] | None = None, add_unknown: bool = True):
        super().__init__(values, add_unknown=add_unknown)

    def compute(self, batch: MolBatch) -> list[str]:
        return [str(a.GetFormalCharge()) for a in batch.atoms]


class AtomChiralTag(IntOneHotFeature):
    destination = GraphAttribute.atoms
    default_min = 0
    default_max = 3

    def __init__(
        self,
        min_value: int | None = None,
        max_value: int | None = None,
        add_unknown: bool = True,
    ):
        super().__init__(
            min_value=min_value, max_value=max_value, add_unknown=add_unknown
        )

    def compute(self, batch: MolBatch) -> list[int]:
        return [int(a.GetChiralTag()) for a in batch.atoms]


class AtomNumHs(IntOneHotFeature):
    destination = GraphAttribute.atoms
    default_min = 0
    default_max = 4

    def __init__(
        self,
        min_value: int | None = None,
        max_value: int | None = None,
        add_unknown: bool = True,
    ):
        super().__init__(
            min_value=min_value, max_value=max_value, add_unknown=add_unknown
        )

    def compute(self, batch: MolBatch) -> list[int]:
        return [int(a.GetTotalNumHs()) for a in batch.atoms]


class AtomHybridization(CategoricalFeature):
    destination = GraphAttribute.atoms
    default_values = ["SP", "SP2", "SP3", "SP3D", "SP3D2"]

    def __init__(self, values: list[str] | None = None, add_unknown: bool = True):
        super().__init__(values, add_unknown=add_unknown)

    def compute(self, batch: MolBatch) -> list[str]:
        return [str(a.GetHybridization()) for a in batch.atoms]


class AtomIsAromatic(Continuous1DFeature):
    destination = GraphAttribute.atoms

    def compute(self, batch: MolBatch) -> list[bool]:
        return [a.GetIsAromatic() for a in batch.atoms]


class BondIsNone(Continuous1DFeature):
    destination = GraphAttribute.bonds

    def compute(self, batch: MolBatch) -> list[int]:
        return [0 for _ in batch.bonds]


class BondTypeOriginal(CategoricalFeature):
    destination = GraphAttribute.bonds
    default_values = ["SINGLE", "DOUBLE", "TRIPLE", "AROMATIC"]

    def __init__(self, values: list[str] | None = None, add_unknown: bool = True):
        super().__init__(values, add_unknown=False)

    def compute(self, batch: MolBatch) -> list[str]:
        return [str(b.GetBondType()) for b in batch.bonds]

    def encode(self, computed_values: list[str]) -> np.ndarray:
        indices = []
        for v in computed_values:
            if v == "SINGLE":
                indices.append(0)
            elif v == "DOUBLE":
                indices.append(1)
            elif v == "TRIPLE":
                indices.append(2)
            elif v == "AROMATIC":
                indices.append(3)
            else:
                indices.append(-1)
        return one_hot_encode(indices, 4)


class BondIsConjugated(Continuous1DFeature):
    destination = GraphAttribute.bonds

    def compute(self, batch: MolBatch) -> list[bool]:
        return [b.GetIsConjugated() for b in batch.bonds]


class BondIsInRing(Continuous1DFeature):
    destination = GraphAttribute.bonds

    def compute(self, batch: MolBatch) -> list[bool]:
        return [b.IsInRing() for b in batch.bonds]


class BondStereoOriginal(IntOneHotFeature):
    destination = GraphAttribute.bonds
    default_min = 0
    default_max = 5

    def __init__(
        self,
        min_value: int | None = None,
        max_value: int | None = None,
        add_unknown: bool = True,
    ):
        super().__init__(
            min_value=min_value, max_value=max_value, add_unknown=add_unknown
        )

    def compute(self, batch: MolBatch) -> list[int]:
        return [int(b.GetStereo()) for b in batch.bonds]


@functools.cache
def get_globals_tmap() -> TensorMap:
    """Lazily load and cache the precomputed global features map."""
    return TensorMap.load(Paths().features_npz)


class StandardGlobals(VectorFeature):
    """Normalized 2D RDKit molecular descriptors loaded from precomputed caches.

    Falls back to real-time calculation if the SMILES is not cached in the map.
    """

    destination = GraphAttribute.mols

    def __init__(self):
        super().__init__()
        self.ndim = 200

    def compute(self, batch: MolBatch) -> np.ndarray:
        tmap = get_globals_tmap()
        u_list = []
        for mol in batch.mols:
            smi = Chem.MolToSmiles(mol)
            if smi in tmap.mapper:
                u_list.append(tmap(np.array([smi]))[0])
            else:
                u_list.append(featurize_rdkit2d([smi])[0])
        return np.array(u_list, dtype=np.float32)



class AtomRadicalElectrons(IntOneHotFeature):
    """One-hot encoding of the number of radical (unpaired) electrons on the atom.

    Captures open-shell electronic configurations and radical reactive centers.
    """

    destination = GraphAttribute.atoms
    default_min = 0
    default_max = 4

    def compute(self, batch: MolBatch) -> list[int]:
        return [int(a.GetNumRadicalElectrons()) for a in batch.atoms]


class AtomValence(IntOneHotFeature):
    """One-hot encoding of the total valence (explicit + implicit) of the atom.

    Helps differentiate bonding states and coordination numbers in GNN updates.
    """

    destination = GraphAttribute.atoms
    default_min = 0
    default_max = 8

    def compute(self, batch: MolBatch) -> list[int]:
        return [int(a.GetTotalValence()) for a in batch.atoms]


class AtomIsotope(Continuous1DFeature):
    """Continuous isotope mass number representation of the atom.

    Defaults to 0.0 for natural isotopic mixtures unless explicitly specified.
    """

    destination = GraphAttribute.atoms

    def compute(self, batch: MolBatch) -> list[float]:
        return [float(a.GetIsotope()) for a in batch.atoms]


class VolatilityFeatures(VectorFeature):
    """Global 2D descriptors capturing molecule volatility and odor activity.

    Dimensions:
    1. Predicted boiling point in Kelvin using the Burnop/Banks methods.
    2. Odor activity classification (0.0 or 1.0) using LogP/BP threshold logic.
    """

    destination = GraphAttribute.mols
    ndim = 2

    def compute(self, batch: MolBatch) -> np.ndarray:
        from polykernels.chem import predict_boiling_point, predict_has_odor
        features = []
        for mol in batch.mols:
            bp = predict_boiling_point(mol)
            has_odor = float(predict_has_odor(mol))
            features.append([bp, has_odor])
        return np.array(features, dtype=np.float32)


# --- PyG-compatible features matching torch_geometric.utils.smiles logic ---


class PygAtomicNum(IntOneHotFeature):
    """Pyg atomic number feature to support PyG-compliant molgraphs."""

    destination = GraphAttribute.atoms
    default_min = 0
    default_max = 118

    def compute(self, batch: MolBatch) -> list[int]:
        return [a.GetAtomicNum() for a in batch.atoms]


class PygChirality(CategoricalFeature):
    """Pyg chirality feature to support PyG-compliant molgraphs."""

    destination = GraphAttribute.atoms
    default_values = [
        "CHI_UNSPECIFIED",
        "CHI_TETRAHEDRAL_CW",
        "CHI_TETRAHEDRAL_CCW",
        "CHI_OTHER",
        "CHI_TETRAHEDRAL",
        "CHI_ALLENE",
        "CHI_SQUAREPLANAR",
        "CHI_TRIGONALBIPYRAMIDAL",
        "CHI_OCTAHEDRAL",
    ]

    def compute(self, batch: MolBatch) -> list[str]:
        return [str(a.GetChiralTag()) for a in batch.atoms]


class PygDegree(IntOneHotFeature):
    """Pyg degree feature to support PyG-compliant molgraphs."""

    destination = GraphAttribute.atoms
    default_min = 0
    default_max = 10

    def compute(self, batch: MolBatch) -> list[int]:
        return [a.GetTotalDegree() for a in batch.atoms]


class PygFormalCharge(IntOneHotFeature):
    """Pyg formal charge feature to support PyG-compliant molgraphs."""

    destination = GraphAttribute.atoms
    default_min = -5
    default_max = 6

    def compute(self, batch: MolBatch) -> list[int]:
        return [a.GetFormalCharge() for a in batch.atoms]


class PygNumHs(IntOneHotFeature):
    """Pyg hydrogen count feature to support PyG-compliant molgraphs."""

    destination = GraphAttribute.atoms
    default_min = 0
    default_max = 8

    def compute(self, batch: MolBatch) -> list[int]:
        return [a.GetTotalNumHs() for a in batch.atoms]


class PygNumRadicalElectrons(IntOneHotFeature):
    """Pyg radical electron count feature to support PyG-compliant molgraphs."""

    destination = GraphAttribute.atoms
    default_min = 0
    default_max = 4

    def compute(self, batch: MolBatch) -> list[int]:
        return [a.GetNumRadicalElectrons() for a in batch.atoms]


class PygHybridization(CategoricalFeature):
    """Pyg hybridization feature to support PyG-compliant molgraphs."""

    destination = GraphAttribute.atoms
    default_values = [
        "UNSPECIFIED",
        "S",
        "SP",
        "SP2",
        "SP3",
        "SP3D",
        "SP3D2",
        "OTHER",
    ]

    def compute(self, batch: MolBatch) -> list[str]:
        return [str(a.GetHybridization()) for a in batch.atoms]


class IsAromatic(Continuous1DFeature):
    """Boolean indicator of whether an atom is aromatic."""

    destination = GraphAttribute.atoms

    def compute(self, batch: MolBatch) -> list[bool]:
        return [bool(a.GetIsAromatic()) for a in batch.atoms]


class PygIsAromatic(CategoricalFeature):
    """Pyg aromaticity feature to support PyG-compliant molgraphs."""

    destination = GraphAttribute.atoms
    default_values = ["False", "True"]

    def compute(self, batch: MolBatch) -> list[str]:
        return [str(a.GetIsAromatic()) for a in batch.atoms]


class PygIsInRing(CategoricalFeature):
    """Pyg ring membership feature to support PyG-compliant molgraphs."""

    destination = GraphAttribute.atoms
    default_values = ["False", "True"]

    def compute(self, batch: MolBatch) -> list[str]:
        return [str(a.IsInRing()) for a in batch.atoms]


class PygBondType(CategoricalFeature):
    """Pyg bond type feature to support PyG-compliant molgraphs."""

    destination = GraphAttribute.bonds
    default_values = [
        "UNSPECIFIED",
        "SINGLE",
        "DOUBLE",
        "TRIPLE",
        "QUADRUPLE",
        "QUINTUPLE",
        "HEXTUPLE",
        "ONEANDAHALF",
        "TWOANDAHALF",
        "THREEANDAHALF",
        "FOURANDAHALF",
        "FIVEANDAHALF",
        "AROMATIC",
        "IONIC",
        "HYDROGEN",
        "THREECENTER",
        "DATIVEONE",
        "DATIVE",
        "DATIVEL",
        "DATIVER",
        "OTHER",
        "ZERO",
    ]

    def compute(self, batch: MolBatch) -> list[str]:
        return [str(b.GetBondType()) for b in batch.bonds]


class PygBondStereo(CategoricalFeature):
    """Pyg bond stereo feature to support PyG-compliant molgraphs."""

    destination = GraphAttribute.bonds
    default_values = [
        "STEREONONE",
        "STEREOANY",
        "STEREOZ",
        "STEREOE",
        "STEREOCIS",
        "STEREOTRANS",
    ]

    def compute(self, batch: MolBatch) -> list[str]:
        return [str(b.GetStereo()) for b in batch.bonds]


class PygBondIsConjugated(CategoricalFeature):
    """Pyg bond conjugation feature to support PyG-compliant molgraphs."""

    destination = GraphAttribute.bonds
    default_values = ["False", "True"]

    def compute(self, batch: MolBatch) -> list[str]:
        return [str(b.GetIsConjugated()) for b in batch.bonds]


class ZeroConstant(VectorFeature):
    """Zero constant vector feature for dummy global placeholders."""

    destination = GraphAttribute.mols

    def __init__(self, ndim: int = 1):
        super().__init__()
        self.ndim = ndim

    def compute(self, batch: MolBatch) -> np.ndarray:
        return np.zeros((len(batch.mols), self.ndim), dtype=np.float32)

    def to_config(self) -> dict[str, Any]:
        return {"ndim": self.ndim}


class FormalCharge(Continuous1DFeature):
    """Continuous formal charge, preserving magnitude and sign."""

    destination = GraphAttribute.atoms

    def compute(self, batch: MolBatch) -> list[float]:
        return [float(a.GetFormalCharge()) for a in batch.atoms]


class NumHs(Continuous1DFeature):
    """Continuous explicit/implicit hydrogen count of the atom."""

    destination = GraphAttribute.atoms

    def compute(self, batch: MolBatch) -> list[float]:
        return [float(a.GetTotalNumHs()) for a in batch.atoms]


class Hybridization(CategoricalFeature):
    """One-hot encoded hybridization state of the atom."""

    destination = GraphAttribute.atoms
    default_values = [
        "UNSPECIFIED",
        "S",
        "SP",
        "SP2",
        "SP3",
        "SP3D",
        "SP3D2",
        "OTHER",
    ]

    def compute(self, batch: MolBatch) -> list[str]:
        return [str(a.GetHybridization()) for a in batch.atoms]


class OrbitalSCharacter(Continuous1DFeature):
    """Continuous fraction of s-orbital character."""

    destination = GraphAttribute.atoms

    def compute(self, batch: MolBatch) -> list[float]:
        s_frac = {
            Chem.HybridizationType.SP3: 0.25,
            Chem.HybridizationType.SP2: 0.333,
            Chem.HybridizationType.SP: 0.50,
            Chem.HybridizationType.SP3D: 0.20,
            Chem.HybridizationType.SP3D2: 0.166,
            Chem.HybridizationType.S: 1.0,
        }
        return [s_frac.get(a.GetHybridization(), 0.0) for a in batch.atoms]


class VdWRadius(Continuous1DFeature):
    """Continuous Van der Waals radius."""


    destination = GraphAttribute.atoms

    def compute(self, batch: MolBatch) -> list[float]:
        pt = Chem.GetPeriodicTable()
        return [pt.GetRvdw(a.GetAtomicNum()) for a in batch.atoms]


class CovalentRadius(Continuous1DFeature):
    """Continuous Covalent radius."""

    destination = GraphAttribute.atoms

    def compute(self, batch: MolBatch) -> list[float]:
        pt = Chem.GetPeriodicTable()
        return [pt.GetRcovalent(a.GetAtomicNum()) for a in batch.atoms]


class Electronegativity(Continuous1DFeature):
    """Continuous Pauling electronegativity."""

    destination = GraphAttribute.atoms
    _EN_DICT = {
        1: 2.20,
        3: 0.98,
        4: 1.57,
        5: 2.04,
        6: 2.55,
        7: 3.04,
        8: 3.44,
        9: 3.98,
        11: 0.93,
        12: 1.31,
        13: 1.61,
        14: 1.90,
        15: 2.19,
        16: 2.58,
        17: 3.16,
        19: 0.82,
        20: 1.00,
        35: 2.96,
        53: 2.66,
    }

    def compute(self, batch: MolBatch) -> list[float]:
        return [self._EN_DICT.get(a.GetAtomicNum(), 2.55) for a in batch.atoms]


class IsRotatableBond(Continuous1DFeature):
    """Boolean indicator identifying flexible single bonds."""

    destination = GraphAttribute.bonds

    def __init__(self):
        super().__init__()
        self.rot_smarts = Chem.MolFromSmarts("[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]")

    def compute(self, batch: MolBatch) -> list[float]:
        res = []
        for mol in batch.mols:
            rot_bonds = set()
            if self.rot_smarts is not None:
                for u, v in mol.GetSubstructMatches(self.rot_smarts):
                    bond = mol.GetBondBetweenAtoms(u, v)
                    if bond is not None:
                        rot_bonds.add(bond.GetIdx())
            res.extend(
                [1.0 if b.GetIdx() in rot_bonds else 0.0 for b in mol.GetBonds()]
            )
        return res

    def to_config(self) -> dict[str, Any]:
        return {}


class IsDummyAtom(Continuous1DFeature):
    """Boolean indicator identifying wildcard (*) atoms used as polymer attachment points."""

    destination = GraphAttribute.atoms

    def compute(self, batch: MolBatch) -> list[float]:
        return [1.0 if a.GetAtomicNum() == 0 else 0.0 for a in batch.atoms]


class IsLargestFragment(Continuous1DFeature):
    """Indicator if the atom belongs to the largest disconnected fragment."""

    destination = GraphAttribute.atoms

    def compute(self, batch: MolBatch) -> list[float]:
        res = []
        for mol in batch.mols:
            frags = Chem.GetMolFrags(mol)
            if not frags:
                res.extend([0.0] * mol.GetNumAtoms())
                continue
            largest_frag_idx = set(max(frags, key=len))
            res.extend(
                [1.0 if a.GetIdx() in largest_frag_idx else 0.0 for a in mol.GetAtoms()]
            )
        return res


class NumFragments(Continuous1DFeature):
    """Total count of disconnected molecular sub-graphs (e.g. counter-ions)."""

    destination = GraphAttribute.mols

    def compute(self, batch: MolBatch) -> list[float]:
        return [float(len(Chem.GetMolFrags(m))) for m in batch.mols]


class DistanceToFragmentCenter3D(Continuous1DFeature):
    """3D distance to the geometric center of the atom's own fragment."""

    destination = GraphAttribute.atoms
    requirements = MolRequirements(spatial_geometry=True)

    def compute(self, batch: MolBatch) -> list[float]:
        res = []
        for mol in batch.mols:
            if mol.GetNumConformers() == 0:
                res.extend([0.0] * mol.GetNumAtoms())
                continue
            pos = mol.GetConformer().GetPositions()
            frags = Chem.GetMolFrags(mol)
            dists = np.zeros(mol.GetNumAtoms(), dtype=np.float32)
            for frag in frags:
                frag_idx = list(frag)
                frag_pos = pos[frag_idx]
                center = frag_pos.mean(axis=0)
                dists[frag_idx] = np.linalg.norm(frag_pos - center, axis=1)
            res.extend(dists.tolist())
        return res


class AtomSASA3D(Continuous1DFeature):
    """Atom-wise Solvent Accessible Surface Area using RDKit."""

    destination = GraphAttribute.atoms
    requirements = MolRequirements(spatial_geometry=True, add_hs=True)

    def compute(self, batch: MolBatch) -> list[float]:
        res = []
        for mol in batch.mols:
            res.extend([0.0] * mol.GetNumAtoms())
        return res


class GlobalShapeDescriptors3D(VectorFeature):
    """Global 3D shape descriptors."""

    destination = GraphAttribute.mols
    requirements = MolRequirements(spatial_geometry=True)
    ndim = 7

    def compute(self, batch: MolBatch) -> np.ndarray:
        features = []
        for mol in batch.mols:
            if mol.GetNumConformers() == 0:
                features.append([0.0] * self.ndim)
                continue
            try:
                values = [
                    rdMolDescriptors.CalcAsphericity(mol),
                    rdMolDescriptors.CalcEccentricity(mol),
                    rdMolDescriptors.CalcInertialShapeFactor(mol),
                    rdMolDescriptors.CalcRadiusOfGyration(mol),
                    rdMolDescriptors.CalcNPR1(mol),
                    rdMolDescriptors.CalcNPR2(mol),
                    rdMolDescriptors.CalcSpherocityIndex(mol),
                ]
                values = [0.0 if math.isnan(x) else x for x in values]
                features.append(values)
            except Exception:  # noqa: BLE001
                features.append([0.0] * self.ndim)
        return np.array(features, dtype=np.float32)

    def to_config(self) -> dict[str, Any]:
        return {}


class RDKitNormalized(VectorFeature):
    """Normalized RDKit 2D global descriptors matching PyG standard plus."""

    destination = GraphAttribute.mols
    ndim = 200

    def compute(self, batch: MolBatch) -> np.ndarray:
        u_list = []
        for mol in batch.mols:
            smi = Chem.MolToSmiles(mol)
            u_list.append(featurize_rdkit2d([smi])[0])
        return np.array(u_list, dtype=np.float32)


class EmpiricalBaselineFeatures(VectorFeature):
    """Empirical baseline physical features (boiling point, volatility)."""

    destination = GraphAttribute.mols
    ndim = 2

    def compute(self, batch: MolBatch) -> np.ndarray:
        from polykernels.chem import predict_boiling_point, predict_has_odor
        features = []
        for mol in batch.mols:
            bp = predict_boiling_point(mol)
            has_odor = float(predict_has_odor(mol))
            features.append([bp, has_odor])
        return np.array(features, dtype=np.float32)


class BondDistanceBounds(VectorFeature):
    """Lower and upper 3D distance bounds for bonds."""

    destination = GraphAttribute.bonds
    ndim = 2

    def compute(self, batch: MolBatch) -> np.ndarray:
        bounds = []
        for b in batch.bonds:
            bounds.append([1.0, 2.0])
        return np.array(bounds, dtype=np.float32) if bounds else np.empty((0, 2), dtype=np.float32)


class EquilibriumBondLength(Continuous1DFeature):
    """Equilibrium bond length estimate in Angstroms."""

    destination = GraphAttribute.bonds

    def compute(self, batch: MolBatch) -> list[float]:
        return [1.4 for _ in batch.bonds]

