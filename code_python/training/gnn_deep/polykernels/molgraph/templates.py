"""Serializable configuration templates for graph featurization."""

import dataclasses
from enum import StrEnum
from typing import Any

import draccus

import polykernels.molgraph.features as features


class GraphType(StrEnum):
    """Supported graph templates for molecular featurization."""

    STANDARD = "standard"
    STANDARD_PLUS = "standard_plus"
    FULL = "full"
    RIGR = "rigr"
    MINIMAL = "minimal"
    POLY = "poly"
    FULL_GRAPH = "full_graph"
    INTRINSIC = "intrinsic"
    CUSTOM = "custom"


draccus.encode.register(GraphType, lambda x: x.value)


def get_template(
    graph_type: GraphType | str | "Template",
    global_dim: int = 64,
) -> "Template":
    """Retrieve template corresponding to the given GraphType enum, string identifier, or Template object."""
    if isinstance(graph_type, Template):
        return graph_type
    match GraphType(graph_type):
        case GraphType.STANDARD:
            return StandardTemplate()
        case GraphType.STANDARD_PLUS:
            return StandardPlusTemplate()
        case GraphType.FULL:
            return FullTemplate()
        case GraphType.RIGR:
            return RigrTemplate()
        case GraphType.MINIMAL:
            return MinimalTemplate()
        case GraphType.POLY:
            return PolyTemplate()
        case GraphType.FULL_GRAPH:
            return FullGraphTemplate()
        case GraphType.INTRINSIC:
            return IntrinsicGraphTemplate()
        case GraphType.CUSTOM:
            return get_custom_full(global_dim=global_dim)
        case _:
            raise ValueError(f"Unknown graph type: {graph_type}")


@dataclasses.dataclass
class Template:
    """A collection of node, edge, and global features defining a graph format."""

    name: str
    node_features: list[features.GraphFeature]
    edge_features: list[features.GraphFeature]
    global_features: list[features.GraphFeature]

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "node_features": [
                {"class": f.__class__.__name__, "config": f.to_config()}
                for f in self.node_features
            ],
            "edge_features": [
                {"class": f.__class__.__name__, "config": f.to_config()}
                for f in self.edge_features
            ],
            "global_features": [
                {"class": f.__class__.__name__, "config": f.to_config()}
                for f in self.global_features
            ],
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Template":
        node_feats = []
        for item in data.get("node_features", []):
            feat_class = getattr(features, item["class"])
            node_feats.append(feat_class(**item.get("config", {})))

        edge_feats = []
        for item in data.get("edge_features", []):
            feat_class = getattr(features, item["class"])
            edge_feats.append(feat_class(**item.get("config", {})))

        global_feats = []
        for item in data.get("global_features", []):
            feat_class = getattr(features, item["class"])
            global_feats.append(feat_class(**item.get("config", {})))

        return cls(
            name=data["name"],
            node_features=node_feats,
            edge_features=edge_feats,
            global_features=global_feats,
        )


def StandardTemplate() -> Template:  # noqa: N802
    """Standard feature configuration with 1-dim zero global (177-dim node, 30-dim edge, 1-dim global) matching PyG SMILES featurizer."""
    return Template(
        name="standard",
        node_features=[
            features.PygAtomicNum(),
            features.PygChirality(),
            features.PygDegree(),
            features.PygFormalCharge(),
            features.PygNumHs(),
            features.PygNumRadicalElectrons(),
            features.PygHybridization(),
            features.PygIsAromatic(),
            features.PygIsInRing(),
        ],
        edge_features=[
            features.PygBondType(),
            features.PygBondStereo(),
            features.PygBondIsConjugated(),
        ],
        global_features=[
            features.ZeroConstant(ndim=1),
        ],
    )




def StandardPlusTemplate() -> Template:  # noqa: N802
    """Standard plus feature configuration with normalized RDKit global descriptors (177-dim node, 30-dim edge, 200-dim global) matching PyG SMILES featurizer."""
    return Template(
        name="standard_plus",
        node_features=[
            features.PygAtomicNum(),
            features.PygChirality(),
            features.PygDegree(),
            features.PygFormalCharge(),
            features.PygNumHs(),
            features.PygNumRadicalElectrons(),
            features.PygHybridization(),
            features.PygIsAromatic(),
            features.PygIsInRing(),
        ],
        edge_features=[
            features.PygBondType(),
            features.PygBondStereo(),
            features.PygBondIsConjugated(),
        ],
        global_features=[
            features.RDKitNormalized(),
        ],
    )


def FullTemplate() -> Template:  # noqa: N802
    """Full modular feature configuration with physical, topological, and 3D descriptors."""
    return Template(
        name="full",
        node_features=[
            features.AtomicSymbolImpute(),
            features.PeriodicGroup(),
            features.AtomIsotope(),
            features.AtomicWeight(),
            features.VdWRadius(),
            features.CovalentRadius(),
            features.Electronegativity(),
            features.Degree(),
            features.AtomValence(),
            features.AtomRadicalElectrons(),
            features.FormalCharge(),
            features.NumHs(),
            features.OrbitalSCharacter(),
            features.PartialCharges(),
            features.IsAromatic(),
            features.HydrogenBonding(),
            features.CrippenFeatures(),
            features.InRingSizes(),
            features.ChiralTag(),
            features.ChiralCenter(),
            features.AtomSASA3D(),
            features.DistanceToFragmentCenter3D(),
            features.GraphwaveEmbeddings(ndim=12, step_size=3.0, scales=[2.0, 10.0]),
        ],
        edge_features=[
            features.BondType(),
            features.BondIsConjugated(),
            features.InRing(),
            features.BondStereo(),
            features.IsRotatableBond(),
            features.BoundDistance(),
            features.SpatialDistance(),
            features.CoulombInteraction(),
        ],
        global_features=[
            features.GlobalShapeDescriptors3D(),
            features.RdkitDescriptors(),
            features.EmpiricalBaselineFeatures(),
        ],
    )


def PolyTemplate() -> Template:  # noqa: N802
    """Template tailored for multi-fragment complexes, salts, and polymers."""
    return Template(
        name="poly",
        node_features=[
            features.AtomicSymbolImpute(),
            features.IsDummyAtom(),
            features.PeriodicGroup(),
            features.AtomIsotope(),
            features.IsLargestFragment(),
            features.Degree(),
            features.AtomValence(),
            features.AtomRadicalElectrons(),
            features.FormalCharge(),
            features.NumHs(),
            features.OrbitalSCharacter(),
            features.PartialCharges(),
            features.IsAromatic(),
            features.HydrogenBonding(),
            features.CrippenFeatures(),
            features.InRingSizes(),
            features.ChiralTag(),
            features.ChiralCenter(),
            features.AtomSASA3D(),
            features.DistanceToFragmentCenter3D(),
            features.GraphwaveEmbeddings(ndim=12, step_size=3.0, scales=[2.0, 10.0]),
        ],
        edge_features=[
            features.BondType(
                values=["SINGLE", "DOUBLE", "TRIPLE", "AROMATIC", "ZERO", "IONIC"],
                add_unknown=True,
            ),
            features.BondIsConjugated(),
            features.InRing(),
            features.BondStereo(),
            features.IsRotatableBond(),
            features.BoundDistance(),
            features.SpatialDistance(),
            features.CoulombInteraction(),
        ],
        global_features=[
            features.NumFragments(),
            features.GlobalShapeDescriptors3D(),
            features.RdkitDescriptors(),
        ],
    )


def RigrTemplate() -> Template:  # noqa: N802
    """RIGR template (Resonance-Invariant Graph Representation) based on FullTemplate."""
    return Template(
        name="rigr",
        node_features=[
            features.AtomicSymbolImpute(),
            features.PeriodicGroup(),
            features.Degree(),
            features.AtomicWeight(),
            features.InRingSizes(),
            features.GraphwaveEmbeddings(ndim=12, step_size=3.0, scales=[2.0, 10.0]),
            features.AtomIsotope(),
        ],
        edge_features=[
            features.InRing(),
            features.BoundDistance(),
            features.SpatialDistance(),
            features.CoulombInteraction(),
        ],
        global_features=[
            features.RdkitDescriptors(),
            features.EmpiricalBaselineFeatures(),
        ],
    )


def MinimalTemplate() -> Template:  # noqa: N802
    """Minimal template (low-dimensional, SMILES-only featurization)."""
    return Template(
        name="minimal",
        node_features=[
            features.AtomicSymbolImpute(),
            features.IsAromatic(),
        ],
        edge_features=[
            features.BondType(add_unknown=True),
        ],
        global_features=[
            features.ZeroConstant(ndim=1),
        ],
    )


def FullGraphTemplate() -> Template:  # noqa: N802
    """Full modular feature configuration without 3D spatial coordinate dependency."""
    return Template(
        name="full_graph",
        node_features=[
            features.AtomicSymbolImpute(),
            features.PeriodicGroup(),
            features.AtomIsotope(),
            features.AtomicWeight(),
            features.VdWRadius(),
            features.CovalentRadius(),
            features.Electronegativity(),
            features.Degree(),
            features.AtomValence(),
            features.AtomRadicalElectrons(),
            features.FormalCharge(),
            features.NumHs(),
            features.OrbitalSCharacter(),
            features.PartialCharges(),
            features.IsAromatic(),
            features.HydrogenBonding(),
            features.CrippenFeatures(),
            features.InRingSizes(),
            features.ChiralTag(),
            features.ChiralCenter(),
            features.GraphwaveEmbeddings(ndim=12, step_size=3.0, scales=[2.0, 10.0]),
        ],
        edge_features=[
            features.BondType(),
            features.BondIsConjugated(),
            features.InRing(),
            features.BondStereo(),
            features.IsRotatableBond(),
            features.BoundDistance(),
        ],
        global_features=[
            features.RdkitDescriptors(),
            features.EmpiricalBaselineFeatures(),
        ],
    )


def IntrinsicGraphTemplate() -> Template:  # noqa: N802
    """Rich molecular graph features without conformer-derived or molecular-global inputs."""
    return Template(
        name="intrinsic",
        node_features=[
            features.AtomicSymbolImpute(),
            features.PeriodicGroup(),
            features.AtomIsotope(),
            features.AtomicWeight(),
            features.VdWRadius(),
            features.CovalentRadius(),
            features.Electronegativity(),
            features.Degree(),
            features.AtomValence(),
            features.AtomRadicalElectrons(),
            features.FormalCharge(),
            features.NumHs(),
            features.OrbitalSCharacter(),
            features.PartialCharges(),
            features.IsAromatic(),
            features.HydrogenBonding(),
            features.CrippenFeatures(),
            features.InRingSizes(),
            features.ChiralTag(),
            features.ChiralCenter(),
            features.GraphwaveEmbeddings(ndim=12, step_size=3.0, scales=[2.0, 10.0]),
        ],
        edge_features=[
            features.BondType(),
            features.BondIsConjugated(),
            features.InRing(),
            features.BondStereo(),
            features.IsRotatableBond(),
            features.BoundDistance(),
            features.EquilibriumBondLength(),
        ],
        global_features=[],
    )


def get_custom_full(global_dim: int = 64) -> Template:  # noqa: N802
    """Custom template using intrinsic node/edge features with custom global placeholder slots."""
    return Template(
        name="custom",
        node_features=[
            features.AtomicSymbolImpute(),
            features.PeriodicGroup(),
            features.AtomIsotope(),
            features.AtomicWeight(),
            features.VdWRadius(),
            features.CovalentRadius(),
            features.Electronegativity(),
            features.Degree(),
            features.AtomValence(),
            features.AtomRadicalElectrons(),
            features.FormalCharge(),
            features.NumHs(),
            features.OrbitalSCharacter(),
            features.PartialCharges(),
            features.IsAromatic(),
            features.HydrogenBonding(),
            features.CrippenFeatures(),
            features.InRingSizes(),
            features.ChiralTag(),
            features.ChiralCenter(),
            features.GraphwaveEmbeddings(ndim=12, step_size=3.0, scales=[2.0, 10.0]),
        ],
        edge_features=[
            features.BondType(),
            features.BondIsConjugated(),
            features.InRing(),
            features.BondStereo(),
            features.IsRotatableBond(),
            features.BoundDistance(),
            features.EquilibriumBondLength(),
        ],
        global_features=[
            features.ZeroConstant(ndim=global_dim),
        ],
    )

