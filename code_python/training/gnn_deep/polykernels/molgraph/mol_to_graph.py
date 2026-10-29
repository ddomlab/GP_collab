"""MolToGraph class for transforming SMILES into molecular graphs."""

from collections.abc import Sequence
import copy

import numpy as np
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import torch
import torch_geometric

from polykernels.chem import get_robust_geometry
import polykernels.molgraph.features as features
from polykernels.molgraph.templates import Template


class MolToGraph:
    """Transform SMILES strings into PyTorch Geometric graph datasets.

    Handles RDKit molecule sanitization, robust 3D/2D conformer generation,
    and batched featurization of nodes, edges, and global descriptors
    according to a specified template.
    """

    def __init__(self, template: Template, *, use_3d: bool = True) -> None:
        self.template = template
        self.use_3d = use_3d
        self.requirements = features.MolRequirements()
        for f in (
            template.node_features
            + template.edge_features
            + template.global_features
        ):
            self.requirements.stereo_tags |= f.requirements.stereo_tags
            self.requirements.spatial_geometry |= f.requirements.spatial_geometry
            self.requirements.add_hs |= f.requirements.add_hs

    def prepare_molecule(self, mol: Chem.Mol) -> Chem.Mol:
        """Apply sanitization, stereochemistry, and coordinate conformation.

        If spatial geometry is required, embeds the molecule using ETKDGv3 and
        MMFF force-field optimization. Raises a ValueError if 3D conformer
        embedding fails.

        Args:
            mol: Input RDKit molecule.

        Returns:
            Chem.Mol: Prepared molecule with computed conformations.
        """
        mol = copy.deepcopy(mol)
        if self.requirements.add_hs:
            mol = Chem.AddHs(mol)
        if self.requirements.stereo_tags:
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        if self.requirements.spatial_geometry:
            mol = get_robust_geometry(mol, use_3d=self.use_3d)
            if mol.GetNumConformers() == 0:
                AllChem.Compute2DCoords(mol)
                if mol.GetNumConformers() == 0:
                    raise ValueError("Molecule conformer generation failed completely.")
            # Only warn on an actual 3D->2D fallback; 2D is expected (not a failure) when
            # self.use_3d is False, e.g. for multi-fragment "graphs" with no real joint
            # geometry (see DatasetSpec.use_3d in polykernels.regression_data).
            if self.use_3d and not mol.GetConformer().Is3D():
                from loguru import logger
                logger.warning("Molecule failed 3D embedding. Fell back to 2D conformer.")
            conf = mol.GetConformer()
            for bond in mol.GetBonds():
                i = bond.GetBeginAtomIdx()
                j = bond.GetEndAtomIdx()
                p1 = conf.GetAtomPosition(i)
                p2 = conf.GetAtomPosition(j)
                dist = p1.Distance(p2)
                bond.SetDoubleProp(features._SPATIAL_KEY, dist)
        return mol

    def __call__(
        self, smiles: str | Sequence[str]
    ) -> torch_geometric.data.Data | list[torch_geometric.data.Data]:
        """Convert SMILES strings into PyTorch Geometric Data objects.

        Parses the inputs, runs robust 3D/2D coordinate embedding, evaluates
        template-specified features in a memory-efficient batched layout,
        and assembles sorted sparse graph objects.

        Args:
            smiles: A single SMILES string or a sequence of SMILES strings.

        Returns:
            A single PyG Data object or a list of PyG Data objects.
        """
        is_single = isinstance(smiles, str)
        smi_list = [smiles] if is_single else list(smiles)

        mols = []
        valid_smiles = []
        for s in smi_list:
            mol = Chem.MolFromSmiles(s)
            if mol is None:
                raise ValueError(f"Invalid SMILES string: {s}")
            prepared = self.prepare_molecule(mol)
            mols.append(prepared)
            valid_smiles.append(Chem.MolToSmiles(prepared))

        batch = features.MolBatch(mols)

        # Compute node features
        node_feats_list = []
        for f in self.template.node_features:
            node_feats_list.append(f(batch))
        x_batch = (
            np.hstack(node_feats_list)
            if node_feats_list
            else np.empty((len(batch.atoms), 0), dtype=np.float32)
        )

        # Compute edge features
        edge_feats_list = []
        for f in self.template.edge_features:
            edge_feats_list.append(f(batch))
        e_batch = (
            np.hstack(edge_feats_list)
            if edge_feats_list
            else np.empty((len(batch.bonds), 0), dtype=np.float32)
        )

        # Compute global features
        global_feats_list = []
        for f in self.template.global_features:
            global_feats_list.append(f(batch))
        u_batch = (
            np.hstack(global_feats_list)
            if global_feats_list
            else np.empty((len(batch.mols), 0), dtype=np.float32)
        )

        atom_offsets = np.cumsum([0] + batch.natoms)
        bond_offsets = np.cumsum([0] + batch.nbonds)

        pyg_data_list = []
        for idx, mol in enumerate(mols):
            n_atoms = batch.natoms[idx]

            x_k = torch.from_numpy(x_batch[atom_offsets[idx] : atom_offsets[idx + 1]]).float()


            b_start = bond_offsets[idx]
            b_end = bond_offsets[idx + 1]
            e_k_raw = e_batch[b_start:b_end]

            if e_k_raw.shape[0] > 0:
                # Fast numpy extraction of bond indices
                bids = np.array([(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()], dtype=np.int64)
                b_index_np = np.empty((2, 2 * bids.shape[0]), dtype=np.int64)
                b_index_np[0, 0::2] = bids[:, 0]
                b_index_np[1, 0::2] = bids[:, 1]
                b_index_np[0, 1::2] = bids[:, 1]
                b_index_np[1, 1::2] = bids[:, 0]

                # Convert to tensor using from_numpy
                b_index = torch.from_numpy(b_index_np)
                # Replicates list extension loop using numpy repeat along rows
                b_feats = torch.from_numpy(np.repeat(e_k_raw, 2, axis=0)).float()

                perm = (b_index[0] * n_atoms + b_index[1]).argsort()

                b_index = b_index[:, perm]
                b_feats = b_feats[perm]
            else:
                b_index = torch.empty((2, 0), dtype=torch.long)
                b_feats = torch.empty((0, e_batch.shape[1]), dtype=torch.float)

            u_k = torch.from_numpy(u_batch[idx : idx + 1]).float()

            data = torch_geometric.data.Data(
                x=x_k,
                edge_index=b_index,
                edge_attr=b_feats,
                u=u_k,
                smiles=valid_smiles[idx],
            )
            pyg_data_list.append(data)

        return pyg_data_list[0] if is_single else pyg_data_list
