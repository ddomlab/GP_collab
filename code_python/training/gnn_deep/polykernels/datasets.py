"""Dataset wrappers and helpers for training models."""

from collections.abc import Sequence

import numpy as np
import torch
from torch.utils.data import Dataset
import torch_geometric.data

from polykernels import molgraph


class GraphDataset(Dataset):
    """Dataset wrapper for PyTorch Geometric Data objects and labels."""

    def __init__(
        self, graphs: list[torch_geometric.data.Data], labels: torch.Tensor
    ) -> None:
        self.graphs = graphs
        self.labels = labels

    def __len__(self) -> int:
        return len(self.graphs)

    def __getitem__(self, idx: int) -> tuple[torch_geometric.data.Data, torch.Tensor]:
        return self.graphs[idx], self.labels[idx]


def prepare_graph_dataset(
    smiles: Sequence[str] | np.ndarray,
    labels: np.ndarray,
    graph_type: str | molgraph.GraphType,
    *,
    use_3d: bool = True,
) -> GraphDataset:
    """Featurize a list of SMILES strings into graphs and wrap them in a GraphDataset.

    Args:
        use_3d: Embed real 3D conformers (default). Set False to use fast, always-
            succeeding 2D coordinates instead -- appropriate when a SMILES string
            doesn't represent one physically coherent molecule (e.g. two unrelated
            fragments joined with "."), since 3D embedding has no meaningful joint
            geometry to solve for and can be extremely slow or fail outright.
    """
    from polykernels.utils import get_logger

    logger = get_logger()
    template = molgraph.get_template(graph_type)
    mol_to_graph = molgraph.MolToGraph(template, use_3d=use_3d)

    logger.info(f"Featurizing {len(smiles)} SMILES into graphs using template: {template.name}")
    graphs = mol_to_graph(smiles)
    if isinstance(graphs, torch_geometric.data.Data):
        graphs = [graphs]

    return GraphDataset(
        graphs=graphs,
        labels=torch.tensor(labels, dtype=torch.float32),
    )


def overwrite_dataset_globals(
    dataset: GraphDataset,
    external_globals: np.ndarray,
) -> GraphDataset:
    """Overwrites the global feature vector `u` in all graphs of the dataset with external data."""
    assert len(dataset.graphs) == len(external_globals), (
        f"Length mismatch: {len(dataset.graphs)} graphs vs {len(external_globals)} rows of external data."
    )
    for idx, graph in enumerate(dataset.graphs):
        graph.u = torch.from_numpy(external_globals[idx]).float().unsqueeze(0)
    return dataset
