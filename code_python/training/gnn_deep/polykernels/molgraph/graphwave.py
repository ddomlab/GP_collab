"""Graphwave diffusion-based structural node features using numpy and scipy."""

import numpy as np
import scipy.sparse.csgraph as spgraph


def _get_graph_spectrum(adj: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lap = spgraph.laplacian(adj)
    eigen_values, eigen_vectors = np.linalg.eigh(lap)
    max_eigen = np.max(eigen_values)
    if max_eigen > 0:
        eigen_values = eigen_values / max_eigen
    return eigen_values, eigen_vectors


def _heat_wavelets(
    eigen_values: np.ndarray, eigen_vectors: np.ndarray, scale: float
) -> np.ndarray:
    kernel = np.diag(np.exp(-scale * eigen_values))
    return eigen_vectors @ kernel @ eigen_vectors.T


def _characteristic_function(
    wavelet_basis: np.ndarray, t: np.ndarray
) -> np.ndarray:
    scaled_waves = np.einsum("ij,k->ijk", wavelet_basis, t)
    return np.mean(np.exp(scaled_waves * 1j), axis=0)


def get_graphwave_embeddings(
    adj: np.ndarray, ndim: int, step_size: float, scale: float
) -> np.ndarray:
    sample_number = ndim // 2
    sample_number = max(sample_number, 1)
    t = np.arange(sample_number) * step_size
    eigen_values, eigen_vectors = _get_graph_spectrum(adj)
    wavelet_basis = _heat_wavelets(eigen_values, eigen_vectors, scale)
    wavelets = _characteristic_function(wavelet_basis, t)
    # wavelets shape: [num_nodes, sample_number]
    # We return real and imag parts stacked horizontally to get [num_nodes, ndim]
    return np.hstack([wavelets.real, wavelets.imag])
