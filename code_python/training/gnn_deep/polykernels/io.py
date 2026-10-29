"""Input/Output and serialization utilities for numpy arrays and datasets."""

from collections.abc import Sequence
import lzma
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.api.types import is_integer_dtype, is_string_dtype


def save_npz(fname: Path, **kwargs: np.ndarray) -> None:
    """Save compressed numpy arrays using lzma serialization."""
    with lzma.open(fname, "wb") as afile:
        np.savez(file=afile, **kwargs)


def load_npz(fname: Path) -> dict[str, np.ndarray]:
    """Load lzma-compressed npz dataset files into a dictionary."""
    with lzma.open(fname, "rb") as afile:
        data = dict(np.load(afile, allow_pickle=False))
    return data


def _is_valid_index(index: pd.Index) -> bool:
    """Check if the index is purely string or integer type and not object."""
    return (is_integer_dtype(index) or is_string_dtype(index)) and index.dtype != object


class TensorMap:
    """Mapper for indexing and retrieving precomputed tensor values using a 1-D index."""

    def __init__(
        self,
        index_data: Sequence[str] | Sequence[int] | np.ndarray,
        values: np.ndarray,
    ) -> None:
        self.values = values
        index = pd.Index(index_data)
        if not _is_valid_index(index):
            raise TypeError(
                f"Index data must be of string or integer type (avoid object), got index with dtype '{index.dtype}'."
            )
        self.mapper = pd.Series(data=np.arange(len(values)), index=index)

    def save(
        self,
        fname: Path,
        index_key: str = "index",
        values_key: str = "values",
    ) -> None:
        """Serialize index and values array to a compressed npz file."""
        dtype = str if is_string_dtype(self.mapper.index) else int
        data = {
            index_key: self.mapper.index.to_numpy(dtype=dtype),
            values_key: self.values,
        }
        save_npz(fname, **data)

    @classmethod
    def load(
        cls,
        fname: Path,
        index_key: str = "index",
        values_key: str = "values",
    ) -> "TensorMap":
        """Load and instantiate a TensorMap from a compressed npz file."""
        data = load_npz(fname)
        if index_key not in data:
            raise KeyError(f"Index key '{index_key}' not found in npz data.")
        if values_key not in data:
            raise KeyError(f"Values key '{values_key}' not found in npz data.")
        return cls(data[index_key], data[values_key])

    def __call__(self, query: str | int | Sequence[str] | Sequence[int] | np.ndarray) -> np.ndarray:
        """Query the mapper to retrieve corresponding values."""
        return self.values[self.mapper[query]]
