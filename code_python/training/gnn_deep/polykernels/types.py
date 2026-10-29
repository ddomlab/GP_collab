"""Type definitions and tensor shape aliases using jaxtyping.

Usage in function signatures::

    from polykernels.types import FloatTensor

    def forward(self, x: FloatTensor["B D"]) -> FloatTensor["B C"]:  # noqa: F722
        ...
"""

from typing import TypeAlias

import jaxtyping as jt
from jaxtyping import Float
import numpy as np
from torch import Tensor

FloatTensor: TypeAlias = Float[Tensor, "..."]


class IntArrayMeta(type):
    """Metaclass to make IntArray generic over jaxtyping shape strings."""

    def __getitem__(cls, shape: str):
        return jt.Int[np.ndarray, shape]


class IntArray(metaclass=IntArrayMeta):
    """Subscriptable type alias for 1-D integer numpy arrays with explicit shapes."""

    pass
