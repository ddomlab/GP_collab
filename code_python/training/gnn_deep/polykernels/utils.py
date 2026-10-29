"""Generic, reusable utility functions for the project."""

import dataclasses
import datetime
import enum
import os
from pathlib import Path
import random
import sys
import typing

from loguru import logger
from loguru._logger import Logger
import numpy as np
import optuna
import torch

from polykernels.paths_and_constants import Paths

_P = Paths()


def round3(val: float | int | np.number | torch.Tensor) -> float:
    """Round a number, numpy scalar, or tensor to 3 decimal places."""
    return round(float(val), 3)


def set_seed(seed: int = 42) -> None:
    """Set random seed for reproducibility across random, numpy, and PyTorch."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    # Full fp32 matmuls on a GPU too (no TF32). TF32 rounds to ~1e-3 relative, which degrades
    # the GP kernel matrices and breaks the SHAP / GP-SHAP / IG sum checks.
    torch.set_float32_matmul_precision("highest")
    # Deterministic GPU kernels, so rerunning a job reproduces its model exactly, as on CPU.
    # Otherwise the atomic adds in message passing differ run to run, and a GNN fold's
    # predictions moved by up to 1.6 (on a target with std 32) between identical runs.
    # cuBLAS reads its workspace setting when CUDA initializes, which comes after this.
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    torch.use_deterministic_algorithms(True, warn_only=True)


def cast_random_state(
    seed: int | np.random.RandomState | None = None,
) -> np.random.RandomState:
    """Cast seed to a numpy.random.RandomState instance."""
    if seed is None or isinstance(seed, (int, np.integer)):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError(
        f"{seed} cannot be used to seed a numpy.random.RandomState instance"
    )


def get_logger() -> Logger:
    """Configure and return a custom loguru logger."""
    logger.remove()
    custom_format = (
        "<green>{time:%y-%m-%d %H:%M:%S}</green> | "
        "<level>{level.name:.1s}</level> | "
        "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - "
        "<level>{message}</level>"
    )
    logger.add(sys.stderr, format=custom_format, colorize=True)
    return logger


def enum_values[E](enum_cls: type[E]) -> list[str]:
    """Return a list of string values for a StrEnum class."""
    return [e.value for e in enum_cls]


def _format_float(val: float) -> str:
    """Format float: scientific notation for very large/small values, else 4 decimals."""
    if val == 0.0:
        return "0.0"
    abs_val = abs(val)
    if abs_val < 1e-3 or abs_val >= 1e4:
        return f"{val:.4e}"
    s = f"{val:.4f}".rstrip("0").rstrip(".")
    return s if "." in s or "e" in s else f"{s}.0"


def _format_dataclass(val: typing.Any, indent: int, spaces: int) -> str:
    """Format a dataclass recursively into an indented representation."""
    prefix = " " * indent
    sub_prefix = " " * (indent + spaces)
    fields = dataclasses.fields(val)
    if not fields:
        return f"{val.__class__.__name__}()"
    lines = []
    for f in fields:
        v = getattr(val, f.name)
        formatted_v = format_config(v, indent + spaces, spaces)
        lines.append(f"{sub_prefix}{f.name}={formatted_v.lstrip()}")
    joined = ",\n".join(lines) + ","
    return f"{val.__class__.__name__}(\n{joined}\n{prefix})"


def _format_dict(val: dict[typing.Any, typing.Any], indent: int, spaces: int) -> str:
    """Format a dictionary recursively into an indented representation."""
    prefix = " " * indent
    sub_prefix = " " * (indent + spaces)
    if not val:
        return "{}"
    lines = []
    for k, v in val.items():
        formatted_v = format_config(v, indent + spaces, spaces)
        lines.append(f"{sub_prefix}{repr(k)}: {formatted_v.lstrip()}")
    joined = ",\n".join(lines) + ","
    return f"{{\n{joined}\n{prefix}}}"


def _format_sequence(
    val: list[typing.Any] | tuple[typing.Any, ...] | set[typing.Any],
    indent: int,
    spaces: int,
) -> str:
    """Format sequence, choosing multiline for complex/long structures."""
    prefix = " " * indent
    sub_prefix = " " * (indent + spaces)
    if not val:
        empty_reprs = {list: "[]", tuple: "()", set: "set()"}
        return empty_reprs[type(val)]

    is_complex = any(
        dataclasses.is_dataclass(x) or isinstance(x, (dict, list, tuple, set))
        for x in val
    )
    if is_complex or len(val) > 10:
        lines = []
        for x in val:
            formatted_x = format_config(x, indent + spaces, spaces)
            lines.append(f"{sub_prefix}{formatted_x.lstrip()}")
        joined = ",\n".join(lines) + ","
        start, end = {list: ("[", "]"), tuple: ("(", ")"), set: ("{", "}")}[type(val)]
        return f"{start}\n{joined}\n{prefix}{end}"

    items = [format_config(x, 0, spaces) for x in val]
    start, end = {list: ("[", "]"), tuple: ("(", ")"), set: ("{", "}")}[type(val)]
    return f"{start}{', '.join(items)}{end}"


def format_config(val: typing.Any, indent: int = 0, spaces: int = 2) -> str:
    """Format configuration objects into a highly readable string.

    Floating point numbers are formatted with scientific notation if they are
    very large (>= 1e4) or very small (< 1e-3), otherwise they are rounded
    to at most 4 decimal places.
    """
    if dataclasses.is_dataclass(val):
        return _format_dataclass(val, indent, spaces)
    if isinstance(val, dict):
        return _format_dict(val, indent, spaces)
    if isinstance(val, (list, tuple, set)):
        return _format_sequence(val, indent, spaces)
    if isinstance(val, bool):
        return str(val)
    if isinstance(val, (int, np.integer)):
        return str(int(val))
    if isinstance(val, (float, np.floating)):
        return _format_float(float(val))
    if isinstance(val, Path):
        return repr(str(val))
    if isinstance(val, enum.Enum):
        return f"{val.__class__.__name__}.{val.name}"
    return repr(val)


def humanize_int(num: int) -> str:
    """Format parameter count into a human-readable string (e.g., 105k, 1.2M)."""
    if num >= 1_000_000:
        val = num / 1_000_000
        return f"{val:.1f}M".replace(".0M", "M")
    if num >= 1_000:
        val = num / 1_000
        return f"{val:.1f}k".replace(".0k", "k")
    return str(num)


class EtaCallback:
    """Callback to calculate and log the estimated time remaining for HPO."""

    def __init__(self, n_trials: int):
        self.n_trials = n_trials

    def __call__(self, study: typing.Any, trial: typing.Any) -> None:
        completed = [
            t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE
        ]
        n_completed = len(completed)
        if n_completed == 0:
            return

        durations = [
            (t.datetime_complete - t.datetime_start).total_seconds()
            for t in completed
            if t.datetime_complete and t.datetime_start
        ]
        if not durations:
            return

        avg_sec = sum(durations) / len(durations)
        remaining = self.n_trials - n_completed

        if remaining <= 0:
            logger.info(
                f"--- HPO Progress: {n_completed}/{self.n_trials} trials complete (Final Trial) ---"
            )
            return

        eta = str(datetime.timedelta(seconds=int(remaining * avg_sec)))
        logger.info(
            f"--- HPO Progress: {n_completed}/{self.n_trials} trials complete --- "
            f"Avg Duration: {avg_sec:.1f}s | "
            f"Estimated Time Remaining: {eta}"
        )
