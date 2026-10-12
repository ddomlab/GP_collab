"""Add label-encoded and Mordred molecular representations to datasets."""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from mordred import Calculator
import mordred.descriptors
from rdkit.Chem import Mol, MolFromSmiles, MolToSmiles
from sklearn.preprocessing import OrdinalEncoder


HERE = Path(__file__).resolve().parent
DATASETS = HERE.parent.parent / "datasets" / "Validation datasets"

# ``filter_data.py`` uses imports relative to the training directory, so make
# that directory importable when this file is run directly.
TRAINING_DIR = HERE.parent / "training"
if str(TRAINING_DIR) not in sys.path:
    sys.path.insert(0, str(TRAINING_DIR))

from filter_data import _get_dataset_features  # noqa: E402


DATASET_GROUPS: tuple[str, ...] = (
    "Beyond molecular structure_ critically assessing machine learning for "
    "designing organic photovoltaic materials and devices|"
    "Beyond molecular structure_seifrid_imputed",
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based "
    "Organic Recovery|flux_data_imputed",
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based "
    "Organic Recovery|separation_data_imputed",
    "Machine Learning-Enabled Prediction and High-Throughput Screening of "
    "Polymer Membranes for Pervaporation Separation|"
    "cleaned_dataset_pervaporation_membranes_wang",
    "Miniaturization of Popular Reactions from the Medicinal Chemists Toolbox "
    "for Ultrahigh_Throughput Experimentation|cleaned_suzuki_synthesis",
    "Understanding and Designing a High-Performance Ultrafiltration Membrane "
    "Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed",
)


def _representation_name(smiles_column: str) -> str:
    """Return the structural-unit name from a ``<unit> SMILES`` column."""
    suffix = " SMILES"
    if smiles_column.endswith(suffix):
        return smiles_column[: -len(suffix)]
    return smiles_column


def canonicalize_columns(data: pd.DataFrame, column: str) -> pd.DataFrame:
    """Return a copy with one SMILES column converted to canonical SMILES."""
    if column not in data.columns:
        raise KeyError(f"SMILES column {column!r} is missing from the dataset.")

    canonical_smiles: list[str] = []
    invalid_rows: list[object] = []
    for index, smiles in data[column].items():
        mol = (
            MolFromSmiles(smiles)
            if isinstance(smiles, str) and smiles.strip()
            else None
        )
        if mol is None:
            invalid_rows.append(index)
            canonical_smiles.append("")
        else:
            canonical_smiles.append(MolToSmiles(mol, canonical=True))

    if invalid_rows:
        preview = ", ".join(map(str, invalid_rows[:10]))
        remainder = " ..." if len(invalid_rows) > 10 else ""
        raise ValueError(
            f"Column {column!r} contains {len(invalid_rows)} missing or invalid "
            f"SMILES value(s) at row index(es): {preview}{remainder}"
        )

    new_data = data.copy()
    new_data[column] = canonical_smiles
    return new_data


def get_mordred_dict(
    mol: Mol, calculator: Calculator | None = None
) -> dict[str, object]:
    """Calculate all two-dimensional Mordred descriptors for an RDKit molecule."""
    if mol is None:
        raise ValueError(
            "Cannot calculate Mordred descriptors for an invalid molecule."
        )
    calc = calculator or Calculator(mordred.descriptors, ignore_3D=True)
    return {str(name): value for name, value in calc(mol).asdict().items()}


class MordredCalculator:
    """Calculate and attach Mordred descriptors for one SMILES column."""

    def __init__(
        self,
        smile_source: pd.DataFrame,
        oligomer_represenation: str = "SMILES",
    ) -> None:
        if oligomer_represenation not in smile_source.columns:
            raise KeyError(
                f"SMILES column {oligomer_represenation!r} is missing from the dataset."
            )
        self.smile_source = smile_source
        self.oligomer_represenation = oligomer_represenation
        self.calculator = Calculator(mordred.descriptors, ignore_3D=True)

    def assign_Mordred(self) -> pd.DataFrame:
        """Add a ``<unit>_Mordred`` column and return the augmented frame."""
        smiles = self.smile_source[self.oligomer_represenation]
        unique_smiles = pd.unique(smiles)

        descriptor_records: list[dict[str, object]] = []
        for value in unique_smiles:
            mol = MolFromSmiles(value)
            if mol is None:  # Defensive: canonicalize_columns normally catches this.
                raise ValueError(
                    f"Invalid SMILES {value!r} in {self.oligomer_represenation!r}."
                )
            descriptor_records.append(get_mordred_dict(mol, self.calculator))

        # Expand the cached results back to the original rows before filtering.
        # This preserves the exact behavior of the original implementation while
        # avoiding recalculation when a molecule occurs more than once.
        raw_descriptor_lookup = dict(zip(unique_smiles, descriptor_records))
        descriptors = pd.DataFrame.from_records(
            smiles.map(raw_descriptor_lookup), index=smiles.index
        )

        # Keep the original Mordred filtering sequence: calculation failures are
        # object-valued, descriptors containing NaN are removed, and descriptors
        # whose sample variance is exactly zero are removed.
        descriptors = descriptors.infer_objects()
        descriptors = descriptors.select_dtypes(exclude=["object"])
        descriptors = descriptors.dropna(axis=1, how="any")
        descriptor_variances = descriptors.var(numeric_only=True)
        invariant_descriptors = descriptor_variances[
            descriptor_variances.eq(0)
        ].index
        descriptors = descriptors.drop(columns=invariant_descriptors)
        if descriptors.empty:
            raise ValueError(
                f"No usable Mordred descriptors were generated for "
                f"{self.oligomer_represenation!r}."
            )

        output_column = f"{_representation_name(self.oligomer_represenation)}_Mordred"
        result = self.smile_source.copy()
        result[output_column] = descriptors.to_dict(orient="records")
        print(
            f"Assigned {len(descriptors.columns)} Mordred descriptors to "
            f"{output_column} ({len(unique_smiles)} unique molecules)."
        )
        return result


def label_encod(data: pd.DataFrame, structure_columns: list[str]) -> pd.DataFrame:
    """Add a separate integer label encoding for each SMILES column."""
    missing = [column for column in structure_columns if column not in data.columns]
    if missing:
        raise KeyError(f"Missing structure column(s): {missing}")
    if not structure_columns:
        return data.copy()

    df = data.copy()
    # Each input feature gets its own vocabulary. This is important for data
    # sets containing distinct molecule roles such as donor and acceptor.
    encoder = OrdinalEncoder(dtype=np.int64)
    encoded = encoder.fit_transform(df[structure_columns])
    for position, column in enumerate(structure_columns):
        df[f"{_representation_name(column)}_LE"] = encoded[:, position]
    return df


def _split_dataset_group(group: str) -> tuple[str, str]:
    try:
        paper, dataset_name = group.split("|", maxsplit=1)
    except ValueError as exc:
        raise ValueError(
            f"Invalid dataset group {group!r}; expected '<paper>|<dataset_name>'."
        ) from exc
    if not paper or not dataset_name:
        raise ValueError(
            f"Invalid dataset group {group!r}; paper and dataset name are required."
        )
    return paper, dataset_name


def _load_dataset_features(
    paper: str, dataset_name: str
) -> tuple[pd.DataFrame, list[str], list[str], list[str]]:
    """Load through the training helper, including NumPy 1/2 pickle compatibility."""
    try:
        return _get_dataset_features(DATASETS, paper, dataset_name)
    except ModuleNotFoundError as exc:
        # Pickles created with NumPy 2 refer to ``numpy._core``. NumPy 1 uses
        # ``numpy.core`` for the same modules, so alias them and retry without
        # degrading ndarray columns by round-tripping through CSV.
        if exc.name not in {"numpy._core", "numpy._core.numeric"}:
            raise
        import numpy.core as numpy_core
        import numpy.core.numeric as numpy_core_numeric

        sys.modules.setdefault("numpy._core", numpy_core)
        sys.modules.setdefault("numpy._core.numeric", numpy_core_numeric)
        return _get_dataset_features(DATASETS, paper, dataset_name)


def _encode_dataset_group(group: str) -> tuple[str, Path]:
    """Encode and save one dataset; designed to run in a worker process."""
    paper, dataset_name = _split_dataset_group(group)
    dataset, _, _, molecular_units = _load_dataset_features(paper, dataset_name)
    structure_columns = [f"{unit} SMILES" for unit in molecular_units]

    for column in structure_columns:
        dataset = canonicalize_columns(dataset, column)
    dataset = label_encod(dataset, structure_columns)
    for column in structure_columns:
        dataset = MordredCalculator(dataset, column).assign_Mordred()

    destination = DATASETS / paper / dataset_name
    dataset.to_pickle(destination.with_suffix(".pkl"))
    dataset.to_csv(destination.with_suffix(".csv"), index=False)
    return group, destination


def encode_structures(
    dataset_groups: Iterable[str] = DATASET_GROUPS,
    max_workers: int | None = None,
) -> dict[str, Path]:
    """Encode datasets in parallel and return their output locations.

    Three workers are used by default to balance speed against the substantial
    memory needed by expanded Mordred dictionaries. Pass ``max_workers=1`` for
    sequential execution or a larger value on a machine with sufficient RAM.
    """
    groups = tuple(dataset_groups)
    if not groups:
        return {}

    if max_workers is None:
        max_workers = min(3, len(groups), os.cpu_count() or 1)
    if max_workers == 0:
        raise ValueError("max_workers cannot be 0; use -1 for all CPU cores.")
    if max_workers > 0:
        max_workers = min(max_workers, len(groups))

    encoded_datasets: dict[str, Path] = {}
    if max_workers == 1:
        for group in groups:
            completed_group, destination = _encode_dataset_group(group)
            encoded_datasets[completed_group] = destination
            print(f"Saved encoded dataset: {destination}")
        return encoded_datasets

    completed_datasets = Parallel(n_jobs=max_workers, prefer="processes")(
        delayed(_encode_dataset_group)(group) for group in groups
    )
    for completed_group, destination in completed_datasets:
        encoded_datasets[completed_group] = destination
        print(f"Saved encoded dataset: {destination}")

    return encoded_datasets


if __name__ == "__main__":
    # Import this file by its real module name before starting joblib workers.
    # Functions defined in ``__main__`` are not reliably pickleable by process
    # backends when the file is launched directly from an IDE.
    import molecular_encoder

    molecular_encoder.encode_structures(max_workers=-1)
