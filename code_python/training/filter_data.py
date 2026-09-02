from pathlib import Path
from typing import Optional, Union, Dict, Tuple

import torch
from all_factories import radius_to_bits, imputer_factory
import pandas as pd
from unrolling_utils import unrolling_factory
import numpy as np
import os
from Gpytorch_sskkernel import pad, build_one_hot, encode_string


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
tkwargs = {"dtype": torch.float, "device": device}

def cutoff_filteration(data: pd.DataFrame, lower_cutoff: Optional[float], upper_cutoff: Optional[float], target_feature: str) -> pd.DataFrame:
    if lower_cutoff is not None:
        data = data.drop(data.index[data[target_feature] < lower_cutoff])
    if upper_cutoff is not None:
        data = data.drop(data.index[data[target_feature] > upper_cutoff])
    return data.reset_index(drop=True)  # Reset index after filtering


def apply_cutoff(data: pd.DataFrame, cutoffs: Dict[str, Tuple[Optional[float], Optional[float]]]) -> pd.DataFrame:
    df = data.copy()
    for feature, (lower_cutoff, upper_cutoff) in cutoffs.items():
        df = cutoff_filteration(data=df, lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff, target_feature=feature)
    return df


def ensure_long_path(path: Path) -> Path:
    """Ensures Windows handles long paths by adding '\\?\' if needed."""
    path_str = str(path)
    if os.name == 'nt' and len(path_str) > 250:
        return Path(f"\\\\?\\{path_str}")
    return path


DATASET_CONFIG = {
    "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices": {
        "features": [
            "HOMO_D (eV)", "LUMO_D (eV)", "Eg_D (eV)", "Ehl_D (eV)",
            "HOMO_A (eV)", "LUMO_A (eV)", "Eg_A (eV)", "Ehl_A (eV)",
            "D:A ratio (m/m)", "solvent additive conc. (% v/v)",
            "temperature of thermal annealing",
            "HTL energy level (eV)", "ETL energy level (eV)"
        ],
        "polymer_unit": ["Donor", "Acceptor"],
        "target": ["calculated PCE (%)"],
    },

    "Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation": {
        "features": [
            "Xn", "Mw (g/mol)", "PDI", "Concentration (mg/ml)",
            "Temperature SANS/SLS/DLS/SEC (K)",
            "polymer dP", "polymer dD", "polymer dH",
            "solvent dP", "solvent dD", "solvent dH"
        ],
        "polymer_unit": ["Monomer"],
        "target": ["log Rg (nm)"],
    },

    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery": {
        "features": [
            'Selective_layer_thickness',
            'Filler_size',
            'Filler_concentration',
            'Crosslinker_concentration',
            'Organic_Hildebrand_solubility_parameter',
            'Organic_molecular_volume',
            'Effective_area',
            'Mass_ratio(A/B)',
            'Experimental_temperature',
            'Downstream_pressure',
        ],
        "polymer_unit": ["polymers"],
        "target": {
            "separation_data_imputed": ["log (Separation factor)"],
            "flux_data_imputed": ["log (Total flux)"]
        }
    },

    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning": {
        "features": [
            'P_MW','surface tension (mN/m)','P_wt%',
            'pore maker molecular weight (Da)',
            'pore maker wt% (pore maker to total weight)',
            'additive wt% (additive to polymer)',
            'solubility parameter (MPa1/2)',
            'TMP (bar)',
            'oraganic compound concentration (mg/L)',
            'organic compound size (Da)',
            'foulant concentration (mg/L)',
        ],
        "polymer_unit": ["polymer"],
        "target": [
            r'water permeability (LMH/bar)',
            'organic compound removal (%)',
            'flux decline ratio (%)',
            'flux recovery ratio (%)',
            'reversible fouling ratio (%)',
            'irreversible fouling ratio(%)',
        ],
    },

    "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation": {
        "features": [
            'Contact angle θ （°）',
            'Thickness ℓ (um)',
            'xw (wt%)',
            'Temperature (°C)',
            'Permeate pressure (mbar)',
            'Solvent solubility parameter （MPa1/2）',
        ],
        "polymer_unit": ["polymer"],
        "target": ['log (Total flux)', 'log (Separation factor)'],
    },
        
    "Miniaturization of Popular Reactions from the Medicinal Chemists Toolbox for Ultrahigh_Throughput Experimentation": {
        "features": [
            "Equiv AcOH", "Equiv TTIP", "Solvent_Encoded", "Reaction Conc (mM)"
        ],
        "polymer_unit": ["Substrate"],
        "target": ["Approx Conv (%)"],
    },

}


def _get_dataset_features(dataset_basket: Path, paper: str, dataset_name: str):
    dataset_path = ensure_long_path(dataset_basket / paper / f"{dataset_name}.pkl")

    if not dataset_path.exists():
        raise FileNotFoundError(f"Dataset '{dataset_name}' not found in {dataset_basket}")

    dataset: pd.DataFrame = pd.read_pickle(dataset_path)

    if paper not in DATASET_CONFIG:
        raise ValueError(f"No predefined features found for paper '{paper}'")

    entry = DATASET_CONFIG[paper]
    features: list[str] = entry["features"]

    target: list[str] = (
        entry["target"].get(dataset_name)
        if isinstance(entry["target"], dict)
        else entry["target"]
    )
    polymer_unit: list[str] = entry["polymer_unit"]

    return dataset, features, target, polymer_unit



def get_structural_info(fp:str,poly_unit_name:list[str],radius:int=None,vector:str=None)->Tuple:
       
        if fp == "Mordred" or fp == "MACCS":
            fp_features = [f"{unit}_{fp}" for unit in poly_unit_name]
            unrolling_featurs = {"representation": fp,
                                "unit_name":poly_unit_name,
                                "col_names": fp_features}
            return fp_features, unrolling_featurs
        
        if fp == "ECFP":
            n_bits = radius_to_bits[radius]
            fp_features = [
                        f"{unit}_{fp}{2 * radius}_{vector}_{n_bits}bits"
                        for unit in poly_unit_name
                        ]
            unrolling_featurs = {
                                "representation": fp,
                                "radius": radius,
                                "n_bits": n_bits,
                                "vector_type": vector,
                                "unit_name": poly_unit_name,
                                "col_names": fp_features,
                                }
            return fp_features, unrolling_featurs
        if fp == "HDF":
            n_bits = 64
            fp_features = [f"{unit} SMILES" for unit in poly_unit_name]
            unrolling_featurs = {
                                "representation": fp,
                                "radius": radius,
                                "n_bits": n_bits,
                                "unit_name": poly_unit_name,
                                "SMILES_columns": fp_features,
                                }
            return fp_features, unrolling_featurs
        
        if fp in {"SSK", "MG"}:
            fp_features = [f"{unit} SMILES" for unit in poly_unit_name]
            unrolling_featurs = {
                                "representation": fp,
                                "unit_name": poly_unit_name,
            }
            return fp_features, unrolling_featurs
        else:
            unrolling_featurs = {
                                "representation": None,
                                "unit_name": None,
                                }
            return None, unrolling_featurs


def sanitize_dataset(
    df: pd.DataFrame, target_feat:list, dropna: bool, **kwargs
) -> pd.DataFrame:
    """
    Sanitize the training features and targets in case the target features contain NaN values.

    Args:
        training_features: Training features.
        targets: Targets.
        dropna: Whether to drop NaN values.
        **kwargs: Keyword arguments to pass to filter_dataset.

    Returns:
        Sanitized training features and targets.
    """
    if dropna:
        df: pd.DataFrame = df.dropna(subset=target_feat)
        return df
    else:
        return df


def _ssk_emb(smiles):

    maxlen = np.max([len(x) for x in smiles])
    # get alphabet of characters used in candidate set (to init SSK)
    alphabet = list({l for word in smiles for l in word})
    print(f'alphabet \n {alphabet} \n length of alphabet {len(alphabet)}')
    print(f'maxlen {maxlen}')
    embds, index = build_one_hot(alphabet)
    embds = embds.to(**tkwargs)
    all_encoded_smiles = torch.cat([pad(encode_string(x, index), maxlen).unsqueeze(0) for x in smiles], dim=0)

    parameters = {
        "maxlen": maxlen,
        "alphabet": alphabet,
        "embds": embds,
        "index": index
    }
    return all_encoded_smiles, parameters


def filter_dataset(
    raw_dataset: pd.DataFrame,
    structure_feats: Optional[list[str]], # can be None
    scalar_feats: Optional[list[str]], # like conc, temp,
    target_feats: list[str], 
    cutoff: Dict[str, Tuple[Optional[float], Optional[float]]]=False,
    dropna: bool = True,
    unroll: Union[dict, list, None] = None,
    cluster_type:Optional[str]= None,
    feat_to_impute: Optional[list[str]]=None,
    imputer: Optional[str]=None,
    # **kwargs,
) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    """
    Filter the dataset.

    Args:
        raw_dataset: Raw dataset.
        structure_feats: Structure features.
        scalar_feats: Scalar features.
        target_feats: Target features.

    Returns:
        Input features and targets.
    """
    structure_feats = list(structure_feats or [])
    scalar_feats = list(scalar_feats or [])
    target_feats = list(target_feats or [])

    if not target_feats:
        raise ValueError("target_feats must contain at least one target column.")
    if not structure_feats and not scalar_feats:
        raise ValueError(
            "At least one structural or continuous feature is required."
        )

    # Add multiple lists together as long as they are non-empty.
    kernel_parameters = {}
    all_feats: list[str] = [
        feat
        for feat_list in [structure_feats, scalar_feats, target_feats]
        if feat_list
        for feat in feat_list
    ]
    if cluster_type:
        all_feats.append(cluster_type)
        if cluster_type=="substructure cluster":
            all_feats.append("Side Chain Cluster")

    # Avoid duplicate DataFrame columns if a caller accidentally repeats a
    # feature name while preserving the requested order.
    all_feats = list(dict.fromkeys(all_feats))
    dataset: pd.DataFrame = raw_dataset.loc[:, all_feats].copy()
    dataset = sanitize_dataset(dataset,
        target_feats, dropna=dropna
        )
    if cutoff:
        dataset = apply_cutoff(dataset,cutoff)

    if unroll is None:
        unroll_configs = []
    elif isinstance(unroll, dict):
        unroll_configs = [unroll]
    elif isinstance(unroll, list) and all(
        isinstance(config, dict) for config in unroll
    ):
        unroll_configs = unroll
    else:
        raise TypeError(
            "unroll must be a dictionary, a list of dictionaries, or None."
        )

    structural_frames: list[pd.DataFrame] = []

    def _aligned_frame(frame: pd.DataFrame) -> pd.DataFrame:
        if not isinstance(frame, pd.DataFrame):
            frame = pd.DataFrame(frame)
        if len(frame) != len(dataset):
            raise ValueError(
                "Unrolled structural features have a different number of rows "
                "than the filtered dataset."
            )
        if not frame.index.equals(dataset.index):
            frame = frame.copy()
            frame.index = dataset.index
        return frame

    for config in unroll_configs:
        config = config.copy()
        if isinstance(config.get("unit_name"), str):
            config["unit_name"] = [config["unit_name"]]
        representation = config.get("representation")

        if representation in {"Mordred", "MACCS", "ECFP", "HDF"}:
            if not structure_feats:
                raise ValueError(
                    f"structure_feats is required for {representation} features."
                )
            if representation in {"Mordred", "MACCS"}:
                config_units = config.get("unit_name")
                if isinstance(config_units, str):
                    config_units = [config_units]
                source_columns = config.get("col_names")
                if not config_units or not source_columns:
                    raise ValueError(
                        "unit_name and col_names are required for "
                        f"{representation}."
                    )
                if len(config_units) != len(source_columns):
                    raise ValueError(
                        "unit_name and col_names must contain the same number "
                        f"of entries for {representation}."
                    )

                for unit, source_column in zip(config_units, source_columns):
                    unit_config = {
                        **config,
                        "unit_name": unit,
                        "col_names": [source_column],
                    }
                    structural_frames.append(_aligned_frame(
                        unrolling_factory[representation](
                            dataset[[source_column]].reset_index(drop=True),
                            **unit_config,
                        )
                    ))
            else:
                structural_frames.append(_aligned_frame(
                    unrolling_factory[representation](
                        dataset[structure_feats], **config
                    )
                ))
            continue

        if representation is None:
            continue

        config_units = config.get("unit_name")
        if isinstance(config_units, str):
            config_units = [config_units]
        elif config_units is None:
            raise ValueError(
                f"unit_name is required for {representation} features."
            )

        for unit in config_units:
            feat_name = f"{unit} SMILES"
            if feat_name not in dataset.columns:
                raise ValueError(
                    f"Structural column {feat_name!r} is missing from the "
                    "dataset."
                )

            if representation.lower() == "ssk":
                all_encoded_smiles, parameters = _ssk_emb(
                    dataset[feat_name].values
                )
                kernel_parameters[f"fp_{unit}"] = parameters
                structural_frames.append(_aligned_frame(pd.DataFrame(
                    all_encoded_smiles,
                    index=dataset.index,
                    columns=[
                        f"{unit}_ssk_emb_{i}"
                        for i in range(all_encoded_smiles.shape[1])
                    ],
                )))
            else:
                structural_frames.append(dataset[[feat_name]])

    structure_features: pd.DataFrame = (
        pd.concat(structural_frames, axis=1)
        if structural_frames
        else dataset[[]]
    )
    if structure_feats and structure_features.shape[1] == 0:
        raise ValueError(
            "Structural features were requested, but no structural columns "
            "were produced. Check the unroll configuration."
        )

    if scalar_feats:
        scalar_features: pd.DataFrame = dataset.loc[:, scalar_feats].copy()
        if feat_to_impute:
            missing_impute_cols = [
                feature for feature in feat_to_impute
                if feature not in scalar_features.columns
            ]
            if missing_impute_cols:
                raise ValueError(
                    "Columns requested for imputation are not continuous "
                    f"features: {missing_impute_cols}."
                )
            if imputer not in imputer_factory:
                raise ValueError(f"Unknown imputer: {imputer!r}.")
            scalar_features[feat_to_impute] = imputer_factory[
                imputer
            ].fit_transform(scalar_features[feat_to_impute])

    else:
        scalar_features: pd.DataFrame = dataset[[]]

    training_features: pd.DataFrame = pd.concat(
        [structure_features, scalar_features], axis=1
    )

    if not training_features.index.equals(dataset.index):
        raise RuntimeError(
            "Feature rows became misaligned while constructing the training "
            "dataset."
        )
    if training_features.columns.has_duplicates:
        duplicate_cols = training_features.columns[
            training_features.columns.duplicated()
        ].tolist()
        raise ValueError(f"Duplicate training feature columns: {duplicate_cols}.")

    targets = dataset.loc[:, target_feats].to_numpy()

    # if not (scalars_available and struct_available):
    new_struct_feats: list[str] = structure_features.columns.tolist()

    if cluster_type:
        if cluster_type == "substructure cluster":
            substructure_labels = dataset[cluster_type].squeeze().to_numpy()
            side_chain_labels = dataset["Side Chain Cluster"].squeeze().to_numpy()
            c_labels = {"substructure cluster": substructure_labels,
                        "Side Chain Cluster": side_chain_labels}
        else:
            c_labels = dataset[cluster_type].squeeze().to_numpy()
        return training_features, targets, new_struct_feats, kernel_parameters, c_labels

    return training_features, targets, new_struct_feats, kernel_parameters 
