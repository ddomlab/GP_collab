import json
from pathlib import Path
import sys
from typing import Optional, Union, Dict, Tuple
from all_factories import radius_to_bits,cutoffs
import pandas as pd
from unrolling_utils import unrolling_factory
import numpy as np
import os


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
    }
}


def _get_dataset_features(dataset_basket: Path, paper: str, dataset_name: str):
    dataset_path = ensure_long_path(dataset_basket / paper / f"{dataset_name}.pkl")

    if not dataset_path.exists():
        raise FileNotFoundError(f"Dataset '{dataset_name}' not found in {dataset_basket}")

    dataset = pd.read_pickle(dataset_path)

    if paper not in DATASET_CONFIG:
        raise ValueError(f"No predefined features found for paper '{paper}'")

    entry = DATASET_CONFIG[paper]
    features = entry["features"]

    target = (
        entry["target"].get(dataset_name)
        if isinstance(entry["target"], dict)
        else entry["target"]
    )
    polymer_unit = entry["polymer_unit"]

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
        if fp == "SMILES":
            return ["Monomer SMILES"], None
        else:
              return None, None


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


def filter_dataset(
    raw_dataset: pd.DataFrame,
    structure_feats: Optional[list[str]], # can be None
    scalar_feats: Optional[list[str]], # like conc, temp,
    target_feats: list[str], 
    cutoff: Dict[str, Tuple[Optional[float], Optional[float]]]=False,
    dropna: bool = True,
    unroll: Union[dict, list, None] = None,
    cluster_type:Optional[str]= None,
    **kwargs,
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
    # Add multiple lists together as long as they are not NoneType
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

    dataset: pd.DataFrame = raw_dataset[all_feats]
    dataset = sanitize_dataset(dataset,
        target_feats, dropna=dropna, **kwargs
        )
    if cutoff:
        dataset = apply_cutoff(dataset,cutoff)

    if unroll:
        if isinstance(unroll, dict):
            structure_features: pd.DataFrame = unrolling_factory[
                unroll["representation"]](dataset[structure_feats], **unroll)
        elif isinstance(unroll, list):
            multiple_unrolled_structure_feats: list[pd.DataFrame] = []
            for unroll_dict in unroll:
                single_structure_feat: pd.DataFrame = filter_dataset(
                    dataset,
                    # structure_feats=unroll_dict["columns"],
                    structure_feats=unroll_dict["col_names"],
                    scalar_feats=[],
                    target_feats=[],
                    # dropna=dropna,
                    dropna=False,
                    unroll=unroll_dict,
                )[0]
                multiple_unrolled_structure_feats.append(single_structure_feat)
            structure_features: pd.DataFrame = pd.concat(
                multiple_unrolled_structure_feats, axis=1
            )
        else:
            raise ValueError(f"Unroll must be a dict or list, not {type(unroll)}")
    

    
    elif structure_feats:
        structure_features: pd.DataFrame = dataset[structure_feats]
    else:
        structure_features: pd.DataFrame = dataset[[]]

    if scalar_feats:
        scalar_features: pd.DataFrame = dataset[scalar_feats]
    else:
        scalar_features: pd.DataFrame = dataset[[]]

    training_features: pd.DataFrame = pd.concat(
        [structure_features, scalar_features], axis=1
    )

    targets = dataset[target_feats].squeeze()
    targets = np.vstack(targets.values)

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
        return training_features, targets, new_struct_feats, c_labels

    return training_features, targets, new_struct_feats