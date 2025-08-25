import json
from pathlib import Path
from typing import Optional, Union, Dict, Tuple

import pandas as pd
from unrolling_utils import unrolling_factory
import numpy as np
# HERE: Path = Path(__file__).resolve().parent
# DATASETS: Path = HERE.parent.parent / "datasets"

# with open(HERE / "filters.json", "r") as f:
#     FILTERS: dict[str, list[str]] = json.load(f)

# with open(HERE / "subsets.json", "r") as f:
#     SUBSETS: dict[str, list[str]] = json.load(f)



# Filtration with optional lower and upper cutoff and resetting index
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
    new_num_feats: list[str] = scalar_features.columns.tolist()
    training_test_shape: Dict ={
                                "targets_shape": targets.shape,
                                "training_features_shape": training_features.shape
                                }
    
    if cluster_type:
        if cluster_type == "substructure cluster":
            substructure_labels = dataset[cluster_type].squeeze().to_numpy()
            side_chian_labels = dataset["Side Chain Cluster"].squeeze().to_numpy()
            c_labels = {"substructure cluster": substructure_labels,
                        "Side Chain Cluster": side_chian_labels}
        else:    
            c_labels = dataset[cluster_type].squeeze().to_numpy()
        return training_features, targets, new_struct_feats, c_labels,training_test_shape
    
    return training_features, targets, new_struct_feats,training_test_shape