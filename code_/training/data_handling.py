import json
import os
from pathlib import Path
from types import NoneType
from typing import Dict, Optional,Tuple

import numpy as np
import pandas as pd
from all_factories import radius_to_bits

HERE: Path = Path(__file__).resolve().parent
ROOT: Path = HERE.parent.parent

feature_abbrev: Dict[str, str] = {
    "Lp (nm)":          "Lp",
    "Rg1 (nm)":         "Rg",
    "Rh (IW avg log)":  "Rh",
    "Concentration (mg/ml)":            "concentration",
    "Temperature SANS/SLS/DLS/SEC (K)":         "temperature",
    "Mn (g/mol)":         "Mn", 
    "Mw (g/mol)":         "Mw", 
    "First Peak":         "Rh First Peak",
    "Second Peak":         "Rh Second Peak",
    "First Peak":         "Rh First Peak",
    "Third Peak":         "Rh Third Peak",
    "canonical_name": "Polymers cluster",
}


def remove_unserializable_keys(d: Dict) -> Dict:
    """Remove unserializable keys from a dictionary."""
    new_d: dict = {k: v for k, v in d.items() if
                   isinstance(v, (str, int, float, bool, NoneType, np.floating, np.integer))}
    return new_d


class NumpyArrayEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, tuple):
            return list(obj)
        else:
            return super(NumpyArrayEncoder, self).default(obj)

def get_cv_splits(score_for_indices):
    indices = {}
    for seed, values in score_for_indices.items():
        if isinstance(seed, int) and "indices" in values:
            indices[seed] = values["indices"]
    return indices

def ensure_long_path(path):
        """Ensures Windows handles long paths by adding '\\?\' if needed."""
        path_str = str(path)
        if os.name == 'nt' and len(path_str) > 250:  
            return Path(f"\\\\?\\{path_str}")
        return path



def _save(scores: Optional[Dict[int, Dict[str, float]]],
          predictions: Optional[pd.DataFrame],
          ground_truth: Optional[Dict],
          feat_length_scale: Optional[Dict],
          results_dir: Path,
          regressor_type: str,
          imputer: Optional[str],
          representation: str,
          pu_type : Optional[str],
          radius : Optional[int],
          vector : Optional[str],
          numerical_feats: Optional[list[str]],
          hypop: bool=True,
          transform_type:Optional[str]=None,
          target_transformer:Optional[str]=None,
          learning_curve:bool=False,
          special_file_name:Optional[str]=None,
          ) -> None:
    
    results_dir.mkdir(parents=True, exist_ok=True)
    
    # just scaler
    if numerical_feats and pu_type==None:
        short_num_feats = 'all_num'
        fname_root = f"({short_num_feats})_{regressor_type}"
        fname_root = f"{fname_root}_{imputer}" if imputer else fname_root
    
    # scaler and structural mixed
    if numerical_feats and pu_type:
        short_num_feats = 'all_num'
        if radius:
            fname_root = f"({representation}{radius}_{vector}_{radius_to_bits[radius]}-{short_num_feats})_{regressor_type}"
        
        else:
            fname_root = f"({representation}-{short_num_feats})_{regressor_type}"
            fname_root = f"{fname_root}_{imputer}" if imputer else fname_root

    #just structural
    if numerical_feats==None and pu_type: 
        if radius:
            fname_root = f"({representation}{radius}_{vector}_{radius_to_bits[radius]})_{regressor_type}"

        else:
            fname_root = f"({representation})_{regressor_type}"
    
    fname_root =f"{fname_root}_hypOFF" if hypop==False else fname_root
    fname_root =f"{fname_root}_{transform_type}" if transform_type else f"{fname_root}_feat_transformerOFF"
    fname_root = f"{fname_root}_{target_transformer}" if target_transformer else f"{fname_root}_target_transformerOFF"
    fname_root = f"{fname_root}_lc" if learning_curve else fname_root
    fname_root = f"{fname_root}_{special_file_name}" if special_file_name else fname_root
    print("Filename:", fname_root)

    if scores:
        # cv_indices = get_cv_splits(scores)
        # print("CV Indices:", cv_indices)
        scores_file: Path = results_dir / f"{fname_root}_scores.json"
        with open(scores_file, "w") as f:
            json.dump(scores, f, cls=NumpyArrayEncoder, indent=2)
        
        # indices_file: Path = results_dir / f"{fname_root}_indices.json"
        # with open(indices_file, "w") as f:
        #     json.dump(cv_indices, f, cls=NumpyArrayEncoder, indent=2)

    if predictions is not None:
        if isinstance(predictions, pd.DataFrame):
            predictions_file: Path = results_dir / f"{fname_root}_predictions.csv"
            predictions.to_csv(predictions_file, index=False)
        elif isinstance(predictions, dict):
            predictions_file: Path = results_dir / f"{fname_root}_predictions.json"
            
            with open(predictions_file, "w") as f:
                json.dump(predictions, f, cls=NumpyArrayEncoder, indent=2)


    if ground_truth:
        # print(ground_truth)
        cluster_ground_truth:Path = results_dir / f"{fname_root}_ClusterTruth.json"
        with open(cluster_ground_truth, "w") as f:
            json.dump(ground_truth, f, cls=NumpyArrayEncoder, indent=2)
    
    if feat_length_scale:
        data_shape_file:Path = results_dir / f"{fname_root}_length_scale.json"
        with open(data_shape_file, "w") as f:
            json.dump(feat_length_scale, f, cls=NumpyArrayEncoder, indent=2)

    
    print('Done Saving scores!')


def save_results(scores:Optional[Dict[int, Dict[str, float]]]=None,
                 predictions: Optional[pd.DataFrame]=None,
                 ground_truth: Optional[pd.DataFrame]=None,
                 feat_length_scale:Optional[Dict]=None,
                 target_features: list=None,
                 regressor_type: str=None,
                 kernel: Optional[str]=None,
                 TEST : bool =True,
                 representation: str=None,
                 pu_type : Optional[str]=None,
                 radius : Optional[int]=None,
                 vector : Optional[str]=None,
                 numerical_feats: Optional[list[str]]=None,
                 imputer: Optional[str] = None,
                 output_dir_name: str = "results",
                 cutoff: Optional[Dict[str, Tuple[Optional[float], Optional[float]]]] =None,
                 hypop: Optional[bool]=True,
                 transform_type:Optional[str]=None,
                 target_transformer:Optional[str]=None,
                 classification:bool=False,
                 clustering_method:str=None,
                 learning_curve:bool=False,
                 special_folder_name:Optional[str]=None,
                 special_file_name:Optional[str]=None,
                 ) -> None:
    
    targets_dir: str = "-".join([feature_abbrev.get(target, target) for target in target_features])
    feature_ids = []
    
    regressor_type = f"{regressor_type}({kernel})" if kernel is not None else regressor_type
    
    if pu_type:
        feature_ids.append(pu_type)
    if numerical_feats:
        feature_ids.append('scaler')
    features_dir: str = "_".join(feature_ids)
    if cutoff:
        cutoff_parameter = "-".join(feature_abbrev.get(key,key) for key in cutoff)
    f_root_dir = f"classification_target_{targets_dir}" if classification else  f"target_{targets_dir}"
    f_root_dir = f"OOD_{f_root_dir}" if clustering_method else f_root_dir
    # f_root_dir = f"{f_root_dir}_{target_transformer}FT" if target_transformer else f_root_dir
    f_root_dir = f"{f_root_dir}_filter_({cutoff_parameter})" if cutoff else f_root_dir
    f_root_dir = f"{f_root_dir}_{special_folder_name}" if special_folder_name else f_root_dir

    results_dir: Path = ROOT / output_dir_name / f_root_dir
    clustering_method= feature_abbrev.get(clustering_method, clustering_method) if clustering_method else None
    results_dir: Path = results_dir / clustering_method if clustering_method else results_dir
    results_dir: Path = results_dir / "test" if TEST else results_dir
    results_dir: Path = results_dir / features_dir
    # if subspace_filter:
    #     results_dir = results_dir / f"subspace_{subspace_filter}"

    _save(scores=scores,
          predictions=predictions,
          ground_truth=ground_truth,
          results_dir=results_dir,
          feat_length_scale=feat_length_scale,
          regressor_type=regressor_type,
          imputer=imputer,
          representation=representation,
          pu_type=pu_type,
          radius=radius,
          vector=vector,
          numerical_feats=numerical_feats,
          hypop=hypop,
          transform_type=transform_type,
          target_transformer=target_transformer,
          learning_curve=learning_curve,
          special_file_name=special_file_name,
          )
    return results_dir

