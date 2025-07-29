import json
from itertools import product
from pathlib import Path
from typing import List, Optional, Dict, Tuple
import os 
# import cmcrameri.cm as cmc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from math import ceil
from visualization_setting import set_plot_style, save_img_path




def ensure_long_path(path):
        """Ensures Windows handles long paths by adding '\\?\' if needed."""
        path_str = str(path)
        if os.name == 'nt' and len(path_str) > 250:  
            return Path(f"\\\\?\\{path_str}")
        return path

HERE: Path = Path(__file__).resolve().parent
results_path = HERE.parent.parent / 'results'/ 'target_log Rg (nm)'


scores_folder_path = results_path /'scaler'
score_file_lc = scores_folder_path / f'(Xn-Mw-PDI-concentration-temperature-polymer dP-polymer dD-polymer dH-solvent dP-solvent dD-solvent dH)_sklearn-GPR.rbf_hypOFF_Standard_length_scale.json'

ls = ensure_long_path(score_file_lc)


with open(ls, 'r') as f:
    ls_data = json.load(f)



flattened = []

for seed, folds in ls_data.items():
    for fold, features in folds.items():
        for feature, value in features.items():
            flattened.append({'feature': feature, 'length_scale': value})

df = pd.DataFrame(flattened)
result = df.groupby('feature')['length_scale'].agg(['mean', 'std']).sort_values('mean', ascending=False)

print(result)