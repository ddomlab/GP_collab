import pandas as pd
from pathlib import Path
from training_utils import train_regressor
from all_factories import radius_to_bits,cutoffs
import json
import numpy as np
import sys
sys.path.append("code_/cleaning")
from clean_dataset import open_json
from argparse import ArgumentParser
from data_handling import save_results

HERE: Path = Path(__file__).resolve().parent
DATASETS: Path = HERE.parent.parent / "datasets"
RESULTS = Path = HERE.parent.parent / "results"

training_df_dir: Path = DATASETS/ "training_dataset"/ "dataset_wo_block_cp_(fp-hsp)_added_additive_dropped.pkl"

w_data = pd.read_pickle(training_df_dir)

print(w_data[['solvent dH','solvent dD','solvent dP']].isna().sum())