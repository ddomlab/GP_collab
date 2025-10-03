import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

HERE: Path = Path(__file__).resolve().parent
VISUALIZATION = HERE.parent/ "visualization"
DATASETS: Path = HERE.parent.parent / "datasets"
RESULTS = Path = HERE.parent.parent / "results"
training_df_dir: Path = DATASETS/ "training_dataset"/ "dataset_wo_block_cp_(fp-hsp)_added_additive_dropped_polyHSP_dropped_peaks_appended_multimodal (40-1000 nm)_added.pkl"
w_data = pd.read_pickle(training_df_dir)


columns_to_keep =  [ 
    "canonical_name", 
    "Mn (g/mol)", "PDI", "Mw (g/mol)", "Concentration (mg/ml)", "Temperature SANS/SLS/DLS/SEC (K)",
    "polymer dP", "polymer dD", "polymer dH", "solvent dP", "solvent dD", "solvent dH", "Ra",
    "Monomer SMILES",	"Dimer SMILES",	"Trimer SMILES", "RRU Monomer SMILES",
    "RRU Dimer SMILES",	"RRU Trimer SMILES",
    "Monomer_Mordred", "Monomer_MACCS", "Monomer_ECFP6_count_512bits", "Monomer_ECFP8_count_1024bits",
    "Monomer_ECFP10_count_2048bits", "Monomer_ECFP12_count_4096bits", "Monomer_ECFP6_binary_512bits",
    "Monomer_ECFP8_binary_1024bits", "Monomer_ECFP10_binary_2048bits", "Monomer_ECFP12_binary_4096bits",
    "Dimer_Mordred", "Dimer_MACCS", "Dimer_ECFP6_count_512bits", "Dimer_ECFP8_count_1024bits",
    "Dimer_ECFP10_count_2048bits", "Dimer_ECFP12_count_4096bits", "Dimer_ECFP6_binary_512bits",
    "Dimer_ECFP8_binary_1024bits", "Dimer_ECFP10_binary_2048bits", "Dimer_ECFP12_binary_4096bits",
    "Trimer_Mordred", "Trimer_MACCS", "Trimer_ECFP6_count_512bits", "Trimer_ECFP8_count_1024bits",
    "Trimer_ECFP10_count_2048bits", "Trimer_ECFP12_count_4096bits", "Trimer_ECFP6_binary_512bits",
    "Trimer_ECFP8_binary_1024bits", "Trimer_ECFP10_binary_2048bits", "Trimer_ECFP12_binary_4096bits",
    "RRU Monomer_Mordred", "RRU Monomer_MACCS", "RRU Monomer_ECFP6_count_512bits", "RRU Monomer_ECFP8_count_1024bits",
    "RRU Monomer_ECFP10_count_2048bits", "RRU Monomer_ECFP12_count_4096bits", "RRU Monomer_ECFP6_binary_512bits",
    "RRU Monomer_ECFP8_binary_1024bits", "RRU Monomer_ECFP10_binary_2048bits", "RRU Monomer_ECFP12_binary_4096bits",
    "RRU Dimer_Mordred", "RRU Dimer_MACCS", "RRU Dimer_ECFP6_count_512bits", "RRU Dimer_ECFP8_count_1024bits",
    "RRU Dimer_ECFP10_count_2048bits", "RRU Dimer_ECFP12_count_4096bits", "RRU Dimer_ECFP6_binary_512bits",
    "RRU Dimer_ECFP8_binary_1024bits", "RRU Dimer_ECFP10_binary_2048bits", "RRU Dimer_ECFP12_binary_4096bits",
    "RRU Trimer_Mordred", "RRU Trimer_MACCS", "RRU Trimer_ECFP6_count_512bits", "RRU Trimer_ECFP8_count_1024bits",
    "RRU Trimer_ECFP10_count_2048bits", "RRU Trimer_ECFP12_count_4096bits", "RRU Trimer_ECFP6_binary_512bits",
    "RRU Trimer_ECFP8_binary_1024bits", "RRU Trimer_ECFP10_binary_2048bits", "RRU Trimer_ECFP12_binary_4096bits",
    "Rg1 (nm)", "multimodal Rh", 
]

out_put = w_data[columns_to_keep]
out_put = out_put[out_put["Rg1 (nm)"].notna() | out_put["multimodal Rh"].notna()]

out_put.to_pickle(DATASETS/"training_dataset"/"cleaned_dataset_gp_colab.pkl")
