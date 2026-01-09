import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import pandas as pd
import json
from visualization_setting import save_img_path
from matplotlib.ticker import MaxNLocator
import os

HERE = Path(__file__).resolve().parent
DATASETS = HERE.parent.parent/'datasets'
JSONS =  DATASETS/'json_resources'
VISUALIZATION = HERE.parent/ "visualization"
training_df_dir: Path = DATASETS/ "training_dataset"/"Rg data with clusters aging imputed.pkl"
w_data = pd.read_pickle(training_df_dir)


def plot_data_reduction_steps()->None:

    with open(JSONS/'data_summary_monitor.json', 'r') as f:
        data = json.load(f)
    targets = data["Target"]
    initial_cleaned = np.array(data["Initial Cleaned"]).squeeze()
    block_copolymer = initial_cleaned - np.array(data["After dropping block copolymer"]).squeeze()
    solid_additives = np.array(data["After dropping block copolymer"]).squeeze() - np.array(data["After dropping solid additives"]).squeeze()
    polymer_hsps = np.array(data["After dropping solid additives"]).squeeze() - np.array(data["After dropping unknown polymer HSPs"]).squeeze()


    def clean_target_labels(targets):
        # Replace specific labels with cleaned versions
        replacements = {
            "Rh (IW avg log)": "Rh (nm)",
            "Rg1 (nm)": "Rg (nm)"
        }
        return [replacements.get(target, target) for target in targets]

    targets = clean_target_labels(targets)

    fig, ax = plt.subplots(figsize=(7, 6))
    bar_width = 0.5
    x = np.arange(len(targets)) 

    bars = ax.bar(x, initial_cleaned, bar_width, label='Final dataset', color='#66C2A5')
    ax.bar(x, block_copolymer, bar_width, bottom=initial_cleaned - block_copolymer, label='Dropped Block Copolymer', color='#EA6F00')
    ax.bar(x, solid_additives, bar_width, bottom=initial_cleaned - block_copolymer - solid_additives, label='Dropped Solid Additives', color='#FFAA5D')
    ax.bar(x, polymer_hsps, bar_width, bottom=initial_cleaned - block_copolymer - solid_additives - polymer_hsps, label='Dropped unknown Polymer HSPs', color='#FFC997')

    remaining_data_points = initial_cleaned - block_copolymer - solid_additives - polymer_hsps
    for i, bar in enumerate(bars):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height()/2, str(remaining_data_points[i]), 
                    ha='center', va='bottom', fontsize=15, color='black')



    ax.set_ylabel('Number of Data Points', fontsize=20, fontweight='bold')
    # ax.set_title('Data Reduction at Each Step of Cleaning', fontsize=21, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(targets, rotation=45, ha='right')
    ax.tick_params(axis='x', labelsize=21)
    ax.tick_params(axis='y', labelsize=18)

    ax.legend(fontsize=12)
    plt.tight_layout()
    save_img_path(HERE/'analysis and test','Data Reduction at Each Step of cleaning.png')
    plt.show()
    plt.close()


def plot_hanson_space(df, hsp_material:str):

    dP = df[f"{hsp_material} dP"]
    dD = df[f"{hsp_material} dD"]
    dH = df[f"{hsp_material} dH"]

    dP_freq = dP.map(dP.value_counts())
    # dD_freq = dD.map(dD.value_counts())
    # dH_freq = dH.map(dH.value_counts())

    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    maxlocator = {
        "integer":True,
        "prune":"lower",
        "nbins":6,
                  }
    
    sc = axes[0].scatter(dP, dD, c=dP_freq, cmap="viridis",edgecolors='none')
    axes[0].set_xlabel(r'$\delta P$ (MPa$^{1/2}$)')
    axes[0].set_ylabel(r'$\delta D$ (MPa$^{1/2}$)')
    axes[0].xaxis.set_major_locator(MaxNLocator(**maxlocator))
    axes[0].yaxis.set_major_locator(MaxNLocator(**maxlocator))


    # Scatter plot for dP vs dH
    sc = axes[1].scatter(dP, dH, c=dP_freq, cmap="viridis",edgecolors='none')
    axes[1].set_xlabel(r'$\delta P$ (MPa$^{1/2}$)')
    axes[1].set_ylabel(r'$\delta H$ (MPa$^{1/2}$)')
    axes[1].xaxis.set_major_locator(MaxNLocator(**maxlocator))
    axes[1].yaxis.set_major_locator(MaxNLocator(**maxlocator))


    # Scatter plot for dD vs dH
    sc = axes[2].scatter(dD, dH, c=dP_freq, cmap='viridis',edgecolors='none')
    axes[2].set_xlabel(r'$\delta D$ (MPa$^{1/2}$)')
    axes[2].set_ylabel(r'$\delta H$ (MPa$^{1/2}$)')
    axes[2].xaxis.set_major_locator(MaxNLocator(**maxlocator))
    axes[2].yaxis.set_major_locator(MaxNLocator(**maxlocator))
    cbar3 = fig.colorbar(sc, ax=axes[2], orientation="vertical", shrink=0.7, aspect=25)
    cbar3.set_label("Frequency")
    cbar3.set_ticks(np.linspace(cbar3.vmin, cbar3.vmax, 3))  # <-- set 3 ticks
    fig.suptitle(f"{hsp_material} Hansen Solubility Parameters", fontsize=14, fontweight="bold")

    plt.tight_layout()
    # plt.show()
    fname = f"Scattering of {hsp_material} Hansen Solubility Parameters"
    save_img_path(VISUALIZATION/"analysis and test",f"{fname}.png")
    plt.close()






# if __name__ == "__main__":
    # plot_data_reduction_steps()
    # plot_hanson_space(w_data,'solvent')