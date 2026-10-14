#!/bin/bash

# Submit the missing regular tree-model runs found in
# Tree_and_GP_structure_and_count_GPU.pkl.

DATE=$(date +%Y%m%d)
missing_jobs=(
    "RF|LE|Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices|Beyond molecular structure_seifrid_imputed"
    "RF|Mordred|Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices|Beyond molecular structure_seifrid_imputed"
    "XGBR|LE|Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices|Beyond molecular structure_seifrid_imputed"
    "XGBR|Mordred|Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices|Beyond molecular structure_seifrid_imputed"
    "NGB|LE|Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices|Beyond molecular structure_seifrid_imputed"
    "NGB|Mordred|Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices|Beyond molecular structure_seifrid_imputed"
    "RF|Mordred|Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation|Rg data with clusters aging imputed"
    "XGBR|Mordred|Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation|Rg data with clusters aging imputed"
    "NGB|Mordred|Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation|Rg data with clusters aging imputed"
)

output_root="/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}"
job_index=0

for job in "${missing_jobs[@]}"; do
    IFS='|' read -r model representation paper dataset <<< "$job"
    output_dir="${output_root}/${paper}"
    dataset_tag=${dataset// /_}
    mkdir -p "$output_dir"

    job_index=$((job_index + 1))
    bsub <<EOT

#BSUB -n 6
#BSUB -W 30:00
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=32GB]"
#BSUB -J "tree_missing_${representation}_${DATE}_${job_index}"
#BSUB -o "${output_dir}/${model}_${representation}_${dataset_tag}_CPU.out"
#BSUB -e "${output_dir}/${model}_${representation}_${dataset_tag}_CPU.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12

python ../train_structure_numerical.py --regressor_type "$model" \
                                        --representation "$representation" \
                                        --dataset "$dataset" \
                                        --paper "$paper"

EOT
done

echo "Submitted ${job_index} missing tree-model jobs."
