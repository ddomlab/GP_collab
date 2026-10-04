#!/bin/bash

# Submit the missing tree-model structure- and continuous-cluster OOD runs.
# Each entry explicitly specifies its model, dataset, and clustering method.

DATE=$(date +%Y%m%d)
missing_jobs=(
    # structure_cluster: log (Separation factor)
    "NGB|Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed|structure_cluster"

    # continuous_cluster: log Rg (nm)
    "RF|Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation|Rg data with clusters aging imputed|continuous_cluster"
    "XGBR|Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation|Rg data with clusters aging imputed|continuous_cluster"
    "NGB|Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation|Rg data with clusters aging imputed|continuous_cluster"
)

output_root="/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}"
job_index=0

for job in "${missing_jobs[@]}"; do
    IFS='|' read -r model paper dataset clustering_method <<< "$job"
    output_dir="${output_root}/${paper}"
    dataset_tag=${dataset// /_}
    mkdir -p "$output_dir"

    job_index=$((job_index + 1))
    bsub <<EOT

#BSUB -n 6
#BSUB -W 5:00
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=32GB]"
#BSUB -J "tree_ood_missing_${DATE}_${job_index}"
#BSUB -o "${output_dir}/${model}_${dataset_tag}_${clustering_method}_CPU.out"
#BSUB -e "${output_dir}/${model}_${dataset_tag}_${clustering_method}_CPU.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12

python ../train_structure_numerical.py --regressor_type "$model" \
                                        --dataset "$dataset" \
                                        --paper "$paper" \
                                        --clustering_method "$clustering_method"

EOT
done

echo "Submitted ${job_index} tree-model OOD jobs."
