#!/bin/bash

# Submit the missing HDF tree-model result identified in the master-data check.
# train_structure_numerical.py trains every target configured for each dataset
# in filter_data.py.

DATE=$(date +%Y%m%d)
models=("XGBR")
hdf_depths=(2)
hdf_bits=(512)
# This dataset produces all six missing ultrafiltration target results.
selected_training_sets=(
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed"
)

output_root="/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}"
job_index=0

for training_set in "${selected_training_sets[@]}"; do
    IFS='|' read -r paper dataset <<< "$training_set"
    output_dir="${output_root}/${paper}"
    dataset_tag=${dataset// /_}
    mkdir -p "$output_dir"

    for depth in "${hdf_depths[@]}"; do
        for n_bits in "${hdf_bits[@]}"; do
            for model in "${models[@]}"; do
                job_index=$((job_index + 1))
                bsub <<EOT

#BSUB -n 6
#BSUB -W 10:59
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=32GB]"
#BSUB -J "tree_HDF_d${depth}_b${n_bits}_${DATE}_${job_index}"
#BSUB -o "${output_dir}/${model}_HDF_depth${depth}_${n_bits}bit_${dataset_tag}_CPU.out"
#BSUB -e "${output_dir}/${model}_HDF_depth${depth}_${n_bits}bit_${dataset_tag}_CPU.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12

python ../train_structure_numerical.py --regressor_type "$model" \
                                        --representation "HDF" \
                                        --radius "$depth" \
                                        --n_bits "$n_bits" \
                                        --dataset "$dataset" \
                                        --paper "$paper"

EOT
            done
        done
    done
done
