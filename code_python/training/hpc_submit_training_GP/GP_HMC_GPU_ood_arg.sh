#!/bin/bash

# Train GpyroHMC OOD models for every dataset group and target. The training
# entrypoint iterates all targets in each dataset; the seven groups below cover
# all 13 configured targets. This submits 56 jobs:
# 7 dataset groups x 4 mixing methods x 2 clustering methods.

DATE=$(date +%Y%m%d)
model="GpyroHMC"
k_fps=("TanimotoMatern32")
k_counts=("Matern32")
k_mixing_methods=("sum" "product" "(count:+)x(fp:+)" "(count:x)+(fp:x)")
clustering_methods=("continuous_cluster" "structure_cluster")
output_root="/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}"
job_index=0

dataset_groups=(
    "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices|Beyond molecular structure_seifrid_imputed"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|flux_data_imputed"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed"
    "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation|cleaned_dataset_pervaporation_membranes_wang"
    "Miniaturization of Popular Reactions from the Medicinal Chemists Toolbox for Ultrahigh_Throughput Experimentation|cleaned_suzuki_synthesis"
    "Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation|Rg data with clusters aging imputed"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed"
)

submit_job() {
    local paper="$1"
    local dataset="$2"
    local fp_kernel="$3"
    local count_kernel="$4"
    local mixing_method="$5"
    local clustering_method="$6"
    local output_dir="${output_root}/${paper}"
    local dataset_tag=${dataset// /_}

    job_index=$((job_index + 1))
    mkdir -p "$output_dir"

    bsub <<EOT
#BSUB -n 1
#BSUB -W 55:10
#BSUB -q gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -R "rusage[mem=32GB]"
#BSUB -R "select[a10 || a30 || a100 || l40 || h100]"
#BSUB -J "gpyrohmc_ood_${DATE}_${job_index}"
#BSUB -o "${output_dir}/${model}_${dataset_tag}_${fp_kernel}_${count_kernel}_${mixing_method}_${clustering_method}_GPU.out"
#BSUB -e "${output_dir}/${model}_${dataset_tag}_${fp_kernel}_${count_kernel}_${mixing_method}_${clustering_method}_GPU.err"

source ~/.bashrc
module load cuda/12.1
module load gcc/9.3.0
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12

python ../train_structure_numerical.py --K_fp "$fp_kernel" \
                                      --K_count "$count_kernel" \
                                      --Kernel_mixing_method "$mixing_method" \
                                      --paper "$paper" \
                                      --dataset "$dataset" \
                                      --regressor_type "$model" \
                                      --clustering_method "$clustering_method"
EOT
}

for dataset_group in "${dataset_groups[@]}"; do
    IFS='|' read -r paper dataset <<< "$dataset_group"
    for fp_kernel in "${k_fps[@]}"; do
        for count_kernel in "${k_counts[@]}"; do
            for mixing_method in "${k_mixing_methods[@]}"; do
                for clustering_method in "${clustering_methods[@]}"; do
                    submit_job \
                        "$paper" \
                        "$dataset" \
                        "$fp_kernel" \
                        "$count_kernel" \
                        "$mixing_method" \
                        "$clustering_method"
                done
            done
        done
    done
done

echo "Submitted ${job_index} GpyroHMC OOD GPU jobs."
