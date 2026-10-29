#!/bin/bash

# Submit the missing IID GpyroHMC GPU results for the selected Tanimoto
# fingerprint kernels, count kernels, and hybridization methods. The missing
# configuration IDs below were generated from
# Tree_and_GP_count_and_fingerprint_GPU.pkl using OOF_R2 as the completeness
# criterion on 2026-10-29.

DATE=$(date +%Y%m%d)
model="GpyroHMC"
output_root="/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}"
job_index=0

k_fps=("TanimotoRBF" "TanimotoMatern32" "TanimotoMatern52")
k_counts=("RBF" "Matern32" "Matern52")
k_mixing_methods=(
    "sum"
    "product"
    "(count:+)x(fp:x)"
    "(count:x)+(fp:x)"
)

submit_job() {
    local paper="$1"
    local dataset="$2"
    local target_feature="$3"
    local fp_kernel="$4"
    local count_kernel="$5"
    local mixing_method="$6"
    local output_dir="${output_root}/${paper}"
    local dataset_tag=${dataset// /_}
    local target_tag=${target_feature// /_}
    target_tag=${target_tag//\//_}

    job_index=$((job_index + 1))
    mkdir -p "$output_dir"

    bsub <<EOT
#BSUB -n 1
#BSUB -W 40:59
#BSUB -q gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -R "rusage[mem=32GB]"
#BSUB -R "select[a10 || a30 || a100 || l40 || h100]"
#BSUB -J "gpyrohmc_missing_${DATE}_${job_index}"
#BSUB -o "${output_dir}/${model}_${dataset_tag}_${target_tag}_${fp_kernel}_${count_kernel}_${mixing_method}_GPU.out"
#BSUB -e "${output_dir}/${model}_${dataset_tag}_${target_tag}_${fp_kernel}_${count_kernel}_${mixing_method}_GPU.err"

source ~/.bashrc
module load cuda/12.1
module load gcc/9.3.0
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12

python ../train_structure_numerical.py --K_fp "$fp_kernel" \
                                      --K_count "$count_kernel" \
                                      --Kernel_mixing_method "$mixing_method" \
                                      --paper "$paper" \
                                      --dataset "$dataset" \
                                      --target_features "$target_feature" \
                                      --regressor_type "$model"
EOT
}

# A configuration ID is its zero-based position in the Cartesian product:
# k_fps x k_counts x k_mixing_methods, with mixing method varying fastest.
# Each entry is paper|dataset|target|space-separated missing configuration IDs.
missing_target_configs=(
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed|log (Separation factor)|1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 20 21 22 23 25 26 27 28 30 31 32 33 34 35"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|flux_data_imputed|log (Total flux)|1 2 3 5 6 7 8 9 10 11 14 15 20 21 22 23 26 27 29 30 31 32 33 34 35"
    "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation|cleaned_dataset_pervaporation_membranes_wang|log (Separation factor)|0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 20 21 22 23 24 26 27 28 29 30 31 32 33 34 35"
    "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation|cleaned_dataset_pervaporation_membranes_wang|log (Total flux)|0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 20 21 22 23 24 26 27 28 29 30 31 32 33 34 35"
    "Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation|Rg data with clusters aging imputed|log Rg (nm)|0 8 20 28"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|flux decline ratio (%)|2 14 22 26 30 34"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|flux recovery ratio (%)|2 14 22 26 30 34"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|irreversible fouling ratio(%)|2 14 22 26 30 34"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|organic compound removal (%)|2 14 22 26 30 34"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|reversible fouling ratio (%)|2 14 22 26 30 34"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|water permeability (LMH/bar)|14 22 26 30 34"
)

mixing_count=${#k_mixing_methods[@]}
count_kernel_count=${#k_counts[@]}
configs_per_fp=$((count_kernel_count * mixing_count))

for target_config in "${missing_target_configs[@]}"; do
    IFS='|' read -r paper dataset target_feature missing_config_ids <<< "$target_config"

    for config_id in $missing_config_ids; do
        fp_index=$((config_id / configs_per_fp))
        config_remainder=$((config_id % configs_per_fp))
        count_index=$((config_remainder / mixing_count))
        mixing_index=$((config_remainder % mixing_count))

        submit_job \
            "$paper" \
            "$dataset" \
            "$target_feature" \
            "${k_fps[$fp_index]}" \
            "${k_counts[$count_index]}" \
            "${k_mixing_methods[$mixing_index]}"
    done
done

echo "Submitted ${job_index} missing GpyroHMC GPU jobs."
