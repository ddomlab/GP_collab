#!/bin/bash

# Submit only GpyroHMC GPU runs missing from
# Tree_and_GP_GPU_count_and_fingerprint.pkl for these FP/count pairs:
#   - TanimotoMatern32 / Matern32
#   - Matern32 / Matern32
#
# The 24 missing dataset-target results below require 14 jobs: the Wang
# pervaporation and ultrafiltration datasets train multiple targets per run.

DATE=$(date +%Y%m%d)
model="GpyroHMC"
output_root="/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}"
job_index=0

submit_job() {
    local paper="$1"
    local dataset="$2"
    local fp_kernel="$3"
    local count_kernel="$4"
    local mixing_method="$5"
    local output_dir="${output_root}/${paper}"

    job_index=$((job_index + 1))
    mkdir -p "$output_dir"

    bsub <<EOT
#BSUB -n 1
#BSUB -W 65:10
#BSUB -q gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -R "rusage[mem=32GB]"
#BSUB -R "select[a10 || a30 || a100 || l40 || h100]"
#BSUB -J "gpyrohmc_missing_${DATE}_${job_index}"
#BSUB -o "${output_dir}/${model}_${dataset}_${fp_kernel}_${count_kernel}_${mixing_method}_GPU.out"
#BSUB -e "${output_dir}/${model}_${dataset}_${fp_kernel}_${count_kernel}_${mixing_method}_GPU.err"

source ~/.bashrc
module load cuda/12.1
module load gcc/9.3.0
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12

python ../train_structure_numerical.py --K_fp "$fp_kernel" \
                                      --K_count "$count_kernel" \
                                      --Kernel_mixing_method "$mixing_method" \
                                      --paper "$paper" \
                                      --dataset "$dataset" \
                                      --regressor_type "$model"
EOT
}

missing_jobs=(
    # TanimotoMatern32 / Matern32
    "Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation|Rg data with clusters aging imputed|TanimotoMatern32|Matern32|sum"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|flux_data_imputed|TanimotoMatern32|Matern32|sum"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed|TanimotoMatern32|Matern32|sum"
    # Both Wang pervaporation targets are trained by each invocation.
    "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation|cleaned_dataset_pervaporation_membranes_wang|TanimotoMatern32|Matern32|sum"
    "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation|cleaned_dataset_pervaporation_membranes_wang|TanimotoMatern32|Matern32|(count:+)x(fp:x)"

    # Matern32 / Matern32
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|flux_data_imputed|Matern32|Matern32|sum"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed|Matern32|Matern32|sum"
    "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation|cleaned_dataset_pervaporation_membranes_wang|Matern32|Matern32|product"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|flux_data_imputed|Matern32|Matern32|(count:+)x(fp:x)"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed|Matern32|Matern32|(count:+)x(fp:x)"
    "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation|cleaned_dataset_pervaporation_membranes_wang|Matern32|Matern32|(count:+)x(fp:x)"
    # This invocation trains all six ultrafiltration targets.
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|Matern32|Matern32|(count:+)x(fp:x)"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed|Matern32|Matern32|(count:x)+(fp:x)"
    "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation|cleaned_dataset_pervaporation_membranes_wang|Matern32|Matern32|(count:x)+(fp:x)"
)

for job in "${missing_jobs[@]}"; do
    IFS='|' read -r paper dataset fp_kernel count_kernel mixing_method <<< "$job"
    submit_job \
        "$paper" \
        "$dataset" \
        "$fp_kernel" \
        "$count_kernel" \
        "$mixing_method"
done

echo "Submitted ${job_index} GpyroHMC GPU jobs."
