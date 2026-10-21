#!/bin/bash

# Submit the four missing unique-fingerprint-permutation GPytorchMAP runs for
# the irreversible fouling ratio target only.

DATE=$(date +%Y%m%d)
model="GPytorchMAP"
fp_kernel="TanimotoMatern32"
count_kernel="Matern32"
output_root="/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}"
job_index=0

submit_job() {
    local paper="$1"
    local dataset="$2"
    local fp_kernel="$3"
    local count_kernel="$4"
    local target_feature="$5"
    local mixing_method="$6"
    local output_dir="${output_root}/${paper}"
    local target_tag=${target_feature// /_}

    job_index=$((job_index + 1))
    mkdir -p "$output_dir"

    bsub <<EOT
#BSUB -n 1
#BSUB -W 30:20
#BSUB -q gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -R "rusage[mem=32GB]"
#BSUB -R "select[a10 || a30 || a100 || l40 || h100]"
#BSUB -J "gpmap_missing_${DATE}_${job_index}"
#BSUB -o "${output_dir}/${model}_${dataset}_${target_tag}_${fp_kernel}_${count_kernel}_${mixing_method}_GPU.out"
#BSUB -e "${output_dir}/${model}_${dataset}_${target_tag}_${fp_kernel}_${count_kernel}_${mixing_method}_GPU.err"

source ~/.bashrc
module load cuda/12.1
module load gcc/9.3.0
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12

python ../train_structure_numerical_v2.py --K_fp "$fp_kernel" \
                                        --K_count "$count_kernel" \
                                        --Kernel_mixing_method "$mixing_method" \
                                        --paper "$paper" \
                                        --dataset "$dataset" \
                                        --target_features "$target_feature" \
                                        --regressor_type "$model"
EOT
}

missing_jobs=(
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|irreversible fouling ratio(%)|sum"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|irreversible fouling ratio(%)|product"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|irreversible fouling ratio(%)|(count:x)+(fp:x)"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|irreversible fouling ratio(%)|(count:+)x(fp:+)"
)

for job in "${missing_jobs[@]}"; do
    IFS='|' read -r paper dataset target_feature mixing_method <<< "$job"
    submit_job \
        "$paper" \
        "$dataset" \
        "$fp_kernel" \
        "$count_kernel" \
        "$target_feature" \
        "$mixing_method"
done

echo "Submitted ${job_index} GPytorchMAP GPU jobs."
