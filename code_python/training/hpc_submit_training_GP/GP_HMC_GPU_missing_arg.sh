#!/bin/bash

# Submit the single GpyroHMC OOD GPU result missing from the continuous-cluster
# master data: the Wang pervaporation separation-factor target.

DATE=$(date +%Y%m%d)
model="GpyroHMC"
output_root="/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}"
job_index=0

submit_job() {
    local paper="$1"
    local dataset="$2"
    local target_feature="$3"
    local fp_kernel="$4"
    local count_kernel="$5"
    local mixing_method="$6"
    local clustering_method="$7"
    local output_dir="${output_root}/${paper}"
    local dataset_tag=${dataset// /_}
    local target_tag=${target_feature// /_}
    target_tag=${target_tag//\//_}

    job_index=$((job_index + 1))
    mkdir -p "$output_dir"

    bsub <<EOT
#BSUB -n 1
#BSUB -W 55:59
#BSUB -q gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -R "rusage[mem=32GB]"
#BSUB -R "select[a10 || a30 || a100 || l40 || h100]"
#BSUB -J "gpyrohmc_missing_${DATE}_${job_index}"
#BSUB -o "${output_dir}/${model}_${dataset_tag}_${target_tag}_${fp_kernel}_${count_kernel}_${mixing_method}_${clustering_method}_GPU.out"
#BSUB -e "${output_dir}/${model}_${dataset_tag}_${target_tag}_${fp_kernel}_${count_kernel}_${mixing_method}_${clustering_method}_GPU.err"

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
                                      --regressor_type "$model" \
                                      --clustering_method "$clustering_method"
EOT
}

missing_jobs=(
    "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation|cleaned_dataset_pervaporation_membranes_wang|log (Separation factor)|Matern32|Matern32|(count:x)+(fp:x)|continuous_cluster"
)

for job in "${missing_jobs[@]}"; do
    IFS='|' read -r paper dataset target_feature fp_kernel count_kernel mixing_method clustering_method <<< "$job"
    submit_job \
        "$paper" \
        "$dataset" \
        "$target_feature" \
        "$fp_kernel" \
        "$count_kernel" \
        "$mixing_method" \
        "$clustering_method"
done

echo "Submitted ${job_index} GpyroHMC GPU jobs."
