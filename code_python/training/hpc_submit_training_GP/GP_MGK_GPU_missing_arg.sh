#!/bin/bash

# Submit only the missing MGK GPU configurations for the multi-target datasets.

DATE=$(date +%Y%m%d)
model="MGK"
fp_kernel="Graph"
count_kernel="Matern32"
feature_mode="per_feature"
output_root="/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}"
job_index=0

submit_job() {
    local paper="$1"
    local dataset="$2"
    local mixing_method="$3"
    local output_dir="${output_root}/${paper}"

    job_index=$((job_index + 1))
    mkdir -p "$output_dir"

    bsub <<EOT
#BSUB -n 2
#BSUB -W 20:55
#BSUB -q gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -R "rusage[mem=32GB]"
#BSUB -R "select[a10 || a30 || a100 || l40 || h100]"
#BSUB -J "mgk_missing_${DATE}_${job_index}"
#BSUB -o "${output_dir}/${model}_${dataset}_${mixing_method}_${feature_mode}_GPU.out"
#BSUB -e "${output_dir}/${model}_${dataset}_${mixing_method}_${feature_mode}_GPU.err"

source ~/.bashrc
module load cuda/12.1
module load gcc/9.3.0
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12

python ../train_structure_numerical.py --K_fp "$fp_kernel" \
                                      --K_count "$count_kernel" \
                                      --Kernel_mixing_method "$mixing_method" \
                                      --paper "$paper" \
                                      --dataset "$dataset" \
                                      --kernel_feature_mode "$feature_mode" \
                                      --regressor_type "$model"
EOT
}

paper="Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery"
for dataset in "separation_data_imputed" "flux_data_imputed"; do
    submit_job "$paper" "$dataset" "sum"
    submit_job "$paper" "$dataset" "(count:+)x(graph:x)"
    submit_job "$paper" "$dataset" "(count:x)+(graph:x)"
done

paper="Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation"
dataset="cleaned_dataset_pervaporation_membranes_wang"
submit_job "$paper" "$dataset" "(count:+)x(graph:x)"

paper="Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning"
dataset="cleaned_dataset_Ultrafiltration Membrane_imputed"
submit_job "$paper" "$dataset" "sum"

echo "Submitted ${job_index} MGK GPU jobs."
