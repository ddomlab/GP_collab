#!/bin/bash

# Submit only GpyroHMC GPU runs missing from
# Tree_and_GP_count_and_fingerprint_GPU.pkl. Four result rows remain missing:
# Matern32 / Matern32 with the count:+ × fp:+ hybridization. They require
# three jobs because the Wang pervaporation run produces both targets.

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
#BSUB -W 71:59
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
    # Polymer Design: one target per dataset.
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|flux_data_imputed|Matern32|Matern32|(count:+)x(fp:+)"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed|Matern32|Matern32|(count:+)x(fp:+)"
    # This run produces both Wang pervaporation targets.
    "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation|cleaned_dataset_pervaporation_membranes_wang|Matern32|Matern32|(count:+)x(fp:+)"
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
