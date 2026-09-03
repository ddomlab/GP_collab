#!/bin/bash

# Submit only the GPytorchMAP GPU configurations that are currently missing
# from the multi-target pervaporation result sets.

DATE=$(date +%Y%m%d)
model="GPytorchMAP"
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
#BSUB -W 1:59
#BSUB -q short_gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -R "rusage[mem=8GB]"
#BSUB -R "select[a10 || a30 || a100 || l40 || h100]"
#BSUB -J "gpmap_missing_${DATE}_${job_index}"
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

# This dataset trains both "log (Total flux)" and
# "log (Separation factor)", so each call fills two target-level gaps.
paper="Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation"
dataset="cleaned_dataset_pervaporation_membranes_wang"

submit_job "$paper" "$dataset" "Matern32"          "RBF"      "(count:x)+(fp:x)"
submit_job "$paper" "$dataset" "TanimotoMatern52" "Matern32" "(count:+)x(fp:x)"
submit_job "$paper" "$dataset" "TanimotoMatern52" "Matern32" "(count:x)+(fp:x)"
submit_job "$paper" "$dataset" "TanimotoRBF"      "Matern32" "(count:x)+(fp:x)"

# flux_data_imputed contains only the "log (Total flux)" target.
paper="Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery"
dataset="flux_data_imputed"

submit_job "$paper" "$dataset" "Tanimoto"    "Matern52" "(count:+)x(fp:x)"
submit_job "$paper" "$dataset" "TanimotoRBF" "Matern32" "(count:x)+(fp:x)"
