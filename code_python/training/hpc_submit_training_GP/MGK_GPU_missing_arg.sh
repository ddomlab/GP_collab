#!/bin/bash

# Submit only the MGK GPU runs whose results are missing from
# Tree_and_GP_count_and_fingerprint_GPU.pkl. There are 12 missing result rows,
# represented by 10 jobs because the Wang pervaporation dataset trains both
# targets in each invocation.

DATE=$(date +%Y%m%d)
model="MGK"
fp_kernel="Graph"
feature_mode="per_feature"
output_root="/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}"
job_index=0

submit_job() {
    local paper="$1"
    local dataset="$2"
    local count_kernel="$3"
    local mixing_method="$4"
    local output_dir="${output_root}/${paper}"

    job_index=$((job_index + 1))
    mkdir -p "$output_dir"

    bsub <<EOT
#BSUB -n 2
#BSUB -W 55:55
#BSUB -q gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -R "rusage[mem=32GB]"
#BSUB -R "select[a10 || a30 || a100 || l40 || h100]"
#BSUB -J "mgk_missing_${DATE}_${job_index}"
#BSUB -o "${output_dir}/${model}_${dataset}_${count_kernel}_${mixing_method}_${feature_mode}_GPU.out"
#BSUB -e "${output_dir}/${model}_${dataset}_${count_kernel}_${mixing_method}_${feature_mode}_GPU.err"

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

# calculated PCE (%) -- 558 datapoints
paper="Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices"
dataset="Beyond molecular structure_seifrid_imputed"

submit_job "$paper" "$dataset" "RBF"      "sum"
submit_job "$paper" "$dataset" "Matern52" "product"
submit_job "$paper" "$dataset" "Matern52" "sum"

# log (Total flux)
paper="Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery"
dataset="flux_data_imputed"

submit_job "$paper" "$dataset" "Matern52" "product"
submit_job "$paper" "$dataset" "Matern52" "sum"

# log (Separation factor)
paper="Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery"
dataset="separation_data_imputed"

submit_job "$paper" "$dataset" "RBF"      "sum"

# Both pervaporation targets are trained by each invocation.
paper="Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation"
dataset="cleaned_dataset_pervaporation_membranes_wang"

submit_job "$paper" "$dataset" "RBF"      "(count:x)+(graph:x)"
submit_job "$paper" "$dataset" "Matern52" "(count:x)+(graph:x)"

# Approx Conv (%) -- 768 datapoints
paper="Miniaturization of Popular Reactions from the Medicinal Chemists Toolbox for Ultrahigh_Throughput Experimentation"
dataset="cleaned_suzuki_synthesis"

submit_job "$paper" "$dataset" "RBF"      "(count:x)+(graph:x)"
submit_job "$paper" "$dataset" "Matern52" "(count:+)x(graph:x)"

echo "Submitted ${job_index} MGK GPU jobs."
