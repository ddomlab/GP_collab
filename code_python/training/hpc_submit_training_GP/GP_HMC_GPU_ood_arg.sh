#!/bin/bash
set -euo pipefail

# Rerun GpyroHMC for the 13 structure-cluster OOD tasks whose 29 cluster
# scores for RUSC are NaN in the GPU master CSV. Each job handles one target.

DATE=$(date +%Y%m%d)
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
model="GpyroHMC"
fp_kernel="Matern32"
count_kernel="Matern32"
mixing_method="(count:+)x(fp:x)"
clustering_method="structure_cluster"
output_root="/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}"
job_count=0

# One entry per paper|dataset|target.
dataset_targets=(
    "Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation|Rg data with clusters aging imputed|log Rg (nm)"
    "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices|Beyond molecular structure_seifrid_imputed|calculated PCE (%)"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed|log (Separation factor)"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|flux_data_imputed|log (Total flux)"
    "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation|cleaned_dataset_pervaporation_membranes_wang|log (Separation factor)"
    "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation|cleaned_dataset_pervaporation_membranes_wang|log (Total flux)"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|flux decline ratio (%)"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|flux recovery ratio (%)"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|irreversible fouling ratio(%)"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|organic compound removal (%)"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|reversible fouling ratio (%)"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|water permeability (LMH/bar)"
    "Miniaturization of Popular Reactions from the Medicinal Chemists Toolbox for Ultrahigh_Throughput Experimentation|cleaned_suzuki_synthesis|Approx Conv (%)"
)

submit_job() {
    local paper="$1"
    local dataset="$2"
    local target_feature="$3"
    local output_dir="${output_root}/${paper}"
    local dataset_tag=${dataset// /_}
    local target_tag=${target_feature// /_}
    target_tag=${target_tag//\//_}

    job_count=$((job_count + 1))
    mkdir -p "$output_dir"

    bsub <<EOT
#BSUB -n 1
#BSUB -W 1:59
#BSUB -q short_gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -R "rusage[mem=32GB]"
#BSUB -R "select[a10 || a30 || a100 || l40 || h100]"
#BSUB -J "gpyrohmc_structure_rusc_${DATE}_${job_count}"
#BSUB -o "${output_dir}/${model}_${dataset_tag}_${target_tag}_${fp_kernel}_${count_kernel}_${mixing_method}_${clustering_method}_GPU.out"
#BSUB -e "${output_dir}/${model}_${dataset_tag}_${target_tag}_${fp_kernel}_${count_kernel}_${mixing_method}_${clustering_method}_GPU.err"

source ~/.bashrc
module load cuda/12.1
module load gcc/9.3.0
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12

python "${script_dir}/../train_structure_numerical.py" \
    --K_fp "$fp_kernel" \
    --K_count "$count_kernel" \
    --Kernel_mixing_method "$mixing_method" \
    --paper "$paper" \
    --dataset "$dataset" \
    --target_features "$target_feature" \
    --regressor_type "$model" \
    --clustering_method "$clustering_method"
EOT
}

for dataset_target in "${dataset_targets[@]}"; do
    IFS='|' read -r paper dataset target_feature <<< "$dataset_target"
    submit_job "$paper" "$dataset" "$target_feature"
done

echo "Submitted ${job_count} GpyroHMC structure-cluster OOD GPU jobs."
