#!/bin/bash

# Submit the currently missing MGK Graph/Matern32 OOD results for these
# hybridizations: sum, product, (count:+)x(graph:x), and (count:x)+(graph:x).
#
# Continuous-cluster OOD has no MGK results for these configurations (28 jobs).
# Structure-cluster OOD needs 8 more jobs; its additive configuration is
# already complete and is intentionally omitted. Total: 36 jobs.

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
    local clustering_method="$4"
    local output_dir="${output_root}/${paper}"
    local dataset_tag=${dataset// /_}

    job_index=$((job_index + 1))
    mkdir -p "$output_dir"

    bsub <<EOT
#BSUB -n 1
#BSUB -W 72:50
#BSUB -q gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -R "rusage[mem=32GB]"
#BSUB -R "select[a10 || a30 || a100 || l40 || h100]"
#BSUB -J "mgk_ood_missing_${DATE}_${job_index}"
#BSUB -o "${output_dir}/${model}_${dataset_tag}_${count_kernel}_${mixing_method}_${clustering_method}_GPU.out"
#BSUB -e "${output_dir}/${model}_${dataset_tag}_${count_kernel}_${mixing_method}_${clustering_method}_GPU.err"

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
                                      --regressor_type "$model" \
                                      --clustering_method "$clustering_method"
EOT
}

continuous_jobs=(
    "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices|Beyond molecular structure_seifrid_imputed"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|flux_data_imputed"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed"
    "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation|cleaned_dataset_pervaporation_membranes_wang"
    "Miniaturization of Popular Reactions from the Medicinal Chemists Toolbox for Ultrahigh_Throughput Experimentation|cleaned_suzuki_synthesis"
    "Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation|Rg data with clusters aging imputed"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed"
)
continuous_mixing_methods=(
    "sum"
    "product"
    "(count:+)x(graph:x)"
    "(count:x)+(graph:x)"
)

for job in "${continuous_jobs[@]}"; do
    IFS='|' read -r paper dataset <<< "$job"
    for mixing_method in "${continuous_mixing_methods[@]}"; do
        submit_job "$paper" "$dataset" "$mixing_method" "continuous_cluster"
    done
done

# Structure-cluster gaps only. Shared Wang and ultrafiltration invocations
# retrain their other targets as part of producing the missing target(s).
structure_jobs=(
    # sum: calculated PCE, separation factor, Wang separation factor,
    # five ultrafiltration targets, and Approx Conv.
    "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices|Beyond molecular structure_seifrid_imputed|sum"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed|sum"
    "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation|cleaned_dataset_pervaporation_membranes_wang|sum"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|sum"
    "Miniaturization of Popular Reactions from the Medicinal Chemists Toolbox for Ultrahigh_Throughput Experimentation|cleaned_suzuki_synthesis|sum"
    # product: calculated PCE and log Rg.
    "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices|Beyond molecular structure_seifrid_imputed|product"
    "Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation|Rg data with clusters aging imputed|product"
    # (count:+)x(graph:x): calculated PCE.
    "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices|Beyond molecular structure_seifrid_imputed|(count:+)x(graph:x)"
)

for job in "${structure_jobs[@]}"; do
    IFS='|' read -r paper dataset mixing_method <<< "$job"
    submit_job "$paper" "$dataset" "$mixing_method" "structure_cluster"
done

echo "Submitted ${job_index} MGK OOD GPU jobs."
