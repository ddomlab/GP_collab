#!/bin/bash

# Submit missing GPytorchMAP structure- and continuous-cluster OOD runs for
# Matern32/TanimotoMatern32 FP kernels, the Matern32 count kernel, and the four
# selected mixing methods. Multi-target datasets are submitted once per
# configuration even when only one target result is missing.

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
    local clustering_method="$6"
    local output_dir="${output_root}/${paper}"
    local dataset_tag=${dataset// /_}

    job_index=$((job_index + 1))
    mkdir -p "$output_dir"

    bsub <<EOT
#BSUB -n 1
#BSUB -W 30:20
#BSUB -q gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -R "rusage[mem=32GB]"
#BSUB -R "select[a10 || a30 || a100 || l40 || h100]"
#BSUB -J "gpmap_ood_missing_${DATE}_${job_index}"
#BSUB -o "${output_dir}/${model}_${dataset_tag}_${fp_kernel}_${count_kernel}_${mixing_method}_${clustering_method}_GPU.out"
#BSUB -e "${output_dir}/${model}_${dataset_tag}_${fp_kernel}_${count_kernel}_${mixing_method}_${clustering_method}_GPU.err"

source ~/.bashrc
module load cuda/12.1
module load gcc/9.3.0
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12

python ../train_structure_numerical.py --K_fp "$fp_kernel" \
                                        --K_count "$count_kernel" \
                                        --Kernel_mixing_method "$mixing_method" \
                                        --paper "$paper" \
                                        --dataset "$dataset" \
                                        --regressor_type "$model" \
                                        --clustering_method "$clustering_method"
EOT
}

missing_jobs=(
    # structure_cluster: 13 missing target rows represented by 7 jobs.
    "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices|Beyond molecular structure_seifrid_imputed|Matern32|Matern32|(count:+)x(fp:+)|structure_cluster"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed|Matern32|Matern32|(count:+)x(fp:+)|structure_cluster"
    "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation|cleaned_dataset_pervaporation_membranes_wang|Matern32|Matern32|(count:+)x(fp:+)|structure_cluster"
    "Miniaturization of Popular Reactions from the Medicinal Chemists Toolbox for Ultrahigh_Throughput Experimentation|cleaned_suzuki_synthesis|Matern32|Matern32|(count:+)x(fp:+)|structure_cluster"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|Matern32|Matern32|(count:x)+(fp:x)|structure_cluster"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed|TanimotoMatern32|Matern32|product|structure_cluster"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|flux_data_imputed|TanimotoMatern32|Matern32|(count:+)x(fp:+)|structure_cluster"

    # continuous_cluster: 31 missing target rows represented by 24 jobs.
    # Calculated PCE (%).
    "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices|Beyond molecular structure_seifrid_imputed|TanimotoMatern32|Matern32|sum|continuous_cluster"
    "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices|Beyond molecular structure_seifrid_imputed|TanimotoMatern32|Matern32|product|continuous_cluster"
    "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices|Beyond molecular structure_seifrid_imputed|TanimotoMatern32|Matern32|(count:x)+(fp:x)|continuous_cluster"
    "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices|Beyond molecular structure_seifrid_imputed|TanimotoMatern32|Matern32|(count:+)x(fp:+)|continuous_cluster"
    "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices|Beyond molecular structure_seifrid_imputed|Matern32|Matern32|sum|continuous_cluster"

    # Organic-recovery separation factor.
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed|TanimotoMatern32|Matern32|sum|continuous_cluster"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed|TanimotoMatern32|Matern32|product|continuous_cluster"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed|TanimotoMatern32|Matern32|(count:x)+(fp:x)|continuous_cluster"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed|TanimotoMatern32|Matern32|(count:+)x(fp:+)|continuous_cluster"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed|Matern32|Matern32|sum|continuous_cluster"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed|Matern32|Matern32|(count:x)+(fp:x)|continuous_cluster"

    # Organic-recovery total flux.
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|flux_data_imputed|TanimotoMatern32|Matern32|sum|continuous_cluster"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|flux_data_imputed|TanimotoMatern32|Matern32|product|continuous_cluster"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|flux_data_imputed|TanimotoMatern32|Matern32|(count:x)+(fp:x)|continuous_cluster"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|flux_data_imputed|TanimotoMatern32|Matern32|(count:+)x(fp:+)|continuous_cluster"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|flux_data_imputed|Matern32|Matern32|product|continuous_cluster"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|flux_data_imputed|Matern32|Matern32|(count:x)+(fp:x)|continuous_cluster"

    # Both Wang pervaporation targets are trained by each invocation.
    "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation|cleaned_dataset_pervaporation_membranes_wang|Matern32|Matern32|sum|continuous_cluster"
    "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation|cleaned_dataset_pervaporation_membranes_wang|Matern32|Matern32|(count:x)+(fp:x)|continuous_cluster"

    # Ultrafiltration: the first job supplies all six missing targets; the
    # remaining three rerun the dataset to supply irreversible fouling.
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|TanimotoMatern32|Matern32|product|continuous_cluster"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|Matern32|Matern32|product|continuous_cluster"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|Matern32|Matern32|(count:x)+(fp:x)|continuous_cluster"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed|Matern32|Matern32|(count:+)x(fp:+)|continuous_cluster"

    # Approx Conv (%).
    "Miniaturization of Popular Reactions from the Medicinal Chemists Toolbox for Ultrahigh_Throughput Experimentation|cleaned_suzuki_synthesis|Matern32|Matern32|product|continuous_cluster"
)

for job in "${missing_jobs[@]}"; do
    IFS='|' read -r paper dataset fp_kernel count_kernel mixing_method clustering_method <<< "$job"
    submit_job \
        "$paper" \
        "$dataset" \
        "$fp_kernel" \
        "$count_kernel" \
        "$mixing_method" \
        "$clustering_method"
done

echo "Submitted ${job_index} GPytorchMAP OOD GPU jobs."
