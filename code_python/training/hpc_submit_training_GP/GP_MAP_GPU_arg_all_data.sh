#!/bin/bash

# Submit all GPytorchMAP kernel configurations for the explicitly selected
# datasets. Each invocation trains every target configured for the corresponding
# paper/dataset pair in filter_data.py.

DATE=$(date +%Y%m%d)
model="GPytorchMAP"
k_fps=("TanimotoMatern32")
k_counts=("Matern32")
k_mixing_methods=("sum" "product" "(count:x)+(fp:x)" "(count:+)x(fp:+)")

selected_training_sets=(
    # "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices|Beyond molecular structure_seifrid_imputed"
    # "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|flux_data_imputed"
    # "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed"
    # "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation|cleaned_dataset_pervaporation_membranes_wang"
    "Miniaturization of Popular Reactions from the Medicinal Chemists Toolbox for Ultrahigh_Throughput Experimentation|cleaned_suzuki_synthesis"
    # "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed"
    "Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation|Rg data with clusters aging imputed"
)

output_root="/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}"
job_index=0

for training_set in "${selected_training_sets[@]}"; do
    IFS='|' read -r paper dataset <<< "$training_set"
    output_dir="${output_root}/${paper}"
    dataset_tag=${dataset// /_}
    mkdir -p "$output_dir"

    for mixing_method in "${k_mixing_methods[@]}"; do
        for fp_kernel in "${k_fps[@]}"; do
            for count_kernel in "${k_counts[@]}"; do
                job_index=$((job_index + 1))

                bsub <<EOT
#BSUB -n 1
#BSUB -W 1:59
#BSUB -q short_gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -R "rusage[mem=32GB]"
#BSUB -R "select[a10 || a30 || a100 || l40 || h100]"
#BSUB -J "gpmap_all_${DATE}_${job_index}"
#BSUB -o "${output_dir}/${model}_${dataset_tag}_${fp_kernel}_${count_kernel}_${mixing_method}_GPU.out"
#BSUB -e "${output_dir}/${model}_${dataset_tag}_${fp_kernel}_${count_kernel}_${mixing_method}_GPU.err"

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
            done
        done
    done
done

echo "Submitted ${job_index} GPytorchMAP GPU jobs."
