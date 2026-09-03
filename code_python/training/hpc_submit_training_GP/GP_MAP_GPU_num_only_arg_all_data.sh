#!/bin/bash

DATE=$(date +%Y%m%d)
model="GPytorchMAP"
k_counts=("Matern32" "Matern52" "RBF")
# With no fingerprint kernels, the other grouped mixing expressions collapse
# to one of these two count-kernel compositions.
k_mixing_methods=("sum" "product")

selected_training_sets=(
    "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices|Beyond molecular structure_seifrid_imputed"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|flux_data_imputed"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed"
    "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation|cleaned_dataset_pervaporation_membranes_wang"
    "Miniaturization of Popular Reactions from the Medicinal Chemists Toolbox for Ultrahigh_Throughput Experimentation|cleaned_suzuki_synthesis"
    "Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation|Rg data with clusters aging imputed"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed"
)

output_root="/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}"
job_index=0

for training_set in "${selected_training_sets[@]}"; do
    IFS='|' read -r paper dataset <<< "$training_set"
    output_dir="${output_root}/${paper}"
    dataset_tag=${dataset// /_}
    mkdir -p "$output_dir"

    for mixing_method in "${k_mixing_methods[@]}"; do
        for count_kernel in "${k_counts[@]}"; do
            job_index=$((job_index + 1))

            sbatch <<EOT
#!/bin/bash


#SBATCH --ntasks=1
#SBATCH --time=01:59:00
#SBATCH --partition=gpu_partners
#SBATCH --qos=short_gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --job-name="numerical_only_${DATE}_${job_index}"
#SBATCH --output="${output_dir}/${model}_${dataset_tag}_${count_kernel}_${mixing_method}_GPU.out"
#SBATCH --error="${output_dir}/${model}_${dataset_tag}_${count_kernel}_${mixing_method}_GPU.err"

source ~/.bashrc
module load cuda/12.1
module load gcc/9.3.0
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12

python ../train_numerical.py --K_count "$count_kernel" \
                             --Kernel_mixing_method "$mixing_method" \
                             --paper "$paper" \
                             --dataset "$dataset" \
                             --regressor_type "$model"

EOT
        done
    done
done
