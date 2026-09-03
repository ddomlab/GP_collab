#!/bin/bash

DATE=$(date +%Y%m%d)
model="MGK"
paper="Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning"
datasets=("cleaned_dataset_Ultrafiltration Membrane_imputed")
output_dir=/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}/${paper}
mkdir -p "$output_dir"


k_mixing_methods=("(count:+)x(graph:x)" "(count:+)x(graph:+)" "(count:x)+(graph:x)" "product" "sum")
k_feature_modes=("per_feature")
k_fps=("Graph")
k_counts=("Matern32")


for mixing_method in "${k_mixing_methods[@]}"; do
    for feature_mode in "${k_feature_modes[@]}"; do
        for dataset in "${datasets[@]}"; do
            for fp_kernel in "${k_fps[@]}"; do
                for count_kernel in "${k_counts[@]}"; do
                    sbatch <<EOT
#!/bin/bash



#SBATCH --ntasks=1
#SBATCH --time=01:55:00
#SBATCH --partition=gpu_partners
#SBATCH --qos=short_gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --job-name="structure_numerical_${DATE}"
#SBATCH --output="${output_dir}/${model}_${mixing_method}_${feature_mode}_GPU.out"
#SBATCH --error="${output_dir}/${model}_${mixing_method}_${feature_mode}_GPU.err"

source ~/.bashrc
module load cuda/12.1
module load gcc/9.3.0
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12

python ../train_structure_numerical.py  --K_fp $fp_kernel \
                                        --K_count $count_kernel \
                                        --Kernel_mixing_method "$mixing_method" \
                                        --paper "$paper" \
                                        --dataset "$dataset" \
                                        --kernel_feature_mode $feature_mode \
                                        --regressor_type "$model" 


EOT
    
                done
            done
        done
    done
done
