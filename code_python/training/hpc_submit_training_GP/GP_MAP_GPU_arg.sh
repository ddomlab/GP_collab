#!/bin/bash

DATE=$(date +%Y%m%d)
paper="Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning"
datasets=("cleaned_dataset_Ultrafiltration Membrane_imputed")
##flux_data_imputed
output_dir=/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}/${paper}
mkdir -p "$output_dir"

model="GPytorchMAP"
k_fps=("Tanimoto" "TanimotoRBF" "TanimotoMatern32" "TanimotoMatern52" "Matern32" "Matern52" "RBF")
k_counts=("Matern32" "Matern52" "RBF")
k_mixing_methods=("sum" "product" "averageProduct" "(count:+)x(fp:x)" "(count:x)+(fp:x)")


for mixing_method in "${k_mixing_methods[@]}"; do
    for dataset in "${datasets[@]}"; do
        for fp_kernel in "${k_fps[@]}"; do
            for count_kernel in "${k_counts[@]}"; do
                sbatch <<EOT
#!/bin/bash


#SBATCH --ntasks=1
#SBATCH --time=01:59:00
#SBATCH --partition=gpu_partners
#SBATCH --qos=short_gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=8G
#SBATCH --job-name="structure_numerical_${DATE}_${model}"
#SBATCH --output="${output_dir}/${model}_${fp_kernel}_${count_kernel}_${mixing_method}_GPU.out"
#SBATCH --error="${output_dir}/${model}_${fp_kernel}_${count_kernel}_${mixing_method}_GPU.err"

source ~/.bashrc
module load cuda/12.1
module load gcc/9.3.0
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12

python ../train_structure_numerical.py  --K_fp $fp_kernel \
                                        --K_count $count_kernel \
                                        --Kernel_mixing_method "$mixing_method" \
                                        --paper "$paper" \
                                        --dataset "$dataset" \
                                        --regressor_type "$model" 

EOT
            done
        done
    done
done
