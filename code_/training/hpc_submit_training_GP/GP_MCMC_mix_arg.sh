#!/bin/bash

DATE=$(date +%Y%m%d)
model="GPMixMCMC"
output_dir=/share/ddomlab/sdehgha2/working-space/main/_colab_GP/GP_collab/results/HPC_history/hpc_${DATE}
mkdir -p "$output_dir"

k_fps=("RBF" "Matern32" "Matern52")
k_counts=("Matern32" "Matern52" "RBF")
k_mixing_methods=("sum" "product" "averageProduct") 

for mixing_method in "${k_mixing_methods[@]}"; do
    for fp_kernel in "${k_fps[@]}"; do
        for count_kernel in "${k_counts[@]}"; do
            bsub <<EOT



#BSUB -n 6
#BSUB -W 33:00
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=32GB]"
#BSUB -J "structure_numerical_${DATE}"
#BSUB -o "${output_dir}/GpyroMCMC_${fp_kernel}_${count_kernel}_${mixing_method}.out"
#BSUB -e "${output_dir}/GpyroMCMC_${fp_kernel}_${count_kernel}_${mixing_method}.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/torch_cpu
python ../train_structure_numerical.py --K_fp $fp_kernel \
                                        --K_count $count_kernel \
                                        --Kernel_mixing_method $mixing_method 

EOT

        done
    done
done