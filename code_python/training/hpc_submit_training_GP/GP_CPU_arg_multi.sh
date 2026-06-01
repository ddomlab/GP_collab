#!/bin/bash

DATE=$(date +%Y%m%d)
model="GpyroHMC"
paper="Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation"
dataset="Rg data with clusters aging imputed"
k_fps=("TanimotoMatern32")
k_counts=("Matern32")           
k_mixing_methods=("averageProduct") 
##flux_data_imputed
output_dir=/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}/${paper}
mkdir -p "$output_dir"

for mixing_method in "${k_mixing_methods[@]}"; do
    for fp_kernel in "${k_fps[@]}"; do
        for count_kernel in "${k_counts[@]}"; do
            bsub <<EOT

#BSUB -n 6
#BSUB -W 50
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=32GB]"
#BSUB -J "structure_numerical_${DATE}"
#BSUB -o "${output_dir}/${model}_${fp_kernel}_${count_kernel}_${mixing_method}.out"
#BSUB -e "${output_dir}/${model}_${fp_kernel}_${count_kernel}_${mixing_method}.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12
python ../train_structure_numerical.py --K_fp $fp_kernel \
                                        --K_count $count_kernel \
                                        --Kernel_mixing_method $mixing_method \
                                        --dataset "$dataset" \
                                        --paper "$paper" \
                                        --regressor_type "$model" 


EOT

        done
    done
done