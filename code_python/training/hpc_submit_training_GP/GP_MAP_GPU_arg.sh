#!/bin/bash

DATE=$(date +%Y%m%d)
# model="GPMixMCMC"
paper="Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation"
dataset="Rg data with clusters aging imputed"
##flux_data_imputed
output_dir=/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}/${paper}
mkdir -p "$output_dir"

k_fps=("SSK")
k_counts=("RBF")
k_mixing_methods=("product")


for mixing_method in "${k_mixing_methods[@]}"; do
    for fp_kernel in "${k_fps[@]}"; do
        for count_kernel in "${k_counts[@]}"; do
            bsub <<EOT


#BSUB -n 2
#BSUB -W 40
#BSUB -q gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -R "rusage[mem=32GB]"
#BSUB -J "structure_numerical_${DATE}"
#BSUB -o "${output_dir}/GPytorchMAP_${fp_kernel}_${count_kernel}_${mixing_method}_GPU.out"
#BSUB -e "${output_dir}/GPytorchMAP_${fp_kernel}_${count_kernel}_${mixing_method}_GPU.err"

source ~/.bashrc
module load cuda/12.6
conda activate /usr/local/usrapps/ddomlab/sdehgha2/gpu_please

python ../train_structure_numerical.py --K_fp $fp_kernel \
                                        --K_count $count_kernel \
                                        --Kernel_mixing_method $mixing_method \
                                        --paper "$paper" \
                                        --dataset "$dataset"

EOT

        done
    done
done