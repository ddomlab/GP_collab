#!/bin/bash

DATE=$(date +%Y%m%d)
model="GpyroHMC"  #"GpyroHMC", "GPytorchMAP"
paper="Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery"
dataset="flux_data_imputed"
k_fps=("TanimotoMatern32")
k_counts=("Matern32")           
k_mixing_methods=("averageProduct" "product") 

##flux_data_imputed
output_dir=/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}/${paper}
mkdir -p "$output_dir"

for mixing_method in "${k_mixing_methods[@]}"; do
    for fp_kernel in "${k_fps[@]}"; do
        for count_kernel in "${k_counts[@]}"; do
            bsub <<EOT

#BSUB -n 1
#BSUB -W 35:35
#BSUB -q gpu
#BSUB -R "select[a10]"
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -R "rusage[mem=16GB]"
#BSUB -J "structure_numerical_${DATE}"
#BSUB -o "${output_dir}/${model}_${fp_kernel}_${count_kernel}_${mixing_method}_GPU.out"
#BSUB -e "${output_dir}/${model}_${fp_kernel}_${count_kernel}_${mixing_method}_GPU.err"

source ~/.bashrc
module load cuda/12.1
module load gcc/9.3.0
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