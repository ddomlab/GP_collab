#!/bin/bash

DATE=$(date +%Y%m%d)
model="MGK"  #"GpyroHMC", "GPytorchMAP"
paper="Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery"
dataset="separation_data_imputed"
k_fps=("Graph")
k_counts=("Matern32")           
k_mixing_methods=("product") 
k_feature_mode="per_feature" #for MGK

##flux_data_imputed
output_dir=/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}/${paper}
mkdir -p "$output_dir"

for mixing_method in "${k_mixing_methods[@]}"; do
    for fp_kernel in "${k_fps[@]}"; do
        for count_kernel in "${k_counts[@]}"; do
            bsub <<EOT

#BSUB -n 1
#BSUB -W 6:10
#BSUB -q gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -R "rusage[mem=4GB]"
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
                                        --regressor_type "$model" \
                                        --kernel_feature_mode "$k_feature_mode"




EOT

        done
    done
done

#BSUB -R "select[a10]"