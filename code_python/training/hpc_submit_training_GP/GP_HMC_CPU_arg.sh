#!/bin/bash

DATE=$(date +%Y%m%d)
paper="Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery"
datasets=("flux_data_imputed" "separation_data_imputed")
##flux_data_imputed
output_dir=/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}/${paper}
mkdir -p "$output_dir"

model="GpyroHMC"
k_fps=("Tanimoto" "TanimotoRBF" "TanimotoMatern32" "TanimotoMatern52" "Matern32" "Matern52" "RBF")
k_counts=("Matern32" "Matern52" "RBF")
k_mixing_methods=("sum" "product" "averageProduct" "(count:+)x(fp:x)" "(count:+)x(fp:+)" "(count:x)+(fp:x)")


for mixing_method in "${k_mixing_methods[@]}"; do
    for dataset in "${datasets[@]}"; do
        for fp_kernel in "${k_fps[@]}"; do
            for count_kernel in "${k_counts[@]}"; do
                bsub <<EOT


#BSUB -n 4
#BSUB -W 40:50
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=32GB]"
#BSUB -J "structure_numerical_${DATE}_${model}_CPU"
#BSUB -o "${output_dir}/${model}_${fp_kernel}_${count_kernel}_${mixing_method}_CPU.out"
#BSUB -e "${output_dir}/${model}_${fp_kernel}_${count_kernel}_${mixing_method}_CPU.err"

source ~/.bashrc
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
