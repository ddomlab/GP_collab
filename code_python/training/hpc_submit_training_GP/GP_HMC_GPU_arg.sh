#!/bin/bash

DATE=$(date +%Y%m%d)
paper="Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery"
datasets=("flux_data_imputed" "separation_data_imputed")
##flux_data_imputed
output_dir=/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}/${paper}
mkdir -p "$output_dir"

model="GpyroHMC"
k_fps=("Tanimoto" "TanimotoRBF" "TanimotoMatern32" "TanimotoMatern52")
k_counts=("RBF")
k_mixing_methods=("sum" "product" "(count:+)x(fp:x)" "(count:x)+(fp:x)")


for mixing_method in "${k_mixing_methods[@]}"; do
    for dataset in "${datasets[@]}"; do
        for fp_kernel in "${k_fps[@]}"; do
            for count_kernel in "${k_counts[@]}"; do
                bsub <<EOT


#BSUB -n 1
#BSUB -W 60:10
#BSUB -q gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -R "rusage[mem=32GB]"
#BSUB -R "select[a10 || a30 || a100 || l40 || h100]"
#BSUB -J "structure_numerical_${DATE}_${model}"
#BSUB -o "${output_dir}/${model}_${fp_kernel}_${count_kernel}_${mixing_method}_GPU.out"
#BSUB -e "${output_dir}/${model}_${fp_kernel}_${count_kernel}_${mixing_method}_GPU.err"

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
