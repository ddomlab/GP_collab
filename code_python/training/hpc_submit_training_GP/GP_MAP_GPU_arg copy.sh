#!/bin/bash

DATE=$(date +%Y%m%d)
paper="Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation"
datasets=("cleaned_dataset_pervaporation_membranes_wang")
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


#BSUB -n 1
#BSUB -W 23:20
#BSUB -q gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -R "rusage[mem=8GB]"
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
