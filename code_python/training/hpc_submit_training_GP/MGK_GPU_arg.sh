#!/bin/bash

DATE=$(date +%Y%m%d)
# model="GPMixMCMC"
# paper="Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation"
# dataset="Rg data with clusters aging imputed"
paper="Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation"
datasets=("Rg data with clusters aging imputed")
##flux_data_imputed
output_dir=/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}/${paper}
mkdir -p "$output_dir"


k_mixing_methods=("product")
k_feature_modes=("per_feature")

for mixing_method in "${k_mixing_methods[@]}"; do
    for feature_mode in "${k_feature_modes[@]}"; do
        for dataset in "${datasets[@]}"; do
            bsub <<EOT


#BSUB -n 1
#BSUB -W 00:45
#BSUB -q gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -R "rusage[mem=32GB]"
#BSUB -J "structure_numerical_${DATE}"
#BSUB -o "${output_dir}/MGK_${mixing_method}_${feature_mode}_GPU.out"
#BSUB -e "${output_dir}/MGK_${mixing_method}_${feature_mode}_GPU.err"

source ~/.bashrc
module load cuda/12.1
module load gcc/9.3.0
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12

python ../train_structure_numerical.py  --Kernel_mixing_method $mixing_method \
                                        --paper "$paper" \
                                        --dataset "$dataset" \
                                        --kernel_feature_mode $feature_mode

EOT
        done
    done
done