#!/bin/bash

DATE=$(date +%Y%m%d)
model="GPytorchMAPRegressor"
paper="Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation"
dataset="Rg data with clusters aging imputed"
mixing_method="product"

##flux_data_imputed
output_dir=/share/ddomlab/sdehgha2/working-space/main/_colab_GP/GP_collab/results/HPC_history/hpc_${DATE}/${paper}
mkdir -p "$output_dir"

bsub <<EOT

#BSUB -n 6
#BSUB -W 2:00
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=32GB]"
#BSUB -J "structure_numerical_${DATE}"
#BSUB -o "${output_dir}/MAP_${mixing_method}_${feature_mode}_CPU.out"
#BSUB -e "${output_dir}/MAP_${mixing_method}_${feature_mode}_CPU.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/torch_cpu
python ../train_structure_numerical.py --dataset "$dataset" \
                                        --paper "$paper" \
                                        --Kernel_mixing_method $mixing_method \

EOT