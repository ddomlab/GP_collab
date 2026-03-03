#!/bin/bash

DATE=$(date +%Y%m%d)
# model="GPMixMCMC"
paper="Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation"
dataset="Rg data with clusters aging imputed"
##flux_data_imputed
output_dir=/share/ddomlab/sdehgha2/working-space/main/_colab_GP/GP_collab/results/HPC_history/hpc_${DATE}/${paper}
mkdir -p "$output_dir"

fp_kernel="SSK"
count_kernel="RBF"
mixing_method="product"


bsub <<EOT


#BSUB -n 1
#BSUB -W 10
#BSUB -q gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -R "rusage[mem=32GB]"
#BSUB -J "structure_numerical_${DATE}"
#BSUB -o "${output_dir}/GPytorchMAP_${fp_kernel}_${count_kernel}_${mixing_method}_GPU.out"
#BSUB -e "${output_dir}/GPytorchMAP_${fp_kernel}_${count_kernel}_${mixing_method}_GPU.err"

source ~/.bashrc
module load cuda/12.1
conda activate /usr/local/usrapps/ddomlab/sdehgha2/torch_gpu

python ../train_structure_numerical.py


EOT

