#!/bin/bash

# DATE=$(date +%Y%m%d)

output_dir=/share/ddomlab/sdehgha2/working-space/main/P1_pls-dataset/pls-dataset-space/PLS-Dataset/results/HPC_history/hpc_2025
mkdir -p "$output_dir"

bsub <<EOT

#BSUB -n 6
#BSUB -W 04:00
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=16GB]"
#BSUB -J "structural_numerical_2025"
#BSUB -o "${output_dir}/structural_numerical_2025.out"
#BSUB -e "${output_dir}/structural_numerical_2025.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/torch_cpu
python ../train_numerical_only.py

EOT