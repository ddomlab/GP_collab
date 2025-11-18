#!/bin/bash

DATE=$(date +%Y%m%d)

output_dir=/share/ddomlab/sdehgha2/working-space/main/P1_pls-dataset/pls-dataset-space/PLS-Dataset/results/HPC_history/hpc_${DATE}
mkdir -p "$output_dir"

bsub <<EOT

#BSUB -n 6
#BSUB -W 04:00
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=16GB]"
#BSUB -J "structural_numerical_${DATE}"
#BSUB -o "${output_dir}/structural_numerical_${DATE}.out"
#BSUB -e "${output_dir}/structural_numerical_${DATE}.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/pls-dataset-env
python ../train_structure_numerical.py

EOT