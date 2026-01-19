#!/bin/bash

DATE=$(date +%Y%m%d)
model="RF-XGB-NGB"
paper="Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery"
dataset="separation_data_imputed"
##flux_data_imputed
output_dir=/share/ddomlab/sdehgha2/working-space/main/_colab_GP/GP_collab/results/HPC_history/hpc_${DATE}/${paper}
mkdir -p "$output_dir"

bsub <<EOT

#BSUB -n 6
#BSUB -W 2:00
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=32GB]"
#BSUB -J "structure_numerical_${DATE}"
#BSUB -o "${output_dir}/structure_numerical${model}_${DATE}.out"
#BSUB -e "${output_dir}/structure_numerical${model}_${DATE}.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/torch_cpu
python ../train_structure_numerical.py --dataset "$dataset" \
                                        --paper "$paper" \
                                         

EOT