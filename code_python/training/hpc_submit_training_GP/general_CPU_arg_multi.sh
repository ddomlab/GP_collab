#!/bin/bash

DATE=$(date +%Y%m%d)
# model="GPytorchMAPRegressor"
paper="Miniaturization of Popular Reactions from the Medicinal Chemists Toolbox for Ultrahigh_Throughput Experimentation"
dataset="cleaned_suzuki_synthesis"
models=("XGBR" "RF" "NGB")
# k_fps=("TanimotoRBF")
# k_counts=("RBF")
# k_mixing_methods=("sum") 

##flux_data_imputed
output_dir=/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}/${paper}
mkdir -p "$output_dir"

for model in "${models[@]}"; do
    bsub <<EOT

#BSUB -n 5
#BSUB -W 6
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=32GB]"
#BSUB -J "structure_numerical_${DATE}"
#BSUB -o "${output_dir}/GPytorchMAP_${fp_kernel}_${count_kernel}_${mixing_method}.out"
#BSUB -e "${output_dir}/GPytorchMAP_${fp_kernel}_${count_kernel}_${mixing_method}.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12
python ../train_structure_numerical.py --regressor_type $model \
                                        --dataset "$dataset" \
                                        --paper "$paper"


EOT
done