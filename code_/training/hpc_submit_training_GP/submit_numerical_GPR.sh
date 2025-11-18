#!/bin/bash
output_dir=/share/ddomlab/sdehgha2/working-space/main/P1_pls-dataset/pls-dataset-space/PLS-Dataset/results
# Define arrays for regressor types, targets, and models
target_to_assess=('log Rg (nm)')
models_to_run=('sklearn-GPR')
scalers=('Standard' 'Robust Scaler')
kernels=('matern' 'rbf')
# Loop through each combination of regressor, target, and model
for target in "${target_to_assess[@]}"; do
    for model in "${models_to_run[@]}"; do
        for kernel in "${kernels[@]}"; do
            for scaler in "${scalers[@]}"; do

              bsub <<EOT

#BSUB -n 8
#BSUB -W 20:01
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=16GB]"
#BSUB -J "${model}" 
#BSUB -o "${output_dir}/numerical_feats_${model}_${kernel}_${scaler}.out"
#BSUB -e "${output_dir}/numerical_feats_${model}_${kernel}_${scaler}.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/pls-dataset-env
python ../train_numerical_only.py --target_features "${target}" \
                                  --regressor_type "${model}" \
                                  --kernel "${kernel}" \
                                  --transform_type "${scaler}" \
                                  --numerical_feats 'Xn' 'Mw (g/mol)' 'PDI' 'Concentration (mg/ml)' 'Temperature SANS/SLS/DLS/SEC (K)' "polymer dP" "polymer dD" "polymer dH" 'solvent dP' 'solvent dD' 'solvent dH'

conda deactivate

EOT
            done
        done
    done
done