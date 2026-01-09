#!/bin/bash
output_dir=/share/ddomlab/sdehgha2/working-space/main/P1_pls-dataset/pls-dataset-space/PLS-Dataset/results

# Correctly define models and numerical features
target_to_asses=("Rh (IW avg log)")
models_to_run=("GPR")
kernels=("matern" "rbf")
scaler_types=('Standard' 'Robust Scaler')

for target in "${target_to_asses[@]}"; do
    for model in "${models_to_run[@]}"; do
        for kernel in "${kernels[@]}"; do
            for scaler in "${scaler_types[@]}"; do
                bsub <<EOT

#BSUB -n 8
#BSUB -W 48:02
#BSUB -R span[ptile=4]
#BSUB -R "rusage[mem=32GB]"
#BSUB -J "poly_HSP_with_${model}_on_${target}_on_polysize_solvent"
#BSUB -o "${output_dir}/polysize_solvent_with_${model}_on_${target}.out"
#BSUB -e "${output_dir}/polysize_solvent_with_${model}_on_${target}.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/pls-dataset-env
python ../train_numerical_only.py --target_features "${target}" \
                                  --regressor_type "${model}" \
                                  --kernel "${kernel}" \
                                  --transform_type "${scaler}" \
                                  --numerical_feats "Concentration (mg/ml)" "Temperature SANS/SLS/DLS/SEC (K)" "solvent dP" "solvent dD" "solvent dH"\
                                  --columns_to_impute "Concentration (mg/ml)" "Temperature SANS/SLS/DLS/SEC (K)" \
                                  --imputer mean
                                  
conda deactivate

EOT 
            done
        done
    done
done