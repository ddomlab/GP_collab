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

#BSUB -n 6
#BSUB -W 48:02
#BSUB -R span[ptile=2]
#BSUB -R "rusage[mem=16GB]"
#BSUB -J "poly_HSP_with_${model}_on_${target}_no_Ra_poly"
#BSUB -o "${output_dir}/solvent_HSP_with_${model}_on_${target}_no_Ra_poly.out"
#BSUB -e "${output_dir}/solvent_HSP_with_${model}_on_${target}_no_Ra_poly.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/pls-dataset-env
python ../train_numerical_only.py --target_features "${target}" \
                                  --regressor_type "${model}" \
                                  --kernel "${kernel}" \
                                  --transform_type "${scaler}" \
                                  --numerical_feats "polymer dP" "polymer dD" "polymer dH" "solvent dP" "solvent dD" "solvent dH" "Ra"
                                  
conda deactivate

EOT     
            done
        done
    done
done

