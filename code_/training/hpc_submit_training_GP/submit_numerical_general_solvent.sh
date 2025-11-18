#!/bin/bash
output_dir=/share/ddomlab/sdehgha2/working-space/main/P1_pls-dataset/pls-dataset-space/PLS-Dataset/results/hpc_out

# Correctly define models and numerical features
target_to_asses=("Rg1 (nm)" "Rh (IW avg log)")
models_to_run=("GPR")
numerical_feats=("polymer dP" "polymer dD" "polymer dH")
kernels=("matern" "rbf")


for target in "${target_to_asses[@]}"; do
        for model in "${models_to_run[@]}"; do
            for feats in "${numerical_feats[@]}"; do
                for kernel in "${kernels[@]}"; do
                    bsub <<EOT

#BSUB -n 8
#BSUB -W 10:01
#BSUB -R span[ptile=4]
#BSUB -R "rusage[mem=32GB]"
#BSUB -J "poly_HSP_with_${model}_on_${target}_signle_polyHSP"
#BSUB -o "${output_dir}/poly_HSP_with_${model}_on_${target}_signle_polyHSP.out"
#BSUB -e "${output_dir}/poly_HSP_with_${model}_on_${target}_signle_polyHSP.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/pls-dataset-env
python ../train_numerical_only.py --target_features "${target}" \
                                  --regressor_type "${model}" \
                                  --kernel "${kernel}" \
                                  --numerical_feats "${feats}" 





conda deactivate

EOT
            done
        done
    done
done