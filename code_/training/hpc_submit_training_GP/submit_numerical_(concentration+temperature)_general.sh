#!/bin/bash
output_dir=/share/ddomlab/sdehgha2/working-space/main/P1_pls-dataset/pls-dataset-space/PLS-Dataset/results

# Correctly define models and numerical features
target_to_asses=("Rh (IW avg log)" "Rg1 (nm)")
models_to_run=("RF" "MLR" "DT")
numerical_feat_list=("abs(solvent dD - polymer dD)" "abs(solvent dP - polymer dP)" "abs(solvent dH - polymer dH)")

for target in "${target_to_asses[@]}"; do
        for model in "${models_to_run[@]}"; do
            for feats in "${numerical_feat_list[@]}"; do
                
                bsub <<EOT

#BSUB -n 8
#BSUB -W 15:01
#BSUB -R span[ptile=4]
##BSUB -x
#BSUB -R "rusage[mem=32GB]"
#BSUB -J "${model}_with_${feats}_on_${target}"
#BSUB -o "${output_dir}/${feats}_with_${model}_with_${feats}_on_${target}.out"
#BSUB -e "${output_dir}/${feats}_with_${model}_with_${feats}_on_${target}.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/pls-dataset-env

python ../train_numerical_only.py --target_features "${target}" \
                                  --regressor_type "${model}" \
                                  --numerical_feats "${feats}" \

conda deactivate

EOT
        done
    done
done


