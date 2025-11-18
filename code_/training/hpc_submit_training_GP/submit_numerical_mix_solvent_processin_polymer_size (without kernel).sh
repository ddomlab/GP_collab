#!/bin/bash
output_dir=/share/ddomlab/sdehgha2/working-space/main/P1_pls-dataset/pls-dataset-space/PLS-Dataset/results

# Correctly define models and numerical features
target_to_asses=('Rh (1_1000 nm) (highest intensity)')
models_to_run=("RF" "XGBR")
# scaler_types=("Robust Scaler")

for target in "${target_to_asses[@]}"; do
    for model in "${models_to_run[@]}"; do
        # for scaler in "${scaler_types[@]}"; do
            bsub <<EOT

#BSUB -n 6
#BSUB -W 25:01
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=16GB]"
#BSUB -J "numerical_${model}_polymer_size_feats_on_${target}_all_num_20250220"
#BSUB -o "${output_dir}/numerical_${model}__${target}_20250220.out"
#BSUB -e "${output_dir}/numerical_${model}__${target}_20250220.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/pls-dataset-env
python ../train_numerical_only.py --target_features "${target}" \
                                    --regressor_type "${model}" \
                                    --numerical_feats 'PDI' 'Mn (g/mol)' 'Mw (g/mol)' 'Concentration (mg/ml)' "Temperature SANS/SLS/DLS/SEC (K)" "polymer dP" "polymer dD" "polymer dH" "solvent dP" "solvent dD" "solvent dH" \
                                    --columns_to_impute "PDI" "Temperature SANS/SLS/DLS/SEC (K)" "Concentration (mg/ml)" \
                                    --special_impute 'Mw (g/mol)' \
                                    --imputer mean


conda deactivate

EOT
        # done
    done
done
                                    # --transform_type "${scaler}" \

# "PDI" "Mw"
# "concentration" "temperature"
# "PDI" "Mw" "concentration" "temperature"
# "PDI" "Mw" "concentration" "temperature" "Ra"
# "PDI" "Mw" "concentration" "temperature" "solvent dP" "solvent dD" "solvent dH"
# "PDI" "Mw" "concentration" "temperature" "polymer dP" "polymer dD" "polymer dH" "solvent dP" "solvent dD" "solvent dH"
# "solvent dP" "solvent dD" "solvent dH"
# "polymer dP" "polymer dD" "polymer dH"
# "polymer dP" "polymer dD" "polymer dH" "solvent dP" "solvent dD" "solvent dH"

# "polymer dP" "polymer dD" "polymer dH" "solvent dP" "solvent dD" "solvent dH" "Ra"
# 'Mn (g/mol)' 'PDI' 'Mw (g/mol)' "Concentration (mg/ml)" "Temperature SANS/SLS/DLS/SEC (K)" "solvent dP" "solvent dD" "solvent dH"
# --special_impute 'Mw (g/mol)' \
# --columns_to_impute 'Concentration (mg/ml)' 'Temperature SANS/SLS/DLS/SEC (K)' \
# --imputer mean

#                                     --columns_to_impute "PDI" "Temperature SANS/SLS/DLS/SEC (K)" "Concentration (mg/ml)" \
    #  --special_impute 'Mw (g/mol)' \
    #                                 --imputer mean

    #--special_impute 'Mw (g/mol)' \