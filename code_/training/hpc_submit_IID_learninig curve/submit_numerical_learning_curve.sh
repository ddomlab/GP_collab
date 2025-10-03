#!/bin/bash
output_dir=/share/ddomlab/sdehgha2/working-space/main/P1_pls-dataset/pls-dataset-space/PLS-Dataset/results/hpc_20250430
mkdir -p "$output_dir"

regressors=('XGBR' 'NGB' 'RF')
targets=('log Rg (nm)')



for regressor in "${regressors[@]}"; do
  for target in "${targets[@]}"; do
      bsub <<EOT

#BSUB -n 6
#BSUB -W 10:01
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=8GB]"
#BSUB -J "${regressor}_${target}_20250430"  
#BSUB -o "${output_dir}/${regressor}_${target}_20250430.out"
#BSUB -e "${output_dir}/${regressor}_${target}_20250430.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/pls-dataset-env

python ../train_structure_numerical_leaning_curve.py --target_features '${target}' \
                                                    --regressor_type '${regressor}' \
                                                    --numerical_feats 'Xn' 'Mw (g/mol)' 'PDI'  'Concentration (mg/ml)' 'Temperature SANS/SLS/DLS/SEC (K)' "polymer dP" "polymer dD" "polymer dH" 'solvent dP' 'solvent dD' 'solvent dH' "Dark/light" "Aging time (hour)" "To Aging Temperature (K)" "Sonication/Stirring/heating Temperature (K)" "Merged Stirring /sonication/heating time(min)" 

EOT

  done
done
