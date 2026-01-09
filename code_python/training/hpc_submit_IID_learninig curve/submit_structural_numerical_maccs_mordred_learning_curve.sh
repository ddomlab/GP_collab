#!/bin/bash
output_dir=/share/ddomlab/sdehgha2/working-space/main/P1_pls-dataset/pls-dataset-space/PLS-Dataset/results/hpc_20250407
mkdir -p "$output_dir"

regressors=("NGB")
targets=('log Rg (nm)')
models=("Mordred")
poly_representations=('Trimer')


for regressor in "${regressors[@]}"; do
  for target in "${targets[@]}"; do
    for fp in "${models[@]}"; do
      for oligo_rep in "${poly_representations[@]}"; do
          bsub <<EOT

#BSUB -n 6
#BSUB -W 7:01
#BSUB -R span[ptile=1]
#BSUB -R "rusage[mem=16GB]"
#BSUB -J "structure_numerical_mordred_NGB_generalizibility_wo_hypo"  
#BSUB -o "${output_dir}/structure_numerical_mordred_generalizability_${regressor}_wo_hypo.out"
#BSUB -e "${output_dir}/structure_numerical_mordred_generalizability_${regressor}_wo_hypo.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/pls-dataset-env

python ../train_structure_numerical_leaning_curve.py --target_features "${target}" \
                                                    --representation "${fp}" \
                                                    --regressor_type "${regressor}" \
                                                    --oligomer_representation "${oligo_rep}" \
                                                    --numerical_feats 'Mw (g/mol)' 'PDI' 'Concentration (mg/ml)' 'Temperature SANS/SLS/DLS/SEC (K)' "polymer dP" "polymer dD" "polymer dH" 'solvent dP' 'solvent dD' 'solvent dH' \

EOT
      done
    done
  done
done
