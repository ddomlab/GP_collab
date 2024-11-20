#!/bin/bash
output_dir=/share/ddomlab/sdehgha2/working-space/main/P1_pls-dataset/pls-dataset-space/PLS-Dataset/results
# Define arrays for regressor types, targets, and models
regressors=("NGB" "XGBR" "RF")
targets=("Rh (IW avg log)")
models=("mordred" "maccs")
poly_representations=('Monomer' 'Dimer' 'Trimer' 'RRU Monomer' 'RRU Dimer' 'RRU Trimer')
scaler_types=('Standard')




# Loop through each combination of regressor, target, and model
for regressor in "${regressors[@]}"; do
  for target in "${targets[@]}"; do
    for model in "${models[@]}"; do
      for oligo_rep in "${poly_representations[@]}"; do
        for scaler in "${scaler_types[@]}"; do

      bsub <<EOT

#BSUB -n 8
#BSUB -W 60:01
#BSUB -R span[ptile=4]
#BSUB -R "rusage[mem=32GB]"
#BSUB -J "structure_only_mordred_${regressor}_with_${scaler}_on_${targets}"  
#BSUB -o "${output_dir}/structure_only_mordred_${regressor}_with_${scaler}_on_${targets}.out"
#BSUB -e "${output_dir}/structure_only_mordred_${regressor}_with_${scaler}_on_${targets}.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/pls-dataset-env
python ../train_structure_only.py $model --regressor_type $regressor --target "$target" --oligo_type "$oligo_rep" --transform_type "$scaler"

EOT   
        done
      done
    done
  done
done