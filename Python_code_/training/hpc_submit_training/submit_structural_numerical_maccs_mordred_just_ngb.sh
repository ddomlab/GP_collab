#!/bin/bash
output_dir=/share/ddomlab/sdehgha2/working-space/main/P1_pls-dataset/pls-dataset-space/PLS-Dataset/results
# Define arrays for regressor types, targets, and models
regressors=("NGB")
targets=("Rg1 (nm)" "Lp (nm)")
models=("maccs" "mordred")

# Define the output directory (make sure to create it beforehand)
output_dir="/path/to/output"  # Change this to the appropriate directory

# Loop through each combination of regressor, target, and model
for regressor in "${regressors[@]}"; do
  for target in "${targets[@]}"; do
    for model in "${models[@]}"; do
      bsub <<EOT
#BSUB -n 8
#BSUB -W 60:05
#BSUB -R span[ptile=4]
#BSUB -R "rusage[mem=32GB]"
#BSUB -J "${model}_${regressor}"  # Job name
#BSUB -o ${output_dir}/${model}_${regressor}.out
#BSUB -e ${output_dir}/${model}_${regressor}.err

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/pls-dataset-env
python train_structure_numerical.py $model --regressor_type $regressor --target "$target"
EOT
    done
  done
done
