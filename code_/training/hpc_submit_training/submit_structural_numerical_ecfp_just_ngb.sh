#!/bin/bash
output_dir=/share/ddomlab/sdehgha2/working-space/main/P1_pls-dataset/pls-dataset-space/PLS-Dataset/results
# Define arrays for regressor types, targets, and models
regressors=("NGB")
targets=("Rg1 (nm)" "Lp (nm)")


# Loop through each combination of regressor, target, and model
for regressor in "${regressors[@]}"; do
  for target in "${targets[@]}"; do
    for radius in 3 4 5 6; do
      for vector in "count" "binary"; do
        # Submitting the job using bsub
        bsub <<EOT
#BSUB -n 8
#BSUB -W 60:05
#BSUB -R span[ptile=4]
#BSUB -R "rusage[mem=32GB]"
#BSUB -J "ecfp_radius${radius}_vector${vector}"  # Job name
#BSUB -o ${output_dir}/ecfp_run_radius${radius}_vector${vector}.out
#BSUB -e ${output_dir}/ecfp_err_radius${radius}_vector${vector}.out

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/pls-dataset-env
python train_structure_numerical.py ecfp --regressor_type $regressor --radius $radius --vector $vector --target "$target"
EOT
      done
    done
  done
done
