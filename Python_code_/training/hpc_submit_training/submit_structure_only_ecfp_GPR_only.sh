#!/bin/bash
output_dir=/share/ddomlab/sdehgha2/working-space/main/P1_pls-dataset/pls-dataset-space/PLS-Dataset/results
# Define arrays for regressor types, targets, and models
regressors=("GPR")
targets=("Rg1 (nm)")
radii=(6) 
poly_representations=('Trimer' 'RRU Monomer')
vectors=("binary" "count")

# Loop through each combination of regressor, target, and model
for regressor in "${regressors[@]}"; do
  for target in "${targets[@]}"; do
    for radius in "${radii[@]}"; do
      for vector in "${vectors[@]}"; do
        for oligo_rep in "${poly_representations[@]}"; do
      
        bsub <<EOT
#BSUB -n 8
#BSUB -W 30:01
#BSUB -R span[ptile=4]
#BSUB -R "rusage[mem=32GB]"
#BSUB -J "ecfp_radius_tructure_only_${regressor}" 
#BSUB -o "${output_dir}/structure_only_ecfp_${regressor}_with_RScaler.out"
#BSUB -e "${output_dir}/structure_only_ecfp_${regressor}_with_RScaler.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/pls-dataset-env
python ../train_structure_only.py ecfp \
    --regressor_type "$regressor" \
    --radius "$radius" \
    --vector "$vector" \
    --target "$target" \
    --oligo_type "$oligo_rep" \
    --kernel "tanimoto"

EOT
        done
      done
    done
  done
done