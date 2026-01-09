#!/bin/bash
output_dir=/share/ddomlab/sdehgha2/working-space/main/P1_pls-dataset/pls-dataset-space/PLS-Dataset/results/hpc_20250321
mkdir -p "$output_dir"

regressors=("RF") 
targets=('log Rg (nm)')
vectors=('count')
radii=(3 4 5 6) 
models=("ECFP")
poly_representations=('Monomer' 'Dimer' 'Trimer' 'RRU Monomer' 'RRU Trimer' 'RRU Dimer')


for regressor in "${regressors[@]}"; do
  for target in "${targets[@]}"; do
    for fp in "${models[@]}"; do
      for oligo_rep in "${poly_representations[@]}"; do
        for radius in "${radii[@]}"; do
          for vector in "${vectors[@]}"; do
            bsub <<EOT
          
#BSUB -n 4
#BSUB -W 30:01
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=8GB]"
#BSUB -J "${fp}_${regressor}_${oligo_rep}_${target}_${radius}_${vector}_20250321"  
#BSUB -o "${output_dir}/${fp}_${regressor}_${oligo_rep}_${target}_${radius}_${vector}_20250321.out"
#BSUB -e "${output_dir}/${fp}_${regressor}_${oligo_rep}_${target}_${radius}_${vector}_20250321.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/pls-dataset-env
python ../train_structure_only.py --target_features "${target}" \
                                  --representation "${fp}" \
                                  --regressor_type "${regressor}" \
                                  --oligomer_representation "${oligo_rep}" \
                                  --radius "${radius}" \
                                  --vector "${vector}" 
EOT
          done
        done
      done
    done
  done
done

#'log First Peak (e-5 place holder)' 'log Second Peak (e-5 place holder)' 'log Third Peak (e-5 place holder)'
#'log First Peak wo placeholder' 'log Second Peak wo placeholder' 'log Third Peak wo placeholder'