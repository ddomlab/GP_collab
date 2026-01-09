#!/bin/bash
output_dir=/share/ddomlab/sdehgha2/working-space/main/P1_pls-dataset/pls-dataset-space/PLS-Dataset/results/hpc_20250531
mkdir -p "$output_dir"

regressors=("XGBR")
targets=('log Rg (nm)')
models=("ECFP")
radii=(3) 
vectors=("count")
poly_representations=('Trimer')
group_out=('substructure cluster' 'KM4 polymer_solvent HSP cluster') 

for regressor in "${regressors[@]}"; do
  for target in "${targets[@]}"; do
    for fp in "${models[@]}"; do
      for oligo_rep in "${poly_representations[@]}"; do
        for radius in "${radii[@]}"; do
          for vector in "${vectors[@]}"; do
            for group in "${group_out[@]}"; do
              bsub <<EOT



#BSUB -n 6
#BSUB -W 10:05
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=16GB]"
#BSUB -J "${regressor}_${target}_${fp}_${oligo_rep}_${group}_20250531"  
#BSUB -o "${output_dir}/${regressor}_${target}_${fp}_${oligo_rep}_${radius}_${vector}_${group}_20250531.out"
#BSUB -e "${output_dir}/${regressor}_${target}_${fp}_${oligo_rep}_${radius}_${vector}_${group}_20250531.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/pls-dataset-env
python ../make_ood_learning_curve.py --target_features "${target}" \
                                      --representation "${fp}" \
                                      --regressor_type "${regressor}" \
                                      --radius "${radius}" \
                                      --vector "${vector}" \
                                      --oligomer_representation "${oligo_rep}" \
                                      --numerical_feats 'Xn' 'Mw (g/mol)' 'PDI' \
                                      --clustering_method "${group}" \



EOT
            done
          done
        done
      done
    done
  done
done

# 'Monomer' 'Dimer' 'Trimer' 'RRU Monomer' 'RRU Dimer'