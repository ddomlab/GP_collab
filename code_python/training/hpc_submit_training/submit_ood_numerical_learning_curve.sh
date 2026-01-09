#!/bin/bash
output_dir=/share/ddomlab/sdehgha2/working-space/main/P1_pls-dataset/pls-dataset-space/PLS-Dataset/results/hpc_20250624
mkdir -p "$output_dir"

regressors=("RF" "XGBR")
targets=('log Rg (nm)')
group_out=('substructure cluster' 'KM4 polymer_solvent HSP cluster') 
# 'KM4 Mordred_Polysize cluster' 'KM4 ECFP6_Count_512bit cluster' 'KM4 Mordred_Polysize cluster'
# 'KM4 polymer_solvent HSP cluster' 'KM4 Mordred_Polysize cluster'
# 'KM4 ECFP6_Count_512bit cluster' 'HBD3 MACCS cluster' 'KM3 Mordred cluster' 'substructure cluster' 'KM4 polymer_solvent HSP and polysize cluster' 'KM5 polymer_solvent HSP and polysize cluster' 'KM4 polymer_solvent HSP cluster' 'KM4 Mordred_Polysize cluster'

for regressor in "${regressors[@]}"; do
  for target in "${targets[@]}"; do
      for group in "${group_out[@]}"; do
          bsub <<EOT



#BSUB -n 18
#BSUB -W 7:30
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=16GB]"
#BSUB -x
#BSUB -J "${regressor}_${target}_${fp}_${oligo_rep}_${group}_20250624_lc"  
#BSUB -o "${output_dir}/${regressor}_${target}_${fp}_${oligo_rep}_${group}_numerical_20250624_lc.out"
#BSUB -e "${output_dir}/${regressor}_${target}_${fp}_${oligo_rep}_${group}_numerical_20250624_lc.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/pls-dataset-env
python ../make_ood_learning_curve.py --target_features "${target}" \
                                      --regressor_type "${regressor}" \
                                      --numerical_feats 'Xn' 'Mw (g/mol)' 'PDI' 'Concentration (mg/ml)' 'Temperature SANS/SLS/DLS/SEC (K)' "polymer dP" "polymer dD" "polymer dH" 'solvent dP' 'solvent dD' 'solvent dH' "Dark/light" "Aging time (hour)" "To Aging Temperature (K)" "Sonication/Stirring/heating Temperature (K)" "Merged Stirring /sonication/heating time(min)" \
                                      --clustering_method "${group}" \



EOT
    done
  done
done

# 'Monomer' 'Dimer' 'Trimer' 'RRU Monomer' 'RRU Dimer'