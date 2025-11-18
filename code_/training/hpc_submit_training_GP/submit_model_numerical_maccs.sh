#!/bin/bash
output_dir=/share/ddomlab/sdehgha2/working-space/main/P1_pls-dataset/pls-dataset-space/PLS-Dataset/results

models_to_run=("RF" "MLR")
functions_to_run=("numerical" "numerical_maccs")

for model in "${models_to_run[@]}"; do
    for function in "${functions_to_run[@]}"; do
        bsub <<EOT




#BSUB -n 8
#BSUB -W 30:01
#BSUB -R span[ptile=4]
##BSUB -x
#BSUB -R "rusage[mem=32GB]"
#BSUB -J finger_numerical_maccs
#BSUB -o ${output_dir}/mordred_run_${model}_${function}.out
#BSUB -e ${output_dir}/mordred_err_${model}_${function}.err


source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/pls-dataset-env
python train_structure_numeric.py --model $model --function $function

conda deactivate

EOT
    done
done