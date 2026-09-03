#!/bin/bash

DATE=$(date +%Y%m%d)
# model="GPytorchMAPRegressor"
paper="Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery"
dataset="flux_data_imputed"
models=("RF")
# k_fps=("TanimotoRBF")
# k_counts=("RBF")
# k_mixing_methods=("sum") 

##flux_data_imputed
output_dir=/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}/${paper}
mkdir -p "$output_dir"

for model in "${models[@]}"; do
    sbatch <<EOT
#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1
#SBATCH --time=01:30:00
#SBATCH --mem=32G
#SBATCH --job-name="structure_numerical_${DATE}"
#SBATCH --output="${output_dir}/${model}_${fp_kernel}_${count_kernel}_${mixing_method}.out"
#SBATCH --error="${output_dir}/${model}_${fp_kernel}_${count_kernel}_${mixing_method}.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12
python ../train_structure_numerical.py --regressor_type $model \
                                        --dataset "$dataset" \
                                        --paper "$paper"


EOT
done
