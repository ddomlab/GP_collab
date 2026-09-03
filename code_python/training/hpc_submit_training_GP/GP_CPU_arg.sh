#!/bin/bash

DATE=$(date +%Y%m%d)
model="GPytorchMAPRegressor"
paper="Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation"
dataset="Rg data with clusters aging imputed"
mixing_method="product"

##flux_data_imputed
output_dir=/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}/${paper}
mkdir -p "$output_dir"

sbatch <<EOT
#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --job-name="structure_numerical_${DATE}"
#SBATCH --output="${output_dir}/MAP_${mixing_method}_CPU.out"
#SBATCH --error="${output_dir}/MAP_${mixing_method}_CPU.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12
python ../train_structure_numerical.py --dataset "$dataset" \
                                        --paper "$paper" \
                                        --Kernel_mixing_method "$mixing_method"

EOT
