#!/bin/bash

DATE=$(date +%Y%m%d)
model="GpyroHMC"
paper="Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices"
dataset="Beyond molecular structure_seifrid_imputed"
k_fps=("TanimotoMatern32")
k_counts=("Matern32")           
k_mixing_methods=("averageProduct") 
##flux_data_imputed
output_dir=/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}/${paper}
mkdir -p "$output_dir"

for mixing_method in "${k_mixing_methods[@]}"; do
    for fp_kernel in "${k_fps[@]}"; do
        for count_kernel in "${k_counts[@]}"; do
            sbatch <<EOT
#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --time=2-12:20:00
#SBATCH --mem=32G
#SBATCH --job-name="structure_numerical_${DATE}"
#SBATCH --output="${output_dir}/${model}_${fp_kernel}_${count_kernel}_${mixing_method}.out"
#SBATCH --error="${output_dir}/${model}_${fp_kernel}_${count_kernel}_${mixing_method}.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12
python ../train_structure_numerical.py --K_fp $fp_kernel \
                                        --K_count $count_kernel \
                                        --Kernel_mixing_method "$mixing_method" \
                                        --dataset "$dataset" \
                                        --paper "$paper" \
                                        --regressor_type "$model" 


EOT

        done
    done
done
