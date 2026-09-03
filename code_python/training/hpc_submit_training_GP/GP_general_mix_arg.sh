#!/bin/bash

DATE=$(date +%Y%m%d)
# model="GPMixMCMC"
paper="Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery"
dataset="flux_data_imputed"
##flux_data_imputed
output_dir=/share/ddomlab/sdehgha2/working-space/main/_colab_GP/GP_collab/results/HPC_history/hpc_${DATE}/${paper}
mkdir -p "$output_dir"

k_fps=("TanimotoRBF" "TanimotoMatern32" "TanimotoMatern52" "Tanimoto" "RBF" "Matern32" "Matern52")
k_counts=("RBF" "Matern32" "Matern52")
k_mixing_methods=("sum" "product" "(count:+)x(fp:x)" "(count:+)x(fp:+)" "(count:x)+(fp:x)") 

for mixing_method in "${k_mixing_methods[@]}"; do
    for fp_kernel in "${k_fps[@]}"; do
        for count_kernel in "${k_counts[@]}"; do
            sbatch <<EOT
#!/bin/bash



#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --job-name="structure_numerical_${DATE}"
#SBATCH --output="${output_dir}/GPytorchMAP_${fp_kernel}_${count_kernel}_${mixing_method}.out"
#SBATCH --error="${output_dir}/GPytorchMAP_${fp_kernel}_${count_kernel}_${mixing_method}.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/torch_cpu
python ../train_structure_numerical.py --K_fp $fp_kernel \
                                        --K_count $count_kernel \
                                        --Kernel_mixing_method "$mixing_method" \
                                        --paper "$paper" \
                                        --dataset "$dataset" \

EOT

        done
    done
done
