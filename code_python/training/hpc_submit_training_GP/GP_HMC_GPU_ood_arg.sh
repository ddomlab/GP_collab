#!/bin/bash

# Train the two missing continuous-cluster GpyroHMC OOD models with one
# independently scheduled job per target.

DATE=$(date +%Y%m%d)
model="GpyroHMC"
k_fps=("TanimotoMatern32")
k_counts=("Matern32")
k_mixing_methods=("(count:x)+(fp:x)")
clustering_methods=("continuous_cluster")
output_root="/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}"
job_index=0

dataset_targets=(
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|flux_data_imputed|log (Total flux)"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed|log (Separation factor)"
)

submit_job() {
    local paper="$1"
    local dataset="$2"
    local target_feature="$3"
    local fp_kernel="$4"
    local count_kernel="$5"
    local mixing_method="$6"
    local clustering_method="$7"
    local output_dir="${output_root}/${paper}"
    local dataset_tag=${dataset// /_}
    local target_tag=${target_feature// /_}
    target_tag=${target_tag//\//_}

    job_index=$((job_index + 1))
    mkdir -p "$output_dir"

    bsub <<EOT
#BSUB -n 1
#BSUB -W 60:10
#BSUB -q gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -R "rusage[mem=32GB]"
#BSUB -R "select[a10 || a30 || a100 || l40 || h100]"
#BSUB -J "gpyrohmc_ood_${DATE}_${job_index}"
#BSUB -o "${output_dir}/${model}_${dataset_tag}_${target_tag}_${fp_kernel}_${count_kernel}_${mixing_method}_${clustering_method}_GPU.out"
#BSUB -e "${output_dir}/${model}_${dataset_tag}_${target_tag}_${fp_kernel}_${count_kernel}_${mixing_method}_${clustering_method}_GPU.err"

source ~/.bashrc
module load cuda/12.1
module load gcc/9.3.0
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12

python ../train_structure_numerical.py --K_fp "$fp_kernel" \
                                      --K_count "$count_kernel" \
                                      --Kernel_mixing_method "$mixing_method" \
                                      --paper "$paper" \
                                      --dataset "$dataset" \
                                      --target_features "$target_feature" \
                                      --regressor_type "$model" \
                                      --clustering_method "$clustering_method"
EOT
}

for dataset_target in "${dataset_targets[@]}"; do
    IFS='|' read -r paper dataset target_feature <<< "$dataset_target"
    for fp_kernel in "${k_fps[@]}"; do
        for count_kernel in "${k_counts[@]}"; do
            for mixing_method in "${k_mixing_methods[@]}"; do
                for clustering_method in "${clustering_methods[@]}"; do
                    submit_job \
                        "$paper" \
                        "$dataset" \
                        "$target_feature" \
                        "$fp_kernel" \
                        "$count_kernel" \
                        "$mixing_method" \
                        "$clustering_method"
                done
            done
        done
    done
done

echo "Submitted ${job_index} GpyroHMC OOD GPU jobs."
