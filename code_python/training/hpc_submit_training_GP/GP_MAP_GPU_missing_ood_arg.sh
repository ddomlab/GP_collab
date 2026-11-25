#!/bin/bash
set -euo pipefail

# Rerun the structure-cluster log Rg task for the NaN RUSC score in PPV.
# The training command evaluates all structure clusters for this target.

DATE=$(date +%Y%m%d)
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
model="GPytorchMAP"
paper="Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation"
dataset="Rg data with clusters aging imputed"
target_feature="log Rg (nm)"
fp_kernel="Matern32"
count_kernel="Matern32"
mixing_method="(count:x)+(fp:x)"
clustering_method="structure_cluster"
output_root="/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}"
output_dir="${output_root}/${paper}"
dataset_tag=${dataset// /_}
target_tag=${target_feature// /_}

mkdir -p "$output_dir"
bsub <<EOT
#BSUB -n 1
#BSUB -W 1:59
#BSUB -q short_gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -R "rusage[mem=32GB]"
#BSUB -R "select[a10 || a30 || a100 || l40 || h100]"
#BSUB -J "gpmap_structure_rusc_${DATE}_1"
#BSUB -o "${output_dir}/${model}_${dataset_tag}_${target_tag}_${fp_kernel}_${count_kernel}_${mixing_method}_${clustering_method}_GPU.out"
#BSUB -e "${output_dir}/${model}_${dataset_tag}_${target_tag}_${fp_kernel}_${count_kernel}_${mixing_method}_${clustering_method}_GPU.err"

source ~/.bashrc
module load cuda/12.1
module load gcc/9.3.0
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12

python "${script_dir}/../train_structure_numerical.py" \
    --K_fp "$fp_kernel" \
    --K_count "$count_kernel" \
    --Kernel_mixing_method "$mixing_method" \
    --paper "$paper" \
    --dataset "$dataset" \
    --target_features "$target_feature" \
    --regressor_type "$model" \
    --clustering_method "$clustering_method"
EOT

echo "Submitted 1 GPytorchMAP structure-cluster OOD GPU job."
