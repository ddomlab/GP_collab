#!/bin/bash

# Submit the GPytorchMAP GPU configurations currently missing from the
# non-test research datasets. These 29 jobs produce 34 target-level results.

DATE=$(date +%Y%m%d)
model="GPytorchMAP"
output_root="/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}"
job_index=0

submit_job() {
    local paper="$1"
    local dataset="$2"
    local fp_kernel="$3"
    local count_kernel="$4"
    local mixing_method="$5"
    local output_dir="${output_root}/${paper}"

    job_index=$((job_index + 1))
    mkdir -p "$output_dir"

    bsub <<EOT
#BSUB -n 1
#BSUB -W 20:20
#BSUB -q gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -R "rusage[mem=8GB]"
#BSUB -R "select[a10 || a30 || a100 || l40 || h100]"
#BSUB -J "gpmap_missing_${DATE}_${job_index}"
#BSUB -o "${output_dir}/${model}_${dataset}_${fp_kernel}_${count_kernel}_${mixing_method}_GPU.out"
#BSUB -e "${output_dir}/${model}_${dataset}_${fp_kernel}_${count_kernel}_${mixing_method}_GPU.err"

source ~/.bashrc
module load cuda/12.1
module load gcc/9.3.0
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12

python ../train_structure_numerical.py --K_fp "$fp_kernel" \
                                        --K_count "$count_kernel" \
                                        --Kernel_mixing_method "$mixing_method" \
                                        --paper "$paper" \
                                        --dataset "$dataset" \
                                        --regressor_type "$model"
EOT
}

# calculated PCE (%) -- 558 datapoints
paper="Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices"
dataset="Beyond molecular structure_seifrid_imputed"

submit_job "$paper" "$dataset" "Tanimoto"          "Matern32" "(count:+)x(fp:+)"
submit_job "$paper" "$dataset" "TanimotoRBF"       "RBF"      "(count:+)x(fp:+)"
submit_job "$paper" "$dataset" "TanimotoRBF"       "RBF"      "(count:+)x(fp:x)"
submit_job "$paper" "$dataset" "TanimotoMatern32" "RBF"      "(count:+)x(fp:+)"
submit_job "$paper" "$dataset" "TanimotoMatern32" "RBF"      "(count:x)+(fp:x)"
submit_job "$paper" "$dataset" "TanimotoMatern52" "RBF"      "(count:x)+(fp:x)"
submit_job "$paper" "$dataset" "Matern32"         "Matern52" "(count:+)x(fp:+)"
submit_job "$paper" "$dataset" "Matern32"         "RBF"      "(count:+)x(fp:x)"
submit_job "$paper" "$dataset" "Matern52"         "Matern32" "(count:x)+(fp:x)"
submit_job "$paper" "$dataset" "Matern52"         "Matern32" "(count:+)x(fp:x)"
submit_job "$paper" "$dataset" "Matern52"         "Matern52" "(count:x)+(fp:x)"
submit_job "$paper" "$dataset" "RBF"              "Matern52" "(count:x)+(fp:x)"
submit_job "$paper" "$dataset" "RBF"              "RBF"      "(count:+)x(fp:+)"
submit_job "$paper" "$dataset" "RBF"              "RBF"      "(count:+)x(fp:x)"

# log Rg (nm) -- 256 datapoints
paper="Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation"
dataset="Rg data with clusters aging imputed"

submit_job "$paper" "$dataset" "TanimotoMatern52" "Matern52" "averageProduct"
submit_job "$paper" "$dataset" "TanimotoMatern52" "RBF"      "(count:+)x(fp:+)"
submit_job "$paper" "$dataset" "Matern52"         "Matern32" "(count:x)+(fp:x)"
submit_job "$paper" "$dataset" "Matern52"         "Matern32" "(count:+)x(fp:x)"
submit_job "$paper" "$dataset" "Matern52"         "RBF"      "(count:x)+(fp:x)"
submit_job "$paper" "$dataset" "RBF"              "RBF"      "(count:+)x(fp:x)"

# log (Separation factor); the Total flux configurations are complete.
paper="Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery"
dataset="separation_data_imputed"

submit_job "$paper" "$dataset" "TanimotoMatern32" "RBF" "sum"

# Each call trains both log (Separation factor) and log (Total flux), so these
# five jobs produce ten missing target-level results.
paper="Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation"
dataset="cleaned_dataset_pervaporation_membranes_wang"

submit_job "$paper" "$dataset" "Tanimoto"          "RBF"      "(count:+)x(fp:+)"
submit_job "$paper" "$dataset" "TanimotoRBF"       "Matern52" "(count:+)x(fp:+)"
submit_job "$paper" "$dataset" "TanimotoMatern32" "RBF"      "(count:+)x(fp:+)"
submit_job "$paper" "$dataset" "Matern32"         "Matern52" "averageProduct"
submit_job "$paper" "$dataset" "Matern52"         "RBF"      "averageProduct"

# Approx Conv (%) -- 768 datapoints
paper="Miniaturization of Popular Reactions from the Medicinal Chemists Toolbox for Ultrahigh_Throughput Experimentation"
dataset="cleaned_suzuki_synthesis"

submit_job "$paper" "$dataset" "TanimotoRBF" "Matern32" "(count:+)x(fp:x)"
submit_job "$paper" "$dataset" "TanimotoRBF" "Matern52" "(count:+)x(fp:x)"
submit_job "$paper" "$dataset" "TanimotoRBF" "RBF"      "(count:+)x(fp:x)"

echo "Submitted ${job_index} GPytorchMAP GPU jobs."
