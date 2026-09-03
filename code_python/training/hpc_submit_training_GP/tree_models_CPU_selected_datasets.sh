#!/bin/bash

# Submit RF, XGBR, and NGB only for the explicitly selected datasets below.
# No datasets are discovered automatically. For each selected dataset,
# train_structure_numerical.py trains its target(s) configured in filter_data.py.

DATE=$(date +%Y%m%d)
models=("RF" "XGBR" "NGB")

# Selected datasets and their configured targets:
# - Beyond molecular structure_seifrid_imputed: calculated PCE (%)
# - flux_data_imputed: log (Total flux)
# - separation_data_imputed: log (Separation factor)
# - cleaned_dataset_pervaporation_membranes_wang: log (Total flux),
#   log (Separation factor)
# - cleaned_suzuki_synthesis: Approx Conv (%)
# - Rg data with clusters aging imputed: log Rg (nm)
# - cleaned_dataset_Ultrafiltration Membrane_imputed: all six configured
#   ultrafiltration targets
selected_training_sets=(
    "Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic materials and devices|Beyond molecular structure_seifrid_imputed"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|flux_data_imputed"
    "Machine Learning for Polymer Design to Enhance Pervaporation-Based Organic Recovery|separation_data_imputed"
    "Machine Learning-Enabled Prediction and High-Throughput Screening of Polymer Membranes for Pervaporation Separation|cleaned_dataset_pervaporation_membranes_wang"
    "Miniaturization of Popular Reactions from the Medicinal Chemists Toolbox for Ultrahigh_Throughput Experimentation|cleaned_suzuki_synthesis"
    "Robust Learning from Literature Data_Model Generalizability and Uncertainty for Predicting Conjugated Polymer Solution Conformation|Rg data with clusters aging imputed"
    "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning|cleaned_dataset_Ultrafiltration Membrane_imputed"
)

output_root="/share/ddomlab/sdehgha2/working_space/GP_collab/results/HPC_history/hpc_${DATE}"
job_index=0

for training_set in "${selected_training_sets[@]}"; do
    IFS='|' read -r paper dataset <<< "$training_set"
    output_dir="${output_root}/${paper}"
    dataset_tag=${dataset// /_}
    mkdir -p "$output_dir"

    for model in "${models[@]}"; do
        job_index=$((job_index + 1))

        sbatch <<EOT
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --job-name="tree_selected_${DATE}_${job_index}"
#SBATCH --output="${output_dir}/${model}_${dataset_tag}_CPU.out"
#SBATCH --error="${output_dir}/${model}_${dataset_tag}_CPU.err"

source ~/.bashrc
conda activate /usr/local/usrapps/ddomlab/sdehgha2/env12

python ../train_structure_numerical.py --regressor_type "$model" \
                                        --dataset "$dataset" \
                                        --paper "$paper"
EOT
    done
done
