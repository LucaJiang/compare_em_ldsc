#!/bin/bash
#SBATCH -A pa_bios_department
#SBATCH -p special_bios
#SBATCH -J UKBHeight_calculate_ldscore
#SBATCH -o ./log/slurm_out_%j.log
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128GB
#SBATCH -t 2:00:00

# Get the global settings
source GlobalSettings.sh
# Create the directory if it does not exist
if [ ! -d $output_path ]; then
    mkdir -p $output_path
fi

if [ ! -d $log_path ]; then
    mkdir -p $log_path
fi
# Calculate ldscore
echo "Start calculating ldscore at $(date)."
source activate $r_env
Rscript UKBheight/cal_ld.R -d $data_path_name -N $N_SAMPLE -P $N_SNP -o $data_path_name >> $log_path/cal_ld.log 2>&1
conda deactivate
echo "Finish calculating ldscore at $(date)."
# END