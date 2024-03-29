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
source shell/GlobalSettings.sh
# Calculate ldscore
echo "Start calculating ldscore at $(date)."
source activate $r_env
Rscript $code_path/cal_ld.R -d $data_path_name -N $N_SAMPLE -P $N_SNP >> $log_path/cal_ld.log 2>&1
conda deactivate
echo "Finish calculating ldscore at $(date)."
# END