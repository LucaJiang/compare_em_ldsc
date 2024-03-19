#!/bin/bash
#SBATCH -A pa_bios_department
#SBATCH -p special_bios
#SBATCH -J UKBHeight_unnamed
#SBATCH -o /home/wjiang49/UKBheight/log/slurm_out_%j.log
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH -t 10:00

# Get the global settings
source GlobalSettings.sh
echo "Start visualizing the results at $(date)."
source activate $python_env
python3 $code_path/visualize.py -d $output_path -o $output_path >> $log_path/visualize.log 2>&1
conda deactivate
echo "Finish visualizing the results at $(date)."
# END