#!/bin/bash
#SBATCH -A pa_bios_department
#SBATCH -p special_bios
#SBATCH -J UKBHeight_unnamed
#SBATCH -o /home/wjiang49/UKBheight/log/slurm_out_%j.log
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH -t 1:00:00

# Get the global settings
source GlobalSettings.sh
# Run irwls algorithm to estimate the heritability
echo "Start running irwls algorithm at $(date), Running h = ${h_list[*]}."
source activate $python_env
for h_value in ${h_list[@]}
do
    aim_folder=$output_path/h$h_value
    if ls $aim_folder/summary_data_h* 1> /dev/null 2>&1; then
        python3 $code_path/irwls.py -d $aim_folder -o $output_path -N $N_SAMPLE >> $log_path/irwls.log 2>&1
    else
        echo "Summary data file dose not exist in ${aim_folder}!"
    fi
done
conda deactivate
echo "Finish running irwls algorithm at $(date)."
# END