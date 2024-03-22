#!/bin/bash
#SBATCH -A pa_bios_department
#SBATCH -p special_bios
#SBATCH -J UKBHeight_unnamed
#SBATCH -o /home/wjiang49/UKBheight/log/slurm_out_%A%a.log
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH -t 12:00:00
#SBATCH --array=1-3

# Get the global settings
source GlobalSettings.sh
# Calculate task load for each task array
if [ -n "$SLURM_ARRAY_TASK_ID" ]
then
    h_per_task=$(((${#h_list[@]} + $SLURM_ARRAY_TASK_MAX - 1) / $SLURM_ARRAY_TASK_MAX ))
    start_index=$((($SLURM_ARRAY_TASK_ID - 1) * $h_per_task))
    end_index=$(($start_index + $h_per_task - 1))
    # Ensure end_index is not greater than the length of h_list
    if [ $end_index -ge ${#h_list[@]} ]
    then
        end_index=$((${#h_list[@]} - 1))
    fi
    h_list=("${h_list[@]:$start_index:$((end_index - start_index + 1))}")
fi
# Run EM algorithm
echo "Start running EM algorithm at $(date), Running h = ${h_list[*]}."
source activate $python_env
for h_value in ${h_list[@]}
do
    aim_folder=$output_path/h$h_value
    if ls $aim_folder/phenotypeh* 1> /dev/null 2>&1; then
        python3 $code_path/lmm_em.py -pd $aim_folder -gd $output_path/genotypes.csv -o $output_path >> $log_path/em.log 2>&1
    else
        echo "Phenotype file dose not exist in ${aim_folder}!"
    fi
done
conda deactivate
echo "Finish running EM algorithm at $(date)."
# END