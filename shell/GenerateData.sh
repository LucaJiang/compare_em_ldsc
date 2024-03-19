#!/bin/bash
#SBATCH -A pa_bios_department
#SBATCH -p special_bios
#SBATCH -J UKBHeight_gendata
#SBATCH -o /home/wjiang49/UKBheight/log/slurm_out_%A%a.log
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB
#SBATCH -t 6:00:00
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
# Generate Simulation Data
echo "Generate Simulation Data: Start time: $(date), Running h = ${h_list[*]}."
source activate $r_env
for i in $(seq 1 $repeat_times)
do
    echo "This is loop $i at $(date)."
    for h_value in ${h_list[@]}
    do
        if [ ! -d $output_path/h$h_value ]; then
            mkdir $output_path/h$h_value
        fi
        for r_value in ${r_list[@]}
        do
            for s_value in ${s_list[@]}
            do
                Rscript $code_path/gen_simul_data.R -d $data_path_name -N $N_SAMPLE -P $N_SNP -o $output_path -H $h_value -r $r_value -s $s_value >> $log_path/gendata.log 2>&1
            done
        done
    done
done
conda deactivate
echo "Finish generating data at $(date)."