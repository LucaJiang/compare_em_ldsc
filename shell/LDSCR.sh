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
source shell/GlobalSettings.sh
# Run ldsc.r algorithm to estimate the heritability
echo "Start running ldsc in R algorithm at $(date), Running h = ${h_list[*]}."

source activate $r_env
for h_value in ${h_list[@]}
do
    aim_folder=$output_path/h$h_value
    summary_data_list=$(ls $aim_folder/summary_data_h*)
    if [ -n "$summary_data_list" ]; then
        for summary_data in $summary_data_list
        do
            Rscript $code_path/ldsc.R -d $summary_data -N $N_SAMPLE -o $output_path >> $log_path/ldsc.log 2>&1
        done
    else
        echo "Summary data file dose not exist in ${aim_folder}!"
    fi
done
conda deactivate
echo "Finish running ldsc in R algorithm at $(date)."
# END