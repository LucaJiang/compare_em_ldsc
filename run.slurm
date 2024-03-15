#!/bin/bash
#SBATCH -A pa_bios_department
#SBATCH -p special_bios
#SBATCH -J UKBHeight_unnamed
#SBATCH -o /home/wjiang49/UKBheight/log/slurm_out_%j.log
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB
#SBATCH -t 12:00:00


# ------------- GLOBAL SETTINGS ------------#
N_SNP=40000 # Set the number of SNPs
N_SAMPLE=2000 # Set the number of samples
repeat_times=10 # Set repeat times of generate data
h_list=($(seq 0.1 0.1 0.9)) # Set heritability list
r_list=(1e-3 1e-2 5e-1) # Set Causal rate list
s_list=(0.1 1 10) # Set sigma_beta
# ! If enable job array, set to # of array 
# ! and add #SBATCH --array=1-3
# SLURM_ARRAY_TASK_MAX=1

# Set data and paths
data_path_name=/home/wjiang49/scratch/height_ukb_50k_chr22
output_path=/home/wjiang49/scratch/UKBsimdata
code_path=/home/wjiang49/UKBheight
log_path=/home/wjiang49/scratch/UKBsimdata/log

# Set conda environment name
r_env="r4"
python_env="nsf"

# ------------ Initialization -------------- #
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

# sleep random time to avoid the conflict
sleep $((RANDOM % 10))
echo "$SLURM_JOB_NAME: Start time: $(date), Running h = ${h_list[*]}."

# Create the directory if it does not exist
if [ ! -d $output_path ]; then
    mkdir -p $output_path
fi

if [ ! -d $log_path ]; then
    mkdir -p $log_path
fi

# -------------- PART ONE ----------------- #
#            Calculate ldscore
# ----------------------------------------- #
# # ! Run this part only once, when the ldscore is not calculated
# echo "Start calculating ldscore."
# source activate $r_env
# Rscript UKBheight/cal_ld.R -d $data_path_name -N $N_SAMPLE -P $N_SNP -o $output_path >> $log_path/ldscore.log 2>&1
# conda deactivate

# -------------- PART TWO ------------------#
#             Generate data
# ----------------------------------------- #
# echo "Start generating data."
# source activate $r_env
# for i in $(seq 1 $repeat_times)
# do
#     echo "This is loop $i at $(date)."
#     for h_value in ${h_list[@]}
#     do
#         if [ ! -d $output_path/h$h_value ]; then
#             mkdir $output_path/h$h_value
#         fi
#         for r_value in ${r_list[@]}
#         do
#             for s_value in ${s_list[@]}
#             do
#                 Rscript $code_path/gen_simul_data.R -d $data_path_name -N $N_SAMPLE -P $N_SNP -o $output_path -H $h_value -r $r_value -s $s_value >> $log_path/gendata.log 2>&1
#             done
#         done
#     done
# done
# conda deactivate

# -------------- PART THREE --------------- #
#            Run the EM algorithm
# ----------------------------------------- #
# echo "Start running EM algorithm."
# source activate $python_env
# for h_value in ${h_list[@]}
# do
#     aim_folder=$output_path/h$h_value
#     if ls $aim_folder/phenotypeh* 1> /dev/null 2>&1; then
#         python3 $code_path/lmm_em.py -pd $aim_folder -gd $output_path/genotypes.csv -o $output_path >> $log_path/em.log 2>&1
#     else
#         echo "Phenotype file dose not exist in ${aim_folder}!"
#     fi
# done
# conda deactivate

# -------------- PART FOUR ---------------- #
#        RUN the ldscore regression
# ----------------------------------------- #
# echo "Start running ldscore regression."
# source activate $r_env
# for h_value in ${h_list[@]}
# do
#     aim_folder=$output_path/h$h_value
#     summary_data_list=$(ls $aim_folder/summary_data_h*)
#     if [ -n "$summary_data_list" ]; then
#         for summary_data in $summary_data_list
#         do
#             Rscript $code_path/ldsc.R -d $summary_data -N $N_SAMPLE -o $output_path >> $log_path/ldsc.log 2>&1
#         done
#     else
#         echo "Summary data file dose not exist in ${aim_folder}!"
#     fi
# done
# conda deactivate

# -------------- PART FIVE ---------------- #
#         visualizing the results
# ----------------------------------------- #
# # ! Run this part only once
# echo "Start visualizing the results."
# source activate $python_env
# python3 $code_path/visualize.py -d $output_path -o $output_path >> $log_path/visualize.log 2>&1
# conda deactivate

# -------------- Finally ------------------#
echo "Finish all the jobs at $(date)."
