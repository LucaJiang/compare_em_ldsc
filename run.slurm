#!/bin/bash
#SBATCH -A pa_bios_department
#SBATCH -p special_bios
#SBATCH -J UKBHeight_gendata
#SBATCH -o /home/wjiang49/UKBheight/slurm_out_%j.log
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=512GB
#SBATCH -t 12:00:00


echo "Start time: $(date)"
echo "Job Name: $SLURM_JOB_NAME"

# -------------- PART ONE ------------------#
## Calculate ldscore
# source activate r4
# Rscript UKBheight/cal_ld.R 
# conda deactivate

# -------------- PART TWO ------------------#
## Generate data
# # Set the parameters
# repeat_times=10

# h_list=($(seq 0.1 0.1 0.9))
# # ! Set the ntasks-per-node=2
# # h_list1=($(seq 0.6 0.1 0.9))
# # h_list2=($(seq 0.1 0.1 0.5))
# # if [ "$SLURM_ARRAY_TASK_ID" -eq "1" ]
# # then
# #     h_list=("${h_list1[@]}")
# # else
# #     h_list=("${h_list2[@]}")
# # fi

# r_list=(1e-3 1e-2 1e-1 5e-1)
# s_list=(0.1 0.5 1 2 5 10)

# TEST
repeat_times=1
h_list=($(seq 0.1 0.1 0.2))
r_list=(1e-3 1e-2)
s_list=(0.1 0.5)

# Set data path
data_path=/home/wjiang49/scratch/height_ukb_50k_chr22
output_path=/home/wjiang49/scratch/UKBsimdata
code_path=/home/wjiang49/UKBheight
log_path=/home/wjiang49/scratch

# # echo "Job: generate data & em. Repeat $repeat_times."

# Create the directory if it does not exist
if [ ! -d $output_path ]; then
    mkdir $output_path
fi

source activate r4
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
                Rscript $code_path/gen_simul_data.R -d $data_path -o $output_path -h $h_value -r $r_value -s $s_value >> $log_path/gendata.log 2>&1
            done
        done
    done
done
conda deactivate


# -------------- PART THREE ------------------#
# # Run the EM algorithm
# source activate nsf
# for h_value in ${h_list[@]}
# do
#     aim_folder=$output_path/h$h_value
#     if ls $aim_folder/phenotypeh* 1> /dev/null 2>&1; then
#         python3 $code_path/lmm_em.py -pd $aim_folder -g $output_path/genotypes.csv -o $output_path >> $log_path/em.log 2>&1
#     else
#         echo "Phenotype file dose not exist!"
#     fi
# done
# conda deactivate


# -------------- PART FOUR ------------------#
## RUN the Ldscore regression
# source activate r4
# for h_value in ${h_list[@]}
# do
#     aim_folder=$output_path/h$h_value
#     if ls $aim_folder/summary_data_h* 1> /dev/null 2>&1; then
#         Rscript $code_path/runLDSC.r -d $aim_folder -o $output_path >> $log_path/ldsc.log 2>&1
#     else
#         echo "Summary data file dose not exist in $(aim_folder)!"
#     fi
# done
# conda deactivate

# -------------- Finally ------------------#
echo "Finish all the jobs at $(date)."
