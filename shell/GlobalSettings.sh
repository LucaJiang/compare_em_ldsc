# ------------ GLOBAL SETTINGS ------------ #
N_SNP=1217 # Set the number of SNPs
N_SAMPLE=1000 # Set the number of samples
repeat_times=50 # Set repeat times of generate data
h_list=($(seq 0.1 0.1 0.9)) # Set heritability list
r_list=(0.1 0.3 0.5 0.7 1) # Set Causal rate list
s_list=(1) # Set sigma_beta
# ! If enable job array, set to # of array 
# ! and add #SBATCH --array=1-3
SLURM_ARRAY_TASK_MAX=3

# Set data and paths
data_path=/home/wjiang49/scratch
data_name=height_ukb_50k_TNNI1
output_path=/home/wjiang49/scratch/UKBTNNI_rep50h9r5
code_path=/home/wjiang49/UKBheight
log_path=/home/wjiang49/scratch/UKBTNNI_rep50h9r5/log
data_path_name=${data_path}/${data_name}

# Set conda environment name
r_env="r4"
python_env="nsf"

# Random sleep to avoid the conflict 
sleep $((RANDOM % 10))

# Create the directory if it does not exist
if [ ! -d $output_path ]; then
    mkdir -p $output_path
fi

if [ ! -d $log_path ]; then
    mkdir -p $log_path
fi
# END
