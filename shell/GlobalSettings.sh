# ------------ GLOBAL SETTINGS ------------ #
N_SNP=1217 # Set the number of SNPs
N_SAMPLE=1000 # Set the number of samples
repeat_times=50 # Set repeat times of generate data
h_list=($(seq 0.1 0.1 0.9)) # Set heritability list
r_list=(0.1, 0.3, 0.5, 0.7, 1) # Set Causal rate list
s_list=(1) # Set sigma_beta
# ! If enable job array, set to # of array 
# ! and add #SBATCH --array=1-3
SLURM_ARRAY_TASK_MAX=3

# Set data and paths
data_path_name=/home/wjiang49/scratch/height_ukb_50k_TNNI1
output_path=/home/wjiang49/scratch/UKBTNNI
code_path=/home/wjiang49/UKBheight
log_path=/home/wjiang49/scratch/UKBTNNI/log

# Set conda environment name
r_env="r4"
python_env="nsf"
# END