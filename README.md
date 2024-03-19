# Compare the Result of EM and LDSC

## Introduction

...

## Structure

* [Calculate LDScore](cal_ld.R): Calculate the LD score
* [Generate_Simulation_Data](gen_simul_data.R): Generate the simulation data
* [Run EM](lmm_em.py): Run EM algorithm
* [Run LDSC](irwls.py): Run LDSC algorithm
* [Visualize](visualize.py): Visualize the result
* [Run All in Server](run.slurm): Run all the above in server

## Usage

Set parameters in the GLOBAL SETTINGS part of `run.slurm` 

The Usage of each script is described in the header of each script.

### UKBTNNI Data

1. Upload the UKBTNNI plink data (.bed, .bim, .fam) to server and change the setting in [GlobalSettings.sh](shell/GlobalSettings.sh) accordingly.

2. cd to the directory of the scripts, for example:

```{bash}
cd /home/wjiang49/UKBheight/
```

The log file will be saved in the above directory.

3. Calculate LD score and record the job id for the following steps:

```{bash}
calLd_id=${sbatch --parsable shell/CalculateLDScore.sh}
```

This will calculate the LD score and save the result in the same directory as the plink data.

4. Generate the simulation data:

```{bash}
gendata_id=${sbatch --parsable --dependency=afterok:${calLd_id} shell/GenerateData.sh}
```

5.1 Run EM:

```{bash}
em_id=${sbatch --parsable --dependency=afterok:${gendata_id} shell/RunEM.sh}
```

5.2 Run LDSC:

```{bash}
ldsc_id=${sbatch --parsable --dependency=afterok:${gendata_id} shell/RunLDSC.sh}
```

6. Visualize the result:

```{bash}
sbatch --dependency=afterok:${em_id}:${ldsc_id} shell/Visualize.sh
```

Here's a script to run all the above:

```{bash}
cd /home/wjiang49/UKBheight/
calLd_id=${sbatch --parsable shell/CalculateLDScore.sh}
gendata_id=${sbatch --parsable --dependency=afterok:${calLd_id} shell/GenerateData.sh}
em_id=${sbatch --parsable --dependency=afterok:${gendata_id} shell/RunEM.sh}
ldsc_id=${sbatch --parsable --dependency=afterok:${gendata_id} shell/RunLDSC.sh}
sbatch --dependency=afterok:${em_id}:${ldsc_id} shell/Visualize.sh
```
