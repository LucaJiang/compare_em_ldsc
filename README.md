# Compare the Result of EM and LDSC

- [Compare the Result of EM and LDSC](#compare-the-result-of-em-and-ldsc)
  - [Introduction](#introduction)
  - [Structure](#structure)
  - [Usage](#usage)
    - [UKBTNNI Data: Run All](#ukbtnni-data-run-all)
    - [UKBTNNI Data: Step by Step](#ukbtnni-data-step-by-step)

## Introduction

This repository provides a pipeline to compare the result of EM algorithm and LDSC algorithm in the context of GWAS. The input data is the UK Biobank data for height and the output is the comparison of the heritability estimation by EM and LDSC. The pipeline includes the following steps:

## Structure

* [Calculate LDScore](cal_ld.R): Calculate the LD score
* [Generate_Simulation_Data](gen_simul_data.R): Generate the simulation data
* [Run EM](lmm_em.py): Run EM algorithm
* [Run LDSC in py](irwls.py): Run LDSC algorithm in python
  * [Run LDSC in R](ldsc.R): Run LDSC algorithm in R
* [Visualize](visualize.py): Visualize the result
* [Run All in Server](run.slurm): Run all the above in server
* [shell](shell): Shell scripts to run the above sequentially (or separately) in server

## Usage

The Usage of each script is described in the header of each script. You can set parameters in the GLOBAL SETTINGS part of `run.slurm` to run the pipeline. Or follow the following steps to run the pipeline step by step.

### UKBTNNI Data: Run All

Here's a script to run all the above. Remember to change the settings in [GlobalSettings.sh](shell/GlobalSettings.sh) accordingly.

`cd` to the directory of the code scripts, then run the following script:

```{bash}
calLd_id=$(sbatch --parsable shell/CalculateLDScore.sh)
gendata_id=$(sbatch --parsable --dependency=afterok:$calLd_id shell/GenerateData.sh)
em_id=$(sbatch --parsable --dependency=afterok:$gendata_id shell/EM.sh)
ldsc_id=$(sbatch --parsable --dependency=afterok:$gendata_id shell/LDSC.sh)
sbatch --dependency=afterok:$em_id:$ldsc_id shell/Visual.sh
```

### UKBTNNI Data: Step by Step

1. Upload the UKBTNNI plink data (.bed, .bim, .fam) to server and change the setting in [GlobalSettings.sh](shell/GlobalSettings.sh) accordingly.

The structure of the files should be look like this:

```{bash}
> tree ${data_path} -L 1
.
├── height_ukb_50k_TNNI1.bed
├── height_ukb_50k_TNNI1.bim
├── height_ukb_50k_TNNI1.fam
├── ...
```

2. cd to the directory of the code scripts, for example:

```{bash}
cd /home/wjiang49/UKBheight/
```

The slurm log file will be saved in the above directory like this:

```{bash}
log
├── slurm_out_614489.log
├── slurm_out_614490.log
├── slurm_out_614491.log
├── slurm_out_615299_1.log
├── slurm_out_615299_2.log
├── slurm_out_615299_3.log
```

3. Calculate LD score and record the job id for the following steps:

```{bash}
calLd_id=$(sbatch --parsable shell/CalculateLDScore.sh)
```

The .ldscore file will be saved in the same directory as the plink data.

This will calculate the LD score and save the result in the same directory as the plink data.

4. Generate the simulation data:

```{bash}
gendata_id=${sbatch --parsable --dependency=afterok:${calLd_id} shell/GenerateData.sh}
```

The results will be saved in the `${output_path}` directory. If such directory does not exist, it will be created. The generated data shall look like this:

```{bash}
.
├── genotypes.csv
├── h0.1
...
├── h0.9
└── log
```

The summary data and phenotype data will be saved in the `h0.*` directory. The log of all the process will be saved in the `log` directory.

5. Run EM and LDSC. The job will begin after the generation of the simulation data. The log will be saved in the `${output_path}/log` directory and results will be saved in the `${output_path}` directory.

i. Run EM:

```{bash}
em_id=${sbatch --parsable --dependency=afterok:${gendata_id} shell/EM.sh}
```

ii. Run LDSC:

```{bash}
ldsc_id=${sbatch --parsable --dependency=afterok:${gendata_id} shell/LDSC.sh}
```

After that, the `${output_path}` directory should look like this:

```{bash}
.
├── em_results.csv
├── genotypes.csv
├── h0.1
...
├── h0.9
├── ldsc_results.csv
└── log
```

6. Visualize the result:

```{bash}
sbatch --dependency=afterok:${em_id}:${ldsc_id} shell/Visual.sh
```
Also, the figures will be saved in the `${output_path}` directory.

END