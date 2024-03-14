# Codes for generate simulated data with command line arguments
# Usage:
# If doesn't have genotype data:
# Rscript /home/wjiang49/UKBheight/gen_simul_data.R -d /home/wjiang49/scratch/height_ukb_50k_chr22 -o /home/wjiang49/scratch/UKBsimdata/ -h 0.6 -r 0.01 -s 2 >> /home/wjiang49/scratch/gendata.log 2>&1
# If already have genotype data:
# Rscript /home/wjiang49/UKBheight/gen_simul_data.R -d /home/wjiang49/scratch/height_ukb_50k_chr22 -g /home/wjiang49/scratch/UKBsimdata/genotypes.csv -o /home/wjiang49/scratch/UKBsimdata/ -h 0.6 -r 0.01 -s 2 >> /home/wjiang49/scratch/gendata.log 2>&1

# File structure of the generated data:
# /BASE/genotypes.csv
# /BASE/h${heritability}/phenotypeh${heritability}_s${sigma_beta}_r${causal_rate}_%m%d%H%M%S.csv
# /BASE/h${heritability}/simul_data_h${heritability}_s${sigma_beta}_r${causal_rate}_%m%d%H%M%S.csv


paste0("Begin at ", Sys.time())

## Set parameters
MAX_GENE <- 40000
MAX_SAMPLE <- 2000

## parse command line arguments
library(optparse, quietly = TRUE)

option_list <- list(
    make_option(c("-d", "--data_path"),
        type = "character", default = NULL,
        help = "Data path", metavar = "character"
    ),
    make_option(c("-h", "--heritability"),
        type = "numeric", default = 0.5,
        help = "Heritability [default= %default]", metavar = "numeric"
    ),
    make_option(c("-o", "--output_path"),
        type = "character", default = "./",
        help = "Output file path [default= %default]", metavar = "character"
    ),
    make_option(c("-r", "--causal_rate"),
        type = "numeric", default = 0.01, # 1e-3 -> 0.5
        help = "Causal rate [default= %default]", metavar = "numeric"
    ),
    make_option(c("-s", "--sigma_beta"),
        type = "numeric", default = 2,
        help = "Standard deviation of effect size [default= %default]", metavar = "numeric"
    ),
    make_option(c("-g", "--genotype_data"),
        type = "character", default = NULL,
        help = "Genotype data file path", metavar = "character"
    )
)

opt_parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
opt <- parse_args(opt_parser)

data.path <- opt$data_path
ldsc.path <- paste0(data.path, ".ldscore")
genotype.data <- opt$genotype_data
heritability <- opt$heritability
output.path <- opt$output_path
causal.rate <- opt$causal_rate
sigma_beta <- opt$sigma_beta

# # TEST
# # For Remote
# data.path <- "/home/wjiang49/scratch/height_ukb_50k_chr22"
# genotype.data <- "/home/wjiang49/scratch/UKBsimdata/genotypes.csv"
# output.path <- "/home/wjiang49/scratch/UKBsimdata/"
# # For Local
# # data.path <- "height_ukb_50k_chr22"
# # genotype.data <- "genotypes.csv"
# # output.path <- "."
# ---commom---
# ldsc.path <- paste0(data.path, ".ldscore")
# # genotype.data <- NULL
# heritability <- 0.6
# causal.rate <- 0.035
# sigma_beta <- 2

## Set output file name for identification
date.time <- gsub("\\.", "", format(Sys.time(), "%m%d%H%M%OS4"))
summary.data.name <- paste0("summary_data_h", as.character(heritability), "_s", as.character(sigma_beta), "_r", as.character(causal.rate), "_", date.time, ".csv")
phenotype.name <- paste0("phenotypeh", as.character(heritability), "_s", as.character(sigma_beta), "_r", as.character(causal.rate), "_", date.time, ".csv")

## Check if the output path exists
if (!dir.exists(output.path)) {
    dir.create(output.path)
}

paste0("Set heritability as ", heritability, ". Set causal rate as ", causal.rate, ". Set sigma_beta as ", sigma_beta, ".")

## Load data if genotype data is not provided
if (is.null(genotype.data)) {
    suppressPackageStartupMessages(library(snpStats))
    bed.file <- paste0(data.path, ".bed")
    fam.file <- paste0(data.path, ".fam")
    bim.file <- paste0(data.path, ".bim")
    genotype.data <- read.plink(bed.file, bim.file, fam.file)
    # print("Data loaded.")

    ## Get genotypes
    suppressPackageStartupMessages(library(dplyr))
    genotypes <- genotype.data$genotypes[1:MAX_SAMPLE, 1:MAX_GENE] %>%
        as("numeric") %>%
        as.data.frame(row.names = rownames(genotype.data$genotypes), col.names = colnames(genotype.data$genotypes)) %>%
        apply(., 2, function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))) %>% # impute missing values with mean
        scale() %>% # standardize genotypes (column)
        as.matrix()
    # be careful with sd=0 when n is small
    if (any(is.na(genotypes))) {
        stop("SD is 0. May try to increase the sample size.")
    }
    # print("Genotypes processed.")
} else {
    genotypes <- as.matrix(read.csv(genotype.data, header = TRUE))
}

## Set Ground Truth for Causal SNPs
# set.seed(123)
causal.snps <- sample(1:MAX_GENE, round(MAX_GENE * causal.rate))
# randomly set the effect size of causal SNPs
beta <- as.matrix(rnorm(length(causal.snps), sd = sigma_beta), ncol=1)
phenotype <- genotypes[, causal.snps] %*% beta
# phenotype <- as.numeric(phenotype)

if (any(is.na(phenotype))) {
    stop("Phenotype contains NA values.")
}

## Set heritability of phenotype
pheno.var <- var(phenotype)
residual.var <- pheno.var * (1 / heritability - 1)
residual <- rnorm(MAX_SAMPLE, 0, sqrt(residual.var))
phenotype <- scale(phenotype + residual)

## Save Genotypes and Phenotype for EM
if (is.null(opt$genotype_data)) {
    write.csv(genotypes, file.path(output.path, "genotypes.csv"), row.names = FALSE)
}
## check path exists
if (!dir.exists(file.path(output.path, paste0("h", as.character(heritability))))) {
    dir.create(file.path(output.path, paste0("h", as.character(heritability))))
}
write.csv(data.frame(phenotype = phenotype), file.path(output.path, paste0("h", as.character(heritability)), phenotype.name), row.names = FALSE)

## Calculate z-scores in Summary-Level Data
## Use lm
# y <- phenotype - mean(phenotype)
# x <- as.matrix(genotypes)
# cal_z <- function(xj) {
#     model.coef <- summary(lm(y ~ xj - 1))$coef
#     return(model.coef[1, 1] / model.coef[1, 2])
# }
# z_scores <- apply(x, 2, cal_z)

## Use matrix multiplication, equiv to lm
y <- as.matrix(phenotype - mean(phenotype), ncol = 1)
x <- as.matrix(genotypes)
x2 <- colSums(x^2)
betas <- (t(x) %*% y) / x2
ses <- sqrt(colSums((y %*% matrix(1, 1, MAX_GENE) - x * (matrix(1, MAX_SAMPLE, 1) %*% t(betas)))^2) / ((MAX_SAMPLE - 2) * x2))
z_scores <- betas / ses

## Concatenate LDSCORE
summary_data <- data.frame(SNP = colnames(genotypes), Z = z_scores, row.names = NULL)
ldsc.data <- read.csv(ldsc.path, header = TRUE, stringsAsFactors = FALSE)
summary_data <- cbind(summary_data, LDSCORE = ldsc.data[, "ldscore"])
# save to file
write.csv(summary_data, file.path(output.path, paste0("h", as.character(heritability)), summary.data.name), row.names = FALSE)

# paste0("Summary data saved to ", output.file)
paste0("End at ", Sys.time())
paste0(" ")
# End of file
