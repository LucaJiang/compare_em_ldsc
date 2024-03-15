# Calculate r square and ldscore with snpStats
# Usage: Rscript UKBheight/cal_ld.R -d /home/wjiang49/scratch/height_ukb_50k_chr22 -N 2000 -P 40000 -o /home/wjiang49/scratch/UKBsimdata

paste0("Begin at ", Sys.time())
# Parse command line arguments
library(optparse, quietly = TRUE)
option_list <- list(
    make_option(c("-d", "--data_path"),
        type = "character", default = NULL,
        help = "Data path", metavar = "character"
    ),
    make_option(c("-N", "--num_samples"),
        type = "numeric", default = 2000,
        help = "Sample size [default= %default]", metavar = "numeric"
    ),
    make_option(c("-P", "--num_snps"),
        type = "numeric", default = 40000,
        help = "Number of SNPs [default= # row ldsc]", metavar = "numeric"
    ),
    make_option(c("-o", "--output_path"),
        type = "character", default = "./",
        help = "Output file path [default= %default]", metavar = "character"
    )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
data.path <- opt$data_path
output.path <- opt$output_path
n_sample <- opt$num_samples
n_snp <- opt$num_snps

suppressPackageStartupMessages(library(snpStats))
bed.file <- paste0(data.path, ".bed")
fam.file <- paste0(data.path, ".fam")
bim.file <- paste0(data.path, ".bim")
# check file existence
if (!file.exists(bed.file) || !file.exists(fam.file) || !file.exists(bim.file)) {
    stop("Genotype data not found")
}

# Read data and impute NA with mean
genotypes.data <- read.plink(bed.file, bim.file, fam.file)
suppressPackageStartupMessages(library(dplyr))
genotypes <- genotypes.data$genotypes[1:n_sample, 1:n_snp] %>%
    as("numeric") %>%
    as.data.frame(row.names = rownames(genotypes.data$genotypes), col.names = colnames(genotypes.data$genotypes)) %>%
    apply(., 2, function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))) %>% # impute missing values with mean
    scale() %>% # standardize genotypes (column)
    as.matrix()

snp.ldscore <- rowSums(cor(genotypes, genotypes)^2)
snp.position <- genotypes.data$map$position[1:n_snp]
snp.ldsc.data <- data.frame(chromosome = 22, position = snp.position, ldscore = snp.ldscore)
suppressPackageStartupMessages(library(data.table))
fwrite(snp.ldsc.data, output.path, row.names = TRUE)
