# Calculate r square and ldscore with snpStats
library(optparse, quietly = TRUE)
MAX_GENE <- 40000
MAX_SAMPLE <- 2000


option_list <- list(
    make_option(c("-d", "--data_path"),
        type = "character", default = NULL,
        help = "Data path", metavar = "character"
    ),
    make_option(c("-h", "--heritability"),
        type = "numeric", default = 0.5,
        help = "Heritability [default= %default]", metavar = "numeric"
    ),

suppressPackageStartupMessages(library(snpStats))
suppressPackageStartupMessages(library(data.table))
data.path <- "/home/wjiang49/scratch/height_ukb_50k_chr22"
bed.file <- paste0(data.path, ".bed")
fam.file <- paste0(data.path, ".fam")
bim.file <- paste0(data.path, ".bim")

genotypes.data <- read.plink(bed.file, bim.file, fam.file)

# ! recalcaulate at 2024-03-12 20:20
# snp.r2 <- ld(genotypes.data$genotypes[1:MAX_SAMPLE, 1:MAX_GENE], depth = MAX_GENE - 1, stats = "R.squared", symmetric = TRUE)
# snp.ldscore <- rowSums(snp.r2, na.rm = TRUE)

# ! recalcaulate at 2024-03-13 14:26
suppressPackageStartupMessages(library(dplyr))
genotypes <- genotypes.data$genotypes[1:MAX_SAMPLE, 1:MAX_GENE] %>%
    as("numeric") %>%
    as.data.frame(row.names = rownames(genotypes.data$genotypes), col.names = colnames(genotypes.data$genotypes)) %>%
    # apply(c(1, 2), function(x) 2 - x) %>%
    # mapping rules: geno.code <- Dict$new("A/A" = 2L, "A/B" = 1L, "B/B" = 0L, "NA" = NA)
    apply(., 2, function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))) %>% # impute missing values with mean
    scale() %>% # standardize genotypes (column)
    as.matrix()

snp.ldscore <- rowSums(cor(genotypes, genotypes)^2)

snp.position <- genotypes.data$map$position[1:MAX_GENE]
snp.ldsc.data <- data.frame(chromosome = 22, position = snp.position, ldscore = snp.ldscore)
fwrite(snp.ldsc.data, "/home/wjiang49/scratch/height_ukb_50k_chr22.ldscore", row.names = TRUE)