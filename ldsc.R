# Rscript to run ldsc.r with given data through command line
# Usage:
# Rscript /home/wjiang49/UKBheight/ldsc.R -d /home/wjiang49/scratch/UKBsimdata/h0.6/summary_data_h0.6_s0.1_r0.001_03150017435666.csv -N 2000 -o /home/wjiang49/scratch/UKBsimdata
# Save the log file in output_path/ldsc.log and the result in output_path/ldsc_results.csv

paste0("Begin at ", Sys.time())

## parse command line arguments
library(optparse, quietly = TRUE)

option_list <- list(
    make_option(c("-d", "--data_path"),
        type = "character", default = "summary_data_h0.3_s5_r0.2_0312212952.csv",
        help = "Data path", metavar = "character"
    ),
    make_option(c("-N", "--num_samples"),
        type = "numeric", default = 2000,
        help = "Sample size", metavar = "numeric"
    ),
    make_option(c("-o", "--output_path"),
        type = "character", default = "/Users/lucajiang/learn/CityU/UKBheight/",
        help = "Output file path", metavar = "character"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

data.path <- opt$data_path
output.path <- opt$output_path
n_sample <- opt$num_samples

# test
# data.path <- "/Users/lucajiang/learn/CityU/UKBheight/compare_em_ldsc/test/summary_data_h0.8_s1_r1_03151639444055.csv"
# output.path <- "/Users/lucajiang/learn/CityU/UKBheight/"
# n_sample <- 2000

# calculate h2 ---------------------
get_coef_raw <- function(ldscore_mat, sumstat1, sumstat2, max_int) {
    #########
    # obtain raw coefficients for h2 / genetic covariance estimates
    #########
    p <- nrow(ldscore_mat)
    ncoef <- ncol(ldscore_mat)

    nprod <- sqrt(sumstat1$N) * sqrt(sumstat2$N)
    zprod <- sumstat1$Z * sumstat2$Z
    score_sum <- rowSums(matrix(ldscore_mat[, -1], nrow = p))

    # get averages
    mean_ldscore <- mean(score_sum)
    mean_zprod <- mean(zprod)
    mean_nprod <- mean(nprod)

    # estimate intercept
    zprod_sorted <- sort(zprod, decreasing = F)
    idx <- 0.95 * p
    intercept <- mean(zprod_sorted[1:idx])
    intercept <- ifelse(intercept > max_int, max_int, intercept)

    # get raw coefs
    coef <- (mean_zprod - intercept) / (mean_nprod * mean_ldscore)

    return(list(intercept = intercept, coef = coef))
}

get_pred <- function(coef, x, n1, n2, intercept) {
    #########
    # obtain predicted chi-squared for constructing regression weight
    #########
    p <- nrow(x)

    nprod <- sqrt(n1) * sqrt(n2)
    score_sum <- rowSums(matrix(x[, -1], nrow = p))
    pred <- coef * score_sum * nprod + intercept

    return(pred)
}

update_weight <- function(w, pred) {
    #########
    # update weight using predicted chi-squared for one trait
    #########
    var <- 2 * pred^2

    w <- w / var

    return(w)
}

regression <- function(x, y, constrain_intercept, subtract, nblocks = 200, jknife = F) {
    #########
    # perform least squared regression to get coefs
    #########

    ncoef <- ncol(x)

    xtx <- t(x) %*% x
    xty <- t(x) %*% y

    # obtain coefs
    if (constrain_intercept) {
        coefs <- solve(xtx[-1, -1], xty[-1])
        coefs <- c(subtract, coefs)
    } else {
        coefs <- solve(xtx, xty)
    }

    if (jknife) {
        seperator <- floor(seq(from = 1, to = length(y), length.out = nblocks + 1))
        from <- seperator[-length(seperator)]
        to <- c(seperator[2:length(y)] - 1, length(y))

        coefs_jk <- matrix(0, nblocks, ncoef)
        for (i in 1:nblocks) {
            xtx_blk <- t(x[from[i]:to[i], ]) %*% x[from[i]:to[i], ]
            xty_blk <- t(x[from[i]:to[i], ]) %*% y[from[i]:to[i]]

            xtx_loo <- xtx - xtx_blk
            xty_loo <- xty - xty_blk

            if (constrain_intercept) {
                coefs_jk[i, -1] <- solve(xtx_loo[-1, -1], xty_loo[-1])
            } else {
                coefs_jk[i, ] <- solve(xtx_loo, xty_loo)
            }
        }
        jk_cov <- cov(coefs_jk) * (nblocks - 1) # cov(t(nblocks*c(coefs)-t((nblocks-1)*coefs_jk))) / nblocks
        jk_se <- sqrt(diag(jk_cov))
    } else {
        jk_se <- NA
    }

    return(list(
        coefs = coefs,
        coefs_se = jk_se
    ))
}

get_coef <- function(score, sumstat1, sumstat2, w, constrain_intercept, subtract, nblocks = 200, jknife = F) {
    #########
    # wrapper function to get coefs
    #########
    p <- nrow(score)
    nprod <- sqrt(sumstat1$N) * sqrt(sumstat2$N)
    zprod <- sumstat1$Z * sumstat2$Z
    if (constrain_intercept) {
        zprod <- zprod - subtract
    }
    score[, -1] <- score[, -1] * nprod

    # scale the matrix to improve matrix condition
    nbar <- mean(nprod)
    score[, -1] <- score[, -1] / nbar

    # apply weight to data
    score_w <- score * c(sqrt(w))
    zprod_w <- zprod * c(sqrt(w))

    # perform regression
    reg <- regression(score_w, zprod_w, constrain_intercept, subtract, nblocks = nblocks, jknife = jknife)

    # re-scale coefs
    reg$coefs[-1] <- reg$coefs[-1] / nbar

    if (jknife) {
        reg$coefs_se[-1] <- reg$coefs_se[-1] / nbar
    }

    return(reg)
}

estimate_h2 <- function(sumstat, ldscore, reg_w = 1, constrain_intercept = F, int = 1) {
    # add an additional column for intercept
    ldscore <- cbind(1, ldscore)

    # get initial estimate of coefs for constructing weights
    raw_coef <- get_coef_raw(ldscore, sumstat, sumstat, 1)

    # get predicted chi-squared for constructing weights
    n <- sumstat$N
    pred <- get_pred(raw_coef$coef, ldscore, n, n, raw_coef$intercept)

    # update weight
    reg_w <- update_weight(reg_w, pred)

    # conduct step 1 regression to obtain coefs
    coef <- get_coef(ldscore, sumstat, sumstat, reg_w, constrain_intercept, int)

    return(coef)
}

read_summary_data_with_ldsc <- function(file_path) {
    # read file
    if (!file.exists(file_path)) {
        stop("File not found at ", file_path)
    }
    data <- read.csv(file_path, header = TRUE)
    sumstats <- data.frame(ldscore = data$LDSCORE, Z = data$Z)
    return(sumstats)
}

sumstats <- read_summary_data_with_ldsc(data.path)
fit_ldsc <- estimate_h2(data.frame(Z = sumstats$Z, N = n_sample), sumstats$ldscore, constrain_intercept = T, int = 1)
p <- length(sumstats$ldscore)
h2.est <- fit_ldsc$coefs[2] * p
h2.est <- round(h2.est, 6)
if (is.null(h2.est)) {
    stop("h2 estimation failed")
}

# get real h2 from file name
param <- strsplit(data.path, "summary_data_h")[[1]][2]
param <- strsplit(param, ".csv")[[1]][1]
param <- strsplit(param, "_")[[1]]
h2.real <- as.numeric(param[1])
s <- as.numeric(substring(param[2], 2))
r <- as.numeric(substring(param[3], 2))
data.name <- as.numeric(param[4])

# save result
paste0("est: ", h2.est, " real: ", h2.real)
result <- data.frame(data = data.name, h = h2.real, sigma_beta = s, causal_rate = r, hest = h2.est)
result.file <- file.path(output.path, "ldsc_results.csv")
if (!file.exists(result.file)) {
    write.table(result, file = result.file, sep = ",", col.names = TRUE, row.names = FALSE)
} else {
    write.table(result, file = result.file, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
}

paste0("End at ", Sys.time())
paste0(" ")
# End of file
