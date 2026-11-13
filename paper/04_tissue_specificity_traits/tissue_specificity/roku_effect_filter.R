## =============================================================================
## roku_effect_filter.R -- the post-filter of the ROKU tissue-specific calls is changed to an effect-size threshold (author 2026-10-07)
## Input (read only): ../v11_tissue_specific_ROKU_results.rds (TCC::ROKU, upper.limit = 0.25; produced by ROKU.R)
##              ../v11_raw_fit_fit_results.txt (overall_mean of each tissue = p1 + (1 - p0 - p1) * mu, the ROKU input)
## Rule: a raw ROKU outlier call (outlier = +1 / -1) is kept if and only if
##        |overall_mean of the tissue - median overall_mean of the non-outlier tissues of the promoter| >= DELTA (0.10 = 10 percentage points of the gene's promoter usage)
##      It is not kept when there are fewer than 2 non-outlier tissues (same as the old z test).
## Why the change: the old z test used the mean/SD of all tissues (including the outlier tissues themselves) as reference, and the outlier tissues inflate each other's SD (masking):
##   the fraction kept falls with the number of simultaneously outlying tissues (1 tissue 81% -> 12+ tissues 23%); with >= 12 equally elevated tissues none passes, however large the increase;
##   and neither the normality assumption of z nor the uncorrected 0.05 holds.
## Output: v11_tissue_specific_ROKU_results_delta010.rds (same structure as the original file: outlier unchanged, true_outliers replaced by the new filter,
##       plus delta_matrix = difference from the non-outlier median, and a filter description); summary_delta_filter.tsv
## =============================================================================
suppressMessages(library(data.table))
DELTA <- 0.10
FIT <- "/home/e1327989/project/GTEx/fit"
OUT <- file.path(FIT, "roku_20261007")
r <- readRDS(file.path(FIT, "v11_tissue_specific_ROKU_results.rds"))
t <- fread(file.path(FIT, "v11_raw_fit_fit_results.txt"), select = c("tissue", "promoter", "overall_mean"))
W <- dcast(t, promoter ~ tissue, value.var = "overall_mean")
X <- as.matrix(W[, -1]); rownames(X) <- W$promoter
O <- r$outlier; X <- X[rownames(O), colnames(O)]
K <- O; K[!is.na(K)] <- 0                        # New true_outliers (NA positions stay NA, as in the original file)
Dm <- matrix(NA_real_, nrow(O), ncol(O), dimnames = dimnames(O))
for (i in seq_len(nrow(O))) {
  x <- X[i, ]; o <- O[i, ]
  out <- which(!is.na(o) & o != 0 & !is.na(x)); if (!length(out)) next
  nor <- which(!is.na(o) & o == 0 & !is.na(x)); if (length(nor) < 2) next
  d <- x[out] - median(x[nor]); Dm[i, out] <- d
  K[i, out] <- ifelse(abs(d) >= DELTA, o[out], 0)
}
res <- r; res$true_outliers <- K; res$delta_matrix <- Dm; res$p_values <- NULL
res$filter <- sprintf("keep ROKU call if |overall_mean - median(non-outlier tissues)| >= %.2f (2026-10-07; replaces the z-test)", DELTA)
saveRDS(res, file.path(OUT, "v11_tissue_specific_ROKU_results_delta010.rds"))
cnt <- function(M) { M[is.na(M)] <- 0; c(up_calls = sum(M == 1), down_calls = sum(M == -1), promoters_any = sum(rowSums(M != 0) > 0),
  promoters_up = sum(rowSums(M == 1) > 0), promoters_down = sum(rowSums(M == -1) > 0), promoters_both = sum(rowSums(M == 1) > 0 & rowSums(M == -1) > 0)) }
S <- rbind(roku_raw = cnt(O), ztest_old = cnt(r$true_outliers), delta010 = cnt(K))
S <- data.table(set = rownames(S), S); fwrite(S, file.path(OUT, "summary_delta_filter.tsv"), sep = "\t"); print(S)
cat("promoters assessable:", nrow(O), "| tissues:", ncol(O), "\n")
