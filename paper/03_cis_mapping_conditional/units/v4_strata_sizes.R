## v4_strata_sizes.R -- phenotype only (no association result is read): how many units the v4 delta adds per type and
## per pre-specified stratum, using prepare_y of the final rules (SD < 0.005 skip) and the S2 estimability rule
## (units_v3_estimable.R: mu n_mid >= 80; p0/p1 boundary share 10-90%). Sets the size of every calibration stratum
## before the permuted delta scans are opened.
## Input: the delta freeze objects (lab staging rerun_v4_delta), core/v4/delta_promoters.tsv.
## Output: core/v4/strata_sizes.tsv; core/v4/delta_promoter_features.tsv (per promoter: n, n0, n1, n_mid, SD and mean on the scan samples)
## Usage: Rscript v4_strata_sizes.R   (lab)
suppressPackageStartupMessages(library(data.table))
W <- "/projects/group/bliu/e1327989/GTEx/nofloor_full"
source(file.path(W, "run_code/lab_v4/code/bohr2_data.R"))
D <- fread(file.path(W, "core/v4/delta_promoters.tsv"))
fs <- list.files(file.path(W, "GTEx/QTL_figure/rerun_v4_delta/work/final/freeze"), "[.]rds$", full.names = TRUE)
R <- rbindlist(lapply(fs, function(f) { x <- readRDS(f); pc <- 1 + ncol(x$pcs) + x$K
  rbindlist(lapply(colnames(x$act), function(p) { y0 <- x$act[, p]; py <- prepare_y(y0, pc)
    if (!is.na(py$skip)) return(data.table(tissue = x$tissue, promoter = p, skip = py$skip))
    y <- py$y; n <- length(y); n0 <- sum(y == 0); n1 <- sum(y == 1); ym <- y[y > 0 & y < 1]
    data.table(tissue = x$tissue, promoter = p, skip = NA_character_, n = n, n0 = n0, n1 = n1, n_mid = n - n0 - n1,
               sd_scan = sd(y0[!is.na(y0)]), mean_scan = mean(y0, na.rm = TRUE)) }), fill = TRUE) }), fill = TRUE)
R <- merge(R, unique(D[, .(tissue, promoter, reason)]), by = c("tissue", "promoter"), all.x = TRUE)
U <- rbind(R[is.na(skip), .(tissue, promoter, reason, type = "mu", present = TRUE, estimable = n_mid >= 80, sd_scan, mean_scan)],
           R[is.na(skip), .(tissue, promoter, reason, type = "p0", present = n0 > 0, estimable = n0 / n >= .1 & n0 / n <= .9, sd_scan, mean_scan)],
           R[is.na(skip), .(tissue, promoter, reason, type = "p1", present = n1 > 0, estimable = n1 / n >= .1 & n1 / n <= .9, sd_scan, mean_scan)])[present == TRUE]
fwrite(R[, .(tissue, promoter, reason, skip, n, n0, n1, n_mid, sd_scan = signif(sd_scan, 4), mean_scan = signif(mean_scan, 4))],
       file.path(W, "core/v4/delta_promoter_features.tsv"), sep = "\t")
U[, sd_band := cut(sd_scan, c(0, .02, .05, Inf), right = FALSE, labels = c("SD 0.005-0.02", "SD 0.02-0.05", "SD >= 0.05 (scan samples)"))]
U[, mean_band := cut(mean_scan, c(-Inf, .05, .95, Inf), labels = c("mean <= 0.05", "mean 0.05-0.95", "mean > 0.95"))]
cat("delta promoters:", nrow(R), "; skipped by prepare_y:", R[!is.na(skip), .N], "\n"); print(R[, .N, by = skip])
S <- rbind(U[, .(stratum = "all", units = .N, estimable = sum(estimable)), by = type],
           U[, .(units = .N, estimable = sum(estimable)), by = .(type, stratum = as.character(sd_band))],
           U[, .(units = .N, estimable = sum(estimable)), by = .(type, stratum = as.character(mean_band))],
           U[, .(units = .N, estimable = sum(estimable)), by = .(type, stratum = reason)])
S[, `:=`(E_perm0_1e3 = estimable * 1e-3, E_pooled3_1e3 = 3 * estimable * 1e-3, E_perm0_1e2 = estimable * 1e-2, E_pooled3_1e2 = 3 * estimable * 1e-2)]
setorder(S, type, stratum); fwrite(S, file.path(W, "core/v4/strata_sizes.tsv"), sep = "\t"); print(S, nrows = 100)
