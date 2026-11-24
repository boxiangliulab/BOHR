## diag_v4_calibration.R -- decision on the units the v4 delta adds (pre-specified rules, fixed before
## any permuted delta result existed). Decision stratum = added units of one type (mu, p0, p1) that enter BH
## (estimable, dedup_drop FALSE). Admitted when all three hold:
##   1. fit diagnostics (real arm; diag_v4_fit.R on NSCC): >= 99% of SNP tests converged and >= 99% of null models;
##   2. perm0-2 pooled, at the most stringent alpha in {1e-4, 1e-3, 1e-2, 5e-2} with pooled expectation >= 100: upper
##      limit of the exact Poisson 95% interval of observed / expected <= 1.5;
##   3. pooled count at P <= 1e-4 not above the 99.9% Poisson quantile of its expectation.
## Reported, not decisive: observed / expected at 1e-3 and 1e-4 for perm0 alone and pooled; the same for the v2-set
## units (reference); secondary strata (SD band, mean band, reason added); permuted units reaching q <= 0.05.
## After the decision is written: real-arm counts of added units (information only).
## Inputs: units of rerun_v4 S2 (results/final_nogate_v3m/freeze/<T>.<arm>.units.rds, four arms, fetched from NSCC),
##   core/v4/delta_promoter_features.tsv (v4_strata_sizes.R), core/v4/fit_diagnostics.tsv (diag_v4_fit.R).
## Outputs: core/v4/calibration_decision.tsv, core/v4/calibration_all.tsv, core/v4/real_added_counts.tsv
## Usage: Rscript diag_v4_calibration.R <units folder>   (lab)
suppressPackageStartupMessages(library(data.table))
W <- "/projects/group/bliu/e1327989/GTEx/nofloor_full"; IN <- commandArgs(TRUE)[1]; OUT <- Sys.getenv("V4DIAG_OUT", file.path(W, "core/v4"))   # V4DIAG_OUT: tests only
ALPHA <- c(1e-4, 1e-3, 1e-2, 5e-2); EMIN <- 100; BOUND <- 1.5
F <- fread(file.path(W, "core/v4/delta_promoter_features.tsv"))
fs <- list.files(IN, "[.]units[.]rds$", full.names = TRUE)
U <- rbindlist(lapply(fs, function(f) as.data.table(readRDS(f))[, .(tissue, arm, promoter, type, p, q, estimable, dedup_drop, n, n0, n1)]))
cat("units files:", length(fs), "; arms:", paste(sort(unique(U$arm)), collapse = " "), "; tissues:", uniqueN(U$tissue), "\n")
stopifnot(setequal(unique(U$arm), c("real", "perm0", "perm1", "perm2")))
U <- U[is.finite(p) & estimable %in% TRUE & dedup_drop %in% FALSE]
U <- merge(U, F[, .(tissue, promoter, reason, sd_scan, mean_scan)], by = c("tissue", "promoter"), all.x = TRUE)
U[, set := fifelse(is.na(reason), "v2 set (reference)", "added")]
U[, sd_band := fifelse(set == "added", as.character(cut(sd_scan, c(0, .02, .05, Inf), right = FALSE, labels = c("SD 0.005-0.02", "SD 0.02-0.05", "SD >= 0.05 on scan samples"))), NA_character_)]
U[, mean_band := fifelse(set == "added", as.character(cut(mean_scan, c(-Inf, .05, .95, Inf), labels = c("mean <= 0.05", "mean 0.05-0.95", "mean > 0.95"))), NA_character_)]
pois_hi <- function(k, e) qgamma(0.975, k + 1) / e
pois_lo <- function(k, e) ifelse(k == 0, 0, qgamma(0.025, k) / e)
tails <- function(x) rbindlist(lapply(ALPHA, function(a) data.table(alpha = a, units = nrow(x), obs = sum(x$p <= a))))
## every stratum x arm group x tail
grp <- list(perm0 = "perm0", pooled = c("perm0", "perm1", "perm2"))
strat <- rbind(U[, .(tissue, arm, type, p, q, set, stratum = "all")],
               U[set == "added", .(tissue, arm, type, p, q, set, stratum = sd_band)],
               U[set == "added", .(tissue, arm, type, p, q, set, stratum = mean_band)],
               U[set == "added", .(tissue, arm, type, p, q, set, stratum = reason)])
A <- rbindlist(lapply(names(grp), function(g) strat[arm %in% grp[[g]], tails(.SD), by = .(set, type, stratum)][, arms := g]))
A[, `:=`(expected = units * alpha)][, `:=`(ratio = round(obs / expected, 3), lo95 = round(pois_lo(obs, expected), 3), hi95 = round(pois_hi(obs, expected), 3))]
setcolorder(A, c("set", "type", "stratum", "arms", "alpha", "units", "obs", "expected", "ratio", "lo95", "hi95"))
setorder(A, set, type, stratum, arms, alpha); fwrite(A, file.path(OUT, "calibration_all.tsv"), sep = "\t")
## decision
FD <- fread(file.path(OUT, "fit_diagnostics.tsv"))[set == "added" & estimable %in% TRUE, .(type, fit_snp_conv = snp_conv, fit_null_conv = null_conv)]
D <- rbindlist(lapply(c("mu", "p0", "p1"), function(ty) {
  x <- A[set == "added" & type == ty & stratum == "all" & arms == "pooled"]
  a_star <- suppressWarnings(min(x[expected >= EMIN, alpha])); if (!is.finite(a_star)) a_star <- NA_real_
  r2 <- if (is.na(a_star)) x[alpha == max(ALPHA)] else x[alpha == a_star]
  r3 <- x[alpha == 1e-4]; p0 <- A[set == "added" & type == ty & stratum == "all" & arms == "perm0"]
  ref <- A[set == "v2 set (reference)" & type == ty & stratum == "all" & arms == "pooled"]
  data.table(type = ty, units_pooled = r2$units, alpha_decision = r2$alpha, obs = r2$obs, expected = r2$expected, ratio = r2$ratio, hi95 = r2$hi95,
             rule2_pass = !is.na(a_star) & r2$hi95 <= BOUND,
             obs_1e4 = r3$obs, exp_1e4 = r3$expected, q999_1e4 = qpois(0.999, r3$expected), rule3_pass = r3$obs <= qpois(0.999, r3$expected),
             perm0_ratio_1e3 = p0[alpha == 1e-3, ratio], perm0_hi95_1e3 = p0[alpha == 1e-3, hi95], pooled_ratio_1e3 = x[alpha == 1e-3, ratio], pooled_hi95_1e3 = x[alpha == 1e-3, hi95],
             pooled_ratio_1e4 = r3$ratio, ref_ratio_1e3 = ref[alpha == 1e-3, ratio], ref_ratio_at_alpha = ref[alpha == r2$alpha, ratio]) }))
D <- merge(D, FD, by = "type", all.x = TRUE)
D[, rule1_pass := fit_snp_conv >= 0.99 & fit_null_conv >= 0.99]
D[, admit := rule1_pass %in% TRUE & rule2_pass %in% TRUE & rule3_pass %in% TRUE]
fwrite(D, file.path(OUT, "calibration_decision.tsv"), sep = "\t")
cat("\n== decision (rules fixed before any permuted delta result) ==\n"); print(D[, .(type, units_pooled, alpha_decision, obs, expected, ratio, hi95, rule2_pass, obs_1e4, q999_1e4, rule3_pass, fit_snp_conv, fit_null_conv, rule1_pass, admit)])
cat("\n== reported: observed / expected at 1e-3 and 1e-4 ==\n")
print(A[stratum == "all" & alpha %in% c(1e-3, 1e-4), .(set, type, arms, alpha, units, obs, expected, ratio, hi95)], nrows = 100)
cat("\n== reported: secondary strata of the added units, pooled arms ==\n")
print(A[set == "added" & stratum != "all" & arms == "pooled" & alpha %in% c(1e-3, 1e-2), .(type, stratum, alpha, units, obs, expected, ratio, hi95)], nrows = 100)
## permuted units at q <= 0.05 (BH within tissue x arm, on the whole primary set)
PQ <- U[arm != "real", .(sig_q05 = sum(q <= 0.05, na.rm = TRUE)), by = .(arm, set, type)]
cat("\n== permuted units at q <= 0.05 ==\n"); print(dcast(PQ, set + type ~ arm, value.var = "sig_q05"))
## information only, after the decision: real-arm counts of the added units
RC <- U[arm == "real", .(units = .N, sig_q05 = sum(q <= 0.05, na.rm = TRUE)), by = .(set, type)]
fwrite(RC, file.path(OUT, "real_added_counts.tsv"), sep = "\t"); cat("\n== real arm (information only) ==\n"); print(RC)
