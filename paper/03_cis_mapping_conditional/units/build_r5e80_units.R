#!/usr/bin/env Rscript
## Units v3m of the final rule set nofloor_reads5_E80 ("r5e80": no read floor, gene reads >= 5, endpoint skip 80%) for a
## fresh S3-S7 run in its own root GTEx/QTL_figure/rerun_r5e80 (the v4 root is not changed).
## Source: the v4 units v3m (rerun_v4/results/final_nogate_v3m/freeze/<T>.<arm>.units.rds, all arms).
## Filters (figures_v4/inputs/nofloor_reads5_E80/promoters_excluded.tsv.gz, make_version_inputs.R):
##   endpoint skip: a promoter is not tested in a tissue when its share of 0 or of 1 among the tissue's analysed samples is >= 0.80;
##   gene reads: every promoter of a gene whose median first-intron junction reads over the tissue's analysed samples is < 5.
## An excluded row is marked as not tested (estimable = FALSE, r5e80_excluded = TRUE). q is recomputed as in units v3m: BH within
## the tissue over all types, on rows with estimable & !dedup_drop & finite p (checked: this reproduces the stored v4 q exactly).
## Other columns are kept as they are (q_full and the v1/v2 history columns describe the v4 run).
## Output: rerun_r5e80/results/final_nogate_v3m/freeze/<T>.<arm>.units.rds and rerun_r5e80/build_r5e80_units.tsv
suppressPackageStartupMessages({ library(data.table); library(parallel) })
G <- "/projects/group/bliu/e1327989/GTEx/nofloor_full"
SRC <- file.path(G, "GTEx/QTL_figure/rerun_v4/results/final_nogate_v3m/freeze")
DST <- file.path(G, "GTEx/QTL_figure/rerun_r5e80/results/final_nogate_v3m/freeze")
EX <- fread(cmd = paste("zcat", file.path(G, "figures_v4/inputs/nofloor_reads5_E80/promoters_excluded.tsv.gz")))
ek <- paste(EX$tissue, EX$promoter)
fs <- list.files(SRC, "[.]units[.]rds$")
one <- function(f) {
  x <- as.data.table(readRDS(file.path(SRC, f))); t <- x$tissue[1]; arm <- x$arm[1]
  stopifnot(uniqueN(x$tissue) == 1, uniqueN(x$arm) == 1)
  i0 <- x$estimable %in% TRUE & x$dedup_drop %in% FALSE & is.finite(x$p)
  same <- isTRUE(all.equal(x$q[i0], p.adjust(x$p[i0], "BH")))                      # the stored q is this BH
  sig0 <- sum(i0 & x$q <= 0.05, na.rm = TRUE)
  x[, r5e80_excluded := paste(tissue, promoter) %in% ek]
  x[r5e80_excluded == TRUE, estimable := FALSE]
  x[, q := NA_real_]
  i1 <- x$estimable %in% TRUE & x$dedup_drop %in% FALSE & is.finite(x$p)
  x[i1, q := p.adjust(p, "BH")]
  saveRDS(x, file.path(DST, f))
  data.table(tissue = t, arm = arm, rows = nrow(x), tested_v4 = sum(i0), tested_r5e80 = sum(i1), excluded_rows = sum(x$r5e80_excluded),
             sig_v4 = sig0, sig_r5e80 = sum(x$q <= 0.05, na.rm = TRUE), sig_mu = x[q <= 0.05 & type == "mu", .N], sig_p0 = x[q <= 0.05 & type == "p0", .N],
             sig_p1 = x[q <= 0.05 & type == "p1", .N], thr_v4 = if (sig0) max(x$p[i0][x$q[i0] <= 0.05 & FALSE | p.adjust(x$p[i0], "BH") <= 0.05]) else NA_real_,
             thr_r5e80 = if (any(x$q <= 0.05, na.rm = TRUE)) max(x$p[which(x$q <= 0.05)]) else NA_real_, q_reproduced = same)
}
R <- rbindlist(mclapply(fs, one, mc.cores = 16))
fwrite(R, file.path(G, "GTEx/QTL_figure/rerun_r5e80/build_r5e80_units.tsv"), sep = "\t")
cat("files", nrow(R), "; stored q reproduced in all:", all(R$q_reproduced), "\n")
print(R[, .(files = .N, tested_v4 = sum(tested_v4), tested_r5e80 = sum(tested_r5e80), sig_v4 = sum(sig_v4), sig_r5e80 = sum(sig_r5e80),
            mu = sum(sig_mu), p0 = sum(sig_p0), p1 = sum(sig_p1)), by = arm])
R2 <- R[arm == "real"]; cat("real-arm threshold ratio r5e80 / v4: median", round(median(R2$thr_r5e80 / R2$thr_v4, na.rm = TRUE), 3), "\n")
## cross-check with make_version_inputs.R (units v3c, same P and q as v3m): significant after = units_keep + units_new_topup
K <- fread(cmd = paste("zcat", file.path(G, "figures_v4/inputs/nofloor_reads5_E80/units_keep.tsv.gz")))
N <- fread(cmd = paste("zcat", file.path(G, "figures_v4/inputs/nofloor_reads5_E80/units_new_topup.tsv.gz")))
cat("make_version_inputs.R significant after:", nrow(K) + nrow(N), "; here (real):", sum(R2$sig_r5e80), "\n")
