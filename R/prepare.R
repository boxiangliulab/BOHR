## BOHR data preparation: endpoint rule and family choice, genotype dosage and
## modal imputation, cis window, covariate standardisation and the per-promoter
## orthonormal design.
##
## Ported from the analysis engine of the BOHR paper,
##   Figure4_puQTL_landscape_eQTL_coloc/rerun_v2/code/bohr2_data.R
##   (md5 2c190bfb6c3f5e79b85f27613a849abc).
## Verbatim (byte-identical bodies): gt_dosage (lines 55-61), impute_modal
## (62-68), orth_design (137-149).
## prepare_y (37-52): the four thresholds of the original body (20, 0.95, 0.05,
## 0.05) are arguments whose defaults are those values (author decision
## 2026-11-10: every threshold is an argument); nothing else changed.
## Changed signature, same arithmetic: cis_index (81-82) takes the window and
## the variant start/end positions as arguments instead of the global FLANK
## and a VCF object; cov_matrix (130-136) takes a samples x k covariate matrix
## instead of the tissue object, the PC/PEER split and K.
##
## All functions in this file are internal (not exported).

## ---- old-pipeline promoter rules ---------------------------------------------------------------
## gamlss_functions_optimized.R:548-567 (skips, endpoint deletion) and :34-47 (select_model).
## y: named numeric over the analysis samples (NA allowed); p_cov = number of covariates.
prepare_y <- function(y, p_cov, min_n = 20, max_boundary = 0.95, min_sd = 0.05, endpoint = 0.05) {
  y <- y[!is.na(y)]
  if (length(y) < min_n || length(y) < p_cov) return(list(skip = "n"))
  p0 <- mean(y == 0); p1 <- mean(y == 1)
  if (p0 >= max_boundary || p1 >= max_boundary || sd(y) <= min_sd) return(list(skip = "prop_or_sd"))
  keep <- rep(TRUE, length(y))
  if (p0 <= endpoint) keep <- keep & y > 0
  if (p1 <= endpoint) keep <- keep & y < 1
  y <- y[keep]
  z <- any(y == 0); o <- any(y == 1)
  list(skip = NA_character_, y = y,
       family = if (z && o) "BEINF" else if (z) "BEZI" else if (o) "BEOI" else "BE")
}

## ---- genotypes ---------------------------------------------------------------------------------
gt_dosage <- function(gt) {
  D <- matrix(NA_integer_, nrow(gt), ncol(gt), dimnames = dimnames(gt))
  D[gt == "0|0" | gt == "0/0"] <- 0L
  D[gt == "0|1" | gt == "1|0" | gt == "0/1" | gt == "1/0"] <- 1L
  D[gt == "1|1" | gt == "1/1"] <- 2L
  D
}
impute_modal <- function(D) {                          # rows = SNPs, cols = samples
  for (i in which(rowSums(is.na(D)) > 0)) {
    tb <- tabulate(D[i, ] + 1L, nbins = 3L)
    D[i, is.na(D[i, ])] <- which.max(tb) - 1L         # first maximum = smaller dosage
  }
  D
}

## cis window of the old code: flank(GRanges(chr, pos), 500000, both = TRUE) = [pos-5e5, pos+5e5-1]
## (bohr2_data.R:81-82: geno$start <= pos + FLANK - 1L & geno$end >= max(1L, pos - FLANK)).
## start/end: variant start and end positions (end = start for a SNV); window = FLANK.
cis_index <- function(start, end, pos, window) which(start <= pos + window - 1L & end >= max(1L, pos - window))

## design: intercept + standardised covariates (bohr2_data.R:130-136, where Z = cbind(PC1-5, PEER 1..K)).
## Z: samples x k numeric matrix (k may be 0); scaled over all of its rows (the analysis samples).
cov_matrix <- function(Z) {
  Z <- scale(Z)
  if (any(!is.finite(Z))) stop("constant covariate column")
  cbind(`(Intercept)` = 1, Z)
}

## per-promoter orthonormal basis of the design (same column space, so the maximised likelihood and
## every SNP statistic are unchanged); columns that are collinear at the lm() tolerance are dropped.
## First column stays the intercept (exactly 1).
orth_design <- function(C, tol = 1e-7) {
  q <- qr(C, tol = tol); r <- q$rank
  Q <- qr.Q(q)[, seq_len(r), drop = FALSE]
  s <- sign(diag(qr.R(q))[seq_len(r)]); s[s == 0] <- 1
  Q <- sweep(Q, 2, s, "*") * sqrt(nrow(C))
  Q[, 1] <- 1
  dimnames(Q) <- list(rownames(C), c("(Intercept)", if (r > 1) paste0("z", seq_len(r - 1))))
  attr(Q, "dropped") <- colnames(C)[q$pivot[-seq_len(r)]]
  Q
}
