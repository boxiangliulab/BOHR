## BOHR v2 -- data preparation shared by the pilot and the full rerun.
## Mirrors the old pipeline (code_orig/gamlss_QTL_optimized.R, gamlss_functions_optimized.R) except
## where a documented fix applies:
##   * missing genotypes imputed per SNP (modal dosage over the analysis samples; ties -> smaller
##     dosage). The old get_geno() imputed per DONOR across SNPs with a random tie-break (:236-243).
##   * covariates standardised (mean 0, sd 1). The likelihood maximum and every SNP statistic are
##     invariant to this; it only conditions the design (some PEER columns are O(1e-6)).
suppressPackageStartupMessages({
  library(data.table)
  library(VariantAnnotation)
})

## local defaults; on NSCC the job sets BOHR2_GTEX / BOHR2_ROOT / BOHR2_VCF (only cached inputs are used there)
GTEX <- Sys.getenv("BOHR2_GTEX", "/home/e1327989/project/GTEx")
QLC  <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/QTLtools_linear_compare"   # read-only
V2   <- Sys.getenv("BOHR2_ROOT", "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2")
PATHS <- list(
  act   = file.path(GTEX, "v11/filter_5_recal_att_filter.txt"),
  info  = file.path(GTEX, "v11/Promoter_info_filter_5.txt"),
  pcs   = file.path(GTEX, "metadata/v10_cov_pca.txt"),
  peer  = file.path(GTEX, "peer/v11"),                    # v11_filter_5_peer_factors_<T>.txt
  vcf   = Sys.getenv("BOHR2_VCF", file.path(GTEX, "Rawdata/v9_vcf_MAF_0.05")),   # chr<N>.vcf.gz (+ .tbi)
  proms = file.path(QLC, "perm/bohr"),                    # <T>.chr<k>_promoters.txt = old tested set
  work  = file.path(V2, "work"),
  maps  = file.path(V2, "perm_maps")
)
FLANK <- 500000L; N_PC <- 5L; MIN_GENO <- 5L; CHRS <- paste0("chr", 1:5)

## pilot tissues, current (best) PEER and the permutation index used for newly drawn maps
## (idx continues the rep_permutations convention: Hippocampus 0, Adipose 1, Skin 2).
PILOT <- data.table(
  tissue = c("Skin_Sun_Exposed_Lower_leg", "Adipose_Visceral_Omentum", "Brain_Hippocampus",
             "Heart_Left_Ventricle", "Brain_Cerebellar_Hemisphere", "Kidney_Cortex", "Bladder"),
  k_now  = c(7L, 15L, 3L, 52L, 45L, 10L, 7L),
  idx    = c(2L, 1L, 0L, 3L, 4L, 5L, 6L))

## ---- old-pipeline promoter rules ---------------------------------------------------------------
## gamlss_functions_optimized.R:548-567 (skips, endpoint deletion) and :34-47 (select_model).
## y: named numeric over the analysis samples (NA allowed); p_cov = number of covariates.
prepare_y <- function(y, p_cov) {
  y <- y[!is.na(y)]
  if (length(y) < 20 || length(y) < p_cov) return(list(skip = "n"))
  p0 <- mean(y == 0); p1 <- mean(y == 1)
  if (p0 >= 0.95 || p1 >= 0.95 || sd(y) < 0.005) return(list(skip = "prop_or_sd"))   # final rule: skip only SD < 0.005 (the old pipeline skipped SD <= 0.05)
  keep <- rep(TRUE, length(y))
  if (p0 <= 0.05) keep <- keep & y > 0
  if (p1 <= 0.05) keep <- keep & y < 1
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
## One tabix query for a whole block of promoters. Record selection = readVcf overlap semantics,
## as in the old get_geno(). Returns dosage as samples x SNPs (imputed) + positions.
read_geno <- function(chr, from, to, samples) {
  prm <- ScanVcfParam(which = GRanges(chr, IRanges(max(1L, from), to)), samples = samples,
                      fixed = NA, info = NA, geno = "GT")
  v <- readVcf(TabixFile(file.path(PATHS$vcf, paste0(chr, ".vcf.gz"))), "hg38", prm)
  if (nrow(v) == 0) return(NULL)
  gt <- geno(v)$GT[, samples, drop = FALSE]
  D <- gt_dosage(gt); nmiss <- sum(is.na(D)); D <- impute_modal(D)
  rr <- rowRanges(v)
  list(id = rownames(gt), start = start(rr), end = end(rr), G = t(D), n_missing = nmiss)
}
## cis window of the old code: flank(GRanges(chr, pos), 500000, both = TRUE) = [pos-5e5, pos+5e5-1]
cis_index <- function(geno, pos) which(geno$start <= pos + FLANK - 1L & geno$end >= max(1L, pos - FLANK))

## ---- tissue inputs (heavy: run inside PBS; prep_pilot.R caches the result) --------------------
## U (the permutation slot set) = tissue activity samples (first row per donor) with complete PC1-5
## and PEER 1..k_now -- identical to the BED sample set of QTLtools_linear_compare/build_pheno_cov.py.
build_tissue <- function(tissue, k_now, act_all, pcs, info) {
  proms <- unique(unlist(lapply(CHRS, function(ch) {
    f <- file.path(PATHS$proms, sprintf("%s.%s_promoters.txt", tissue, ch))
    if (file.exists(f)) readLines(f) else character(0) })))
  sub <- act_all[act_all$SMTSD == tissue, , drop = FALSE]
  n_dup <- sum(duplicated(sub$sample))
  sub <- sub[!duplicated(sub$sample), , drop = FALSE]
  peer <- as.matrix(read.table(file.path(PATHS$peer, paste0("v11_filter_5_peer_factors_", tissue, ".txt")),
                               header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
  ok <- sub$sample %in% rownames(pcs) & sub$sample %in% rownames(peer)
  U <- sub$sample[ok]
  U <- U[complete.cases(pcs[U, , drop = FALSE]) & complete.cases(peer[U, seq_len(k_now), drop = FALSE])]
  miss_col <- setdiff(proms, names(sub))
  A <- as.matrix(sub[match(U, sub$sample), intersect(proms, names(sub)), drop = FALSE])
  rownames(A) <- U
  list(tissue = tissue, k_now = k_now, U = U, act = A,
       pcs = pcs[U, , drop = FALSE], peer = peer[U, , drop = FALSE],
       info = info[colnames(A), c("seqnames", "start", "geneId"), drop = FALSE],
       n_dup_rows = n_dup, missing_promoter_cols = miss_col)
}

## ---- permutation maps --------------------------------------------------------------------------
## arm "real" = identity; "perm0" = the published Fig 2e map; "perm1"/"perm2" = rep maps.
## Every map actually used is stored in rerun_v2/perm_maps (prep_pilot.R copies or draws them).
map_file <- function(tissue, arm) file.path(PATHS$maps, sprintf("%s.%s.perm_map.tsv", tissue, arm))
arm_view <- function(td, arm) {
  if (arm == "real") return(list(act = td$act, peer = td$peer))
  m <- fread(map_file(td$tissue, arm))
  stopifnot(setequal(m$slot, td$U), !anyDuplicated(m$donor), setequal(m$donor, td$U))
  d <- setNames(m$donor, m$slot)[td$U]              # slot s receives donor pi(s)
  A <- td$act[d, , drop = FALSE]; P <- td$peer[d, , drop = FALSE]
  rownames(A) <- rownames(P) <- td$U
  list(act = A, peer = P)
}

## ---- covariate cap (PILOT_PROTOCOL.md section 3) ----------------------------------------------
peer_k <- function(config, k_now, n) {
  switch(config,
    C1 = k_now,
    C2 = max(0L, min(k_now, floor(n / 20) - N_PC)),
    C3 = max(0L, min(k_now, floor(n / 30) - N_PC)),
    stop("unknown config ", config))
}
## design: intercept + standardised PC1-5 + standardised PEER 1..K (rows = td$U)
cov_matrix <- function(td, peer, K) {
  Z <- cbind(td$pcs, peer[, seq_len(K), drop = FALSE])
  Z <- scale(Z)
  if (any(!is.finite(Z))) stop("constant covariate column in ", td$tissue)
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
