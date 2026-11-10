## Real-data parity of bohr_conditional(), step B (bohrdev R with BOHR installed).
## Mode "stored": inputs from 04_extract_conditional_inputs.R; the `cis` input is built from the stored units (all chromosomes).
## Mode "cis": inputs from 01_extract_inputs.R and the bohr_cis() run of 02_run_bohr_cis.R on one chromosome; its q values and
##   threshold are recomputed over the tissue (stored units of the other chromosomes + this run's units, rule units_v3m_mirror.R:44),
##   which is what bohr_combine() gives when every chromosome is run with the package.
## Usage: Rscript 05_run_conditional.R stored <inputs.rds> <out.rds> <ncores> [lib]
##        Rscript 05_run_conditional.R cis <inputs.rds> <cis.rds> <out.rds> <ncores> [lib]
a <- commandArgs(TRUE); MODE <- a[1]
if (MODE == "stored") { IN <- a[2]; OUT <- a[3]; NC <- as.integer(a[4]); LIB <- if (length(a) > 4) a[5] else "" } else {
  IN <- a[2]; CIS <- a[3]; OUT <- a[4]; NC <- as.integer(a[5]); LIB <- if (length(a) > 5) a[6] else "" }
if (nzchar(LIB)) library(BOHR, lib.loc = LIB) else library(BOHR)
x <- readRDS(IN)
if (MODE == "stored") {
  cis <- structure(list(snp = x$snp, units = x$units, promoters = NULL, mirror = NULL, threshold = x$threshold,
                        settings = list(window = 5e5, min_genotype_count = 5L, min_n_mid = 80L, boundary_share = c(0.1, 0.9), fdr = 0.05,
                                        acat_cap = TRUE, mirror = TRUE, mirror_tol = 1e-6, control = bohr_control())), class = "bohr_cis")
} else {
  cis <- readRDS(CIS)
  V2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"
  U <- as.data.frame(readRDS(file.path(V2, "results/final_nogate_v3m/freeze", paste0(x$tissue, ".real.units.rds"))))
  other <- U[!U$promoter %in% cis$promoters$promoter, c("promoter", "type", "p", "estimable", "dedup_drop")]
  pk <- data.frame(promoter = cis$units$promoter, type = cis$units$component, p = cis$units$p, estimable = cis$units$estimable,
                   dedup_drop = cis$units$dedup_drop, stringsAsFactors = FALSE)
  allu <- rbind(other, pk); prim <- !allu$dedup_drop & allu$estimable
  q <- rep(NA_real_, nrow(allu)); q[prim] <- p.adjust(allu$p[prim], "BH")
  qp <- q[nrow(other) + seq_len(nrow(pk))]
  cis$units$q <- qp; cis$units$significant <- !is.na(qp) & qp <= 0.05
  cis$threshold <- max(allu$p[!is.na(q) & q <= 0.05])
  cat(sprintf("tissue-wide BH: %d units, threshold %.6g; significant units on this chromosome: %d\n", nrow(allu), cis$threshold, sum(cis$units$significant)))
}
t0 <- proc.time()[["elapsed"]]; c0 <- sum(proc.time()[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE)
res <- bohr_conditional(cis, x$phenotype, x$genotype, x$variants, x$promoters, covariates = x$covariates, ncores = NC,
                        absolute = if (!is.null(x$absolute)) x$absolute else NULL)
res$wall <- proc.time()[["elapsed"]] - t0
res$cpu <- sum(proc.time()[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE) - c0
res$cis_threshold <- cis$threshold; res$tissue <- x$tissue
saveRDS(res, OUT)
print(res); cat(sprintf("wall %.1f s on %d cores; CPU %.0f s\n", res$wall, NC, res$cpu))
