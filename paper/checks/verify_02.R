## run_paper.sh 02.1 and 02.4 (authored for the repository): the stored check of the BOHR engine against gamlss, and the
## stored benchmark summaries against Supplementary Table 2.
## Usage: Rscript verify_02.R gamlss    <run root (rerun_v2)>          <Supplementary Tables folder> <status.tsv>
##        Rscript verify_02.R summaries <figures_bohr2/data_v3_nothr_realjunc> <Supplementary Tables folder> <status.tsv>
a <- commandArgs(TRUE)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
if (a[1] == "gamlss") {
  D <- file.path(a[2], "results/stage1_equiv_run2")
  s <- readLines(file.path(D, "equiv_summary.txt"))[1]
  err <- as.integer(sub(".*errors +([0-9]+).*", "\\1", s))
  add("gamlss check ran without errors (equiv_summary.txt)", isTRUE(err == 0), trimws(s))
  M <- fread(file.path(D, "equiv_models.tsv"))
  for (cl in c("ll_engine", "ll_gamlss")) set(M, j = cl, value = suppressWarnings(as.numeric(M[[cl]])))   # read as text when a cell is not a number
  M[, gamlss_conv := as.logical(gamlss_conv)]
  C <- M[gamlss_conv == TRUE & is.finite(ll_engine) & is.finite(ll_gamlss)]
  C[, d := ll_engine - ll_gamlss]
  add("engine log-likelihood never below a converged gamlss fit by more than 1e-3", C[d < -1e-3, .N] == 0,
      sprintf("%d model fits where gamlss converged: %d within 1e-3, %d with the engine higher (separated or gamlss stopped early), %d with the engine lower",
              nrow(C), C[abs(d) <= 1e-3, .N], C[d > 1e-3, .N], C[d < -1e-3, .N]))
  L <- fread(file.path(D, "equiv_lrt.tsv"))
  add("likelihood-ratio statistics of engine and gamlss (stored table)", nrow(L) > 0, sprintf("%d LRT rows stored", nrow(L)))
} else {
  P <- fread(file.path(a[2], "power.csv")); S2a <- supp(a[3], "S2a")
  names(S2a) <- c("framework", "effect_type", "method", "test", "N", "n", "p05", "p001", "t05", "t001")
  hit <- vapply(seq_len(nrow(S2a)), function(i) { r <- S2a[i]
    P[N == r$N & promoters == r$n & abs(power_cal05 - r$p05) < 1e-9 & abs(power_cal001 - r$p001) < 1e-9, .N] > 0 }, TRUE)
  add("Supplementary Table 2a power values found in the stored power.csv (same N, promoters, power at both alphas)", all(hit),
      sprintf("%d of %d rows of Supplementary Table 2a", sum(hit), nrow(S2a)))
  Th <- fread(file.path(a[2], "thresholds.csv"))
  tv <- unlist(Th[, vapply(Th, is.numeric, TRUE), with = FALSE])
  th <- vapply(seq_len(nrow(S2a)), function(i) any(abs(tv - S2a$t05[i]) < 1e-12) && any(abs(tv - S2a$t001[i]) < 1e-12), TRUE)
  add("Supplementary Table 2a calibrated P thresholds found in the stored thresholds.csv", all(th),
      sprintf("%d of %d rows", sum(th), nrow(S2a)))
}
finish(a[4])
