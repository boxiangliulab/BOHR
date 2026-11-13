## run_paper.sh 04.1 and 04.2 (authored for the repository): re-run tissue-specificity filter and donor-characteristic
## tests against the stored outputs.
## Usage: Rscript check_04.R roku   <new folder> <stored folder (fit/roku_20261007)> <status.tsv>
##        Rscript check_04.R traits <new results/real> <stored results/real> <status.tsv> <tissue> [<tissue> ...]
a <- commandArgs(TRUE)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
if (a[1] == "roku") {
  f <- "v11_tissue_specific_ROKU_results_delta010.rds"
  N <- readRDS(file.path(a[2], f)); O <- readRDS(file.path(a[3], f))
  add("ROKU calls kept by the 0.10 effect filter (true_outliers)", identical(N$true_outliers, O$true_outliers),
      sprintf("%d promoters x %d tissues; identical: %s; calls kept %d (stored %d)", nrow(O$true_outliers), ncol(O$true_outliers),
              identical(N$true_outliers, O$true_outliers), sum(N$true_outliers != 0, na.rm = TRUE), sum(O$true_outliers != 0, na.rm = TRUE)))
  d <- max(abs(N$delta_matrix - O$delta_matrix), na.rm = TRUE)
  add("differences from the reference median (delta_matrix)", identical(is.na(N$delta_matrix), is.na(O$delta_matrix)) && d <= 1e-12,
      sprintf("same missing cells: %s; max |difference| %.3g", identical(is.na(N$delta_matrix), is.na(O$delta_matrix)), d))
  r <- compare_tables(fread(file.path(a[2], "summary_delta_filter.tsv")), fread(file.path(a[3], "summary_delta_filter.tsv")), "set", tol = 0)
  add("summary of the filter (summary_delta_filter.tsv)", r$ok, r$detail)
} else {
  for (t in a[-(1:4)]) for (f in list.files(file.path(a[2], t), "\\.rds$")) {
    N <- readRDS(file.path(a[2], t, f)); O <- readRDS(file.path(a[3], t, f))
    r <- compare_tables(as_tab(N), as_tab(O), "promoter", tol = 1e-6)
    sk <- identical(N$meta$skip_trait, O$meta$skip_trait)
    add(sprintf("%s, %s: per-promoter tests (rules, statistics, estimates)", t, sub("\\.rds$", "", f)), r$ok && sk,
        paste0(r$detail, if (!sk) "; the trait-level skip reason differs" else ""))
  }
  ## the published significant tests of these tissues and characteristics: P values and probability-scale effects (for p0
  ## and p1 the effect comes from the unpenalised fit, also where it is separated; it converges to a limit there)
  trg <- dirname(dirname(dirname(a[3])))
  E <- fr(file.path(trg, "trait_regression_20261106_mirror_abs/results/export_B/trait_regression_units.tsv.gz"))[sig_final == TRUE & tissue %in% a[-(1:4)]]
  rows <- rbindlist(lapply(seq_len(nrow(E)), function(i) { e <- E[i]; f <- file.path(a[2], e$tissue, paste0(e$trait, ".rds")); if (!file.exists(f)) return(NULL)
    r <- as_tab(readRDS(f))[promoter == e$promoter]
    pc <- c(mu = "mu_F_logp", p0 = "nu_firth_logp", p1 = "tau_firth_logp")[[e$type]]; ec <- c(mu = "mu_pe", p0 = "p0_pe", p1 = "p1_pe")[[e$type]]
    data.table(dp = abs(r[[pc]] - e$logp), de = abs(r[[ec]] - e$pe)) }))
  if (nrow(rows)) add("published significant tests of these tissues and characteristics (P and probability-scale effect)",
                      all(rows$dp <= 1e-6 * pmax(1, abs(E$logp[seq_len(nrow(rows))])) & rows$de <= 1e-4),
                      sprintf("%d tests; max |difference| of log P %.3g, of the effect %.3g", nrow(rows), max(rows$dp), max(rows$de)))
}
finish(a[4])
