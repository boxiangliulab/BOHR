## Subset test, checks (authored for the repository): compares the outputs of run_subset.sh and of run_paper.sh sub-steps
## 01.2, 03.1 and 03.2 with the stored results of the paper and prints PASS or FAIL per check. Prints aggregate numbers only.
## Usage: Rscript check_subset.R <OUT> <RUN_ROOT> <chr> <tissue> [<tissue> ...]
## Environment: CHECK_PARTS selects parts 1 (promoter filter), 2 (cis scan) and 3 (signals), default all; with CHECK_OUT the
## rows are appended to that file (run_paper.sh), otherwise they go to <OUT>/summary.tsv and are printed.
suppressMessages(library(data.table))
a <- commandArgs(TRUE); OUT <- a[1]; RUN <- a[2]; CHR <- a[3]; TIS <- a[-(1:3)]
v11 <- file.path(Sys.getenv("GTEX", "/home/e1327989/project/GTEx"), "v11")
res <- list()
add <- function(check, ok, detail) res[[length(res) + 1L]] <<- data.table(check = check, result = if (isTRUE(ok)) "PASS" else "FAIL", detail = detail)
TOL <- 1e-8
PARTS <- strsplit(Sys.getenv("CHECK_PARTS", "1 2 3"), " +")[[1]]

## 1. promoter filter
f1 <- file.path(OUT, "step1_relative_activity.rds")
if (!"1" %in% PARTS) {} else if (!file.exists(f1)) add("1 promoter filter", FALSE, "no output (see logs/01.2.log)") else {
  s1 <- readRDS(f1)
  load(file.path(v11, "Promoter_info_filter_5.RData"))      # info: the paper's final promoter set
  want <- rownames(info)[as.character(info$seqnames) == CHR]
  got <- rownames(s1$relative)
  add(sprintf("1a promoter set (%s)", CHR), setequal(got, want),
      sprintf("%d kept by the package, %d in the paper's set; %d extra, %d missing", length(got), length(want),
              length(setdiff(got, want)), length(setdiff(want, got))))
  load(file.path(v11, "filter_5_recal_att_filter.RData"))   # recal_att_f: samples x promoters, the paper's values
  ids <- intersect(got, colnames(recal_att_f)); ss <- intersect(colnames(s1$relative), rownames(recal_att_f))
  X <- s1$relative[ids, ss, drop = FALSE]; Y <- t(as.matrix(recal_att_f[ss, ids, drop = FALSE])); Y[is.nan(Y)] <- NA
  d <- suppressWarnings(max(abs(X - Y), na.rm = TRUE))
  add("1b relative activity values", identical(is.na(X), is.na(Y)) && d <= TOL,
      sprintf("%d promoters x %d samples; same missing values: %s; max |difference| %.3g", length(ids), length(ss),
              identical(is.na(X), is.na(Y)), d))
}

## 2. cis mapping: per-promoter rules and per-SNP statistics
num <- c("mu_lrt_logp", "nu_lrt_logp", "tau_lrt_logp", "nu_firth_logp", "tau_firth_logp", "mu_est", "nu_est", "tau_est",
         "nu_firth_est", "tau_firth_est")
pcols <- c("family", "skip", "n", "n0", "n1", "n_cov", "elig_mu", "elig_p0", "elig_p1", "n_snp_window", "n_snp_tested")
for (t in if ("2" %in% PARTS) TIS else character()) {
  fn <- file.path(OUT, "shadow", "results", "final", "freeze", t, "real", paste0(CHR, ".rds"))
  fo <- file.path(RUN, "results", "final", "freeze", t, "real", paste0(CHR, ".rds"))
  if (!file.exists(fn)) { add(sprintf("2 cis scan %s", t), FALSE, "no output (see logs/03.1.log)"); next }
  N <- readRDS(fn); O <- readRDS(fo)
  pn <- as.data.table(N$prom)[order(promoter)]; po <- as.data.table(O$prom)[order(promoter)]
  pc <- intersect(pcols, intersect(names(pn), names(po)))
  same_p <- identical(pn$promoter, po$promoter) && isTRUE(all.equal(as.data.frame(pn[, ..pc]), as.data.frame(po[, ..pc]), check.attributes = FALSE))
  add(sprintf("2a promoters and rules, %s %s", t, CHR), same_p, sprintf("%d promoters; columns compared: %s", nrow(pn), paste(pc, collapse = ", ")))
  sn <- as.data.table(N$snp); so <- as.data.table(O$snp)
  k <- c("promoter", "snp"); setkeyv(sn, k); setkeyv(so, k)
  same_rows <- nrow(sn) == nrow(so) && identical(sn$promoter, so$promoter) && identical(sn$snp, so$snp)
  cols <- intersect(num, intersect(names(sn), names(so)))
  ## an unpenalised estimate of nu or tau under complete separation (sep TRUE in both runs) diverges, and the value where
  ## the optimiser stops is arbitrary; such SNPs use the Firth fit, whose estimate and P value are compared as usual
  divg <- function(cl) { p <- sub("_est$", "", cl); s <- paste0(p, "_sep")
    if (cl %in% c("nu_est", "tau_est") && s %in% names(sn) && s %in% names(so)) (sn[[s]] %in% TRUE & so[[s]] %in% TRUE) else rep(FALSE, nrow(sn)) }
  dmax <- if (same_rows) max(vapply(cols, function(cl) suppressWarnings(max(abs(sn[[cl]] - so[[cl]])[!divg(cl)], na.rm = TRUE)), 0)) else NA
  n_divg <- if (same_rows) sum(divg("nu_est")) + sum(divg("tau_est")) else NA
  na_same <- same_rows && all(vapply(cols, function(cl) identical(is.na(sn[[cl]]), is.na(so[[cl]])), TRUE))
  add(sprintf("2b per-SNP statistics, %s %s", t, CHR), same_rows && na_same && dmax <= 1e-6,
      sprintf("%d SNP tests (stored %d); same missing values: %s; max |difference| over %d columns %.3g (%d unpenalised estimates under complete separation not compared)",
              nrow(sn), nrow(so), na_same, length(cols), dmax, n_divg))
}

## 3. independent signals of the re-run units, and 4. against the final signal table
if ("3" %in% PARTS) fin <- fread(cmd = paste("gzip -dc", shQuote(file.path(RUN, "results", "conditional_final_v3c", "freeze_real_signals.tsv.gz"))))   # no R.utils needed
for (t in if ("3" %in% PARTS) TIS else character()) {
  fn <- file.path(OUT, "shadow", "results", "conditional_v3m", "freeze", t, "real", paste0(CHR, ".rds"))
  fo <- file.path(RUN, "results", "conditional_v3m", "freeze", t, "real", paste0(CHR, ".rds"))
  if (!file.exists(fn)) { add(sprintf("3 conditional %s", t), FALSE, "no output (see logs/03.2.log)"); next }
  N <- readRDS(fn); O <- readRDS(fo)
  key <- function(s) if (nrow(s)) s[order(promoter, type, k), paste(promoter, type, k, snp)] else character(0)
  sn <- as.data.table(N$signals); so <- as.data.table(O$signals)
  d <- if (nrow(sn) && nrow(so) && identical(key(sn), key(so)))
         max(abs(sn[order(promoter, type, k)]$logp_at_selection - so[order(promoter, type, k)]$logp_at_selection), na.rm = TRUE) else NA
  add(sprintf("3 signals of the re-run units, %s %s", t, CHR), identical(key(sn), key(so)) && (is.na(d) && !nrow(sn) || isTRUE(d <= 1e-6)),
      sprintf("%d units, %d signals (stored %d); max |difference| of log P at selection %.3g", nrow(N$units), nrow(sn), nrow(so), d))
  if (nrow(sn)) {
    u <- unique(sn[, .(promoter, type)])
    ff <- merge(fin[tissue == t], u, by = c("promoter", "type"))
    kept <- if ("backward_keep" %in% names(ff)) ff[backward_keep == TRUE] else ff
    add(sprintf("4 final signal table, %s %s", t, CHR), setequal(kept$snp, sn$snp) || all(kept$snp %in% sn$snp),
        sprintf("%d signals of these units in the paper's final table (14,800 signals), all among the re-run signals: %s",
                nrow(kept), all(kept$snp %in% sn$snp)))
  }
}
R <- rbindlist(res)
if (nzchar(Sys.getenv("CHECK_OUT"))) {
  fwrite(R, Sys.getenv("CHECK_OUT"), sep = "\t", append = file.exists(Sys.getenv("CHECK_OUT"))); print(R); quit(save = "no")
}
fwrite(R, file.path(OUT, "summary.tsv"), sep = "\t")
w <- max(nchar(R$check))
for (i in seq_len(nrow(R))) cat(sprintf("%-*s  %s  %s\n", w, R$check[i], R$result[i], R$detail[i]))
cat(sprintf("\n%d of %d checks passed\n", sum(R$result == "PASS"), nrow(R)))
