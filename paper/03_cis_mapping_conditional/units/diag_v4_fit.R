## diag_v4_fit.R -- fit diagnostics of the units the v4 delta adds (rule 1 of the admission rules), real arm.
## Per unit type (mu, p0, p1; estimable by the S2 rule from the scan's n, n0, n1: mu n_mid >= 80, p0/p1 boundary share
## 10-90%): the share of SNP tests that converged (mu_conv; nu_firth_conv; tau_firth_conv, the Firth columns taken from
## the main run or else the Firth add-on, as discovery_nogate.R combines them), the share with a finite P, and the share
## of units whose null model converged (null_beta_conv for mu; null_cat_conv for p0/p1).
## Added units = every testable promoter of the delta root (rerun_v4_delta). Reference = the v2-set units of
## chromosomes 20-22 of every tissue (rerun_v2 scans of the same run and engine).
## Only convergence flags and the finiteness of P are read, no association statistic.
## Usage: Rscript diag_v4_fit.R <ncores> [tissue ...]   (NSCC, after the real delta; output $N/v4_diag/fit_*.tsv)
suppressPackageStartupMessages({ library(data.table); library(parallel) })
N <- "/scratch/users/nus/e1327989/nofloor_full"; NC <- as.integer(commandArgs(TRUE)[1])
OUT <- file.path(N, "v4_diag"); dir.create(OUT, showWarnings = FALSE)
FC <- c("nu_firth_logp", "nu_firth_conv", "tau_firth_logp", "tau_firth_conv")
TIS <- sub("[.]rds$", "", list.files(file.path(N, "rerun_v4_delta/work/final/freeze"), "[.]rds$"))
if (length(commandArgs(TRUE)) > 1) TIS <- commandArgs(TRUE)[-1]          # a subset, for a test
one <- function(set, t, c) {
  root <- if (set == "added") "rerun_v4_delta" else "rerun_v2"
  f <- function(r) file.path(N, root, "results", r, "freeze", t, "real", sprintf("chr%d.rds", c))
  if (!file.exists(f("final")) || !file.exists(f("final_firth_addon"))) return(data.table(set = set, tissue = t, chr = c, missing = TRUE))
  M <- readRDS(f("final")); A <- readRDS(f("final_firth_addon"))
  P <- as.data.table(M$prom)[is.na(skip) & !is.na(family)]
  if (!nrow(P)) return(NULL)
  S <- as.data.table(M$snp)[, c("promoter", "snp", "mu_conv", "mu_lrt_logp", FC), with = FALSE]
  for (k in c("mu_conv", FC)) set(S, j = k, value = as.numeric(S[[k]]))
  if (NROW(A$snp)) { ad <- as.data.table(A$snp)[, c("promoter", "snp", FC), with = FALSE]; setnames(ad, FC, paste0("a_", FC))
    S <- merge(S, ad, by = c("promoter", "snp"), all.x = TRUE)
    for (k in FC) S[is.na(get(k)) & !is.na(get(paste0("a_", k))), (k) := as.numeric(get(paste0("a_", k)))]
    S[, paste0("a_", FC) := NULL] }
  U <- rbind(P[, .(promoter, type = "mu", estimable = (n - n0 - n1) >= 80, null_conv = as.numeric(null_beta_conv))],
             P[family %in% c("BEZI", "BEINF"), .(promoter, type = "p0", estimable = n0 / n >= .1 & n0 / n <= .9, null_conv = as.numeric(null_cat_conv))],
             P[family %in% c("BEOI", "BEINF"), .(promoter, type = "p1", estimable = n1 / n >= .1 & n1 / n <= .9, null_conv = as.numeric(null_cat_conv))])
  L <- rbind(S[, .(promoter, type = "mu", conv = mu_conv, fin = is.finite(mu_lrt_logp))],
             S[, .(promoter, type = "p0", conv = nu_firth_conv, fin = is.finite(nu_firth_logp))],
             S[, .(promoter, type = "p1", conv = tau_firth_conv, fin = is.finite(tau_firth_logp))])
  L <- L[, .(snp_tests = .N, snp_conv = sum(conv %in% 1), snp_finite = sum(fin)), by = .(promoter, type)]
  U <- merge(U, L, by = c("promoter", "type"), all.x = TRUE)
  for (k in c("snp_tests", "snp_conv", "snp_finite")) U[is.na(get(k)), (k) := 0L]
  U[, `:=`(set = set, tissue = t, chr = c)] }
jobs <- rbind(CJ(set = "added", t = TIS, c = 1:22), CJ(set = "reference", t = TIS, c = 20:22))
R <- rbindlist(mclapply(seq_len(nrow(jobs)), function(i) tryCatch(one(jobs$set[i], jobs$t[i], jobs$c[i]),
  error = function(e) data.table(set = jobs$set[i], tissue = jobs$t[i], chr = jobs$c[i], error = conditionMessage(e))), mc.cores = NC, mc.preschedule = FALSE), fill = TRUE)
if ("missing" %in% names(R)) { cat("missing scan files:", R[missing %in% TRUE, .N], "\n"); print(head(R[missing %in% TRUE])) }
if ("error" %in% names(R)) { cat("read errors:", R[!is.na(error), .N], "\n"); print(head(R[!is.na(error)])) }
R <- R[!is.na(promoter)]
agg <- function(x) x[, .(units = .N, null_conv = round(mean(null_conv %in% 1), 5), units_any_finite = round(mean(snp_finite > 0), 5),
                         snp_tests = sum(snp_tests), snp_conv = round(sum(snp_conv) / sum(snp_tests), 5), snp_finite = round(sum(snp_finite) / sum(snp_tests), 5),
                         units_snp_conv_lt_99 = sum(snp_conv < 0.99 * snp_tests))]
FD <- R[, agg(.SD), by = .(set, type, estimable)]
FD[, rule1_pass := fifelse(set == "added" & estimable %in% TRUE, snp_conv >= 0.99 & null_conv >= 0.99, NA)]
setorder(FD, set, type, -estimable); fwrite(FD, file.path(OUT, "fit_diagnostics.tsv"), sep = "\t")
fwrite(R[set == "added" & estimable %in% TRUE, agg(.SD), by = .(tissue, type)], file.path(OUT, "fit_by_tissue.tsv"), sep = "\t")
print(FD, nrows = 50); cat("FIT DIAGNOSTICS DONE\n")
