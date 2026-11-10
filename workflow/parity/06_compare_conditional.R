## Real-data parity of bohr_conditional(), step C (original R + data.table): compare with the stored v3c tables
## (results/conditional_final_v3c, the published 14,800 signals) and, for the two refitted units, with v3m.
## Prints counts, largest differences and agreement rates only (no sample or donor IDs).
## Usage: Rscript 06_compare_conditional.R <run.rds> [summary.rds]
options(width = 200)
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); R <- readRDS(a[1]); OUTS <- if (length(a) > 1) a[2] else ""
V2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"
z <- function(f) fread(cmd = paste("zcat", f))
TIS <- R$tissue
for (ver in c("v3c", "v3m")) {
  F <- file.path(V2, "results", paste0("conditional_final_", ver))
  SU <- z(file.path(F, "freeze_real_units.tsv.gz"))[tissue == TIS]; SS <- z(file.path(F, "freeze_real_signals.tsv.gz"))[tissue == TIS]
  SR <- z(file.path(F, "freeze_real_rounds.tsv.gz"))[tissue == TIS]
  PU <- as.data.table(R$units); PS <- as.data.table(R$signals); PR <- as.data.table(R$rounds)
  setnames(PU, "component", "type"); setnames(PS, "component", "type"); setnames(PR, "component", "type")
  SU <- SU[paste(promoter, type) %in% paste(PU$promoter, PU$type)]
  SS <- SS[paste(promoter, type) %in% paste(PU$promoter, PU$type)]; SR <- SR[paste(promoter, type) %in% paste(PU$promoter, PU$type)]
  cat(sprintf("\n==== %s vs stored %s: units package %d, stored %d\n", TIS, ver, nrow(PU), nrow(SU)))
  M <- merge(SU[, .(promoter, type, ns = n_signals, nk = n_signals_backward, st = stop, source)], PU[, .(promoter, type, ns_pk = n_signals, nk_pk = n_kept, st_pk = stop)],
             by = c("promoter", "type"), all = TRUE)
  cat(sprintf("units: matched %d; forward signal count equal %d; kept count equal %d; stop reason equal %d\n", M[!is.na(ns) & !is.na(ns_pk), .N],
              M[ns == ns_pk, .N], M[nk == nk_pk, .N], M[st == st_pk, .N]))
  if (M[ns != ns_pk | nk != nk_pk, .N]) { cat("units that differ (forward/kept, stored vs package; source):\n"); print(M[ns != ns_pk | nk != nk_pk, .(type, ns, ns_pk, nk, nk_pk, source)]) }
  K <- merge(SS[, .(promoter, type, k, snp, sr = selected_round, lp = logp_at_selection, lpc = logp_chi_at_selection, ac = acat_at_selection, le = loo_est, ls = loo_se,
                    ll = loo_logp, la = loo_acat, lt = loo_top_snp, bk = backward_keep)],
             PS[, .(promoter, type, k, snp_pk = variant, sr_pk = selected_round, lp_pk = log_p, lpc_pk = log_p_chi, ac_pk = acat_p, le_pk = loo_est, ls_pk = loo_se,
                    ll_pk = loo_log_p, la_pk = loo_acat_p, lt_pk = loo_top_variant, bk_pk = kept)], by = c("promoter", "type", "k"), all = TRUE)
  cat(sprintf("signals: stored %d (kept %d), package %d (kept %d); matched (promoter, type, k) %d; same lead variant %d; same round %d; same backward call %d; same LOO top variant (NA = NA) %d\n",
              nrow(SS), SS[backward_keep %in% TRUE, .N], nrow(PS), PS[kept %in% TRUE, .N], K[!is.na(snp) & !is.na(snp_pk), .N], K[snp == snp_pk, .N],
              K[sr == sr_pk, .N], K[bk == bk_pk, .N], K[(is.na(lt) & is.na(lt_pk)) | lt %in% lt_pk, .N]))
  dd <- function(a, b, lab, log10 = FALSE) { ok <- is.finite(a) & is.finite(b); d <- abs(a[ok] - b[ok]); if (log10) d <- d / log(10)
    data.table(stat = lab, both = sum(ok), na_mismatch = sum(is.finite(a) != is.finite(b)), max_abs_diff = if (any(ok)) max(d) else NA_real_,
               max_rel_diff = if (any(ok)) max(d / pmax(abs(a[ok]), 1e-300)) else NA_real_) }
  K2 <- K[snp == snp_pk]
  print(rbind(dd(K2$lp, K2$lp_pk, "log P at selection"), dd(K2$lpc, K2$lpc_pk, "mu chi2 log P at selection"), dd(K2$ac, K2$ac_pk, "window ACAT P at selection"),
              dd(K2$ll, K2$ll_pk, "LOO log P"), dd(K2$le, K2$le_pk, "LOO estimate"), dd(K2$ls, K2$ls_pk, "LOO SE"), dd(K2$la, K2$la_pk, "LOO window ACAT P")), digits = 3)
  RR <- merge(SR[, .(promoter, type, round, nt = n_tested, ac = acat_p, bs = best_snp, bl = best_logp)], PR[, .(promoter, type, round, nt_pk = n_tested, ac_pk = acat_p,
              bs_pk = best_variant, bl_pk = best_log_p)], by = c("promoter", "type", "round"), all = TRUE)
  cat(sprintf("rounds: stored %d, package %d, matched %d; same variants tested %d; same best variant %d; max rel. diff of window ACAT P %.3g\n",
              nrow(SR), nrow(PR), RR[!is.na(ac) & !is.na(ac_pk), .N], RR[nt == nt_pk, .N], RR[bs == bs_pk, .N], RR[is.finite(ac) & is.finite(ac_pk), max(abs(ac - ac_pk) / ac)]))
  cat(sprintf("decision margins: smallest |log(ACAT/threshold)| over the package's rounds %.3g; over the package's LOO scans %.3g\n",
              PR[is.finite(acat_p), min(abs(log(acat_p / R$cis_threshold)))], if (PS[is.finite(loo_acat_p), .N]) PS[is.finite(loo_acat_p), min(abs(log(loo_acat_p / R$cis_threshold)))] else NA))
  if (nzchar(OUTS) && ver == "v3c") saveRDS(list(units = M, signals = K, rounds = RR), OUTS)
}
cat(sprintf("\npackage run: wall %.0f s, CPU %.0f s\n", R$wall, R$cpu))
