## Real-data effect of bohr_control(firth_null_start = TRUE) (stage 3b-1): compare two package runs on the same inputs, the
## first with the paper's engine (BOHR_FIRTH_NULL_START=FALSE, or a run of a commit before 3b-1) and the second with the null
## start. Prints counts and largest differences of the Firth statistics, and every decision that changes (base R only).
## Usage: Rscript 10_compare_null_start.R cis   <old.rds> <new.rds> [summary.rds]   (scripts 02: per-variant tests and units)
##        Rscript 10_compare_null_start.R cond  <old.rds> <new.rds> [summary.rds]   (script 05; with BOHR_SAVE_SCANS=1 in both runs,
##                                                                                    every conditional scan is compared too)
##        Rscript 10_compare_null_start.R trait <old.rds> <new.rds> [summary.rds]   (script 08)
options(width = 200)
a <- commandArgs(TRUE); KIND <- a[1]; O <- readRDS(a[2]); N <- readRDS(a[3]); OUTS <- if (length(a) > 3) a[4] else ""
## differences of matched values: counts beyond 1e-9, 1e-6, 1e-3 (absolute) and the largest
dd <- function(x, y, lab) {
  ok <- is.finite(x) & is.finite(y); d <- abs(x[ok] - y[ok])
  data.frame(stat = lab, both = sum(ok), na_mismatch = sum(is.finite(x) != is.finite(y)), identical = sum(d == 0),
             gt_1e9 = sum(d > 1e-9), gt_1e6 = sum(d > 1e-6), gt_1e3 = sum(d > 1e-3), max_abs = if (any(ok)) max(d) else NA_real_)
}
key_match <- function(A, B, keys) { ka <- do.call(paste, c(A[keys], sep = "\r")); kb <- do.call(paste, c(B[keys], sep = "\r"))
  list(a = A[ka %in% kb, , drop = FALSE][order(ka[ka %in% kb]), , drop = FALSE], b = B[kb %in% ka, , drop = FALSE][order(kb[kb %in% ka]), , drop = FALSE],
       only_a = sum(!ka %in% kb), only_b = sum(!kb %in% ka)) }
out <- list()
if (KIND == "cis") {
  m <- key_match(O$snp, N$snp, c("promoter", "variant")); A <- m$a; B <- m$b
  cat(sprintf("variant tests: old %d, new %d, matched %d (only old %d, only new %d)\n", nrow(O$snp), nrow(N$snp), nrow(A), m$only_a, m$only_b))
  tab <- rbind(dd(A$log_p_p0, B$log_p_p0, "p0 log P (Firth)"), dd(A$log_p_p1, B$log_p_p1, "p1 log P (Firth)"),
               dd(A$stat_p0, B$stat_p0, "p0 LRT statistic"), dd(A$stat_p1, B$stat_p1, "p1 LRT statistic"),
               dd(A$est_p0, B$est_p0, "p0 estimate"), dd(A$est_p1, B$est_p1, "p1 estimate"), dd(A$se_p0, B$se_p0, "p0 SE"),
               dd(A$log_p_mu, B$log_p_mu, "mu log P (F)"), dd(A$effect_p0, B$effect_p0, "p0 effect"))
  print(tab, digits = 3, row.names = FALSE)
  for (cc in c("p0", "p1")) {
    d <- abs(A[[paste0("log_p_", cc)]] - B[[paste0("log_p_", cc)]]); i <- which(is.finite(d) & d > 1e-6)
    if (length(i)) { i <- i[order(-d[i])]
      cat(sprintf("%s tests whose log P changed by more than 1e-6 (%d), largest first:\n", cc, length(i)))
      print(data.frame(promoter = A$promoter[i], variant = A$variant[i], log10P_old = A[[paste0("log_p_", cc)]][i] / log(10),
                       log10P_new = B[[paste0("log_p_", cc)]][i] / log(10), z_old = A[[paste0("est_", cc)]][i] / A[[paste0("se_", cc)]][i],
                       z_new = B[[paste0("est_", cc)]][i] / B[[paste0("se_", cc)]][i])[seq_len(min(15, length(i))), ], digits = 4, row.names = FALSE) }
  }
  mu <- key_match(O$units, N$units, c("promoter", "component")); UA <- mu$a; UB <- mu$b
  cat(sprintf("units: old %d, new %d, matched %d; same leader %d; same n_variants %d\n", nrow(O$units), nrow(N$units), nrow(UA),
              sum(UA$leader == UB$leader), sum(UA$n_variants == UB$n_variants)))
  print(rbind(dd(log(UA$p), log(UB$p), "unit log P (ACAT)"), dd(UA$q, UB$q, "q (this run only)")), digits = 3, row.names = FALSE)
  ch <- which(UA$leader != UB$leader | abs(log(UA$p) - log(UB$p)) > 1e-6 | UA$significant != UB$significant)
  if (length(ch)) { cat("units whose leader, ACAT P (> 1e-6 in log) or significance changed:\n")
    print(data.frame(promoter = UA$promoter[ch], component = UA$component[ch], p_old = signif(UA$p[ch], 4), p_new = signif(UB$p[ch], 4),
                     leader_same = UA$leader[ch] == UB$leader[ch], sig_old = UA$significant[ch], sig_new = UB$significant[ch]), row.names = FALSE) }
  out <- list(snp = tab, units_changed = if (length(ch)) UA[ch, ] else NULL)
} else if (KIND == "cond") {
  ku <- function(U) paste(U$promoter, U$component)
  cat(sprintf("units: old %d, new %d; forward signals old %d, new %d; kept old %d, new %d\n", nrow(O$units), nrow(N$units), nrow(O$signals), nrow(N$signals),
              sum(O$signals$kept), sum(N$signals$kept)))
  m <- key_match(O$units, N$units, c("promoter", "component")); UA <- m$a; UB <- m$b
  ch <- which(UA$n_signals != UB$n_signals | UA$n_kept != UB$n_kept | UA$stop != UB$stop)
  cat(sprintf("units matched %d: same forward count %d, same kept count %d, same stop reason %d\n", nrow(UA), sum(UA$n_signals == UB$n_signals),
              sum(UA$n_kept == UB$n_kept), sum(UA$stop == UB$stop)))
  if (length(ch)) { cat("units that changed:\n"); print(data.frame(promoter = UA$promoter[ch], component = UA$component[ch], forward_old = UA$n_signals[ch],
    forward_new = UB$n_signals[ch], kept_old = UA$n_kept[ch], kept_new = UB$n_kept[ch], stop_old = UA$stop[ch], stop_new = UB$stop[ch]), row.names = FALSE) }
  s <- key_match(O$signals, N$signals, c("promoter", "component", "k")); SA <- s$a; SB <- s$b
  cat(sprintf("signals matched %d (only old %d, only new %d); same variant %d; same backward call %d\n", nrow(SA), s$only_a, s$only_b,
              sum(SA$variant == SB$variant), sum(SA$kept == SB$kept)))
  print(rbind(dd(SA$log_p, SB$log_p, "log P at selection"), dd(log(SA$acat_p), log(SB$acat_p), "log window ACAT at selection"),
              dd(SA$loo_log_p, SB$loo_log_p, "LOO log P"), dd(log(SA$loo_acat_p), log(SB$loo_acat_p), "log LOO window ACAT")), digits = 3, row.names = FALSE)
  r <- key_match(O$rounds, N$rounds, c("promoter", "component", "round")); RA <- r$a; RB <- r$b
  rr <- abs(RA$acat_p - RB$acat_p) / RA$acat_p
  cat(sprintf("rounds matched %d (only old %d, only new %d); same best variant %d; window ACAT P: relative difference > 1e-6 in %d, max %.3g\n",
              nrow(RA), r$only_a, r$only_b, sum(RA$best_variant == RB$best_variant), sum(rr > 1e-6, na.rm = TRUE), max(rr, na.rm = TRUE)))
  i <- which(rr > 1e-6)
  if (length(i)) print(data.frame(promoter = RA$promoter[i], component = RA$component[i], round = RA$round[i], acat_old = signif(RA$acat_p[i], 4),
                                  acat_new = signif(RB$acat_p[i], 4), rel = signif(rr[i], 3), best_same = RA$best_variant[i] == RB$best_variant[i],
                                  threshold = signif(N$cis_threshold, 4))[order(-rr[i]), ], row.names = FALSE)
  if (!is.null(O$scans) && !is.null(N$scans)) {
    flat <- function(S) do.call(rbind, lapply(seq_along(S), function(pi) do.call(rbind, lapply(S[[pi]], function(j) if (!is.null(j$res) && nrow(j$res))
      data.frame(key = j$key, type = j$type, cs = paste(j$cs, collapse = ","), variant = j$res$variant, est = j$res$est, se = j$res$se,
                 log_p = j$res$log_p, stringsAsFactors = FALSE)))))
    ## a scan is identified by the unit and its conditioning variants; a leave-one-out scan of a two-signal unit repeats the
    ## round-1 scan (conditioning on the leader), so each distinct scan is counted once
    uniq <- function(F) F[!duplicated(paste(F$key, F$type, F$cs, F$variant, sep = "\r")), , drop = FALSE]
    FA <- uniq(flat(O$scans)); FB <- uniq(flat(N$scans))
    q <- key_match(FA, FB, c("key", "type", "cs", "variant")); A <- q$a; B <- q$b
    ns <- function(F) length(unique(paste(F$key, F$type, F$cs, sep = "\r")))
    cat(sprintf("conditional scans (distinct forward-round and leave-one-out scans): old %d scans, %d tests; new %d scans, %d tests; matched %d tests (only old %d, only new %d)\n",
                ns(FA), nrow(FA), ns(FB), nrow(FB), nrow(A), q$only_a, q$only_b))
    for (cc in c("mu", "p0", "p1")) { w <- A$type == cc
      if (any(w)) print(rbind(dd(A$log_p[w], B$log_p[w], paste(cc, "conditional log P")), dd(A$est[w], B$est[w], paste(cc, "conditional estimate"))), digits = 3, row.names = FALSE) }
    d <- abs(A$log_p - B$log_p); i <- which(is.finite(d) & d > 1e-6)
    if (length(i)) { i <- i[order(-d[i])]; cat("conditional tests whose log P changed by more than 1e-6, largest first:\n")
      print(data.frame(promoter = A$key[i], type = A$type[i], conditioned_on = nchar(A$cs[i]) > 0, variant = A$variant[i],
                       log10P_old = A$log_p[i] / log(10), log10P_new = B$log_p[i] / log(10), z_old = A$est[i] / A$se[i], z_new = B$est[i] / B$se[i])[seq_len(min(20, length(i))), ],
            digits = 4, row.names = FALSE) }
    out$scans <- list(n = nrow(A), changed = if (length(i)) cbind(A[i, ], log_p_new = B$log_p[i], est_new = B$est[i], se_new = B$se[i]) else NULL)
  }
  out$units_changed <- if (length(ch)) UA[ch, ] else NULL
} else if (KIND == "trait") {
  m <- key_match(O$results, N$results, intersect(c("tissue", "trait", "promoter", "component"), names(O$results))); A <- m$a; B <- m$b
  cat(sprintf("tests: old %d, new %d, matched %d\n", nrow(O$results), nrow(N$results), nrow(A)))
  for (cc in c("p0", "p1", "mu")) { w <- A$component == cc
    print(rbind(dd(A$log_p[w], B$log_p[w], paste(cc, "log P")), dd(A$est[w], B$est[w], paste(cc, "estimate")), dd(A$q[w], B$q[w], paste(cc, "q"))),
          digits = 3, row.names = FALSE) }
  d <- abs(A$log_p - B$log_p); i <- which(is.finite(d) & d > 1e-6)
  if (length(i)) { i <- i[order(-d[i])]; cat("tests whose log P changed by more than 1e-6, largest first:\n")
    print(data.frame(promoter = A$promoter[i], component = A$component[i], p_old = signif(A$p[i], 4), p_new = signif(B$p[i], 4),
                     sig_old = A$significant[i], sig_new = B$significant[i])[seq_len(min(20, length(i))), ], row.names = FALSE) }
  ch <- which(A$significant != B$significant | A$eligible != B$eligible)
  cat(sprintf("decisions: significant old %d, new %d; eligible old %d, new %d; changed calls %d\n", sum(A$significant), sum(B$significant),
              sum(A$eligible), sum(B$eligible), length(ch)))
  if (length(ch)) print(data.frame(promoter = A$promoter[ch], component = A$component[ch], q_old = signif(A$q[ch], 4), q_new = signif(B$q[ch], 4),
                                   sig_old = A$significant[ch], sig_new = B$significant[ch]), row.names = FALSE)
  out <- list(changed = if (length(i)) A[i, ] else NULL, calls_changed = if (length(ch)) A[ch, ] else NULL)
} else stop("kind must be cis, cond or trait")
if (nzchar(OUTS)) saveRDS(out, OUTS)
