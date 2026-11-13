## run_conditional_v3m.R (2026-11-04): mirror rule v3 (exact boundary events labelled p0 on the silent promoter; author
## 2026-11-04). Forward conditional + LOO, from the leader, for the units listed in work/cond_v3m/units_to_run.tsv
## (code/cond_v3m_select.R: new p0(partner) units, plus units whose tissue threshold moved enough to continue/truncate),
## with the tissue's v3m threshold (results/final_nogate_v3m). Everything else = run_conditional_twin.R = run_conditional_ng3.R.
## Output: results/conditional_v3m/<def>/<tissue>/real/<chr>.rds
## run_conditional_twin.R (2026-10-02 evening): Runs the units-v3 conditional analysis on the
## DROPPED exact mirror twins (dedup_drop units: mu(partner) for interior_exact, p1/p0(partner) for boundary_exact) of a
## selection of surviving canonical units (work/twin_units_v3_sel.tsv: all multi-signal ones + 500 random single-signal),
## with the tissue's v3 threshold, from the twin's own leader. Everything else is run_conditional_ng3.R.
## Output: results/conditional_twins_v3/<def>/<tissue>/real/<chr>.rds
## run_conditional_ng3.R (2026-10-02 evening): SUPPLEMENT for units v3 (results/final_nogate_v3: p0/p1 minor class >= 10%,
## i.e. 0.10 <= boundary share <= 0.90, on top of mu n_mid >= 80; author). The v3 BH thresholds are looser than v2 in 47
## tissues and equal in 3 (none stricter), so the current forward sequences stay valid prefixes. This run handles:
##   new units (significant under v3, not under v2)            -> full forward run from the leader;
##   continuation units (current result stopped on "window ACAT above threshold" with last-round ACAT <= thr_v3)
##                                                              -> start from the current signal set, keep going.
## "Current result" = the v2 final per unit: the v2 supplement (conditional_nogate_supp) if the unit is there, else the v1
## main run (conditional_nogate). Units unchanged keep that result; v2 units outside v3 are dropped at merge time.
## Also fixes the est-name bug of test_one() (unname the Firth estimate). Real arm only.
## Output: results/conditional_nogate_supp3/<def>/<tissue>/real/<chr>.rds (from_v1 = TRUE marks carried-over signals)
## run_conditional_ng.R (2026-10-01): run_conditional.R applied to the UNGATED, mirror-de-duplicated discovery
## (code/discovery_nogate.R -> results/final_nogate/<def>/<tissue>.<arm>.units.rds; author: no proportion gate;
## mirror table v2). Units = primary-set units (dedup_drop == FALSE) with q <= 0.05 (real) or a 25% sample (perm);
## thr = the tissue's primary BH threshold. Everything else is identical to run_conditional.R.
## Output: results/conditional_nogate/<def>/<tissue>/<arm>/<chr>.rds
## Stage 3: forward conditional analysis (independent signals) per significant unit (promoter x type) of one
## (definition, tissue, arm, chromosome). Choices confirmed by the author on 2026-10-01:
##   * own-type test as at discovery: mu -> mu LRT with the small-sample F correction (amendment 7, adopted 2026-10-01:
##     F = (n_mid - p) * expm1(W / n_mid) ~ F(1, n_mid - p), p = rank of the interior-sample mean design incl. the SNP,
##     phi not counted; leaders raise p by one each when linearly independent); p0 -> nu Firth PLR; p1 -> tau Firth PLR;
##   * plan B: leaders are appended to the covariate matrix, so they enter the beta mean AND both the nu and tau
##     equations of the categorical part (a genotype can move P(0) and P(1) at the same time);
##   * full window: every SNP the main run tested, no marginal-p candidate filter;
##   * stop when the window ACAT of the conditional p exceeds the tissue's discovery threshold (the largest unit p with
##     BH q <= 0.05 in that tissue's real arm); safety cap 10 signals;
##   * afterwards, for units with >= 2 signals, each signal's leave-one-out window statistics (all SNPs conditional on
##     the other signals) for per-signal coloc.
## Perm mode (arm != real): pseudo-leader null. A random share COND_PERM_FRAC (default 0.25) of the arm's eligible
## units is conditioned on its marginal top SNP and tested once; the share of round-1 window ACAT p <= the REAL arm's
## threshold estimates the false secondary-signal rate.
## Inputs are rebuilt with the same code as run_final2.R; the main run's files are only read.
## Usage: Rscript run_conditional.R <def> <tissue> <arm> <chr> <ncores>
## Output (atomic): results/conditional/<def>/<tissue>/<arm>/<chr>.rds = list(units, signals, rounds, loo, meta).
## Decisions use `logp` (mu: F-corrected; p0/p1: Firth); for mu the uncorrected chi2 value is kept alongside as `logp_chi`
## (per SNP, per signal, leave-one-out) and `acat_p_chi` per round, so the correction can be revised later.
suppressPackageStartupMessages(library(parallel))
CODE <- Sys.getenv("BOHR2_CODE", "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/code")
source(file.path(CODE, "bohr2_data.R")); source(file.path(CODE, "bohr2_engine.R")); source(file.path(CODE, "bohr2_addon3.R"))
a <- commandArgs(TRUE); def <- a[1]; tis <- a[2]; arm <- a[3]; chr <- a[4]; NC <- as.integer(a[5])
DEF <- def; TIS <- tis; ARM <- arm; CHR <- chr                 # filter values; never equal to a column name
MAXSIG <- 10L; PERM <- arm != "real"; FRAC <- as.numeric(Sys.getenv("COND_PERM_FRAC", "0.25"))
out_dir <- file.path(Sys.getenv("COND_ROOT", file.path(V2, "results", "conditional_v3m")), def, tis, arm)
out <- file.path(out_dir, paste0(chr, ".rds"))
if (file.exists(out)) { cat("exists:", out, "\n"); quit(save = "no") }
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
t_start <- Sys.time()
acat_log <- function(logp) { p <- exp(logp[is.finite(logp)]); if (!length(p)) return(NA_real_); if (any(p == 0)) return(0)
  if (length(p) > 1) p <- pmin(p, 1 - 1 / length(p))
  sm <- p < 1e-15; t <- numeric(length(p)); t[!sm] <- tan((0.5 - p[!sm]) * pi); t[sm] <- 1 / (p[sm] * pi); st <- mean(t)
  if (st > 1e15) return(1 / (st * pi)); pcauchy(st, lower.tail = FALSE) }

## ---- units and threshold ----
td <- readRDS(file.path(V2, "work", "final", def, paste0(tis, ".rds"))); K <- td$K
RULE <- "firthF_nogate_dedup_nmid80_p01minor10_mirrorv3"
stopifnot(ARM == "real")
UR2 <- as.data.table(readRDS(file.path(V2, "results", "final_nogate_v3m", DEF, paste(TIS, "real", "units.rds", sep = "."))))
thr <- UR2[q <= 0.05, if (.N) max(p) else NA_real_]           # the tissue's v3m primary BH threshold
TW <- fread(Sys.getenv("V3M_LIST", file.path(V2, "work", "cond_v3m", "units_to_run.tsv")))[tissue == TIS & chr == CHR]
U <- UR2[paste(promoter, type) %in% paste(TW$promoter, TW$type) & q <= 0.05][, mode := "new"]
stopifnot(nrow(U) == nrow(TW))
s1 <- data.table()
main_f <- file.path(V2, "results", "final", def, tis, arm, paste0(chr, ".rds"))
empty <- function(why) { meta <- list(def = def, tissue = tis, arm = arm, chr = chr, K = K, thr = thr, n_units = 0L, note = why, incomplete = 0L,
    wall_sec = as.numeric(difftime(Sys.time(), t_start, units = "secs")), date = format(Sys.time()))
  tmp <- paste0(out, ".tmp"); saveRDS(list(units = data.table(), signals = data.table(), rounds = data.table(), loo = data.table(), meta = meta), tmp)
  file.rename(tmp, out); cat(def, tis, arm, chr, ":", why, "\n"); quit(save = "no") }
if (!is.finite(thr)) empty("no significant unit in the tissue's real arm")
if (!nrow(U)) empty("no unit on this chromosome")
M <- readRDS(main_f); stopifnot(identical(M$meta$def, def), identical(M$meta$tissue, tis), identical(M$meta$arm, arm), K == M$meta$K)
snp_of <- split(M$snp$snp, M$snp$promoter)
if (arm == "real") { A <- td$act; P <- td$peer } else {
  m <- fread(file.path(V2, "perm_maps", "final", def, paste0(tis, ".", arm, ".perm_map.tsv")))
  stopifnot(setequal(m$slot, td$U), setequal(m$donor, td$U), !anyDuplicated(m$donor))
  d <- setNames(m$donor, m$slot)[td$U]; A <- td$act[d, , drop = FALSE]; P <- td$peer[d, , drop = FALSE]; rownames(A) <- rownames(P) <- td$U }
C_all <- cov_matrix(td, P, K)

## ---- per-promoter inputs (same steps as run_final2.R) ----
prep_promoter <- function(p) {
  z <- prepare_y(A[, p], ncol(C_all) - 1L); if (!is.na(z$skip)) stop("skip differs from the main run: ", z$skip)
  s <- names(z$y); y <- unname(z$y)
  pos <- td$info[p, "start"]; geno <- read_geno(chr, pos - FLANK, pos + FLANK - 1L, td$U)
  idx <- cis_index(geno, pos); G <- geno$G[s, idx, drop = FALSE]
  cnt <- apply(G, 2, function(x) tabulate(x + 1L, nbins = 3L)); keep <- which(apply(cnt, 2, min) >= MIN_GENO)
  if (!setequal(colnames(G)[keep], snp_of[[p]])) stop("tested SNP set differs from the main run")
  list(promoter = p, family = z$family, y = y, C = C_all[s, , drop = FALSE], G = G[, keep, drop = FALSE])
}
prs <- unique(U$promoter)
PP <- mclapply(prs, function(p) tryCatch(prep_promoter(p), error = function(e) list(promoter = p, err = conditionMessage(e))), mc.cores = NC, mc.preschedule = FALSE)
names(PP) <- prs
bad <- prs[!vapply(PP, function(x) is.list(x) && is.null(x$err), TRUE)]

mid_rank <- function(y, Cc) { m <- y > 0 & y < 1; c(n_mid = sum(m), rank = qr(Cc[m, , drop = FALSE])$rank) }
test_one <- function(type, y, Cc, g, nul, mr = NULL) {
  if (type == "mu") { r <- snp_tests(y, Cc, g, nul)
    n <- mr[["n_mid"]]; pp <- mr[["rank"]] + 1L; W <- pmax(r$mu_stat, 0)
    lF <- if (is.finite(W) && n > pp + 1L) pf((n - pp) * expm1(W / n), 1, n - pp, lower.tail = FALSE, log.p = TRUE) else NA_real_
    c(est = unname(r$mu_est), se = unname(r$mu_se), logp = lF, logp_chi = unname(r$mu_lrt_logp)) }
  else { r <- firth_snp_tests(y, Cc, g, nul)
    if (type == "p0") c(est = unname(r$nu_firth_est), se = unname(r$nu_firth_se), logp = unname(r$nu_firth_logp), logp_chi = NA_real_)
    else c(est = unname(r$tau_firth_est), se = unname(r$tau_firth_se), logp = unname(r$tau_firth_logp), logp_chi = NA_real_) } }
## one pass: for each job (unit i, conditioning set cs, SNPs to test), fit the null with cs appended, test the SNPs
run_pass <- function(jobs) {
  nl <- mclapply(jobs, function(jb) { Pp <- PP[[jb$promoter]]
    Cc <- orth_design(cbind(Pp$C, Pp$G[, jb$cs, drop = FALSE]))
    list(Cc = Cc, nul = promoter_null(Pp$y, Cc, Pp$family), dropped = attr(Cc, "dropped"), mr = mid_rank(Pp$y, Cc)) }, mc.cores = NC, mc.preschedule = FALSE)
  tasks <- list()
  for (b in seq_along(jobs)) { sn <- jobs[[b]]$test; if (!length(sn) || !is.list(nl[[b]])) next
    w <- c(mu = 1, p0 = 4, p1 = 4)[[jobs[[b]]$type]] * length(PP[[jobs[[b]]$promoter]]$y)
    for (ch in split(sn, ceiling(seq_along(sn) / 50))) tasks[[length(tasks) + 1]] <- list(b = b, snps = ch, w = w * length(ch)) }
  tasks <- tasks[order(-vapply(tasks, `[[`, 0, "w"))]
  one <- function(t) { jb <- jobs[[t$b]]; Pp <- PP[[jb$promoter]]; N <- nl[[t$b]]
    rbindlist(lapply(t$snps, function(sn) { v <- tryCatch(test_one(jb$type, Pp$y, N$Cc, Pp$G[, sn], N$nul, N$mr), error = function(e) c(est = NA, se = NA, logp = NA, logp_chi = NA))
      data.table(b = t$b, snp = sn, est = unname(v["est"]), se = unname(v["se"]), logp = unname(v["logp"]), logp_chi = unname(v["logp_chi"])) })) }
  res <- mclapply(tasks, function(t) tryCatch(one(t), error = function(e) NULL), mc.cores = NC, mc.preschedule = FALSE)
  dead <- which(!vapply(res, is.data.table, TRUE))
  if (length(dead)) res[dead] <- mclapply(tasks[dead], function(t) tryCatch(one(t), error = function(e) NULL), mc.cores = min(NC, length(dead)), mc.preschedule = FALSE)
  lost <- sum(!vapply(res, is.data.table, TRUE))
  list(res = rbindlist(res[vapply(res, is.data.table, TRUE)]), lost = lost,
       dropped = vapply(nl, function(x) if (is.list(x)) length(x$dropped) else NA_integer_, 0L), null_ok = vapply(nl, is.list, TRUE))
}

## ---- forward rounds (synchronous over units) ----
UU <- U[!promoter %in% bad]
init_sig <- function(i) { if (UU$mode[i] == "continue") s1[promoter == UU$promoter[i] & type == UU$type[i]][order(k)]$snp else UU$leader[i] }
st <- lapply(seq_len(nrow(UU)), function(i) list(promoter = UU$promoter[i], type = UU$type[i], signals = init_sig(i), active = TRUE, stop = NA_character_))
rounds <- list(); incomplete <- length(bad); sig_rows <- list()
for (i in seq_along(st)) {
  if (UU$mode[i] == "continue") { old <- s1[promoter == UU$promoter[i] & type == UU$type[i]][order(k)]
    sig_rows[[length(sig_rows) + 1]] <- old[, .(promoter, type, k, snp, selected_round, logp_at_selection, logp_chi_at_selection, acat_at_selection)][, from_v1 := TRUE] }
  else sig_rows[[length(sig_rows) + 1]] <- data.table(promoter = st[[i]]$promoter, type = st[[i]]$type, k = 1L, snp = st[[i]]$signals[1],
    selected_round = 0L, logp_at_selection = UU$leader_logp[i], logp_chi_at_selection = NA_real_, acat_at_selection = UU$p[i], from_v1 = FALSE) }
maxround <- if (PERM) 1L else MAXSIG
for (r in seq_len(maxround)) {
  act <- which(vapply(st, `[[`, TRUE, "active")); if (!length(act)) break
  jobs <- lapply(act, function(i) { s <- st[[i]]; list(i = i, promoter = s$promoter, type = s$type, cs = s$signals,
    test = setdiff(colnames(PP[[s$promoter]]$G), s$signals)) })
  ps <- run_pass(jobs); incomplete <- incomplete + ps$lost
  for (jb_i in seq_along(jobs)) { b <- jb_i; i <- jobs[[b]]$i; z <- ps$res[b == jb_i]
    acat <- acat_log(z$logp); best <- if (nrow(z) && any(is.finite(z$logp))) z[which.min(logp)] else NULL
    rounds[[length(rounds) + 1]] <- data.table(promoter = st[[i]]$promoter, type = st[[i]]$type, round = r, n_cond = length(st[[i]]$signals),
      n_tested = sum(is.finite(z$logp)), n_failed = sum(!is.finite(z$logp)), acat_p = acat, acat_p_chi = acat_log(z$logp_chi),
      best_snp = if (is.null(best)) NA_character_ else best$snp,
      best_logp = if (is.null(best)) NA_real_ else best$logp, cols_dropped = ps$dropped[b], null_ok = ps$null_ok[b])
    go <- is.finite(acat) && acat <= thr && !is.null(best) && length(st[[i]]$signals) < MAXSIG && !PERM
    if (go) { st[[i]]$signals <- c(st[[i]]$signals, best$snp)
      sig_rows[[length(sig_rows) + 1]] <- data.table(promoter = st[[i]]$promoter, type = st[[i]]$type, k = length(st[[i]]$signals), snp = best$snp,
        selected_round = r, logp_at_selection = best$logp, logp_chi_at_selection = best$logp_chi, acat_at_selection = acat) }
    else { st[[i]]$active <- FALSE
      st[[i]]$stop <- if (PERM) "perm: one round" else if (!is.finite(acat)) "no finite conditional p" else if (acat > thr) "window ACAT above threshold" else "cap" } }
}

## ---- leave-one-out window statistics for units with >= 2 signals ----
loo <- data.table()
multi <- which(vapply(st, function(s) length(s$signals) >= 2, TRUE))
if (length(multi)) {
  jobs <- unlist(lapply(multi, function(i) { s <- st[[i]]; lapply(seq_along(s$signals), function(k) list(i = i, k = k, promoter = s$promoter, type = s$type,
    cs = s$signals[-k], test = setdiff(colnames(PP[[s$promoter]]$G), s$signals[-k]))) }), recursive = FALSE)
  ps <- run_pass(jobs); incomplete <- incomplete + ps$lost
  loo <- ps$res[, `:=`(promoter = vapply(jobs, `[[`, "", "promoter")[b], type = vapply(jobs, `[[`, "", "type")[b], k = vapply(jobs, `[[`, 0L, "k")[b])][, b := NULL]
}
signals <- rbindlist(sig_rows, fill = TRUE)
if (nrow(loo)) { lk <- merge(signals[, .(promoter, type, k, snp)], loo, by = c("promoter", "type", "k", "snp"))
  signals <- merge(signals, lk[, .(promoter, type, k, loo_est = est, loo_se = se, loo_logp = logp, loo_logp_chi = logp_chi)], by = c("promoter", "type", "k"), all.x = TRUE) }
units <- rbindlist(lapply(st, function(s) data.table(promoter = s$promoter, type = s$type, n_signals = length(s$signals), stop = s$stop)))
units <- merge(units, UU[, .(promoter, type, family, n, n0, n1, unit_p = p, unit_q = q, leader, mode)], by = c("promoter", "type"))
if (length(bad)) units <- rbind(units, data.table(promoter = bad, type = NA_character_, n_signals = NA_integer_, stop = paste("prep error:", vapply(PP[bad], function(x) if (is.list(x)) x$err else "worker died", ""))), fill = TRUE)
meta <- list(def = def, tissue = tis, arm = arm, chr = chr, K = K, rule = RULE, thr = thr, perm_frac = if (PERM) FRAC else NA, n_units = nrow(UU), incomplete = incomplete,
             cores = NC, wall_sec = as.numeric(difftime(Sys.time(), t_start, units = "secs")), host = Sys.info()[["nodename"]], date = format(Sys.time()),
             engine = "bohr2_engine.R snp_tests (mu, F-corrected, p = interior rank + 1) / bohr2_addon3.R firth_snp_tests (p0, p1); leaders appended to the covariates (plan B)")
res <- list(units = units, signals = signals, rounds = rbindlist(rounds), loo = loo, meta = meta)
if (incomplete > 0) { saveRDS(res, sub("\\.rds$", ".incomplete.rds", out)); cat("INCOMPLETE:", incomplete, "\n"); quit(save = "no", status = 0) }
tmp <- paste0(out, ".tmp"); saveRDS(res, tmp); file.rename(tmp, out)
cat(sprintf("%s %s %s %s: %d units, signals %s, %.1f min\n", def, tis, arm, chr, nrow(UU), paste(names(table(units$n_signals)), table(units$n_signals), sep = ":", collapse = " "), meta$wall_sec / 60))
