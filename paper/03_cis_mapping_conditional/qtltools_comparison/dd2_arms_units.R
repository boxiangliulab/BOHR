### dd2_arms_units.R: rebuilds the units of Fig. 2E from BOHR's component tables, for the real arm and the
### three permuted arms, exactly as class_calib.R does with UNIT_MODE = merge_all (its tables are read, not written):
###   component tests = rerun_v4/results/final_nogate_v3m/freeze/<Tissue>.<arm>.units.rds, dedup_drop == FALSE, finite P (every computed test);
###   unit0 = mirror_pair_id for pairs exact in both parts, else the promoter; kept if in the unit list of fig2e_r5e80_units.R;
###   a partial pair (0/1 pattern exact, interior not) is merged into one unit; units with at least one estimable test.
### Calls, per tissue and arm:
###   BOHR_B  (rule B): a component passes BH 5% over all component tests of the tissue, or the unit's Simes P passes BH 5% over units
###           while the μ part and the 0/1 part both have P < 0.05. Real arm: class_calib's exact code (checked against K_units B_bh).
###   BOHR_S  : unit P = Simes over the unit's component tests; BH 5% over units.
###   BOHR_F  (exploratory, pooled evidence): part P = smallest P of the part x number of tests of the part (Bonferroni; parts: μ, 0/1);
###           unit P = Fisher's combination of the two part P values (4 d.f.) when both parts were tested, else the one part P; BH 5% over units.
###   QTLtools: unit P = smallest P of the unit's promoters (K_map), BH 5% over units; kinds nofill (donors without a value left out) and
###           pairs (mean fill; all permuted arms exist), modes raw and rint.
### For every permuted arm: discoveries with BH run inside the arm (the quantity shown in Fig. 2F) and the number of units
### that pass the REAL arm's thresholds (the estimate of false discoveries among the real discoveries).
### Output (this folder): dd2_arms.tsv (tissue x arm x method), dd2_units.tsv.gz (one row per unit, real arm), dd2_checks.txt.
### Usage: Rscript dd2_arms_units.R <ncores>
suppressPackageStartupMessages({ library(data.table); library(parallel) })
NC <- as.integer(commandArgs(TRUE)[1]); if (is.na(NC)) NC <- 4
G <- "/projects/group/bliu/e1327989/GTEx"; W <- file.path(G, "nofloor_full"); K <- file.path(W, "consult/fig2e_20261119"); MV <- file.path(K, "matched_variants")
KO <- file.path(W, "consult/why_fold_20261120/out"); V2 <- file.path(W, "GTEx/QTL_figure/rerun_v4"); OUT <- file.path(W, "figures_v4/r5e80/Figure2/fig2e_deepdive_20261120")
ARMS <- c("real", "perm0", "perm1", "perm2"); KD <- c(nofill = "pairs_nofill", pairs = "pairs"); MODES <- c("raw", "rint")
fr <- function(f, ...) fread(cmd = paste("zcat", shQuote(f)), ...)
X <- fr(file.path(K, "main/out/units_paper_r5e80.tsv.gz"), select = c("set", "tissue", "unit"))[set == "r5e80"]; tis <- sort(unique(X$tissue)); stopifnot(length(tis) == 50)
KM <- fr(file.path(KO, "K_map_r5e80_merge_all.tsv.gz"))[set == "r5e80"]
simes <- function(p) { p <- sort(p); k <- length(p); min(1, min(k * p / seq_len(k))) }
thr_bh <- function(p) { q <- p.adjust(p, "BH"); if (any(q <= 0.05)) max(p[q <= 0.05]) else -Inf }
cls <- function(ty) if (!length(ty)) NA_character_ else if (all(ty == "mu")) "mu only" else if (!any(ty == "mu")) "p0/p1 only" else "both"
EX <- c("leader", "mu_est_at_leader", "p0_est_at_leader", "p1_est_at_leader", "mu_logp_at_leader", "p0_logp_at_leader", "p1_logp_at_leader", "n", "n0", "n1", "n_mid", "family")
one <- function(t) { keep <- X[tissue == t, unit]; um <- KM[tissue == t, .(unit, promoter)]
  CA <- lapply(setNames(ARMS, ARMS), function(arm) { U <- as.data.table(readRDS(file.path(V2, "results/final_nogate_v3m/freeze", paste(t, arm, "units.rds", sep = "."))))[dedup_drop == FALSE & is.finite(p)]
    U[, unit0 := fifelse(boundary_exact %in% TRUE & interior_exact %in% TRUE, mirror_pair_id, promoter)]; U <- U[unit0 %in% keep]
    U[, partial := boundary_exact %in% TRUE & !(interior_exact %in% TRUE)]; U[, unit := fifelse(partial, mirror_pair_id, unit0)]
    U[, c("unit", "unit0", "promoter", "type", "p", "gene", "partial", "estimable", EX), with = FALSE] })
  CA <- lapply(CA, function(d) d[unit %in% d[estimable == TRUE, unit]])
  per_unit <- function(d) d[, { im <- type == "mu"; kmu <- sum(im); kb <- sum(!im); pmu <- if (kmu) min(p[im]) else NA_real_; pb <- if (kb) min(p[!im]) else NA_real_
      a <- if (kmu) min(1, kmu * pmu) else NA_real_; b <- if (kb) min(1, kb * pb) else NA_real_
      .(pS = simes(p), pm = min(p), k = .N, kmu = kmu, kb = kb, pmu = pmu, pb = pb,
        pF = if (kmu && kb) pchisq(-2 * (log(a) + log(b)), 4, lower.tail = FALSE) else if (kmu) a else b) }, by = unit]
  PU <- lapply(CA, per_unit); m <- nrow(PU$real)
  ## ---- real arm, class_calib's code for D, S and B ----
  Cr <- CA$real; tS <- thr_bh(PU$real$pS); tD <- thr_bh(Cr$p); tF <- thr_bh(PU$real$pF); R <- sum(PU$real$pS <= tS); lvlS <- 0.05 * R / m
  Ct <- Cr[, .(p = min(p)), by = .(unit, type)]; Ct[, padj := p.adjust(p, "BH"), by = unit]
  u <- copy(PU$real)
  cs <- Ct[unit %in% u[pS <= tS, unit] & padj <= lvlS, .(cl = cls(type)), by = unit]; u[, S_bh := cs$cl[match(unit, cs$unit)]]
  cd <- Cr[p <= tD, .(cl = cls(type)), by = unit]; u[, D_bh := cd$cl[match(unit, cd$unit)]]
  Cr[, qc := p.adjust(p, "BH")]
  nr <- Cr[, .(near = if (type[which.min(p)] == "mu") "mu only" else "p0/p1 only", qD = min(qc)), by = unit]; u <- merge(u, nr, by = "unit")
  u[!is.na(pmu) & !is.na(pb) & pmu < 0.05 & pb < 0.05, near := "both"]
  u[, B_bh := fifelse(!is.na(D_bh), D_bh, fifelse(!is.na(S_bh) & near == "both", "both", NA_character_))]
  u[, `:=`(S_unit = pS <= tS, F_unit = pF <= tF, qF = p.adjust(pF, "BH"), qS = p.adjust(pS, "BH"))]
  ## the unit's most significant μ row and most significant 0/1 row: effects of all components at that row's leader (same promoter, same variant)
  bm <- Cr[type == "mu"][order(p), .SD[1], by = unit][, .(unit, mu_leader = leader, mu_p = p, m_mu = mu_est_at_leader, m_p0 = p0_est_at_leader, m_p1 = p1_est_at_leader,
                                                        m_lp0 = p0_logp_at_leader, m_lp1 = p1_logp_at_leader)]
  bb <- Cr[type != "mu"][order(p), .SD[1], by = unit][, .(unit, b_type = type, b_leader = leader, b_p = p, b_mu = mu_est_at_leader, b_p0 = p0_est_at_leader, b_p1 = p1_est_at_leader, b_lmu = mu_logp_at_leader)]
  ft <- Cr[order(p), .SD[1], by = unit][, .(unit, gene, n, n0, n1, n_mid, family, partial)]
  u <- merge(merge(merge(u, bm, by = "unit", all.x = TRUE), bb, by = "unit", all.x = TRUE), ft, by = "unit")
  ## ---- BOHR calls per arm: BH inside the arm, and at the real arm's thresholds ----
  callB <- function(a, td, ts) a[, pm <= td | (pS <= ts & !is.na(pmu) & !is.na(pb) & pmu < 0.05 & pb < 0.05)]
  rows <- rbindlist(lapply(ARMS, function(arm) { a <- PU[[arm]]; d <- CA[[arm]]; td <- thr_bh(d$p); ts <- thr_bh(a$pS); tf <- thr_bh(a$pF)
    data.table(arm = arm, method = c("BOHR_B", "BOHR_S", "BOHR_F"), units = nrow(a), tests = c(nrow(d), nrow(a), nrow(a)),
               own_bh = c(sum(callB(a, td, ts)), sum(a$pS <= ts), sum(a$pF <= tf)),
               at_real_thr = c(sum(callB(a, tD, tS)), sum(a$pS <= tS), sum(a$pF <= tF)), own_thr = c(td, ts, tf)) }))
  rows[arm == "real" & method == "BOHR_B", own_bh := u[!is.na(B_bh), .N]]             # the real arm's B by class_calib's exact code (S part needs padj)
  ## ---- QTLtools ----
  for (kd in names(KD)) for (md in MODES) { qf <- function(arm) file.path(MV, KD[[kd]], "v4", md, paste(t, arm, "promoters.tsv.gz", sep = "."))
    rd1 <- function(arm) merge(fr(qf(arm), select = c("promoter", "p")), um, by = "promoter")[!is.na(p), .(p = min(p)), by = unit]
    qr <- rd1("real"); tq <- thr_bh(qr$p); nm <- paste("Q", kd, md, sep = "_")
    if (kd == "nofill") { u[, (paste0("p_", nm)) := qr$p[match(unit, qr$unit)]]; u[, (nm) := !is.na(get(paste0("p_", nm))) & get(paste0("p_", nm)) <= tq]; u[, (paste0("thr_", nm)) := tq] }
    rows <- rbind(rows, rbindlist(lapply(ARMS, function(arm) { if (!file.exists(qf(arm))) return(data.table(arm = arm, method = nm, units = NA_integer_, tests = NA_integer_, own_bh = NA_integer_, at_real_thr = NA_integer_, own_thr = NA_real_))
      q <- if (arm == "real") qr else rd1(arm); ta <- thr_bh(q$p); data.table(arm = arm, method = nm, units = nrow(q), tests = nrow(q), own_bh = sum(q$p <= ta), at_real_thr = sum(q$p <= tq), own_thr = ta) }))) }
  u[, `:=`(tissue = t, tD = tD, tS = tS, tF = tF)]; rows[, tissue := t]
  list(u = u, rows = rows) }
res <- mclapply(tis, one, mc.cores = NC, mc.preschedule = FALSE); bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r$u), TRUE)
if (any(bad)) { print(res[bad][[1]]); stop("failed tissues: ", paste(tis[bad], collapse = ", ")) }
UU <- rbindlist(lapply(res, `[[`, "u")); RR <- rbindlist(lapply(res, `[[`, "rows"))
## ---- checks against class_calib's stored table ----
KU <- fr(file.path(KO, "K_units_r5e80_merge_all.tsv.gz"), select = c("tissue", "unit", "B_bh", "S_bh", "D_bh", "near", "pS", "Q_nofill_raw_bh", "Q_nofill_rint_bh"))   # unit ids of single promoters repeat across tissues: key = tissue + unit
stopifnot(!anyDuplicated(UU[, .(tissue, unit)]), !anyDuplicated(KU[, .(tissue, unit)]))
ck <- merge(UU[, .(tissue, unit, B = fifelse(is.na(B_bh), "", B_bh), S = fifelse(is.na(S_bh), "", S_bh), near, pS, Q_nofill_raw, Q_nofill_rint)], KU, by = c("tissue", "unit"))
C <- c(sprintf("[%s] units: %d here, %d in K_units, %d matched", if (nrow(ck) == nrow(KU) && nrow(UU) == nrow(KU)) "PASS" else "FAIL", nrow(UU), nrow(KU), nrow(ck)),
       sprintf("[%s] BOHR class B_bh identical for every unit (mismatches %d)", if (ck[, all(B == fifelse(is.na(B_bh), "", B_bh))]) "PASS" else "FAIL", ck[, sum(B != fifelse(is.na(B_bh), "", B_bh))]),
       sprintf("[%s] S_bh identical (mismatches %d); near identical (mismatches %d); Simes P max rel. diff %.2g", if (ck[, all(S == fifelse(is.na(S_bh), "", S_bh)) && all(near.x == near.y)]) "PASS" else "FAIL",
               ck[, sum(S != fifelse(is.na(S_bh), "", S_bh))], ck[, sum(near.x != near.y)], ck[, max(abs(pS.x - pS.y) / pmax(pS.y, 1e-300))]),
       sprintf("[%s] QTLtools no-fill calls identical: raw mismatches %d, rint mismatches %d", if (ck[, all(Q_nofill_raw == (Q_nofill_raw_bh %in% TRUE)) && all(Q_nofill_rint == (Q_nofill_rint_bh %in% TRUE))]) "PASS" else "FAIL",
               ck[, sum(Q_nofill_raw != (Q_nofill_raw_bh %in% TRUE))], ck[, sum(Q_nofill_rint != (Q_nofill_rint_bh %in% TRUE))]),
       sprintf("totals, real arm: BOHR_B %d (mu only %d, p0/p1 only %d, both %d); QTLtools no-fill raw %d, rint %d", UU[!is.na(B_bh), .N], UU[B_bh == "mu only", .N], UU[B_bh == "p0/p1 only", .N], UU[B_bh == "both", .N],
               UU[, sum(Q_nofill_raw)], UU[, sum(Q_nofill_rint)]))
tot <- RR[, .(tissues = sum(!is.na(own_bh)), units = sum(units), own_bh = sum(own_bh), at_real_thr = sum(at_real_thr)), by = .(method, arm)]
writeLines(C, file.path(OUT, "dd2_checks.txt")); cat(C, sep = "\n"); print(tot)
fwrite(RR, file.path(OUT, "dd2_arms.tsv"), sep = "\t"); fwrite(tot, file.path(OUT, "dd2_arms_totals.tsv"), sep = "\t")
fwrite(UU, file.path(OUT, "dd2_units.tsv"), sep = "\t"); system(paste("gzip -f", shQuote(file.path(OUT, "dd2_units.tsv")))); cat("DD2_DONE\n")
