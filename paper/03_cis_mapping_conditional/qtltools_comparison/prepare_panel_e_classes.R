#!/usr/bin/env Rscript
### prepare_panel_e_classes.R: the Fig. 2E and Fig. 2F results (Source Data) for the final rules
### (r5e80), 50 GTEx tissues, units = class_calib.R UNIT_MODE merge_all.
###   E  BOHR's discoveries by the part of its model that is significant, against QTLtools raw and rank-normal first (no fill), BH 5% within tissue:
###      only μ / only 0/1 = the component's own test (BH 5% over all component tests of the tissue); μ and 0/1 also takes the units with no passing
###      component, both parts P < 0.05, whose combined P passes BH 5% over units: rule B = Simes, rule BF = Fisher (part P = smallest P x tests of the part).
###      A unit only QTLtools finds is placed by BOHR's P values (both parts P < 0.05 -> μ and 0/1, else the part with the smaller P).
###   F  promoters called in the permuted data (phenotypes and PEER factors permuted together), BH 5% inside each permuted arm, summed over 50 tissues,
###      with the range expected by chance: total calls of a calibrated method whose tests are independent and uniform under the null, simulated
###      with each tissue's number of units (one P per unit), 20,000 replicates; the same checked by simulating BOHR's own rules B and BF with each
###      unit's numbers of μ and 0/1 tests (2,000 replicates).
### Inputs (read only): nofloor_full/figures_v4/r5e80/Figure2/fig2e_deepdive_20261120/{dd11_units.tsv.gz, dd11_arms.tsv, dd2_arms.tsv, dd2_units.tsv.gz}.
### Output: data_v3_E80_nofloor_reads5/panel_e_classes_20261120/{classes.csv, permuted_by_arm.csv, null_band.csv, CHECKS.txt}.
### PART: the part P of rule BF. "bonferroni" (the default, used by earlier versions of the figure):
### min(1, k x smallest P); "sidak": 1 - (1 - smallest P)^k (dd11_fisher_rules.R PART=sidak). Inputs dd11_*_sidak, output panel_e_classes_20261120_sidak;
### the chance band of rule BF is simulated with the same part P. Rule B does not use the part P and is the same in both.
suppressMessages(library(data.table))
DD <- "/projects/group/bliu/e1327989/GTEx/nofloor_full/figures_v4/r5e80/Figure2/fig2e_deepdive_20261120"
PART <- Sys.getenv("PART", "bonferroni"); stopifnot(PART %in% c("bonferroni", "sidak")); SFX <- if (PART == "sidak") "_sidak" else ""
partv <- function(p, k) if (PART == "sidak") -expm1(k * log1p(-pmin(p, 1))) else pmin(1, k * p)   # vectorised part P (k = 0 rows are masked by the caller)
OD <- paste0("/projects/group/bliu/e1327989/GTEx/new_simulation/figures_bohr2/data_v3_E80_nofloor_reads5/panel_e_classes_20261120", SFX); dir.create(OD, showWarnings = FALSE)
U <- fread(cmd = paste("zcat", file.path(DD, paste0("dd11_units", SFX, ".tsv.gz")))); stopifnot(nrow(U) == 226093)
C <- character(0); ok <- function(nm, cond) { C <<- c(C, sprintf("[%s] %s", if (isTRUE(cond)) "PASS" else "FAIL", nm)); invisible(cond) }
RULES <- c(B = "cl_B", BF = "cl_BF")
## ---- E ----
CL <- rbindlist(lapply(names(RULES), function(r) { U[, call := get(r)]; U[, place := fifelse(call, get(RULES[[r]]), near)]
  rbindlist(lapply(c("raw", "rint"), function(m) { Q <- paste0("Q_nofill_", m)
    U[, .(BOHR = sum(call), shared = sum(call & get(Q)), BOHR_only = sum(call & !get(Q)), QTLtools_only = sum(!call & get(Q))), keyby = .(class = place)][, `:=`(rule = r, mode = m)] })) }))
CL[, `:=`(QTLtools = shared + QTLtools_only)][, fold := BOHR / QTLtools]; setcolorder(CL, c("rule", "class", "mode"))
K0 <- fread("/projects/group/bliu/e1327989/GTEx/nofloor_full/consult/why_fold_20261120/out/K_summary_r5e80_merge_all.tsv")[scheme == "B" & thr == "bh" & kind == "nofill" & mode %in% c("raw", "rint") & class != "all"]
ck <- merge(CL[rule == "B"], K0[, .(class, mode, b = bohr, a = also, o = only)], by = c("class", "mode"))
ok("rule B = class_calib.R scheme B (K_summary_r5e80_merge_all.tsv), 3 classes x 2 modes: BOHR, shared, QTLtools only", nrow(ck) == 6 && ck[, all(BOHR == b & shared == a & QTLtools_only == o)])
tots <- CL[, .(s = sum(BOHR)), by = .(rule, mode)]; tots[, want := vapply(rule, function(r) U[, sum(get(r))], 0L)]
ok(sprintf("classes add up to each rule's total (B %d; BF %d)", U[, sum(B)], U[, sum(BF)]), tots[, all(s == want)])
ok("only μ and only 0/1 identical under B and BF (only μ and 0/1 changes)", CL[class != "both", .(d = uniqueN(paste(BOHR, shared, QTLtools_only))), by = .(class, mode)][, all(d == 1)])
if (PART == "sidak") {   # the Šidák part P is never larger than the Bonferroni one, and BH is monotone: every Bonferroni BF call stays a call
  U0 <- fread(cmd = paste("zcat", file.path(DD, "dd11_units.tsv.gz")), select = c("tissue", "unit", "B", "BF")); m0 <- merge(U0, U[, .(tissue, unit, B, BF)], by = c("tissue", "unit"), suffixes = c("_bonf", "_sidak"))
  ok(sprintf("Šidák: same units and the same rule B calls as the Bonferroni run (%d units)", nrow(m0)), nrow(m0) == nrow(U) && m0[, all(B_bonf == B_sidak)])
  ok(sprintf("Šidák: every Bonferroni BF call is still called (BF Bonferroni %d, Šidák %d, added %d, lost %d)", m0[, sum(BF_bonf)], m0[, sum(BF_sidak)],
             m0[, sum(BF_sidak & !BF_bonf)], m0[, sum(BF_bonf & !BF_sidak)]), m0[, !any(BF_bonf & !BF_sidak)]) }
fwrite(CL, file.path(OD, "classes.csv"))
## ---- F: observed ----
AB <- fread(file.path(DD, paste0("dd11_arms", SFX, ".tsv")))[arm != "real" & rule %in% names(RULES)]; AQ <- fread(file.path(DD, "dd2_arms.tsv"))[method %in% c("Q_nofill_raw", "Q_nofill_rint") & arm != "real"]
ok("BOHR permuted arms: 3 arms x 50 tissues for rules B and BF", AB[, .(t = uniqueN(tissue)), by = .(rule, arm)][, .N == 6 && all(t == 50)])
qa <- AQ[!is.na(own_bh), .N, by = .(method, arm)][N == 50]; armsQ <- intersect(qa[method == "Q_nofill_raw", arm], qa[method == "Q_nofill_rint", arm])
P <- rbind(AB[, .(n = sum(own_bh), at_real = sum(at_real_thr)), by = .(rule, arm)][, key := "BOHR"],
           rbindlist(lapply(names(RULES), function(r) AQ[arm %in% armsQ, .(n = sum(own_bh), at_real = sum(at_real_thr)), by = .(method, arm)][, `:=`(rule = r, key = sub("Q_nofill_", "", method))][, method := NULL])))
setcolorder(P, c("rule", "key", "arm")); fwrite(P, file.path(OD, "permuted_by_arm.csv"))
C <- c(C, sprintf("  QTLtools no-fill permuted arms complete for all 50 tissues: %s", paste(armsQ, collapse = ", ")))
## ---- F: expected by chance ----
set.seed(20261120); q <- 0.05
m_t <- U[, .N, by = tissue]$N; nT <- length(m_t); REP <- 20000L; Kk <- 40L
gen <- matrix(0L, REP, nT)                                     # one P per unit, independent uniforms: the smallest Kk order statistics are enough (BH under the null rejects few)
UT <- upper.tri(diag(Kk), diag = TRUE) * 1; IDX <- matrix(seq_len(Kk), REP, Kk, byrow = TRUE)
for (j in seq_len(nT)) { E <- matrix(rexp(REP * Kk), REP, Kk); S <- E %*% UT; tot <- S[, Kk] + rgamma(REP, m_t[j] + 1 - Kk); Uo <- S / tot   # U_(k) = S_k / S_(m+1)
  crit <- matrix(rep(seq_len(Kk) * q / m_t[j], each = REP), REP, Kk); gen[, j] <- as.integer(do.call(pmax, as.data.frame((Uo <= crit) * IDX))) }   # BH: largest k with U_(k) <= k q / m
ok("order-statistic shortcut never saturated (no replicate with Kk calls in a tissue)", max(gen) < Kk)
tot1 <- rowSums(gen)
band <- function(x, lab, reps) { m3 <- colMeans(matrix(sample(x, 3 * 20000, replace = TRUE), 3))
  data.table(model = lab, replicates = reps, expected = mean(x), q025 = unname(quantile(x, 0.025, type = 1)), q975 = unname(quantile(x, 0.975, type = 1)),
             p_ge_7 = mean(x >= 7), mean3_q025 = unname(quantile(m3, 0.025, type = 1)), mean3_q975 = unname(quantile(m3, 0.975, type = 1))) }
NB <- band(tot1, "calibrated, one P per unit (BH over units)", REP)
## BOHR's own rules under the same null: each unit's μ and 0/1 tests independent uniforms; component BH; Simes / Fisher step for units with both parts P < 0.05
D2 <- fread(cmd = paste("zcat", file.path(DD, "dd2_units.tsv.gz")), select = c("tissue", "unit", "kmu", "kb")); stopifnot(nrow(D2) == 226093)
simes_row <- function(M) apply(M, 1, function(v) { v <- sort(v[is.finite(v)]); k <- length(v); min(1, min(k * v / seq_len(k))) })
bh_thr <- function(p) { m <- length(p); s <- sort(p); w <- which(s <= seq_len(m) * q / m); if (length(w)) s[max(w)] else -Inf }
REP2 <- 2000L; simB <- matrix(0L, REP2, nT); simBF <- matrix(0L, REP2, nT); tis <- unique(U$tissue)
for (j in seq_along(tis)) { d <- D2[tissue == tis[j]]; nu <- nrow(d); kmu <- d$kmu; kb <- d$kb; Km <- max(kmu); Kb <- max(kb)
  for (i in seq_len(REP2)) {
    Mm <- matrix(Inf, nu, max(1, Km)); Mb <- matrix(Inf, nu, max(1, Kb))
    for (c in seq_len(Km)) { w <- kmu >= c; Mm[w, c] <- runif(sum(w)) }; for (c in seq_len(Kb)) { w <- kb >= c; Mb[w, c] <- runif(sum(w)) }
    allp <- c(Mm[is.finite(Mm)], Mb[is.finite(Mb)]); tD <- bh_thr(allp)
    pmu <- if (Km) do.call(pmin, as.data.frame(Mm)) else rep(Inf, nu); pb <- if (Kb) do.call(pmin, as.data.frame(Mb)) else rep(Inf, nu); pm <- pmin(pmu, pb)
    a <- ifelse(kmu > 0, partv(pmu, kmu), NA_real_); b <- ifelse(kb > 0, partv(pb, kb), NA_real_)
    pF <- ifelse(!is.na(a) & !is.na(b), pchisq(-2 * (log(a) + log(b)), 4, lower.tail = FALSE), ifelse(!is.na(a), a, b))
    pS <- pm * (kmu + kb); low <- which(pm <= 1e-3); if (length(low)) pS[low] <- simes_row(cbind(Mm, Mb)[low, , drop = FALSE])   # Simes exact where it can matter (P <= 1e-3)
    nom <- kmu > 0 & kb > 0 & pmu < 0.05 & pb < 0.05
    simB[i, j] <- sum(pm <= tD | (nom & pS <= bh_thr(pS))); simBF[i, j] <- sum(pm <= tD | (nom & pF <= bh_thr(pF))) } }
NB <- rbind(NB, band(rowSums(simB), "BOHR rule B (components + Simes) under an independent null", REP2), band(rowSums(simBF), if (PART == "sidak") "BOHR rule BF (components + Fisher, Šidák parts) under an independent null" else "BOHR rule BF (components + Fisher) under an independent null", REP2))
NB[, applies_to := c("QTLtools raw and rank-normal first (one P per unit)", "BOHR, rule B", "BOHR, rule BF")]
C <- c(C, sprintf("  note: chance level differs by procedure: one P per unit %.2f (95%% %d-%d), BOHR rule B %.2f (%d-%d), BOHR rule BF %.2f (%d-%d); Fisher's step adds a second BH family",
                  NB$expected[1], NB$q025[1], NB$q975[1], NB$expected[2], NB$q025[2], NB$q975[2], NB$expected[3], NB$q025[3], NB$q975[3]))
ok("rule B's chance level equals the one-P-per-unit level (expected within 0.2)", abs(NB$expected[2] - NB$expected[1]) < 0.2)
fwrite(NB, file.path(OD, "null_band.csv")); obs <- P[, .(mean = mean(n), max = max(n)), by = .(rule, key)]
C <- c(C, "null range (single arm): ", capture.output(print(NB)), "observed (BH inside each permuted arm, summed over 50 tissues):", capture.output(print(P)))
writeLines(c(sprintf("prepare_panel_e_classes.R, %s", format(Sys.time(), "%Y-%m-%d %H:%M")), C), file.path(OD, "CHECKS.txt")); cat(readLines(file.path(OD, "CHECKS.txt")), sep = "\n")
stopifnot(!any(grepl("^\\[FAIL", C))); cat("PANEL_E_CLASSES_DONE\n")
