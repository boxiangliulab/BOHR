### dd11_fisher_rules.R: BOHR's unit calls under several ways of combining a unit's evidence; Fig. 2E uses rule BF with PART=sidak.
### BOHR's unit calls under different ways of combining the μ part and the 0/1 part, real arm and the three permuted arms (units and tests as
### dd2_arms_units.R = class_calib.R UNIT_MODE merge_all). Part P = smallest P of the part x number of tests of the part (Bonferroni inside the part).
###   B   : a component passes BH 5% over all component tests of the tissue, or (no component passes, both parts P < 0.05) the Simes P passes BH over units.
###   BF  : the same, with Fisher's combination of the two part P values (4 d.f.) in place of Simes in the second step.
###   F   : Fisher only: one P per unit (Fisher when both parts were tested, else the one part P), BH 5% over units.
###   S   : Simes only: one P per unit, BH 5% over units.
###   JF  : one BH family per tissue = every component test + the Fisher P of every unit with both parts; a unit is called when any of its members passes.
### Class of a call: by the components that pass; a unit called only through the combined P is "both" (BF, JF); F and S: by the same-rule class (near).
### For each rule and permuted arm: calls with BH run inside the arm, and units passing the real arm's thresholds. QTLtools calls from dd2_units.tsv.gz.
### Output: dd11_fisher_rules.txt, dd11_units.tsv.gz, dd11_arms.tsv (per tissue, arm and rule). Usage: Rscript dd11_fisher_rules.R <ncores>
### PART: how a part's tests become the part P used by Fisher. "bonferroni" (the default, used in an earlier run
### of this script): min(1, k x smallest P). "sidak": 1 - (1 - smallest P)^k, exact for independent tests; outputs carry the suffix _sidak. Simes, the rescue
### condition (both parts' smallest P < 0.05) and the placement (near) do not use the part P and are the same in both.
suppressPackageStartupMessages({ library(data.table); library(parallel) })
NC <- as.integer(commandArgs(TRUE)[1]); if (is.na(NC)) NC <- 2
PART <- Sys.getenv("PART", "bonferroni"); stopifnot(PART %in% c("bonferroni", "sidak")); SFX <- if (PART == "sidak") "_sidak" else ""
partp <- function(p, k) if (PART == "sidak") -expm1(k * log1p(-p)) else min(1, k * p)
G <- "/projects/group/bliu/e1327989/GTEx"; W <- file.path(G, "nofloor_full"); K <- file.path(W, "consult/fig2e_20261119"); V2 <- file.path(W, "GTEx/QTL_figure/rerun_v4")
setwd(file.path(W, "figures_v4/r5e80/Figure2/fig2e_deepdive_20261120")); fr <- function(f, ...) fread(cmd = paste("zcat", shQuote(f)), ...); ARMS <- c("real", "perm0", "perm1", "perm2")
X <- fr(file.path(K, "main/out/units_paper_r5e80.tsv.gz"), select = c("set", "tissue", "unit"))[set == "r5e80"]; tis <- sort(unique(X$tissue))
simes <- function(p) { p <- sort(p); k <- length(p); min(1, min(k * p / seq_len(k))) }
thr_bh <- function(p) { q <- p.adjust(p, "BH"); if (any(q <= 0.05)) max(p[q <= 0.05]) else -Inf }
cls <- function(ty) if (!length(ty)) NA_character_ else if (all(ty == "mu")) "mu only" else if (!any(ty == "mu")) "p0/p1 only" else "both"
one <- function(t) { keep <- X[tissue == t, unit]
  CA <- lapply(setNames(ARMS, ARMS), function(arm) { U <- as.data.table(readRDS(file.path(V2, "results/final_nogate_v3m/freeze", paste(t, arm, "units.rds", sep = "."))))[dedup_drop == FALSE & is.finite(p)]
    U[, unit0 := fifelse(boundary_exact %in% TRUE & interior_exact %in% TRUE, mirror_pair_id, promoter)]; U <- U[unit0 %in% keep]
    U[, partial := boundary_exact %in% TRUE & !(interior_exact %in% TRUE)]; U[, unit := fifelse(partial, mirror_pair_id, unit0)]; U <- U[unit %in% U[estimable == TRUE, unit]]; U[, .(unit, type, p)] })
  per_unit <- function(d) d[, { im <- type == "mu"; kmu <- sum(im); kb <- sum(!im); pmu <- if (kmu) min(p[im]) else NA_real_; pb <- if (kb) min(p[!im]) else NA_real_
      a <- if (kmu) partp(pmu, kmu) else NA_real_; b <- if (kb) partp(pb, kb) else NA_real_
      .(pS = simes(p), pm = min(p), two = kmu > 0 && kb > 0, pmu = pmu, pb = pb, pF = if (kmu && kb) pchisq(-2 * (log(a) + log(b)), 4, lower.tail = FALSE) else if (kmu) a else b) }, by = unit]
  PU <- lapply(CA, per_unit)
  thr <- function(arm) { a <- PU[[arm]]; d <- CA[[arm]]; list(tD = thr_bh(d$p), tS = thr_bh(a$pS), tF = thr_bh(a$pF), tJ = thr_bh(c(d$p, a[two == TRUE, pF]))) }
  call <- function(a, h) { nom <- !is.na(a$pmu) & !is.na(a$pb) & a$pmu < 0.05 & a$pb < 0.05
    data.table(unit = a$unit, D = a$pm <= h$tD, B = a$pm <= h$tD | (nom & a$pS <= h$tS), BF = a$pm <= h$tD | (nom & a$pF <= h$tF), F = a$pF <= h$tF, S = a$pS <= h$tS, JF = a$pm <= h$tJ | (a$two & a$pF <= h$tJ)) }
  hr <- thr("real"); cr <- call(PU$real, hr); RUL <- c("B", "BF", "F", "S", "JF")
  rows <- rbindlist(lapply(ARMS, function(arm) { own <- call(PU[[arm]], thr(arm)); atr <- call(PU[[arm]], hr)
    data.table(arm = arm, rule = RUL, units = nrow(own), own_bh = vapply(RUL, function(r) sum(own[[r]]), 0), at_real_thr = vapply(RUL, function(r) sum(atr[[r]]), 0)) }))
  ## classes in the real arm
  u <- merge(PU$real, cr, by = "unit"); Cr <- CA$real
  cD <- Cr[p <= hr$tD, .(cl = cls(type)), by = unit]; cJ <- Cr[p <= hr$tJ, .(cl = cls(type)), by = unit]
  u[, near := fifelse(!is.na(pmu) & !is.na(pb) & pmu < 0.05 & pb < 0.05, "both", fifelse(is.na(pb) | (!is.na(pmu) & pmu <= pb), "mu only", "p0/p1 only"))]
  u[, `:=`(cl_B = fifelse(D, cD$cl[match(unit, cD$unit)], fifelse(B, "both", NA_character_)), cl_BF = fifelse(D, cD$cl[match(unit, cD$unit)], fifelse(BF, "both", NA_character_)),
           cl_JF = fifelse(JF, fifelse(unit %in% cJ$unit, cJ$cl[match(unit, cJ$unit)], "both"), NA_character_), tissue = t, tD = hr$tD, tS = hr$tS, tF = hr$tF, tJ = hr$tJ)]
  rows[, tissue := t]; list(u = u, rows = rows) }
res <- mclapply(tis, one, mc.cores = NC, mc.preschedule = FALSE); bad <- vapply(res, function(r) !is.list(r) || is.null(r$u), TRUE); if (any(bad)) { print(res[bad][[1]]); stop("failed") }
UU <- rbindlist(lapply(res, `[[`, "u")); RR <- rbindlist(lapply(res, `[[`, "rows"))
D2 <- fread(cmd = "zcat dd2_units.tsv.gz", select = c("tissue", "unit", "B_bh", "near", "Q_nofill_raw", "Q_nofill_rint")); setnames(D2, "near", "near2")
UU <- merge(UU, D2, by = c("tissue", "unit")); stopifnot(nrow(UU) == 226093, UU[, all(near == near2)], UU[, all(fifelse(is.na(cl_B), "", cl_B) == fifelse(is.na(B_bh), "", B_bh))])   # rule B here = rule B of dd2_arms_units.R, every unit
CL <- c("mu only", "p0/p1 only", "both"); out <- list()
tot <- RR[, .(own_bh = sum(own_bh), at_real_thr = sum(at_real_thr)), by = .(rule, arm)]
out[["1 calls per rule: real arm; permuted arms with BH inside the arm; permuted units passing the real thresholds"]] <- merge(merge(
  tot[arm == "real", .(rule, real = own_bh)], dcast(tot[arm != "real"], rule ~ arm, value.var = "own_bh"), by = "rule"),
  tot[arm != "real", .(permuted_at_real_thresholds_mean = round(mean(at_real_thr), 1)), by = rule], by = "rule")[, share_of_calls := round(permuted_at_real_thresholds_mean / real, 4)][match(c("B", "BF", "F", "S", "JF"), rule)]
for (m in c("raw", "rint")) { Q <- paste0("Q_nofill_", m)
  out[[paste0("2 totals against QTLtools ", m)]] <- rbindlist(lapply(c("B", "BF", "F", "S", "JF"), function(r) UU[, .(rule = r, BOHR = sum(get(r)), QTLtools = sum(get(Q)), shared = sum(get(r) & get(Q)), BOHR_only = sum(get(r) & !get(Q)),
     QTLtools_only = sum(!get(r) & get(Q)))]))[, fold := round(BOHR / QTLtools, 3)]
  ## by class, Fig. 2E placement: a BOHR call by its class under the rule, a unit only QTLtools finds by near
  for (r in c("B", "BF", "JF")) { cc <- paste0("cl_", r); UU[, pl := fifelse(get(r), get(cc), near)]
    out[[paste0("3 by class (draft placement), rule ", r, ", QTLtools ", m)]] <- UU[, .(BOHR = sum(get(r)), QTLtools = sum(get(Q)), shared = sum(get(r) & get(Q)), BOHR_only = sum(get(r) & !get(Q)), QTLtools_only = sum(!get(r) & get(Q))),
       keyby = .(class = factor(pl, levels = CL))][, fold := round(BOHR / QTLtools, 3)] }
  ## by class, one rule for every unit (near)
  out[[paste0("4 by class (one rule for every unit), QTLtools ", m)]] <- dcast(rbindlist(lapply(c("B", "BF", "F", "JF"), function(r) UU[, .(rule = r, BOHR = sum(get(r)), QTLtools = sum(get(Q))), keyby = .(class = factor(near, levels = CL))])),
     class + QTLtools ~ rule, value.var = "BOHR")[, `:=`(fold_B = round(B / QTLtools, 3), fold_BF = round(BF / QTLtools, 3), fold_F = round(F / QTLtools, 3), fold_JF = round(JF / QTLtools, 3))] }
out[["5 what BF adds to the draft's rule"]] <- UU[, .(added = sum(BF & !B), of_these_QTLtools_raw_finds = sum(BF & !B & Q_nofill_raw), of_these_QTLtools_rint_finds = sum(BF & !B & Q_nofill_rint), dropped = sum(B & !BF))]
sink(paste0("dd11_fisher_rules", SFX, ".txt")); cat("PART =", PART, "\n"); for (k in names(out)) { cat("\n####", k, "\n"); print(out[[k]]) }; sink()
fwrite(RR, paste0("dd11_arms", SFX, ".tsv"), sep = "\t")   # per tissue x arm x rule (read by prepare_panel_e_classes.R)
fwrite(UU, paste0("dd11_units", SFX, ".tsv"), sep = "\t"); system(paste0("gzip -f dd11_units", SFX, ".tsv")); cat(readLines(paste0("dd11_fisher_rules", SFX, ".txt")), sep = "\n"); cat("DD11_DONE\n")
