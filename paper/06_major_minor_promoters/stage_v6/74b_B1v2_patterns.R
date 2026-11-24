## 74b_B1v2_patterns.R -- patterns from the B1 refit (73b), with the A1 / A2 aligned contrasts recomputed on
## the same signal sets. Changes from 74 / 82:
##  - aligned frame uses the designated partner promoter P (not the sum of all other promoters);
##  - direction: conditional (loo_est) for signals of multi-signal units, marginal otherwise (dir_c); the marginal one is
##    kept as a sensitivity;
##  - S's usage direction for p1 (S = partner) is the mirror of own usage (-dir_c); agreement with the measured partner
##    usage fit (z_uP) is reported, and A2 is recomputed with the measured direction;
##  - "flat" by equivalence: the 90% CI of the standardised effect (b / SD of the response, per ALT allele) inside
##    +-delta (0.1; 0.05 as sensitivity); otherwise "undetermined" unless |z| >= 1.96;
##  - N level (donors with N > 0) for every signal with >= 30 such donors (no zero-fraction gate);
##  - sensitivity sets are real intersections: exact mirrors; single-signal units; no splice-region candidate; one signal
##    per gene and group (strongest); and all of the last three together.
source(file.path(Sys.getenv("F5S", "/home/e1327989/f5_stage_20261106_v6"), "code/00_helpers.R"))
OUT <- file.path(OUTROOT, "B1v2_multipromoter"); setDTthreads(2); B <- 2000L; ZC <- 1.96; DEL <- c(0.1, 0.05)
F <- fr(file.path(OUT, "B1v2_fits.tsv.gz")); F[, sid := paste(tissue, promoter, type, k, sep = "|")]
cat("fits:", nrow(F), "\n"); print(F[, .N, by = status]); F <- F[status == "ok"]
sg <- fread(file.path(OUTROOT, "A2_eqtl_sqtl/A2_signal_eqtl.tsv"), select = c("sid", "gene", "n_prom_gene", "role", "route", "dir", "z_eq"))
s0 <- fr(file.path(IN, "delivery/signals.tsv.gz"), select = c("tissue", "promoter", "type", "k", "loo_est", "backward_keep"))
s0[, sid := paste(tissue, promoter, type, k, sep = "|")]
s0[, unit := paste(tissue, promoter, type, sep = "|")][, n_unit := .N, by = unit]
s0[, dir_loo := fifelse(is.na(loo_est), NA_real_, fifelse(type == "mu", sign(loo_est), fifelse(type == "p0", -sign(loo_est), sign(loo_est))))]
X <- merge(F, sg, by = "sid"); X <- merge(X, s0[, .(sid, n_unit, dir_loo, bk = as.logical(backward_keep))], by = "sid")
X <- X[bk == TRUE]   # v5: only the signals retained after the backward step
fl <- fr(file.path(IN, "mirror/promoter_flags.tsv.gz"))[, .(tissue, promoter, exact_mirror)]
X <- merge(X, fl, by = c("tissue", "promoter"))
X[, dir_c := fifelse(n_unit > 1 & !is.na(dir_loo), dir_loo, dir)]
cat("signals:", nrow(X), "; multi-signal-unit signals:", X[n_unit > 1, .N], "; marginal vs conditional direction differ:", X[dir_c != dir, .N], "\n")
X[, two := !is.na(route) & route != "" & !is.na(partner) & partner != ""]
X[two == TRUE, S_is_own := fifelse(type == "p0", TRUE, fifelse(type == "p1", FALSE, role == "minor"))]
X[two == TRUE, sS := fifelse(S_is_own, dir_c, -dir_c)]
X[two == TRUE, sS_meas := fifelse(S_is_own, dir_c, fifelse(!is.na(z_uP) & abs(z_uP) >= ZC, sign(z_uP), -dir_c))]
pick <- function(own, a, b) fifelse(own, a, b)
X[two == TRUE, `:=`(
  zS = pick(S_is_own, z_t, z_P) * sS, zN = pick(S_is_own, z_P, z_t) * sS, zTS = z_T * sS,
  dN = pick(S_is_own, b_P / sd_P, b_t / sd_t) * sS, seN = pick(S_is_own, se_P / sd_P, se_t / sd_t),
  zlevN = pick(S_is_own, z_levP, z_levt) * sS, dlevN = pick(S_is_own, b_levP / sd_levP, b_levt / sd_levt) * sS,
  selevN = pick(S_is_own, se_levP / sd_levP, se_levt / sd_levt),
  dT = b_T / sd_T * sS, seT = se_T / sd_T,
  pS = pick(S_is_own, zperm_t, zperm_P) * sS, pN = pick(S_is_own, zperm_P, zperm_t) * sS, pT = zperm_T * sS,
  plevN = pick(S_is_own, zperm_levP, zperm_levt) * sS)]
ncat <- function(z, d, se, del) fifelse(is.na(z), NA_character_, fifelse(z <= -ZC, "opposite", fifelse(z >= ZC, "same",
  fifelse(!is.na(d) & abs(d) + 1.645 * se < del, "equivalent_zero", "undetermined"))))
## v5 (Methods §9): code by |z| regardless of direction (z oriented to the allele raising the minor promoter's usage)
patS <- function(s, n) fifelse(abs(s) >= ZC & abs(n) < ZC, "S_only", fifelse(abs(s) < ZC & abs(n) >= ZC, "N_only",
  fifelse(abs(s) >= ZC & abs(n) >= ZC, fifelse(sign(s) == sign(n), "coactivation", "reciprocal_competition"), "neither")))
X[two == TRUE, `:=`(patternS = patS(zS, zN), patternS_perm = patS(pS, pN),
                    N_cat = ncat(zN, dN, seN, DEL[1]), N_cat05 = ncat(zN, dN, seN, DEL[2]), Nlev_cat = ncat(zlevN, dlevN, selevN, DEL[1]),
                    T_cat = ncat(zTS, dT, seT, DEL[1]), N_cat_perm = ncat(pN, NA_real_, NA_real_, DEL[1]))]
## own frame (signed by own usage, conditional direction), all genes
X[, `:=`(zt = z_t * dir_c, zO = z_O * dir_c, zT = z_T * dir_c)]
pat <- function(a, b) fifelse(a >= ZC & abs(b) < ZC, "target_only", fifelse(abs(a) < ZC & b <= -ZC, "other_only", fifelse(a >= ZC & b <= -ZC, "reciprocal",
  fifelse(a >= ZC & b >= ZC, "both_up", fifelse(a <= -ZC & b <= -ZC, "both_down", fifelse(abs(a) < ZC & abs(b) < ZC, "neither", "other"))))))
X[, pattern := pat(zt, zO)]
## A2 recomputed: eQTL concordance with own usage (conditional direction) and with S usage (mirrored / measured)
X[, conc_own_c := fifelse(abs(z_eq) >= ZC, sign(z_eq) == dir_c, NA)]
X[two == TRUE, `:=`(conc_S_c = fifelse(abs(z_eq) >= ZC, sign(z_eq) == sS, NA), conc_S_meas = fifelse(abs(z_eq) >= ZC, sign(z_eq) == sS_meas, NA))]
## A1 shares (unchanged inputs) for the same sets
a1 <- fread(file.path(OUTROOT, "A1_location/A1_signal_location.tsv"), select = c("sid", "f_S", "f_N", "nlp"))
X <- merge(X, a1, by = "sid", all.x = TRUE)
sp <- fread(file.path(OUTROOT, "A2_eqtl_sqtl/A2_splice_region_by_signal.tsv"))
X[, splice := sid %in% sp[any_cand_own > 0 | any_cand_other > 0, sid]]
X[, grp := fifelse(two == TRUE, route, NA_character_)]
X[, keep1 := { o <- order(-nlp); k <- rep(FALSE, .N); k[o[1]] <- TRUE; k }, by = .(gene, fifelse(is.na(grp), type, grp))]
fwrite(X, file.path(OUT, "B1v2_signal_patterns.tsv.gz"), sep = "\t")

tabs <- list()
add <- function(d, groups, lab, seed) { groups <- intersect(groups, unique(d$group)); if (length(groups) < 1) return(invisible())
  r <- ratio_table(d, groups, B, seed, lab)
  tabs[[length(tabs) + 1]] <<- rbind(r$groups[, .(label, group, n_signals, estimate, lo, hi, contrast = NA_character_, p_boot = NA_real_)],
    if (!is.null(r$contrasts)) r$contrasts[, .(label, group = NA_character_, n_signals = NA_real_, estimate = difference, lo, hi, contrast, p_boot)]) }
cls <- c("mu", "p0", "p1"); rts <- c("onoff_S_minor", "onoff_S_major", "mu")
sets <- list(all = X, exact_mirror = X[exact_mirror == TRUE], single_signal_units = X[n_unit == 1], no_splice = X[splice == FALSE],
             one_per_gene_group = X[keep1 == TRUE], strict = X[n_unit == 1 & splice == FALSE & keep1 == TRUE])
for (nm in names(sets)) { d0 <- sets[[nm]]; A <- d0[two == TRUE]; sdd <- 10L * match(nm, names(sets))
  for (p in c("S_only", "reciprocal_competition", "coactivation", "N_only"))
    add(A[, .(gene, group = route, num = as.numeric(patternS == p), den = 1)], rts, paste0(nm, " | aligned pattern: ", p), sdd + 1L)
  if (nm == "all") for (p in c("S_only", "reciprocal_competition", "coactivation", "N_only"))
    add(A[, .(gene, group = route, num = as.numeric(patternS_perm == p), den = 1)], rts, paste0(nm, " | NULL permuted: aligned pattern: ", p), sdd + 2L)
  for (cc in c("opposite", "same", "equivalent_zero", "undetermined")) {
    add(A[!is.na(N_cat), .(gene, group = route, num = as.numeric(N_cat == cc), den = 1)], rts, paste0(nm, " | N (all donors, partner-specific): ", cc), sdd + 3L)
    add(A[!is.na(Nlev_cat), .(gene, group = route, num = as.numeric(Nlev_cat == cc), den = 1)], rts, paste0(nm, " | N level (donors with N > 0): ", cc), sdd + 4L)
    add(A[!is.na(T_cat), .(gene, group = route, num = as.numeric(T_cat == cc), den = 1)], rts, paste0(nm, " | gene output vs S usage: ", cc), sdd + 5L) }
  if (nm == "all") { add(A[!is.na(N_cat05), .(gene, group = route, num = as.numeric(N_cat05 == "equivalent_zero"), den = 1)], rts, "all | N equivalent to zero at delta 0.05", sdd + 6L)
    add(A[!is.na(N_cat_perm), .(gene, group = route, num = as.numeric(N_cat_perm == "opposite"), den = 1)], rts, "all | NULL permuted: N opposite", sdd + 7L) }
  add(A[!is.na(conc_S_c), .(gene, group = route, num = as.numeric(conc_S_c), den = 1)], rts, paste0(nm, " | A2 eQTL concordant with S usage (mirrored)"), sdd + 8L)
  add(A[!is.na(conc_S_meas), .(gene, group = route, num = as.numeric(conc_S_meas), den = 1)], rts, paste0(nm, " | A2 eQTL concordant with S usage (measured partner usage)"), sdd + 9L)
  add(A[!is.na(f_S) & f_S + f_N > 0, .(gene, group = route, num = f_N, den = f_S + f_N)], rts, paste0(nm, " | A1 share at N among promoter-proximal candidates"), sdd + 10L)
  add(d0[!is.na(conc_own_c), .(gene, group = type, num = as.numeric(conc_own_c), den = 1)], cls, paste0(nm, " | A2 eQTL concordant with own usage"), sdd + 11L)
  for (p in c("target_only", "other_only", "reciprocal")) add(d0[, .(gene, group = type, num = as.numeric(pattern == p), den = 1)], cls, paste0(nm, " | own-frame pattern: ", p), sdd + 12L)
}
## measured partner usage agrees with the mirrored direction (S = partner)
chk <- X[two == TRUE & S_is_own == FALSE & abs(z_uP) >= ZC, .(n = .N, agree = mean(sign(z_uP) == -dir_c)), by = .(route, exact_mirror)]
print(chk); writeT(chk, file.path(OUT, "B1v2_partner_usage_direction_check.tsv"))
TAB <- rbindlist(tabs); writeT(TAB, file.path(OUT, "B1v2_summary.tsv"))
fmtrow <- function(x) if (is.na(x$contrast)) sprintf("  - %s: %s [%s, %s] (n = %s)", x$group, pct(x$estimate), pct(x$lo), pct(x$hi), format(round(x$n_signals, 1))) else
  sprintf("  - %s: %+.1f pp [%+.1f, %+.1f], P = %s", x$contrast, 100 * x$estimate, 100 * x$lo, 100 * x$hi, formatC(x$p_boot, format = "g", digits = 2))
L <- c("# B1 v2 -- partner-specific, conditional refit; equivalence-based 'flat'; sensitivity sets as intersections", "",
       sprintf("Signals: %d (two-promoter aligned frame: %d). Direction: conditional for multi-signal units.", nrow(X), X[two == TRUE, .N]), "")
for (lab in unique(TAB$label)) { L <- c(L, paste0("- ", lab)); x <- TAB[label == lab]; for (i in seq_len(nrow(x))) L <- c(L, fmtrow(x[i])) }
writeLines(L, file.path(OUT, "B1v2_report.md")); cat(L[1:min(length(L), 400)], sep = "\n")
