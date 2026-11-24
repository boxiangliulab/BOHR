#!/usr/bin/env Rscript
### Source Data for Fig. 2B and Supplementary Table 2f: power by the per-allele odds ratio of the simulated effect.
### An effect of 0.5 on the probability scale reads easily as a 50% increase; an odds ratio does not. Calls, tests, thresholds and promoters are
### exactly those of the probability-scale version (prepare_power_by_effect.R -> power_by_effect_N<N>/per_promoter.csv, which
### reproduces Fig. 2A at this N); this script only adds the true per-allele odds ratio of every promoter and re-bins by it.
###
### Truth. Both generators run with effect_scale = "prob" (M1 newsim_direct_nothreshold_indep.R arg 6 default, run by
### newmethod_bohr2/scripts/pbs_m1nothr_g2000.sh with args "N 2000 2"; M2 v12_simulation_nothr.R line 471). The scaled effect
### (effect_size after the proportional scale factor) is added once per allele (dosage d = 0, 1, 2) on the probability scale to
### μ (mean of the intermediate values), p0 and/or p1, the other boundary capped at 1 - p_affected - 0.1, then the safety clamps
### and the BEINF conversion (M1 lines 279-309 and 330-334; M2 lines 147-192). par_at() below repeats those lines and returns the
### realised μ, P(y = 0) = nu / (1 + nu + tau) and P(y = 1) = tau / (1 + nu + tau) at each dosage.
### Per-allele OR = odds(param at d = 1) / odds(param at d = 0), odds(x) = x / (1 - x); the d = 2 vs 0 OR is also given (not the
### square, because the effect is additive on the probability scale). |log OR| pools increases and decreases. Mixed effects
### (μ + p0 on the simulated promoter): both ORs are given and the promoter is binned by the larger |log OR|.
### M2: parameters are those of transcript 1 (tx info); a promoter at transcript position 2 carries the mirror (y = 1 - y1), so
### its μ odds invert and p0 <-> p1 (the relabel() rule of prepare_data_v3.R); |log OR| is unchanged by the mirror.
### Bins (per-allele OR, left-closed): 1-1.25, 1.25-1.5, 1.5-2, 2-3, 3-5, >= 5. Points at the bin's geometric-mean OR.
### Output: data_v3_<DATASET>/power_by_OR_N<N>/ {per_promoter_OR.csv, power_by_OR_bin.csv, CHECKS.txt}. Reads only.
suppressMessages(library(data.table))
DATASET <- Sys.getenv("DATASET", "nothr_realjunc"); NN <- as.integer(Sys.getenv("EFFECT_N", "300"))
FIGROOT <- "/projects/group/bliu/e1327989/GTEx/new_simulation/figures_bohr2"
NM <- "/projects/group/bliu/e1327989/GTEx/new_simulation/newmethod_bohr2"
COMBINE <- Sys.getenv("COMBINE", "simes"); stopifnot(COMBINE %in% c("simes", "fisher", "fisher_sidak")); CSFX <- c(simes = "", fisher = "_fisher", fisher_sidak = "_fishersidak")[[COMBINE]]   # as prepare_data_v3.R
DATA <- file.path(FIGROOT, paste0("data_v3_", DATASET, CSFX)); stopifnot(DATASET %in% c("nothr_realjunc", "E80_nofloor_reads5", "E80_nofloor"))   # E80_*: same design tables
## OR_BINS (Fig. 2B: linear x axis, OR <= 5): "main" (default) = bins 1-1.25-1.5-2-3-5-Inf -> power_by_OR_N<N>/ (unchanged);
## "lin5" = the 3-5 bin split into 3-4 and 4-5 (>= 5 kept in the table, not drawn in Fig. 2B) -> power_by_OR_N<N>_lin5/. Same promoters and calls.
OR_BINS <- Sys.getenv("OR_BINS", "main"); stopifnot(OR_BINS %in% c("main", "lin5"))
IND <- file.path(DATA, sprintf("power_by_effect_N%d", NN)); OUTD <- file.path(DATA, sprintf("power_by_OR_N%d%s", NN, if (OR_BINS == "lin5") "_lin5" else "")); dir.create(OUTD, showWarnings = FALSE)
M1GT <- file.path(NM, sprintf("m1_nothr_g2000/sim_N%d_G2000_realfloor", NN))
M2TX <- function(b) file.path(NM, sprintf("m2_nothr_realjunc/%s/sim_%d/test_output/sim_tx_info_extended.txt", b, NN))

## realised BEINF parameters at dosage d (generator prob-scale branch, verbatim logic)
par_at <- function(et, p0, p1, mid, eff, dmu, dnu, dtau, d) {
  ev <- eff * d; cm <- mid; c0 <- p0; c1 <- p1
  if (et == "middle_value") cm <- mid + dmu * ev
  else if (et == "zero_proportion") { c0 <- p0 + dnu * ev; c1 <- max(min(p1, 1 - c0 - 0.1), 0) }
  else if (et == "one_proportion") { c1 <- p1 + dtau * ev; c0 <- max(min(p0, 1 - c1 - 0.1), 0) }
  else { cm <- mid + dmu * ev; c0 <- p0 + dnu * ev; c1 <- max(min(p1, 1 - c0 - 0.1), 0) }   # mixed_effect
  cm <- min(max(cm, 0.01), 0.99); c0 <- min(max(c0, 0.001), 0.8); c1 <- min(max(c1, 0), 0.8)
  pm <- max(1 - c0 - c1, 0.01); nu <- max(c0 / pm, 1e-4); tau <- max(c1 / pm, 1e-4)
  c(mu = cm, P0 = nu / (1 + nu + tau), P1 = tau / (1 + nu + tau))
}
m1 <- fread(file.path(M1GT, "ground_truth.txt"))
G <- rbind(m1[, .(m = "M1", promoter = gene_id, tx_pos = 1L, et = effect_type, eff = effect_size, p0 = p_zero, p1 = p_one, mid = middle_mean,
                  dmu = dir_mu, dnu = dir_nu, dtau = dir_tau)],
           rbindlist(lapply(c("b1", "b2"), function(b) fread(M2TX(b))[!is.na(promoter_id), .(m = "M2", promoter = sprintf("%s:PromoterID_%s", b, promoter_id),
                  tx_pos = transcript_position, et = effect_type, eff = effect_size, p0 = p_zero, p1 = p_one, mid = middle_mean,
                  dmu = sign(true_mu_effect), dnu = sign(true_nu_effect), dtau = sign(true_tau_effect))])))
G <- unique(G, by = c("m", "promoter"))
P <- G[, { a <- rbind(par_at(et, p0, p1, mid, eff, dmu, dnu, dtau, 0), par_at(et, p0, p1, mid, eff, dmu, dnu, dtau, 1), par_at(et, p0, p1, mid, eff, dmu, dnu, dtau, 2))
  if (tx_pos == 2) a <- cbind(mu = 1 - a[, "mu"], P0 = a[, "P1"], P1 = a[, "P0"])                                  # promoter-level mirror
  .(mu_d0 = a[1, 1], mu_d1 = a[2, 1], mu_d2 = a[3, 1], P0_d0 = a[1, 2], P0_d1 = a[2, 2], P0_d2 = a[3, 2], P1_d0 = a[1, 3], P1_d1 = a[2, 3], P1_d2 = a[3, 3]) },
  by = .(m, promoter, tx_pos, et)]
lo <- function(x) log(x / (1 - x))
P[, `:=`(lOR_mu = lo(mu_d1) - lo(mu_d0), lOR_p0 = lo(P0_d1) - lo(P0_d0), lOR_p1 = lo(P1_d1) - lo(P1_d0),
         lOR2_mu = lo(mu_d2) - lo(mu_d0), lOR2_p0 = lo(P0_d2) - lo(P0_d0), lOR2_p1 = lo(P1_d2) - lo(P1_d0))]
## promoter-level effect type (prepare_data_v3.R relabel: position 2 swaps zero <-> one)
SW <- c(middle_value = "mu", zero_proportion = "p0", one_proportion = "p1", mixed_effect = "mixed")
P[, ptype := SW[et]][tx_pos == 2 & ptype %in% c("p0", "p1"), ptype := c(p0 = "p1", p1 = "p0")[ptype]]
P[, mixed_boundary := fifelse(et == "mixed_effect", fifelse(tx_pos == 2, "p1", "p0"), NA_character_)]
## binning |log OR|: the affected parameter; mixed = the larger of μ and its boundary
P[, `:=`(abs_lOR_bin = fcase(ptype == "mu", abs(lOR_mu), ptype == "p0", abs(lOR_p0), ptype == "p1", abs(lOR_p1),
                             ptype == "mixed", pmax(abs(lOR_mu), fifelse(mixed_boundary == "p0", abs(lOR_p0), abs(lOR_p1)))),
         bin_by = fcase(ptype %in% c("mu", "p0", "p1"), ptype,
                        ptype == "mixed", fifelse(abs(lOR_mu) >= fifelse(mixed_boundary == "p0", abs(lOR_p0), abs(lOR_p1)), "mu", mixed_boundary)))]
## other boundary moved by the min() cap (p0 or p1 effects, and the boundary part of mixed effects)
P[, other_boundary := fcase(ptype == "p0", "p1", ptype == "p1", "p0", ptype == "mixed", fifelse(mixed_boundary == "p0", "p1", "p0"), default = NA_character_)]
P[, other_capped_d2 := !is.na(other_boundary) & abs(fifelse(other_boundary == "p0", lOR2_p0, lOR2_p1)) > 1e-9]

pp <- fread(file.path(IND, "per_promoter.csv"))
FWN <- c(M1 = "Simulated activity", M2 = "Simulated RNA-seq")
P[, framework := FWN[m]]
X <- merge(pp, P, by = c("framework", "promoter"), all.x = TRUE)
BR <- c(1, 1.25, 1.5, 2, 3, 5, Inf); BL <- c("1–1.25", "1.25–1.5", "1.5–2", "2–3", "3–5", "≥ 5")
if (OR_BINS == "lin5") { BR <- c(1, 1.25, 1.5, 2, 3, 4, 5, Inf); BL <- c("1–1.25", "1.25–1.5", "1.5–2", "2–3", "3–4", "4–5", "≥ 5") }
X[, OR_bin := factor(BL[findInterval(exp(abs_lOR_bin), BR)], levels = BL)]
e <- function(v) exp(v)
fwrite(X[order(framework, model, promoter), .(framework, N, promoter, effect_type, simulated_effect_0to2, transcript_position_M2 = fifelse(m == "M2", tx_pos, NA_integer_),
    mu_d0, mu_d1, mu_d2, P0_d0, P0_d1, P0_d2, P1_d0, P1_d1, P1_d2,
    OR_mu_per_allele = e(lOR_mu), OR_p0_per_allele = e(lOR_p0), OR_p1_per_allele = e(lOR_p1),
    OR_mu_d2_vs_d0 = e(lOR2_mu), OR_p0_d2_vs_d0 = e(lOR2_p0), OR_p1_d2_vs_d0 = e(lOR2_p1),
    binned_by = bin_by, abs_OR_per_allele = e(abs_lOR_bin), other_boundary_capped = other_capped_d2, OR_bin,
    model, test, p, threshold_alpha05, called)], file.path(OUTD, "per_promoter_OR.csv"))

wilson <- function(k, n, z = qnorm(.975)) { p <- k / n; c0 <- p + z^2 / (2 * n); s <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)); d <- 1 + z^2 / n
  list((c0 - s) / d, (c0 + s) / d) }
S <- X[, .(promoters = .N, called = sum(called), gm_OR = exp(mean(abs_lOR_bin)), median_OR = exp(median(abs_lOR_bin)),
           min_OR = exp(min(abs_lOR_bin)), max_OR = exp(max(abs_lOR_bin)), mean_effect_0to2 = mean(simulated_effect_0to2),
           min_effect_0to2 = min(simulated_effect_0to2), max_effect_0to2 = max(simulated_effect_0to2),
           n_mu = sum(effect_type == "mu"), n_p0 = sum(effect_type == "p0"), n_p1 = sum(effect_type == "p1"), n_mixed = sum(effect_type == "mixed"),
           n_other_capped = sum(other_capped_d2), threshold = unique(threshold_alpha05)), by = .(framework, model, OR_bin)][order(framework, model, OR_bin)]
S[, power := called / promoters][, c("ci_lo", "ci_hi") := wilson(called, promoters)][, share := promoters / sum(promoters), by = .(framework, model)]
fwrite(S[, .(framework, N = NN, model, OR_bin, promoters, share_of_true_effect_promoters_pct = 100 * share, called, power_pct = 100 * power,
             ci95_lo_pct = 100 * ci_lo, ci95_hi_pct = 100 * ci_hi, geometric_mean_OR = gm_OR, median_OR, min_OR, max_OR,
             mean_effect_0to2, min_effect_0to2, max_effect_0to2, promoters_mu = n_mu, promoters_p0 = n_p0, promoters_p1 = n_p1,
             promoters_mixed = n_mixed, promoters_other_boundary_capped = n_other_capped, threshold_alpha05 = threshold)],
       file.path(OUTD, "power_by_OR_bin.csv"))

## ---------------- checks ----------------
C <- character(0); ok <- function(nm, cond) { C <<- c(C, sprintf("[%s] %s", if (isTRUE(cond)) "PASS" else "FAIL", nm)); invisible(cond) }
ok(sprintf("every promoter matched to its generator parameters (%d of %d rows)", X[!is.na(abs_lOR_bin), .N], nrow(X)), X[, !anyNA(abs_lOR_bin) & !anyNA(OR_bin)])
ok("promoter-level effect type from the generator (+ mirror) = Fig. 2B Source Data effect_type, every promoter", X[, all(ptype == effect_type)])
ok("calls, P values and thresholds carried over unchanged from power_by_effect per_promoter.csv", nrow(X) == nrow(pp) && X[, uniqueN(.SD), .SDcols = c("framework", "promoter", "model")] == nrow(pp))
pw <- fread(file.path(DATA, "power.csv"))[set == "all" & N == NN]
ck <- merge(X[, .(scan = mean(called), n = .N), by = .(m, model)], pw[, .(m, model, power_cal05, promoters)], by = c("m", "model"))
ok(sprintf("pooled over OR bins = power.csv (Fig. 2A) at N = %d, %d framework x model cells", NN, nrow(ck)), ck[, all(abs(scan - power_cal05) < 1e-12 & n == promoters)] && nrow(ck) == 8)
## direction: the affected parameter moves in the generator's direction (M1 dirs; promoter level)
if (OR_BINS == "lin5") {   # the split bins add up to the main 3-5 bin; every other bin identical to the main run (power_by_OR_N<N>/power_by_OR_bin.csv)
  M0 <- fread(file.path(DATA, sprintf("power_by_OR_N%d", NN), "power_by_OR_bin.csv"), encoding = "UTF-8")
  S5 <- copy(S)[, OR_bin := as.character(OR_bin)][OR_bin %in% c("3–4", "4–5"), OR_bin := "3–5"][, .(promoters = sum(promoters), called = sum(called)), by = .(framework, model, OR_bin)]
  cm <- merge(S5, M0[, .(framework, model, OR_bin, p0 = promoters, c0 = called)], by = c("framework", "model", "OR_bin"))
  ok(sprintf("lin5: 3-4 + 4-5 = main 3-5 bin and every other bin = main run (%d framework x model x bin cells; promoters, calls)", nrow(cm)),
     nrow(cm) == nrow(M0) && cm[, all(promoters == p0 & called == c0)])
  ok("lin5: no promoter has OR exactly 5 (the 4-5 bin is OR <= 5; the >= 5 bin is OR > 5)", X[, !any(abs(exp(abs_lOR_bin) - 5) < 1e-12)])
}
ok("|log OR| > 0 except where the clamps absorb the whole effect", X[model == "BOHR", sum(abs_lOR_bin == 0)] <= 10)
C <- c(C, sprintf("  note: %d promoters have OR exactly 1 (p0 or p1 baseline below the 0.001 clamp, effect scaled to ~0)", X[model == "BOHR", sum(abs_lOR_bin == 0)]))
## M1 generator replay against the simulated data: zeros, ones and interior means by gene x dosage
e1 <- new.env(); load(file.path(M1GT, "sim_data.RData"), envir = e1); A <- e1$activity_mat; Gm <- e1$geno_mat
V <- rbindlist(lapply(colnames(A), function(g) { y <- A[, g]; d <- Gm[, g]
  data.table(promoter = g, d = 0:2, n = tabulate(d + 1, 3), z = sapply(0:2, function(k) sum(y[d == k] == 0)), o = sapply(0:2, function(k) sum(y[d == k] == 1)),
             mid = sapply(0:2, function(k) { v <- y[d == k & y > 0 & y < 1]; if (length(v)) mean(v) else NA_real_ })) }))
Pm <- melt(P[m == "M1", .(promoter, mu_d0, mu_d1, mu_d2, P0_d0, P0_d1, P0_d2, P1_d0, P1_d1, P1_d2)], id.vars = "promoter")
Pm[, `:=`(par = sub("_d.*", "", variable), d = as.integer(sub(".*_d", "", variable)))]
V <- merge(V, dcast(Pm, promoter + d ~ par, value.var = "value"), by = c("promoter", "d"))
zsum <- V[, .(obs0 = sum(z), exp0 = sum(n * P0), obs1 = sum(o), exp1 = sum(n * P1))]
r_mu <- V[!is.na(mid), cor(mid, mu)]; r0 <- V[n > 0, cor(z / n, P0)]; r1 <- V[n > 0, cor(o / n, P1)]
ok(sprintf("M1 replay, N = %d, 2000 genes x 3 dosages: observed zeros %d vs expected %.0f, ones %d vs expected %.0f (each within 2%%); r(observed share, P0) = %.4f, r(share of ones, P1) = %.4f, r(interior mean, μ) = %.4f",
           NN, zsum$obs0, zsum$exp0, zsum$obs1, zsum$exp1, r0, r1, r_mu),
   abs(zsum$obs0 / zsum$exp0 - 1) < .02 && abs(zsum$obs1 / zsum$exp1 - 1) < .02 && r0 > .97 && r1 > .97 && r_mu > .97)
## sensitivity: bin by the per-allele geometric mean of the dosage 2 vs 0 step, sqrt(OR_2v0)
X[, abs_lOR2h := fcase(ptype == "mu", abs(lOR2_mu), ptype == "p0", abs(lOR2_p0), ptype == "p1", abs(lOR2_p1),
                       ptype == "mixed", pmax(abs(lOR2_mu), fifelse(mixed_boundary == "p0", abs(lOR2_p0), abs(lOR2_p1)))) / 2]
X[, OR_bin2 := factor(BL[findInterval(exp(abs_lOR2h), BR)], levels = BL)]
same <- X[model == "BOHR", mean(OR_bin == OR_bin2)]
s2 <- dcast(X[model %in% c("BOHR", "raw"), .(pw = 100 * mean(called)), by = .(framework, OR_bin2, model)], framework + OR_bin2 ~ model, value.var = "pw")[, gain := BOHR - raw]
C <- c(C, sprintf("  sensitivity: binned by sqrt(OR d2 vs d0) instead, %.1f%% of promoters keep their bin; largest BOHR gain then in bin %s",
                  100 * same, paste(s2[, .SD[which.max(gain), sprintf("%s %s (%.1f vs %.1f)", framework, OR_bin2, BOHR, raw)], by = framework]$V1, collapse = "; ")))
chk <- c(sprintf("prepare_power_by_OR.R, DATASET %s, COMBINE %s, N = %d, OR_BINS %s, %s", DATASET, COMBINE, NN, OR_BINS, format(Sys.time(), "%Y-%m-%d %H:%M")), C)
writeLines(chk, file.path(OUTD, "CHECKS.txt")); cat(chk, sep = "\n"); stopifnot(!any(grepl("^\\[FAIL", C)))
print(dcast(S[model %in% c("BOHR", "raw")], framework + OR_bin ~ model, value.var = c("power", "promoters", "share", "gm_OR"))[,
  .(framework, OR_bin, n = promoters_BOHR, share = round(100 * share_BOHR, 1), gmOR = round(gm_OR_BOHR, 2), BOHR = round(100 * power_BOHR, 1), raw = round(100 * power_raw, 1),
    gain = round(100 * (power_BOHR - power_raw), 1))])
