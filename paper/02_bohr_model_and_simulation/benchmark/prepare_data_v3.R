#!/usr/bin/env Rscript
### Source Data v3 for Figure 2. No filtering in the simulation (no generation thresholds, no
### eligibility/unit rules -- no BH is involved; every existing unit is tested); the linear comparators are the three
### QTLtools modes computed on the SAME samples and the SAME pre-drawn permutations as BOHR (newmethod_bohr2/scripts/baselines_v3.R):
###   RINTfirst (GTEx order, z = qnorm((rank - 0.5)/n)), raw (QTLtools default), normal (QTLtools --normal, max(rank));
###   plus, in calibration only, a covariate demonstration (ancestry-like covariate C, rho = 0.5, permuted with the genotype).
### BOHR tests: mu = F-corrected LRT (F = (n_mid - 2) expm1(W / n_mid) ~ F(1, n_mid - 2)); p0 / p1 = Firth PLR. Panels a-c use the
### promoter-level P of the SNP (COMBINE below: Simes(mu, p0, p1) by default; Fisher with a Sidak 0/1 part in the final runs), judged at a permutation-calibrated
### threshold (no BH). The mu test alone is kept as view "mu_alone" (supplement); min p ("BOHR_any") only for the record.
### DATASET env: nothr (default; M1 = 2000 genes with independent genotypes; M2 = gene batches b1 + b2 pooled) | realfloor | nothr_realjunc (M2 junctions re-filtered with the real rule).
### M2 promoters come in two-promoter genes; exact complements (y_B = 1 - y_A in every sample, never both 0) carry ONE
### phenotype, so the minor promoter (lower mean usage across the simulated samples) is dropped (Figure 4 mirror rule, exact_mirror on both components);
### types and truths are expressed relative to the kept promoter (transcript_position 2 mirrored: zero <-> one, mu sign).
### Every calibration / power / effect statistic is computed on the de-duplicated set. M1 has one promoter per gene.
### Up-front filter, applied to BOTH methods: a promoter x N unit enters only if n_mid >= 5 (interior
### samples among those kept by prepare_y). Below that the mu beta regression (intercept, SNP, precision) is not estimable
### (n_mid = 4 gave perfect fits and null p < 1e-8); the unit is removed for BOHR and for every QTLtools mode alike.
suppressMessages(library(data.table))
DATASET <- Sys.getenv("DATASET", "nothr"); stopifnot(DATASET %in% c("nothr", "realfloor", "nothr_realjunc", "E80_paper", "E80_nofloor", "E80_nofloor_reads5"))
## E80_*: the same nothr_realjunc design (the same simulated genes and reads) with the promoter rules of a version
## (newmethod_bohr2/scripts/pbs_e80_sim.sh -> results/e80_20261119/<paper|nofloor>/):
##   E80_paper          : junction floor (unique reads >= 5), endpoint skip 80%, sd <= 0.05 skip
##   E80_nofloor        : no junction floor (m2_nofloor, all junctions), endpoint skip 80%, no sd <= 0.05 skip (sd < 0.005 skip)
##   E80_nofloor_reads5 : (the final version) E80_nofloor, then M2 genes with median gene reads < 5 (m2_nofloor/<B>/sim_<N>/gene_reads.tsv) dropped
##                        with their promoters (calibration and power alike). M1 has no reads: M1 = E80_nofloor.
E80 <- startsWith(DATASET, "E80_"); E80V <- if (E80) (if (DATASET == "E80_paper") "paper" else "nofloor") else NA_character_
READS_MIN <- if (DATASET == "E80_nofloor_reads5") 5 else NA_real_
## nothr_realjunc: M2 junctions re-filtered with the real-data rule (unique reads >= 5)
NM  <- "/projects/group/bliu/e1327989/GTEx/new_simulation/newmethod_bohr2"
CAL <- "/projects/group/bliu/e1327989/GTEx/new_simulation/calibration_perm1000"
SEEDTAG <- Sys.getenv("SEEDTAG", "")   # "" = original permutations (seed base 20261001); "_seed2" = seed check (20262025), nothr only
stopifnot(SEEDTAG == "" || DATASET %in% c("nothr", "nothr_realjunc", "E80_nofloor", "E80_nofloor_reads5"))
## COMBINE: the promoter-level BOHR P of panels A-C, ED 2B and the
## power-by-effect tables. "simes" (default; the earlier runs): Simes over the SNP's μ, p0 and p1 tests. "fisher": Fisher's combination of the
## μ part and the 0/1 part (0/1 part = the smaller of the p0 / p1 P values x the number of 0/1 tests; the one part when only one exists), the rule of
## the Fisher version of Fig. 2E. The column keeps the name BOHR_simes so every downstream script runs unchanged; output folder data_v3_<DATASET>_fisher.
## "fisher_sidak" (the final runs): the same Fisher combination with the Šidák 0/1 part, 1 - (1 - smaller P)^k, in place of
## k x smaller P (Bonferroni); output folder data_v3_<DATASET>_fishersidak. "simes" and "fisher" are unchanged.
COMBINE <- Sys.getenv("COMBINE", "simes"); stopifnot(COMBINE %in% c("simes", "fisher", "fisher_sidak")); CSFX <- c(simes = "", fisher = "_fisher", fisher_sidak = "_fishersidak")[[COMBINE]]
fisher2 <- function(mu, p0, p1) { kb <- (!is.na(p0)) + (!is.na(p1)); pb <- suppressWarnings(pmin(p0, p1, na.rm = TRUE)); b <- ifelse(kb > 0, if (COMBINE == "fisher_sidak") -expm1(kb * log1p(-pb)) else pmin(1, kb * pb), NA_real_)
  ifelse(!is.na(mu) & !is.na(b), pchisq(-2 * (log(mu) + log(b)), 4, lower.tail = FALSE), ifelse(!is.na(mu), mu, b)) }
OUT <- file.path("/projects/group/bliu/e1327989/GTEx/new_simulation/figures_bohr2", paste0("data_v3_", DATASET, SEEDTAG, CSFX)); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
Ns <- seq(50, 600, 50); PK <- 2; MIN_NMID <- as.integer(Sys.getenv("MIN_NMID", "5"))
METH <- c(M1 = "Method 1 (direct)", M2 = "Method 2 (polyester)")
## sources: one entry per (method, batch)
M2DIR <- if (DATASET == "nothr_realjunc") "m2_nothr_realjunc" else "m2_nothr"; M2RES <- if (DATASET == "nothr_realjunc") "m2nothr_realjunc" else "m2nothr"
M2B3 <- if (DATASET == "nothr_realjunc") "nothr_rj_M2" else "nothr_M2"
SRC <- if (DATASET %in% c("nothr", "nothr_realjunc")) list(
  list(m = "M1", b = "",   cal = "results/m1nothr_g2000/calib", pow = "results/m1nothr_g2000/power", b3 = "results/baselines_v3/nothr_M1"),
  list(m = "M2", b = "b1", cal = sprintf("results/%s/b1/calib", M2RES), pow = sprintf("results/%s/b1/power", M2RES), b3 = sprintf("results/baselines_v3/%sb1", M2B3),
       pheno = function(N) file.path(NM, sprintf("%s/b1/sim_%d", M2DIR, N))),
  list(m = "M2", b = "b2", cal = sprintf("results/%s/b2/calib", M2RES), pow = sprintf("results/%s/b2/power", M2RES), b3 = sprintf("results/baselines_v3/%sb2", M2B3),
       pheno = function(N) file.path(NM, sprintf("%s/b2/sim_%d", M2DIR, N)))) else list(
  list(m = "M1", b = "",   cal = "results/m1floor/L0.005/calib", pow = "results/m1floor/L0.005/power", b3 = "results/baselines_v3/realfloor_M1"),
  list(m = "M2", b = "",   cal = "results/calib/M2",             pow = "results/power_eff/M2",         b3 = "results/baselines_v3/realfloor_M2",
       pheno = function(N) file.path(CAL, sprintf("big_pheno/M2_N%d", N))))
if (E80) { R0 <- file.path("results/e80_20261119", E80V); M2P <- if (E80V == "paper") "m2_nothr_realjunc" else "m2_nofloor"
  SRC <- c(list(list(m = "M1", b = "", cal = file.path(R0, "m1/calib"), pow = file.path(R0, "m1/power"), b3 = file.path(R0, "baselines_M1"))),
           lapply(c("b1", "b2"), function(bb) { force(bb); list(m = "M2", b = bb, cal = file.path(R0, bb, "calib"), pow = file.path(R0, bb, "power"),
             b3 = file.path(R0, paste0("baselines_M2", bb)), pheno = function(N) file.path(NM, sprintf("%s/%s/sim_%d", M2P, bb, N))) })) }
## baselines layout of the E80 runs: <b3>/{calib,power}/N<N>.rds, as baselines_v3/<x>/{calib,power}
## gene-reads filter (E80_nofloor_reads5): promoter ids of M2 genes with median reads >= READS_MIN in that (batch, N) run
reads_keep <- function(src, N) { d <- src$pheno(N); gr <- fread(file.path(d, "gene_reads.tsv"))
  tf <- file.path(d, "test_output", "sim_tx_info_extended.txt"); tx <- unique(fread(tf)[!is.na(promoter_id), .(id = paste0("PromoterID_", promoter_id), gene_id)])
  stopifnot(mean(unique(tx$gene_id) %in% gr$gene_id) > 0.9)          # same gene ids (gencode) in the design table and proActiv
  tx[gene_id %in% gr[median_reads >= READS_MIN, gene_id], id] }      # a gene absent from gene_reads has no counts, hence no tested promoter
SRC <- lapply(SRC, function(src) { src$b3cal <- if (!nzchar(SEEDTAG)) src$b3 else if (E80) paste0(src$b3, SEEDTAG) else sub("baselines_v3", paste0("baselines_v3", SEEDTAG), src$b3)   # E80: <b3>_seedK (pbs_e80_seedK.sh)
  if (nzchar(SEEDTAG)) src$cal <- paste0(src$cal, SEEDTAG); src })   # calibration runs of the chosen permutation seed; power runs are seed-free
FfromLog <- function(logp, n) { W <- qchisq(logp, 1, lower.tail = FALSE, log.p = TRUE); out <- rep(NA_real_, length(W)); ok <- !is.na(W) & n > PK
  out[ok] <- pf((n[ok] - PK) * expm1(W[ok] / n[ok]), 1, n[ok] - PK, lower.tail = FALSE); out }
tag <- function(b, id) if (nzchar(b)) paste0(b, ":", id) else id

## ---- M2 mirror map and promoter-level labels, per (batch, N) ----
mirror_map <- function(src, N) {
  d <- src$pheno(N); e <- new.env(); load(file.path(d, "recal_att_f.RData"), envir = e); X <- as.data.frame(e$recal_att_f)
  tf <- file.path(d, "test_output", "sim_tx_info_extended.txt"); if (!file.exists(tf)) tf <- file.path(d, "sim_tx_info_extended.txt")   # realfloor big_pheno layout
  tx <- fread(tf)[!is.na(promoter_id)]
  tx <- unique(tx[, .(id = paste0("PromoterID_", promoter_id), gene_id, tx_pos = transcript_position)], by = "id")[id %in% names(X)]
  pr <- tx[, if (.N == 2) .(id1 = id[1], id2 = id[2]), by = gene_id]
  pr[, c("exact", "canon") := { a <- X[[id1]]; b <- X[[id2]]; ok <- !is.na(a) & !is.na(b)
    list(all(abs(a[ok] + b[ok] - 1) < 1e-6) && !any(a[ok] == 0 & b[ok] == 0), if (mean(a[ok]) >= mean(b[ok])) id1 else id2) }, by = gene_id]
  mp <- merge(tx, pr[, .(gene_id, exact, canon)], by = "gene_id", all.x = TRUE)
  mp[, drop := !is.na(exact) & exact & id != canon][, .(id, gene_id, tx_pos, exact, canon, drop)]
}
dedup_ids <- function(ids, mp) {                   # to drop = minor promoters of exact pairs whose major twin is in `ids` (canon = major)
  d <- mp[drop == TRUE & id %in% ids & canon %in% ids, id]; setdiff(ids, d) }
relabel <- function(h, mp) {                       # promoter-level type and truths (pos-2 promoters carry the mirrored effect)
  h <- merge(h, mp[, .(id, tx_pos)], by = "id", all.x = TRUE); stopifnot(!anyNA(h$tx_pos))
  sw <- c(zero_proportion = "one_proportion", one_proportion = "zero_proportion", middle_value = "middle_value", mixed_effect = "mixed_effect")
  h[tx_pos == 2, effect_type := sw[effect_type]]
  h[tx_pos == 2, `:=`(true_mu_effect = -true_mu_effect, true_nu_effect = true_tau_effect, true_tau_effect = true_nu_effect)]; h[] }

## ================= nulls =================
NUL <- rbindlist(lapply(SRC, function(src) rbindlist(lapply(Ns, function(NN) {
  x <- readRDS(file.path(NM, src$cal, sprintf("N%d.rds", NN))); b3 <- readRDS(file.path(NM, src$b3cal, "calib", sprintf("N%d.rds", NN)))$res
  r <- merge(merge(x$res, x$header[, .(id, nmid = n - n0 - n1)], by = "id"), b3, by = c("id", "perm"))
  r <- r[nmid >= MIN_NMID]                                                            # up-front filter, both methods
  if (src$m == "M2") { mp <- mirror_map(src, NN); r <- r[id %in% dedup_ids(unique(r$id), mp)] }
  if (src$m == "M2" && !is.na(READS_MIN)) r <- r[id %in% reads_keep(src, NN)]                # E80_nofloor_reads5 gene-reads filter
  data.table(m = src$m, N = NN, id = tag(src$b, r$id), perm = r$perm, nmid = r$nmid,
             BOHR = FfromLog(r$mu_lrt_logp, r$nmid), BOHR_raw = exp(r$mu_lrt_logp), BOHR_p0 = exp(r$nu_firth_logp), BOHR_p1 = exp(r$tau_firth_logp),
             RINTfirst = r$rintfirst, raw = r$raw, normal = r$normal,
             cov_RINTfirst = r$cov_rintfirst, cov_raw = r$cov_raw, cov_normal = r$cov_normal, boundary_share = r$boundary_share) }), fill = TRUE)), fill = TRUE)
NUL[, BOHR_any := pmin(BOHR, BOHR_p0, BOHR_p1, na.rm = TRUE)]
## Promoter-level BOHR p for panels a-c: Simes over the promoter's existing tests for the
## SNP (mu-F, p0 Firth, p1 Firth) = the cross-layer rule of Fig 2e. min p ("BOHR_any") is kept only for the record (not a p-value).
simes3 <- function(a, b, c) { L <- data.table(r = rep(seq_along(a), 3), p = c(a, b, c))[!is.na(p)]; setorder(L, r, p)
  L[, `:=`(i = seq_len(.N), k = .N), by = r]; S <- L[, .(s = min(k * p / i)), by = r]
  out <- rep(NA_real_, length(a)); out[S$r] <- pmin(S$s, 1); out }
NUL[, BOHR_simes := simes3(BOHR, BOHR_p0, BOHR_p1)]
if (COMBINE %in% c("fisher", "fisher_sidak")) NUL[, BOHR_simes := fisher2(BOHR, BOHR_p0, BOHR_p1)]   # COMBINE: Fisher in place of Simes
cat(sprintf("up-front filter: n_mid >= %d (both methods)\n", MIN_NMID)); cat("calibration promoters per method x N (after mirror de-duplication and the n_mid filter):\n"); print(dcast(NUL[, .(k = uniqueN(id)), by = .(m, N)], m ~ N, value.var = "k"))
MODES <- c("RINTfirst", "raw", "normal")
## per-promoter FPR (mu test for BOHR, the whole-phenotype test for each linear mode), plus the covariate demonstration
VIEWS <- rbind(data.table(view = "main", model = c("BOHR", MODES), col = c("BOHR_simes", MODES)),     # panels b, c: BOHR = the promoter-level P (BOHR_simes, see COMBINE)
               data.table(view = "mu_alone", model = "BOHR", col = "BOHR"),                            # supplement: mu-F test alone
               data.table(view = "uncorrected", model = "BOHR", col = "BOHR_raw"),
               data.table(view = "covariate", model = MODES, col = paste0("cov_", MODES)))
FPR <- rbindlist(lapply(seq_len(nrow(VIEWS)), function(i) { v <- VIEWS[i]
  NUL[, .(fpr = mean(get(v$col) < .05, na.rm = TRUE), hits001 = sum(get(v$col) < 1e-3, na.rm = TRUE), tests = sum(!is.na(get(v$col)))), by = .(m, N, id)][tests > 0][
    , `:=`(view = v$view, model = v$model)] }))[, method := METH[m]]
fwrite(FPR, file.path(OUT, "fpr_perpromoter.csv"))
set.seed(20261005)
TAIL <- FPR[, { pooled <- .SD[, .(h = sum(hits001), K = sum(tests), f = sum(fpr * tests)), by = id]
  B <- replicate(2000, { s <- sample(nrow(pooled), replace = TRUE); sum(pooled$h[s]) / sum(pooled$K[s]) / 1e-3 })
  .(promoters = nrow(pooled), fpr05 = 100 * sum(pooled$f) / sum(pooled$K), tail = sum(pooled$h) / sum(pooled$K) / 1e-3,
    lo = quantile(B, .025, names = FALSE), hi = quantile(B, .975, names = FALSE)) }, by = .(method, view, model)]
TAILN <- FPR[, .(fpr05 = 100 * sum(fpr * tests) / sum(tests), tail = sum(hits001) / sum(tests) / 1e-3), by = .(method, N, view, model)]
fwrite(TAIL, file.path(OUT, "tail_pooled.csv")); fwrite(TAILN, file.path(OUT, "tail_byN.csv"))
qq_of <- function(p) { p <- sort(p[!is.na(p)]); n <- length(p)
  d <- data.table(rank = seq_len(n), n = n, observed = -log10(pmax(p, 1e-300)), expected = -log10(seq_len(n) / (n + 1)))
  d[, `:=`(band_lo = -log10(qbeta(.975, rank, n - rank + 1)), band_hi = -log10(qbeta(.025, rank, n - rank + 1)), lambda = median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(.5, 1))]
  d[(rank <= 3000) | (rank %% pmax(1L, as.integer(n / 4000)) == 0L)] }
QQ <- rbindlist(lapply(names(METH), function(mm) { X <- NUL[m == mm]; if (!nrow(X)) return(NULL)
  rbindlist(lapply(seq_len(nrow(VIEWS)), function(i) qq_of(X[[VIEWS$col[i]]])[, `:=`(view = VIEWS$view[i], model = VIEWS$model[i])]))[, method := METH[[mm]]] }))
fwrite(QQ, file.path(OUT, "qq.csv"))

## ================= thresholds: alpha-quantile of each statistic's own null (pooled calibration promoters), per method x N =================
STATS <- c("BOHR_simes", "BOHR_any", "BOHR", "BOHR_p0", "BOHR_p1", MODES)
TH <- melt(NUL[, c("m", "N", STATS), with = FALSE], id.vars = c("m", "N"), variable.name = "stat", value.name = "p")[
  , .(t05 = quantile(p, .05, na.rm = TRUE, type = 1), t001 = quantile(p, 1e-3, na.rm = TRUE, type = 1)), by = .(m, N, stat)][, stat := as.character(stat)]
fwrite(TH, file.path(OUT, "thresholds.csv"))

## ================= power + effects =================
lab <- c(middle_value = "mu", zero_proportion = "p0", one_proportion = "p1", mixed_effect = "mixed")
PWD <- rbindlist(lapply(SRC, function(src) rbindlist(lapply(Ns, function(NN) {
  x <- readRDS(file.path(NM, src$pow, sprintf("N%d.rds", NN))); b3 <- readRDS(file.path(NM, src$b3, "power", sprintf("N%d.rds", NN)))$res[, perm := NULL]
  h <- x$header[is.na(skip) & effect_size > 0 & effect_type %in% names(lab)]
  if (src$m == "M2") { mp <- mirror_map(src, NN); h <- relabel(h, mp); h <- h[(n - n0 - n1) >= MIN_NMID]; h <- h[id %in% dedup_ids(h$id, mp)] }
  if (src$m == "M2" && !is.na(READS_MIN)) h <- h[id %in% reads_keep(src, NN)]                # E80_nofloor_reads5 gene-reads filter
  h[, set := lab[effect_type]]
  r <- merge(merge(x$res, h[, .(id, set, nmid = n - n0 - n1, true_mu_effect, true_nu_effect, true_tau_effect)], by = "id"), b3, by = "id")
  r <- r[nmid >= MIN_NMID]                                                            # up-front filter, both methods
  data.table(m = src$m, N = NN, id = tag(src$b, r$id), set = r$set, true_mu = 2 * r$true_mu_effect, true_p0 = 2 * r$true_nu_effect, true_p1 = 2 * r$true_tau_effect,
             mu_pe = r$mu_pe, p0_pe = r$p0_pe, p1_pe = r$p1_pe, eff_RINTfirst = r$eff_rintfirst, eff_raw = r$eff_raw, eff_normal = r$eff_normal,
             BOHR = FfromLog(r$mu_lrt_logp, r$nmid), BOHR_p0 = exp(r$nu_firth_logp), BOHR_p1 = exp(r$tau_firth_logp),
             RINTfirst = r$rintfirst, raw = r$raw, normal = r$normal) }), fill = TRUE)), fill = TRUE)
PWD[, BOHR_any := pmin(BOHR, BOHR_p0, BOHR_p1, na.rm = TRUE)]
PWD[, BOHR_simes := simes3(BOHR, BOHR_p0, BOHR_p1)]
if (COMBINE %in% c("fisher", "fisher_sidak")) PWD[, BOHR_simes := fisher2(BOHR, BOHR_p0, BOHR_p1)]   # COMBINE: Fisher in place of Simes
own <- c(mu = "BOHR", p0 = "BOHR_p0", p1 = "BOHR_p1", mixed = "BOHR_simes", all = "BOHR_simes")   # all / mixed: the promoter-level P (Simes, or Fisher under COMBINE)
pw <- function(D, setname) { L <- rbind(D[, .(m, N, id, model = "BOHR", stat = own[[setname]], p = get(own[[setname]]))],
  rbindlist(lapply(MODES, function(md) D[, .(m, N, id, model = md, stat = md, p = get(md))])))
  L <- merge(L, TH, by = c("m", "N", "stat"))
  L[, .(promoters = uniqueN(id), power_cal05 = mean(!is.na(p) & p <= t05), power_cal001 = mean(!is.na(p) & p <= t001), power_nom05 = mean(!is.na(p) & p < .05)), by = .(m, N, model)][, set := setname] }
POW <- rbind(pw(PWD, "all"), rbindlist(lapply(c("mu", "p0", "p1", "mixed"), function(s) pw(PWD[set == s], s))))[, method := METH[m]]
fwrite(POW, file.path(OUT, "power.csv"))
EFF <- rbind(
  PWD[set == "mu", .(m, N, id, facet = "Middle (mu)", true = true_mu, BOHR = mu_pe, RINTfirst = eff_RINTfirst, raw = eff_raw, normal = eff_normal)],
  PWD[set == "mixed", .(m, N, id, facet = "Mixed (mu)", true = true_mu, BOHR = mu_pe, RINTfirst = eff_RINTfirst, raw = eff_raw, normal = eff_normal)],
  PWD[set == "p0", .(m, N, id, facet = "Zero (p0)", true = true_p0, BOHR = p0_pe, RINTfirst = -eff_RINTfirst, raw = -eff_raw, normal = -eff_normal)],
  PWD[set == "mixed" & true_p0 != 0, .(m, N, id, facet = "Mixed (p0)", true = true_p0, BOHR = p0_pe, RINTfirst = -eff_RINTfirst, raw = -eff_raw, normal = -eff_normal)],
  PWD[set == "p1", .(m, N, id, facet = "One (p1)", true = true_p1, BOHR = p1_pe, RINTfirst = eff_RINTfirst, raw = eff_raw, normal = eff_normal)],
  PWD[set == "mixed" & true_p1 != 0, .(m, N, id, facet = "Mixed (p1)", true = true_p1, BOHR = p1_pe, RINTfirst = eff_RINTfirst, raw = eff_raw, normal = eff_normal)])
EFF <- melt(EFF, id.vars = c("m", "N", "id", "facet", "true"), variable.name = "model", value.name = "estimated")[!is.na(estimated) & is.finite(estimated)]
EFF <- EFF[, if (uniqueN(model) == 4) .SD, by = .(m, N, id, facet)][, method := METH[m]]      # common records only
fwrite(EFF, file.path(OUT, "effect.csv"))

## ---- digest ----
cat("\n== true-effect promoters per method x N (de-duplicated) ==\n"); print(dcast(POW[set == "all" & model == "BOHR"], m ~ N, value.var = "promoters"))
cat("\n== a: power % (calibrated .05 | 1e-3), all true-effect promoters ==\n")
print(dcast(POW[set == "all" & N %in% c(50, 100, 200, 600)][, v := sprintf("%.1f | %.1f", 100 * power_cal05, 100 * power_cal001)], m + model ~ N, value.var = "v"))
cat("\n== fold BOHR / RINTfirst, calibrated .05, by set (range over 12 N) ==\n")
a <- dcast(POW, m + set + N ~ model, value.var = "power_cal05"); print(a[, .(min = round(min(BOHR / RINTfirst), 2), max = round(max(BOHR / RINTfirst), 2),
  cells_BOHR_lower = sum(BOHR < RINTfirst), cells_lower_than_any_mode = sum(BOHR < pmax(RINTfirst, raw, normal))), by = .(m, set)])
b <- dcast(POW, m + set + N ~ model, value.var = "power_cal001"); cat("cells BOHR < max(modes): .05", a[BOHR < pmax(RINTfirst, raw, normal), .N], "/", nrow(a), "; 1e-3", b[BOHR < pmax(RINTfirst, raw, normal), .N], "/", nrow(b), "\n")
cat("\n== calibration, pooled over N: FPR % | tail [95% CI] ==\n"); print(TAIL[order(method, view, model)][, .(method, view, model, promoters, fpr05 = round(fpr05, 2), tail = round(tail, 2), lo = round(lo, 2), hi = round(hi, 2))])
cat("\n== calibration by N (FPR % | tail) ==\n"); print(dcast(TAILN[N %in% c(50, 100, 200, 600)][, v := sprintf("%.2f | %.2f", fpr05, tail)], method + view + model ~ N, value.var = "v"))
cat("\n== effect r (pooled N, common records) ==\n"); print(dcast(EFF[, .(r = round(cor(true, estimated), 3)), by = .(method, facet, model)], method + facet ~ model, value.var = "r"))
