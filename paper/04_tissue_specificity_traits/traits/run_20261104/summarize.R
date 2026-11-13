## =============================================================================
## summarize.R -- summary of run 20261104 results (= the same fits as run 20261007, with mirror rule v3; see README.md)
## Usage: Rscript summarize.R [arms...]   (default real sim1; only completed tissue x trait are summarised)
##
## 1. Full base tables (author 2026-10-01: store the full table of base statistics first; later changes add columns to it and never change the original table)
##   full_<arm>.rds / .tsv.gz   : all raw rows of run_trait.R (one row per tissue x trait x promoter, including skipped rows and reasons), no filtering
##   units_<arm>.rds / .tsv.gz  : full table expanded into (promoter, type) units (including untestable, non-finite p, mirror duplicates), each unit with all flags
## 2. Tests (one 1-df test per (promoter, type))
##   mu : mu_F_logp (beta-part LRT + small-sample F correction); for binary traits mu_testable is required (>= 5 per class among the mid-value samples)
##   p0 : nu_firth_logp (Firth-penalised LRT); p1 : tau_firth_logp
## 3. Primary set and multiple testing
##   primary = testable (finite p) and not a duplicate unit of a mirror-pair partner (author rule; Figure 4 shared table v2, by component:
##             for interior_exact drop the partner's mu, for boundary_exact drop the partner's p0/p1; mirror table work/mirror.tsv)
##   author: mirroring is the biology of promoter switching itself; de-duplication only affects counts and BH; in plots the two promoters of one gene are drawn together
##   BH: within tissue x trait x type, on primary, q <= 0.05 (author rule: test everything first, no proportion gate)
##   partner units stay in the full table, duplicate_of points to the matching unit of the canonical (mu->mu, p0->p1, p1->p0);
##   the kept canonical unit carries merged_component (interior/boundary), twin, eq_class (<tissue|gene>|interior/canonical_off/canonical_full)
## 4. Display flags (not in BH): gate_prop (v13 proportion gate), n50 (n >= 50), pass_est (v13 link-scale effect gate)
## 5. Final rules (DECISIONS.md): final = primary & mu n_mid >= 80 & p0/p1 boundary proportion 0.10-0.90 (2026-10-05, unified with QTL)
##     & p0/p1 not completely separated & tissue donors >= 50;
##     q_final = BH within tissue x trait x type on final; sig_final is the final significant set (sig_<arm>.tsv); q / sig (all primary) are kept for comparison (sig_all_<arm>.tsv)
## Output (results/summary/): the full tables above, sig_<arm>.tsv, drops_<arm>.tsv, counts_<arm>_<type>.tsv (type = mu/p0/p1/union,
##       promoter counts on primary), calibration.tsv, null_discovery.tsv
## =============================================================================
source("/projects/group/bliu/e1327989/GTEx/tissue_regression/trait_regression_20261104/code/data.R")
a <- commandArgs(TRUE); ARMS <- if (length(a)) a else c("real", "sim1")
OUT <- file.path(PATHS$res, "summary"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
MIR <- fread(file.path(PATHS$work, "mirror_v3.tsv"))   # Mirror rule v3 (2026-11-04): for boundary_exact keep p0 of both members and drop p1 of both

load_arm <- function(arm) {
  fs <- Sys.glob(file.path(PATHS$res, arm, "*", "*.rds"))
  rbindlist(lapply(fs, function(f) { x <- readRDS(f); r <- x$res
    if (!nrow(r)) return(data.table(tissue = x$meta$tissue, trait = x$meta$trait, skip = x$meta$skip_trait, n_donors = x$meta$n_donors))
    r[, n_donors := x$meta$n_donors]; r }), fill = TRUE)
}
save_both <- function(x, stem) { saveRDS(x, file.path(OUT, paste0(stem, ".rds")))
  fwrite(x, file.path(OUT, paste0(stem, ".tsv.gz")), sep = "\t") }

drop_log <- function(R) {
  rbind(
    R[!is.na(skip), .(N = .N), by = .(reason = sub(":.*", "", skip))][, type := "all"],
    R[is.na(skip) & mu_testable %in% FALSE, .(type = "mu", reason = "mu_untestable_rare_mid", N = .N)],
    R[is.na(skip) & !(mu_testable %in% FALSE) & !is.finite(mu_F_logp), .(type = "mu", reason = "mu_nonfinite", N = .N)],
    R[is.na(skip) & !is.na(firth_error), .(type = "p0/p1", reason = "firth_error", N = .N)],
    R[is.na(skip) & family %in% c("BEINF", "BEZI") & !is.finite(nu_firth_logp), .(type = "p0", reason = "p0_nonfinite", N = .N)],
    R[is.na(skip) & family %in% c("BEINF", "BEOI") & !is.finite(tau_firth_logp), .(type = "p1", reason = "p1_nonfinite", N = .N)],
    R[is.na(skip) & cat_complete_sep %in% TRUE, .(type = "p0/p1", reason = "cat_complete_sep_tested_by_firth", N = .N)],
    fill = TRUE)
}

to_units <- function(R) {
  R <- R[is.na(skip)]
  base <- c("tissue", "trait", "promoter", "gene", "family", "n", "n0", "n1", "n_mid", "p0_raw", "p1_raw",
            "n_minor", "n_minor_mid", "n_donors", "cat_complete_sep")
  U <- rbind(
    R[, c(mget(base), list(type = "mu", testable = mu_testable %in% TRUE, logp = mu_F_logp, logp_chi = mu_lrt_logp,
                           est = mu_est, se = mu_se, pe = mu_pe, conv = mu_conv, sep = mu_sep,
                           gate_prop = (p0_raw + p1_raw) <= 0.8, pass_est = abs(mu_est) >= 0.001))],
    R[family %in% c("BEINF", "BEZI"), c(mget(base), list(type = "p0", testable = TRUE, logp = nu_firth_logp, logp_chi = nu_lrt_logp,
                           est = nu_firth_est, se = nu_firth_se, pe = p0_pe, conv = nu_firth_conv, sep = nu_sep,
                           gate_prop = p0_raw >= 0.2, pass_est = abs(nu_firth_est) >= 0.01))],
    R[family %in% c("BEINF", "BEOI"), c(mget(base), list(type = "p1", testable = TRUE, logp = tau_firth_logp, logp_chi = tau_lrt_logp,
                           est = tau_firth_est, se = tau_firth_se, pe = p1_pe, conv = tau_firth_conv, sep = tau_sep,
                           gate_prop = p1_raw >= 0.2, pass_est = abs(tau_firth_est) >= 0.01))])
  U[, testable := testable & is.finite(logp)]
  U[, p := fifelse(testable, exp(logp), NA_real_)]
  U <- merge(U, MIR[, .(tissue, promoter, n_prom_gene, mirror_pair_id, boundary_exact, interior_exact,
                         mirror_role, mirror_partner, drop_mu, drop_boundary_v2, drop_p0, drop_p1)], by = c("tissue", "promoter"), all.x = TRUE)
  for (v in c("boundary_exact", "interior_exact", "drop_mu", "drop_boundary_v2", "drop_p0", "drop_p1")) set(U, which(is.na(U[[v]])), v, FALSE)
  ## Mirror rule v3: mu follows interior_exact (drop the partner's mu, as in v2); p0/p1 follow boundary_exact (drop p1 of both members, keep p0)
  U[, dedup_drop := fifelse(type == "mu", drop_mu, fifelse(type == "p1", drop_p1, drop_p0))]
  U[, dedup_drop_v2 := fifelse(type == "mu", drop_mu, drop_boundary_v2)]          # de-duplication flags under the v2 rule (for comparison)
  swap <- c(mu = "mu", p0 = "p1", p1 = "p0")
  U[, duplicate_of := fifelse(dedup_drop, paste(mirror_partner, swap[type], sep = ":"), NA_character_)]
  ## Kept unit that absorbed its twin unit: mu = mu of the canonical; boundary = p0 of the two members (each absorbs the other's p1)
  ## The eq_class suffix is independent of the label: p0(canonical) = canonical_off, p0(partner) == p1(canonical) = canonical_full
  U[, merged := fifelse(type == "mu", mirror_role %in% "canonical" & interior_exact, type == "p0" & boundary_exact)]
  U[, merged_component := fifelse(merged, fifelse(type == "mu", "interior", "boundary"), NA_character_)]
  U[, twin := fifelse(merged, paste(mirror_partner, swap[type], sep = ":"), NA_character_)]
  U[, eq_class := NA_character_]
  U[merged | dedup_drop, eq_class := paste(mirror_pair_id, fifelse(type == "mu", "interior",
      fifelse((type == "p0") == (mirror_role == "canonical"), "canonical_off", "canonical_full")), sep = "|")]
  U[, merged := NULL]
  U[, primary := testable & !dedup_drop]
  U[, q := NA_real_]
  U[primary == TRUE, q := p.adjust(p, "BH"), by = .(tissue, trait, type)]
  U[, sig := primary & q <= 0.05]
  U[, n50 := n >= 50]
  ## ---- final rules (DECISIONS.md, 2026-10-03; confirmed stratified on sim1 and perm1): only add flag columns, never delete rows ----
  U[, share := fifelse(type == "p0", n0 / n, fifelse(type == "p1", n1 / n, NA_real_))]
  U[, share_10_90 := type == "mu" | (share >= 0.10 & share <= 0.90)]   # p0/p1: boundary proportion 0.10-0.90 (author 2026-10-05: unified with QTL)
  U[, estimable_mu := !(type == "mu" & n_mid < 80)]                       # mu: n_mid >= 80
  U[, assessable := !(type != "mu" & cat_complete_sep %in% TRUE)]         # p0/p1: categorical part completely separated by covariates -> not evaluable
  U[, report_tissue := n_donors >= 50]                                     # the 6 tissues with donors < 50 are not reported
  U[, final := primary & estimable_mu & share_10_90 & assessable & report_tissue]
  U[, q_final := NA_real_]
  U[final == TRUE, q_final := p.adjust(p, "BH"), by = .(tissue, trait, type)]
  U[, sig_final := final & q_final <= 0.05]
  U[, trait_class := fifelse(trait %in% CONT_TRAITS, "continuous",
                     fifelse(!is.na(n_minor) & n_minor / n < 0.10, "rare_binary", "common_binary"))]
  setcolorder(U, c("tissue", "trait", "promoter", "gene", "type", "final", "sig_final", "q_final", "primary", "sig", "q", "p", "logp"))
  U
}
count_mat <- function(S, tis) {
  M <- dcast(unique(S[, .(tissue, trait, promoter)]), tissue ~ trait, value.var = "promoter", fun.aggregate = length)
  M <- merge(data.table(tissue = tis), M, by = "tissue", all.x = TRUE)
  for (t in setdiff(TRAITS, names(M))) M[, (t) := 0L]
  for (t in TRAITS) set(M, which(is.na(M[[t]])), t, 0L)
  setcolorder(M, c("tissue", TRAITS)); M
}

## lambda at p-value quantile q: ratio of observed to expected chi2(1) (q = 0.5 / 0.1 / 0.01 / 0.001 correspond to λ_500 / λ_900 / λ_990 / λ_999)
lam <- function(p, q) qchisq(quantile(p, q, names = FALSE), 1, lower.tail = FALSE) / qchisq(q, 1, lower.tail = FALSE)
calib <- function(U) U[, .(N = .N, lambda_500 = lam(p, 0.5), lambda_900 = lam(p, 0.1), lambda_990 = lam(p, 0.01),
                          lambda_999 = if (.N >= 5000) lam(p, 0.001) else NA_real_,
                          ratio_1e3 = mean(p < 1e-3) / 1e-3, ratio_1e4 = mean(p < 1e-4) / 1e-4,
                          lambda_500_chi = lam(exp(logp_chi[is.finite(logp_chi)]), 0.5))]

UU <- list()
for (arm in ARMS) {
  R <- load_arm(arm); if (is.null(R) || !nrow(R)) next
  save_both(R, paste0("full_", arm))                                   # full raw rows, no filtering
  D <- drop_log(R)[, arm := arm]; fwrite(D, file.path(OUT, paste0("drops_", arm, ".tsv")), sep = "\t")
  U <- to_units(R)
  save_both(U, paste0("units_", arm))                                  # full unit table (with all flags)
  fwrite(U[sig == TRUE], file.path(OUT, paste0("sig_all_", arm, ".tsv")), sep = "\t")   # significant set without the final rules (for comparison)
  S <- U[sig_final == TRUE]; fwrite(S, file.path(OUT, paste0("sig_", arm, ".tsv")), sep = "\t")   # final significant set
  tis <- sort(unique(R$tissue))
  for (ty in c("mu", "p0", "p1")) fwrite(count_mat(S[type == ty], tis), file.path(OUT, sprintf("counts_%s_%s.tsv", arm, ty)), sep = "\t")
  fwrite(count_mat(S, tis), file.path(OUT, sprintf("counts_%s_union.tsv", arm)), sep = "\t")
  cat(sprintf("%s: rows %d | units %d (testable %d, mirror duplicates dropped %d, primary %d, final %d) | significant all-primary %d | FINAL %d (mu %d, p0 %d, p1 %d) | tissue x trait x promoter %d\n",
              arm, nrow(R), nrow(U), sum(U$testable), sum(U$testable & U$dedup_drop), sum(U$primary), sum(U$final),
              sum(U$sig), nrow(S), sum(S$type == "mu"), sum(S$type == "p0"), sum(S$type == "p1"), nrow(unique(S[, .(tissue, trait, promoter)]))))
  cat("  dropped / flagged rows:\n"); print(D[, .(N = sum(N)), by = .(type, reason)])
  ## Memory: the full tables are already written to disk; keep only the columns needed for calibration and null-discovery proportions (each arm has about 6.7 million raw rows; both arms in memory together would exceed 24 GB)
  UU[[arm]] <- U[primary == TRUE, .(tissue, trait, type, p, logp_chi, trait_class, cat_complete_sep, gate_prop, n_donors, sig, final, sig_final)]
  rm(R, U, S); invisible(gc())
}

## ---- calibration: primary units of the null arms; by trait class, trait, cat_complete_sep stratum, proportion gate, tissue donor count ----
cal <- rbindlist(lapply(intersect(names(UU), c("sim1", "perm1")), function(arm) {
  U <- UU[[arm]]                                                     # already contains primary units only
  rbind(U[, calib(.SD), by = .(type, group = trait_class)][, by := "class"],
        U[, calib(.SD), by = .(type, group = trait)][, by := "trait"],
        U[type != "mu", calib(.SD), by = .(type, group = paste0("cat_complete_sep=", cat_complete_sep %in% TRUE))][, by := "sep"],
        U[, calib(.SD), by = .(type, group = paste0("gate_prop=", gate_prop))][, by := "gate"],
        ## Stratified by tissue donor count: to judge whether small-sample tissues are reliable (author 2026-10-01: decide on exclusion after everything has run)
        U[, calib(.SD), by = .(type, group = cut(n_donors, c(0, 49, 99, 199, 399, Inf),
                                                 labels = c("n<50", "50-99", "100-199", "200-399", ">=400")))][, by := "tissue_n"],
        U[final == TRUE, calib(.SD), by = .(type, group = trait_class)][, by := "final_class"],
        fill = TRUE)[, arm := arm]
}), fill = TRUE)
if (nrow(cal)) { fwrite(cal, file.path(OUT, "calibration.tsv"), sep = "\t")
  print(cal[by %in% c("class", "sep", "gate", "tissue_n", "final_class")][order(arm, by, type, group)]) }

## ---- null-discovery proportion: discoveries of the null arm / discoveries of real (same tissue x trait) ----
if ("real" %in% names(UU)) {
  nd <- rbindlist(lapply(intersect(names(UU), c("sim1", "perm1")), function(arm) {
    key <- unique(UU[[arm]][, .(tissue, trait)])
    rbindlist(lapply(c("sig", "sig_final"), function(sc) {
      r <- merge(UU$real[get(sc) == TRUE], key, by = c("tissue", "trait"))[, .(real = .N), by = type]
      z <- UU[[arm]][get(sc) == TRUE, .(null = .N), by = type]
      m <- merge(r, z, by = "type", all = TRUE); m[is.na(null), null := 0L]
      m[, `:=`(ratio = null / real, arm = arm, set = fifelse(sc == "sig", "all_primary", "final"))]; m })) }))
  if (nrow(nd)) { fwrite(nd, file.path(OUT, "null_discovery.tsv"), sep = "\t"); print(nd) }
}
