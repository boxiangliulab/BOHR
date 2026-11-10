## Real-data parity, step C: compare a bohr_cis() run with the stored results of the original pipeline.
options(width = 220)
## Prints aggregate numbers only (no sample or donor IDs). Usage: Rscript compare_stored.R <run.rds> <tissue> <chr> [out.rds]
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); RUN <- a[1]; TIS <- a[2]; CHR <- a[3]; OUTS <- if (length(a) > 3) a[4] else ""
V2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"
f <- function(root) file.path(V2, "results", root, "freeze", TIS, "real", paste0(CHR, ".rds"))
M <- readRDS(f("final")); Fm <- readRDS(f("final_muF")); Ad <- readRDS(f("final_firth_addon"))
R <- readRDS(RUN)
S <- as.data.table(M$snp)
## p0/p1 Firth: main run (old-gate promoters) or add-on (all other promoters with a 0/1 part) = discovery_nogate.R:28-31
FC <- c("nu_firth_stat", "nu_firth_logp", "nu_firth_est", "nu_firth_se", "nu_firth_conv", "tau_firth_stat", "tau_firth_logp", "tau_firth_est", "tau_firth_se", "tau_firth_conv", "firth_rounds")
ad <- as.data.table(Ad$snp)[, c("promoter", "snp", FC), with = FALSE]; setnames(ad, FC, paste0("a_", FC))
S <- merge(S, ad, by = c("promoter", "snp"), all.x = TRUE)
for (k in FC) { i <- which(is.na(S[[k]]) & !is.na(S[[paste0("a_", k)]])); if (length(i)) set(S, i, k, S[[paste0("a_", k)]][i]) }
S[, paste0("a_", FC) := NULL]
S <- merge(S, as.data.table(Fm$snp)[, .(promoter, snp, mu_F_logp)], by = c("promoter", "snp"), all.x = TRUE)
P <- as.data.table(R$snp); setnames(P, "variant", "snp")
if (length(unique(P$promoter)) < length(unique(S$promoter))) S <- S[promoter %in% c(P$promoter, R$promoters$promoter)]
mk <- merge(S[, .(promoter, snp)], P[, .(promoter, snp)], by = c("promoter", "snp"))
cat(sprintf("== %s %s: stored SNP rows %d, package SNP rows %d, matched %d, only stored %d, only package %d\n", TIS, CHR, nrow(S), nrow(P),
            nrow(mk), nrow(S) - nrow(mk), nrow(P) - nrow(mk)))
K <- merge(S, P, by = c("promoter", "snp"), suffixes = c(".st", ".pk"))
cmp <- function(a, b, name) {
  a <- as.numeric(a); b <- as.numeric(b); fa <- is.finite(a); fb <- is.finite(b); both <- fa & fb
  d <- abs(a[both] - b[both]); rel <- d / pmax(abs(a[both]), 1e-300)
  data.table(column = name, both_finite = sum(both), na_mismatch = sum(fa != fb), exact = sum(d == 0),
             max_abs_diff = if (any(both)) max(d) else NA_real_, max_rel_diff = if (any(both)) max(rel[abs(a[both]) > 1e-8], 0) else NA_real_,
             agree_1e8 = if (any(both)) mean(d <= 1e-8 * pmax(1, abs(a[both]))) else NA_real_,
             agree_1e6 = if (any(both)) mean(d <= 1e-6 * pmax(1, abs(a[both]))) else NA_real_)
}
cols <- c("mu_F_logp", "nu_firth_logp", "tau_firth_logp", "mu_lrt_logp", "nu_lrt_logp", "tau_lrt_logp", "mu_stat", "nu_stat", "tau_stat",
          "nu_firth_stat", "tau_firth_stat", "mu_est", "mu_se", "nu_firth_est", "tau_firth_est", "nu_firth_se", "tau_firth_se", "nu_est", "tau_est",
          "mu_pe", "p0_pe", "p1_pe", "joint_logp", "mu_wald_logp", "ll_beta_full", "ll_cat_full", "alt_af", "maf", "min_geno")
TAB <- rbindlist(lapply(cols, function(k) { st <- if (paste0(k, ".st") %in% names(K)) K[[paste0(k, ".st")]] else K[[k]]
  pk <- if (paste0(k, ".pk") %in% names(K)) K[[paste0(k, ".pk")]] else K[[k]]; cmp(st, pk, k) }))
## logical flags
for (k in c("mu_conv", "nu_conv", "tau_conv", "mu_sep", "nu_sep", "tau_sep", "nu_firth_conv", "tau_firth_conv", "cat_complete_sep")) {
  st <- K[[paste0(k, ".st")]]; pk <- K[[paste0(k, ".pk")]]
  TAB <- rbind(TAB, data.table(column = k, both_finite = sum(!is.na(st) & !is.na(pk)), na_mismatch = sum(is.na(st) != is.na(pk)),
                               exact = sum(st == pk, na.rm = TRUE), max_abs_diff = NA_real_, max_rel_diff = NA_real_,
                               agree_1e8 = mean(st == pk, na.rm = TRUE), agree_1e6 = mean(st == pk, na.rm = TRUE)), fill = TRUE) }
print(TAB, digits = 3)
## promoter table
PS <- as.data.table(M$prom); PP <- as.data.table(R$promoters)
PM <- merge(PS[, .(promoter, family, skip, n, n0, n1, n_cov, n_snp_window, n_snp_tested, null_ll_beta, null_ll_cat)],
            PP[, .(promoter, family_pk = family, skip_pk = skip, n_pk = n, n0_pk = n0, n1_pk = n1, n_cov_pk = n_cov, nw_pk = n_variants_window,
                   nt_pk = n_variants_tested, llb_pk = null_ll_beta, llc_pk = null_ll_cat)], by = "promoter")
cat(sprintf("\n== promoters: %d compared; family equal %d; skip equal %d; n/n0/n1 equal %d; n_cov equal %d; window size equal %d; tested equal %d; max|null ll diff| beta %.3g cat %.3g\n",
            nrow(PM), PM[family %in% family_pk | (is.na(family) & is.na(family_pk)), .N], PM[skip %in% skip_pk | (is.na(skip) & is.na(skip_pk)), .N],
            PM[(n == n_pk & n0 == n0_pk & n1 == n1_pk) | (is.na(n) & is.na(n_pk)), .N], PM[(n_cov == n_cov_pk) | (is.na(n_cov) & is.na(n_cov_pk)), .N],
            PM[(n_snp_window == nw_pk) | (is.na(n_snp_window) & is.na(nw_pk)), .N], PM[n_snp_tested == nt_pk, .N],
            PM[, max(abs(null_ll_beta - llb_pk), na.rm = TRUE)], PM[, max(abs(null_ll_cat - llc_pk), na.rm = TRUE)]))
print(table(stored = PS$skip, useNA = "a")); print(table(package = PP$skip, useNA = "a"))
## units
U <- as.data.table(readRDS(file.path(V2, "results/final_nogate_v3m/freeze", paste0(TIS, ".real.units.rds"))))
UP <- as.data.table(R$units); setnames(UP, "component", "type")
Uc <- U[promoter %in% R$promoters$promoter]
UM <- merge(Uc[, .(promoter, type, p, n_snp, leader, leader_logp, leader_est, leader_pe, estimable, dedup_drop, boundary_exact, interior_exact, mirror_role, q)],
            UP[, .(promoter, type, p_pk = p, n_pk = n_variants, leader_pk = leader, llp_pk = leader_log_p, lest_pk = leader_est, lpe_pk = leader_effect,
                   est_pk = estimable, dd_pk = dedup_drop, be_pk = boundary_exact, ie_pk = interior_exact, role_pk = mirror_role)], by = c("promoter", "type"), all = TRUE)
cat(sprintf("\n== units on %s: stored %d, package %d, matched %d\n", CHR, nrow(Uc), nrow(UP), UM[!is.na(p) & !is.na(p_pk), .N]))
UM2 <- UM[!is.na(p) & !is.na(p_pk)]
cat(sprintf("   ACAT p: max|diff| %.3g, max rel diff %.3g, max |log10 p diff| %.3g; n_snp equal %d/%d; leader equal %d/%d; estimable equal %d/%d; dedup_drop equal %d/%d; boundary_exact equal %d/%d; interior_exact equal %d/%d; role equal %d/%d\n",
            UM2[, max(abs(p - p_pk))], UM2[, max(abs(p - p_pk) / p)], UM2[, max(abs(log10(p) - log10(p_pk)))],
            UM2[n_snp == n_pk, .N], nrow(UM2), UM2[leader == leader_pk, .N], nrow(UM2), UM2[estimable == est_pk, .N], nrow(UM2),
            UM2[dedup_drop == dd_pk, .N], nrow(UM2), UM2[boundary_exact == be_pk, .N], nrow(UM2), UM2[interior_exact == ie_pk, .N], nrow(UM2),
            UM2[mirror_role %in% role_pk | (is.na(mirror_role) & is.na(role_pk)), .N], nrow(UM2)))
cat(sprintf("   leader log p: max|diff| %.3g; leader est: max|diff| %.3g; leader effect: max|diff| %.3g\n",
            UM2[, max(abs(leader_logp - llp_pk))], UM2[, max(abs(leader_est - lest_pk), na.rm = TRUE)], UM2[, max(abs(leader_pe - lpe_pk), na.rm = TRUE)]))
## q: BH over the tissue = stored units of the other chromosomes + the package's units of this chromosome
if (nrow(UP) && length(unique(UP$promoter)) == length(unique(R$promoters$promoter[is.na(R$promoters$skip)]))) {
  other <- U[!promoter %in% R$promoters$promoter, .(promoter, type, p, estimable, dedup_drop)]
  allu <- rbind(other, UP[, .(promoter, type, p, estimable, dedup_drop)])
  allu[, q_pk := NA_real_]; allu[dedup_drop == FALSE & estimable == TRUE, q_pk := p.adjust(p, "BH")]
  QQ <- merge(U[, .(promoter, type, q)], allu[, .(promoter, type, q_pk)], by = c("promoter", "type"))
  Qc <- QQ[promoter %in% R$promoters$promoter]
  cat(sprintf("   q (tissue-wide BH with this chromosome recomputed): all units %d, NA pattern equal %d; max|q diff| all %.3g, %s %.3g; significant (q<=0.05) stored %d package %d, agree %d\n",
              nrow(QQ), QQ[is.na(q) == is.na(q_pk), .N], QQ[, max(abs(q - q_pk), na.rm = TRUE)], CHR, Qc[, max(abs(q - q_pk), na.rm = TRUE)],
              QQ[q <= 0.05, .N], QQ[q_pk <= 0.05, .N], QQ[(q <= 0.05) == (q_pk <= 0.05) | (is.na(q) & is.na(q_pk)), .N]))
  cat(sprintf("   %s units significant: stored %d, package %d\n", CHR, Qc[q <= 0.05, .N], Qc[q_pk <= 0.05, .N]))
}
if (nzchar(OUTS)) saveRDS(list(snp = TAB, units = UM, prom = PM), OUTS)
