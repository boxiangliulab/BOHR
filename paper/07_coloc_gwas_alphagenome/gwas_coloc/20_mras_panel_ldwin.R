#!/usr/bin/env Rscript
# MRAS x cardiovascular GWAS across all tissues under the LD-window rule (author 2026-10-12; GWAS lead P < TMAIN, default 5e-8; the jobs set 1e-7).
# = scripts/05_mras_cardiovascular_panel.R with the region replaced: every MRAS analysis row (significant signals incl.
# per-signal LOO, plus the non-significant primary estimable units) x 4 cardiovascular GWAS takes its best tested pair
# (max PP.H4 over the GWAS loci whose LD window meets the signal's; step 14), and that pair's coloc is re-run on its union
# with keep_res = TRUE for the 95% credible set and v784. A row x trait without a tested pair is written with status
# "not tested" and NA posteriors. The v784 switch direction per unit does not depend on the region (as in 05).
# Output: results/mras/ (same file names as 05; the +/- 200 kb versions are in results/other/pre_ldwin_20261012/mras/).
source("/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/scripts/coloc_final_lib.R")
D <- LDU; V784 <- "chr3_138348784_C_T_b38"
OUT <- Sys.getenv("MRAS_OUT", file.path(RES_DIR, "mras")); dir.create(OUT, showWarnings = FALSE); TMAIN <- as.numeric(Sys.getenv("TMAIN", "5e-8"))   # previews of other tiers (2026-10-12)
cat("tier: GWAS lead P <", TMAIN, "; output:", OUT, "\n")
CV <- c("Hypertension_UKB_20002_1065_hg38.txt", "CAD_CARDIoGRAM_C4D_Nikpay_2015_hg38.txt",
        "Coronary_Artery_Disease_Sakaue_2021_Nature_Genetics_hg38.txt", "HeartAttack_UKB_6150_1_hg38.txt")
GROUP <- c(Artery_Aorta = "Vascular", Artery_Coronary = "Vascular", Artery_Tibial = "Vascular",
           Heart_Atrial_Appendage = "Heart", Heart_Left_Ventricle = "Heart", Kidney_Cortex = "Kidney",
           Esophagus_Muscularis = "Visceral smooth muscle", Esophagus_Gastroesophageal_Junction = "Visceral smooth muscle",
           Colon_Sigmoid = "Visceral smooth muscle", Bladder = "Visceral smooth muscle",
           Brain_Hypothalamus = "HPA axis", Pituitary = "HPA axis", Adrenal_Gland = "HPA axis")
GORD <- c("Vascular", "Heart", "Kidney", "Visceral smooth muscle", "HPA axis", "Other brain", "Other")
signal_snps_region <- function(r, marg, loo, start, end) {          # = 14_coloc_tissue.R (library signal_snps over a region)
  if (r$stat == "loo") {
    x <- loo[promoter == r$promoter & type == r$type & k == r$k, .(snp, logp)]
    x <- merge(x, marg[promoter == r$promoter, .(snp, snp_pos, maf)], by = "snp")
  } else x <- marg[promoter == r$promoter, .(snp, snp_pos, maf, logp = get(PCOL[[r$type]]))]
  x <- x[snp_pos >= start & snp_pos <= end & is.finite(logp) & !is.na(maf) & maf > 0 & maf < 1]
  x[, p := exp(logp)]; x <- x[p > 0]; x[maf > 0.5, maf := 1 - maf]
  x[, gwas_id := paste0(sub("^chr", "", r$chr), ":", snp_pos)]
  x[order(p)][!duplicated(gwas_id), .(snp, gwas_id, snp_pos, maf, p)]
}
A <- readRDS(file.path(RES_DIR, "analysis_set.rds"))[promoter %in% MRAS]
P <- readRDS(file.path(D, "results/coloc_ldwin_pairs.rds"))[promoter %in% MRAS & trait %in% CV]
cat("MRAS analysis rows:", nrow(A), " units:", uniqueN(A[, .(tissue, promoter, type)]), " tissues:", uniqueN(A$tissue), "; candidate pairs x CV traits:", nrow(P), "\n")
best <- P[lead_p < TMAIN & status == "tested"][order(-PP.H4, lead_p)]
best <- best[!duplicated(best[, .(tissue, promoter, type, k, trait)])]
st <- P[lead_p < TMAIN, .(n_pairs_candidate = .N, n_pairs_tested = sum(status == "tested")), by = .(tissue, promoter, type, k, trait)]
X <- list(); dirs <- list()
for (t in unique(A$tissue)) {
  At <- A[tissue == t]
  loo <- if (any(At$stat == "loo")) { l <- readRDS(file.path(V2, COND_ROOT, "loo", paste0(t, ".rds"))); setDT(l)
           l[promoter %in% MRAS, .(snp, logp, promoter, type, k)] } else NULL
  marg <- load_marginal(t, "chr3", MRAS)
  for (i in seq_len(nrow(At))) X[[length(X) + 1]] <- list(r = At[i], marg = marg[promoter == At$promoter[i]], loo = loo)
  v <- marg[snp == V784]
  if (nrow(v)) dirs[[t]] <- merge(unique(At[, .(tissue, promoter, type, q, sig)]),
                                  v[, .(promoter, mu_est, nu_firth_est, tau_firth_est, mu_F_logp, nu_firth_logp, tau_firth_logp)], by = "promoter")
}
Dd <- rbindlist(dirs)
Dd[, `:=`(v784_p = exp(fifelse(type == "mu", mu_F_logp, fifelse(type == "p0", nu_firth_logp, tau_firth_logp))),
          T_usage_dir = fifelse(type == "mu", sign(mu_est), fifelse(type == "p0", -sign(nu_firth_est), sign(tau_firth_est))))]
Dd[, expected_dir := fcase(promoter == "PromoterID_45112", 1, promoter == "PromoterID_45113", -1, default = NA_real_)]
Dd[, switch_consistent := T_usage_dir == expected_dir]
fwrite(Dd, file.path(OUT, "mras_v784_switch_direction_by_unit.tsv"), sep = "\t")

rows <- list(); gdir <- list()
for (tr in CV) {
  gw <- readRDS(file.path(D, "gwas_subset", sub("\\.txt$", ".rds", tr))); setkey(gw, gwas_id); Ng <- gwas_N(tr)
  g0 <- readRDS(file.path(NDIR, "gwas_subset", sub("\\.txt$", ".rds", tr)))          # +/- 200 kb subset: holds v784 for every trait
  g7 <- g0[gwas_id == "3:138348784"]
  gdir[[tr]] <- data.table(trait = tr, v784_gwas_beta = if (nrow(g7)) g7$BETA else NA_real_, v784_gwas_p = if (nrow(g7)) g7$Pvalue else NA_real_)
  for (z in X) { r <- z$r
    b <- best[tissue == r$tissue & promoter == r$promoter & type == r$type & k == r$k & trait == tr]
    s0 <- st[tissue == r$tissue & promoter == r$promoter & type == r$type & k == r$k & trait == tr]
    if (!nrow(b)) {
      rows[[length(rows) + 1]] <- data.table(tissue = r$tissue, promoter = r$promoter, type = r$type, k = r$k, n_signals = r$n_signals,
        stat = r$stat, rule = r$rule, signal_snp = r$signal_snp, gene = r$gene, unit_q = r$q, sig_unit = r$sig, trait = tr,
        status = if (nrow(s0)) "not tested: windows share < 5 SNPs" else sprintf("not tested: no GWAS locus (P < %s) in the LD window", sub("e-0", "e-", format(TMAIN))),   # label was hard-coded "5e-8" (wrong at TMAIN 1e-7; fixed 2026-11-04)
        n_pairs_candidate = if (nrow(s0)) s0$n_pairs_candidate else 0L, n_pairs_tested = 0L)
      next }
    x <- signal_snps_region(r, z$marg, z$loo, b$uni_start, b$uni_end)
    g <- gw[J(unique(x$gwas_id)), nomatch = 0L]
    o <- coloc_signal(r, x, g, Ng, keep_res = TRUE); stopifnot(!is.null(o))
    stopifnot(abs(o$row$PP.H4 - b$PP.H4) < 1e-8)                         # the re-run reproduces step 14
    cs <- cred_set(o$res, o$q)
    o$row[, `:=`(trait = tr, status = "tested", cs95_size = nrow(cs), v784_in_cs95 = V784 %in% cs$snp, gwas_locus_lead = b$gwas_locus_lead,
                 gwas_locus_p = b$lead_p, uni_start = b$uni_start, uni_end = b$uni_end, n_int = b$n_int,
                 n_pairs_candidate = s0$n_pairs_candidate, n_pairs_tested = s0$n_pairs_tested)]
    rows[[length(rows) + 1]] <- o$row
  }
  cat(tr, "done\n")
}
res <- rbindlist(rows, fill = TRUE)
res[, group := fifelse(tissue %in% names(GROUP), GROUP[tissue], fifelse(grepl("^Brain_", tissue), "Other brain", "Other"))]
fwrite(res, file.path(OUT, "mras_cardiovascular_all_units_all_tissues.tsv"), sep = "\t")
fwrite(rbindlist(gdir), file.path(OUT, "mras_v784_gwas_direction.tsv"), sep = "\t")
h <- res[trait == CV[1] & status == "tested"]
pt <- h[, .(group = group[1], n_units = uniqueN(paste(promoter, type)), n_sig_units = uniqueN(paste(promoter, type)[sig_unit]),
            n_HC = sum(PP.H4 >= 0.8), best_PPH4 = max(PP.H4),
            best_unit = paste(sub("PromoterID_", "P", promoter), type, paste0("k", k))[which.max(PP.H4)],
            best_sig_PPH4 = if (any(sig_unit)) max(PP.H4[sig_unit]) else NA_real_,
            HC_units = paste(unique(paste(sub("PromoterID_", "P", promoter[PP.H4 >= 0.8]), type[PP.H4 >= 0.8])), collapse = ";"),
            switch_pair_HC = any(promoter == "PromoterID_45112" & PP.H4 >= 0.8) & any(promoter == "PromoterID_45113" & PP.H4 >= 0.8)),
        by = tissue][order(factor(group, GORD), -best_PPH4)]
fwrite(pt, file.path(OUT, "mras_hypertension_by_tissue_grouped.tsv"), sep = "\t")
gs <- pt[, .(tissues = .N, tissues_with_sig_unit = sum(n_sig_units > 0), tissues_with_HC = sum(n_HC > 0),
             median_best_PPH4 = round(median(best_PPH4), 3), tissues_switch_pair_HC = sum(switch_pair_HC)), by = group][order(factor(group, GORD))]
fwrite(gs, file.path(OUT, "mras_hypertension_by_group.tsv"), sep = "\t")
options(width = 220); print(pt); print(gs); print(rbindlist(gdir))
cat("\nrows by trait and status:\n"); print(res[, .N, by = .(trait = sub("_hg38.txt", "", trait), status)])
cat("\nTissues with >= 1 significant MRAS unit (hypertension, tested):", pt[n_sig_units > 0, .N], " best significant PP.H4 >= 0.8:", pt[best_sig_PPH4 >= 0.8, .N],
    " min best significant:", round(pt[n_sig_units > 0, min(best_sig_PPH4)], 3), pt[n_sig_units > 0][which.min(best_sig_PPH4), tissue], "\n")
cat("\nSwitch direction of the T allele (P45112 expected +, P45113 expected -):\n")
print(Dd[!is.na(expected_dir), .(units = .N, consistent = sum(switch_consistent, na.rm = TRUE), sig_units = sum(sig), sig_consistent = sum(switch_consistent & sig, na.rm = TRUE)), by = promoter])
hh <- h[PP.H4 >= 0.8]; cat("\nHypertension HC rows:", nrow(hh), " v784 in cs95:", sum(hh$v784_in_cs95), "\n"); print(hh[, .N, by = gwas_lead]); print(hh[, .N, by = top_snp_H4])
## comparison with the +/- 200 kb panel
old <- fread(file.path(NDIR, "results/other/pre_ldwin_20261012/mras/mras_cardiovascular_all_units_all_tissues.tsv"))
cmpd <- merge(old[, .(tissue, promoter, type, k, trait, old_PPH4 = PP.H4, old_v784_cs = v784_in_cs95)],
              res[, .(tissue, promoter, type, k, trait, status, new_PPH4 = PP.H4, new_v784_cs = v784_in_cs95)], by = c("tissue", "promoter", "type", "k", "trait"), all = TRUE)
cat("\nold vs new, by trait: rows, HC old, HC new, not tested now, |dPP.H4| median among tested\n")
print(cmpd[, .(rows = .N, HC_old = sum(old_PPH4 >= 0.8, na.rm = TRUE), HC_new = sum(new_PPH4 >= 0.8, na.rm = TRUE), not_tested = sum(status != "tested"),
               med_abs_d = round(median(abs(new_PPH4 - old_PPH4), na.rm = TRUE), 4)), by = .(trait = sub("_hg38.txt", "", trait))])
