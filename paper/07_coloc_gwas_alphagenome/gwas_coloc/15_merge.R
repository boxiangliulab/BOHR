#!/usr/bin/env Rscript
# LD-window GWAS coloc (author 2026-10-12), step 15: merge the per-tissue pair rows (step 14) and summarise per GWAS lead
# threshold tier.
# - Pair level (every candidate pair, every tier): results/coloc_ldwin_pairs.rds.
# - Signal x trait level, per tier T (pairs with GWAS lead P < T): the best tested pair (max PP.H4; ties -> smaller lead P),
#   as Figure 4's eQTL side takes the best call over a signal's pairs (kept consistent with Fig 4). The p12 = 5e-6
#   sensitivity column is the max over the same tested pairs. Columns of the main run's master are kept, so the existing
#   Figure 6 / supplementary scripts can read it, plus n_pairs_candidate, n_pairs_tested and the best pair's locus and
#   windows.
# - Status of a signal x trait without a tested pair (reported, never counted as "no colocalization"):
#   "no GWAS locus" (no window overlap; implicit, counted only), "windows share < 5 SNPs", "< 10 shared SNPs".
# - The significant analysis signals only; the non-significant MRAS units (MRAS panel) go to a separate file.
# Tiers: 5e-8, 1e-7 (main candidates), 1e-6 (sensitivity), 1e-5 (exploratory).
source("/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/scripts/coloc_final_lib.R")
D <- LDU; R <- file.path(D, "results"); dir.create(R, showWarnings = FALSE)
fs <- list.files(file.path(D, "per_tissue"), "^coloc_ldwin_.*\\.rds$", full.names = TRUE)
stopifnot(length(fs) == 50)
P <- rbindlist(lapply(fs, readRDS), fill = TRUE)
saveRDS(P, file.path(R, "coloc_ldwin_pairs.rds"))
A <- readRDS(file.path(RES_DIR, "analysis_set.rds"))
cat("pair rows:", nrow(P), "\n"); print(P[, .N, by = status])
TIERS <- c("5e-8" = 5e-8, "1e-7" = 1e-7, "1e-6" = 1e-6, "1e-5" = 1e-5)
best_tier <- function(t) {
  X <- P[sig_unit == TRUE & lead_p < t]
  st <- X[, .(n_pairs_candidate = .N, n_pairs_tested = sum(status == "tested"), n_pairs_disjoint = sum(status == "different locus (LD windows disjoint)"),
              n_pairs_lt10 = sum(grepl("^not assessable", status))), by = .(tissue, promoter, type, k, trait)]
  B <- X[status == "tested"][order(-PP.H4, lead_p)][, `:=`(PP.H4_p12_5e6_max = max(PP.H4_p12_5e6)), by = .(tissue, promoter, type, k, trait)]
  B <- B[!duplicated(B[, .(tissue, promoter, type, k, trait)])]
  B <- merge(B, st, by = c("tissue", "promoter", "type", "k", "trait"))
  B[, `:=`(PP.H4_p12_5e6_bestpair = PP.H4_p12_5e6, PP.H4_p12_5e6 = PP.H4_p12_5e6_max)][, PP.H4_p12_5e6_max := NULL]
  B[, `:=`(HC = PP.H4 >= 0.8, HC_5e6 = PP.H4_p12_5e6 >= 0.8, tier = names(TIERS)[TIERS == t])]
  nt <- st[n_pairs_tested == 0]
  nt[, status := fifelse(n_pairs_disjoint > 0, "not tested: windows share < 5 SNPs", "not tested: < 10 shared SNPs")]
  list(best = B, untested = nt)
}
S <- lapply(TIERS, best_tier)
nsig <- A[sig == TRUE, .N]
tot <- rbindlist(lapply(names(TIERS), function(tn) { B <- S[[tn]]$best; U <- S[[tn]]$untested
  data.table(tier = tn,
    signal_x_trait_total = nsig * length(TRAITS),
    with_candidate_pair = nrow(B) + nrow(U),
    tested = nrow(B), not_tested_disjoint = U[status == "not tested: windows share < 5 SNPs", .N], not_tested_lt10 = U[status == "not tested: < 10 shared SNPs", .N],
    no_gwas_locus = nsig * length(TRAITS) - nrow(B) - nrow(U),
    traits_tested = uniqueN(B$trait), HC_tests = sum(B$HC), HC_mu = B[HC == TRUE & type == "mu", .N], HC_boundary = B[HC == TRUE & type != "mu", .N],
    HC_p12_5e6 = sum(B$HC_5e6), HC_promoter_tissue_trait = uniqueN(B[HC == TRUE, .(tissue, promoter, trait)]),
    HC_genes = uniqueN(sub("\\..*", "", B[HC == TRUE, gene])), HC_files = uniqueN(B[HC == TRUE, trait]), HC_tissues = uniqueN(B[HC == TRUE, tissue]),
    pairs_tested = P[sig_unit == TRUE & lead_p < TIERS[[tn]] & status == "tested", .N]) }))
fwrite(tot, file.path(R, "summary_totals_by_tier.tsv"), sep = "\t")
options(width = 250); print(tot)
for (tn in names(TIERS)) {
  saveRDS(S[[tn]]$best, file.path(R, sprintf("coloc_ldwin_signal_trait_%s.rds", tn)))
  fwrite(S[[tn]]$untested, file.path(R, sprintf("not_tested_signal_trait_%s.tsv.gz", tn)), sep = "\t")
}
## MRAS non-significant units (MRAS panel), best pair per tier
M <- P[sig_unit == FALSE]
saveRDS(M, file.path(R, "coloc_ldwin_pairs_mras_nonsig.rds"))
## comparison with the +/- 200 kb master at promoter x tissue x trait (max PP.H4 over signals and classes)
old <- readRDS(file.path(NDIR, "results/coloc_final_merged.rds")); setDT(old)
o <- old[, .(old_PPH4 = max(PP.H4)), by = .(tissue, promoter, trait)]
cmp <- rbindlist(lapply(names(TIERS), function(tn) { f <- S[[tn]]$best[, .(new_PPH4 = max(PP.H4)), by = .(tissue, promoter, trait)]
  j <- merge(o, f, by = c("tissue", "promoter", "trait"), all = TRUE)
  data.table(tier = tn, old_HC = j[old_PPH4 >= .8, .N], new_HC = j[new_PPH4 >= .8, .N], both_HC = j[old_PPH4 >= .8 & new_PPH4 >= .8, .N],
             old_HC_not_tested_new = j[old_PPH4 >= .8 & is.na(new_PPH4), .N], old_HC_tested_below_new = j[old_PPH4 >= .8 & new_PPH4 < .8, .N],
             new_HC_not_old = j[new_PPH4 >= .8 & (is.na(old_PPH4) | old_PPH4 < .8), .N],
             pearson_shared = round(j[!is.na(old_PPH4) & !is.na(new_PPH4), cor(old_PPH4, new_PPH4)], 3)) }))
fwrite(cmp, file.path(R, "compare_vs_200kb_by_tier.tsv"), sep = "\t"); print(cmp)
## MRAS
for (tn in names(TIERS)) { B <- S[[tn]]$best[promoter %in% MRAS & HC == TRUE]
  cat(sprintf("\nMRAS HC tests, tier %s: %d (tissues %d; files: %s)\n", tn, nrow(B), uniqueN(B$tissue), paste(sort(unique(sub("_hg38.txt", "", B$trait))), collapse = ", "))) }
## gene ranking (old panel-b metric: distinct HC disease colocalizations, tissue x promoter x trait)
source(file.path(NDIR, "figure6/scripts/trait_maps.R"))
DISCATS <- c("Cardiovascular","Lipids","Glycemic & metabolic","Immune & inflammatory","Neuropsychiatric","Endocrine")
sym <- fread("/home/e1327989/project/GTEx/QTL_figure/figure/data/ensg_symbol.tsv"); setnames(sym, c("gid", "symbol"))
rk <- rbindlist(lapply(names(TIERS), function(tn) { d <- S[[tn]]$best[HC == TRUE]; d[, cat := trait_cat[trait_map[trait]]]; d <- d[cat %in% DISCATS]
  d[, gid := sub("\\..*", "", gene)]
  r <- d[, .(n_tissue_promoter_trait = uniqueN(paste(tissue, promoter, trait)), n_tissues = uniqueN(tissue), n_traits = uniqueN(trait_map[trait])), by = gid]
  r <- merge(r, sym, by = "gid", all.x = TRUE)[order(-n_tissue_promoter_trait)]
  r[, `:=`(rank = frank(-n_tissue_promoter_trait, ties.method = "min"), tier = tn)]; r }))
fwrite(rk, file.path(R, "gene_rank_disease_HC_by_tier.tsv"), sep = "\t")
print(rk[symbol == "MRAS"]); print(rk[, .(n_genes = .N), by = tier])
for (tn in names(TIERS)) { cat("\ntop 10 genes, tier", tn, "\n"); print(head(rk[tier == tn, .(rank, symbol, n_tissue_promoter_trait, n_tissues, n_traits)], 10)) }
