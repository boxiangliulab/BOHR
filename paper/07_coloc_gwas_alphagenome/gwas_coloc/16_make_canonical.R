#!/usr/bin/env Rscript
# LD-window GWAS coloc, step 16 (author 2026-10-12: use 5e-8, rebuild everything): the canonical master and the summary
# files of scripts/04_merge_summarise.R, from the signal x trait table of step 15 at the tier in TIER (the jobs set 1e-7).
# - Master: one row per TESTED signal x trait (significant units; the best
#   pair over the signal's GWAS loci with lead P < TIER), with every column of the previous master, so the Figure 6 and
#   supplementary scripts read it unchanged. "tests" now means tested signal x trait units (a GWAS locus within the LD
#   window); untested units are counted in summary_totals.tsv and listed in not_tested_signal_trait_5e-8.tsv.gz.
# - Summaries with the names 04 wrote (summary_totals/by_type/by_tissue/by_trait/by_statistic, coloc_final_HC_rows,
#   gene_rank_disease_HC_final_1e5/_5e6, MRAS_genome_wide_rank) go to results/canonical_<TAG>/ (TAG from TIER).
source("/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/scripts/coloc_final_lib.R")
# Tier (GWAS lead P threshold) from the environment: 5e-8 first (2026-10-12), then 1e-7 (author: lower the threshold to 1e-7).
TIER <- Sys.getenv("TIER", "5e-8"); TAG <- c("5e-8" = "5e8", "1e-7" = "1e7", "1e-6" = "1e6", "1e-5" = "1e5")[[TIER]]
D <- LDU; R <- file.path(D, "results"); O <- file.path(R, paste0("canonical_", TAG)); dir.create(O, showWarnings = FALSE)
m <- readRDS(file.path(R, sprintf("coloc_ldwin_signal_trait_%s.rds", TIER))); setDT(m)
cat("tier: GWAS lead P <", TIER, "\n")
m[, gwas_min_p := as.numeric(gwas_min_p)]
m[, `:=`(HC = PP.H4 >= 0.8, HC_5e6 = PP.H4_p12_5e6 >= 0.8, class2 = fifelse(type == "mu", "mu", "boundary"))]
setcolorder(m, c("tissue", "promoter", "type", "k", "n_signals", "stat", "rule", "signal_snp", "gene", "unit_q", "sig_unit", "n_snps",
                 "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "PP.H4_p12_5e6", "pqtl_min_p", "gwas_min_p", "top_snp_H4", "top_snp_PP",
                 "gwas_lead", "trait"))
stopifnot(!anyDuplicated(m[, .(tissue, promoter, type, k, trait)]), all(m$sig_unit))
saveRDS(m, if (UNITSET == "v3") file.path(NDIR, sprintf("results/coloc_ldwin_merged_%s_20261012.rds", TAG)) else file.path(RES_DIR, sprintf("coloc_ldwin_merged_%s_%s_%s.rds", TAG, UNITSET, if (UNITSET == "v3c") "20261106" else "20261104")))   # 2026-11-06: v3c dated by its build
fwrite(m[HC | HC_5e6], file.path(O, "coloc_final_HC_rows.tsv"), sep = "\t")
A <- readRDS(file.path(RES_DIR, "analysis_set.rds"))[sig == TRUE]
U <- fread(cmd = paste("zcat", file.path(R, sprintf("not_tested_signal_trait_%s.tsv.gz", TIER))))
u <- m[, .(HC = any(HC), HC_5e6 = any(HC_5e6), best_PPH4 = max(PP.H4)), by = .(tissue, promoter, type, gene, trait)]
tot <- data.table(
  metric = c("tissues", "traits (GWAS files)", "analysis signals (significant units)", "signal x trait units",
             sprintf("signal x trait units with a GWAS locus (lead P < %s) whose LD window overlaps the signal's", TIER),
             "  tested (>= 5 shared SNPs in the window intersection, >= 10 in the union)", "  not tested: windows share < 5 SNPs",
             "signal x trait units without a GWAS locus", "traits with >= 1 tested unit", sprintf("candidate pairs (signal x GWAS locus), lead P < %s", TIER),
             "tested pairs", "HC tests (signal x trait), p12=1e-5 (MAIN)", "  mu", "  boundary (p0/p1)", "HC tests, p12=5e-6 (sensitivity)",
             "unit x trait HC (any signal), 1e-5", "units with >= 1 HC, 1e-5", "unique promoter-tissue-trait HC, 1e-5",
             "HC tests from secondary signals (k >= 2), 1e-5", "HC tests from LOO statistics, 1e-5"),
  value = c(uniqueN(A$tissue), length(TRAITS), nrow(A), nrow(A) * length(TRAITS), nrow(m) + nrow(U), nrow(m), nrow(U),
            nrow(A) * length(TRAITS) - nrow(m) - nrow(U), uniqueN(m$trait), sum(m$n_pairs_candidate) + sum(U$n_pairs_candidate),
            sum(m$n_pairs_tested), sum(m$HC), m[HC == TRUE & type == "mu", .N], m[HC == TRUE & type != "mu", .N], sum(m$HC_5e6), sum(u$HC),
            uniqueN(u[HC == TRUE, .(tissue, promoter, type)]), uniqueN(u[HC == TRUE, .(tissue, promoter, trait)]),
            m[HC == TRUE & k >= 2, .N], m[HC == TRUE & stat == "loo", .N]))
fwrite(tot, file.path(O, "summary_totals.tsv"), sep = "\t")
fwrite(rbind(m[, .(tests = .N, HC = sum(HC), HC_5e6 = sum(HC_5e6)), by = .(class = type)][order(class)],
             m[, .(tests = .N, HC = sum(HC), HC_5e6 = sum(HC_5e6)), by = .(class = class2)][order(-tests)]), file.path(O, "summary_by_type.tsv"), sep = "\t")
fwrite(m[, .(tests = .N, HC = sum(HC), HC_5e6 = sum(HC_5e6), per_1000 = round(1000 * mean(HC), 2)), by = tissue][order(-HC)], file.path(O, "summary_by_tissue.tsv"), sep = "\t")
fwrite(m[, .(tests = .N, HC = sum(HC), HC_5e6 = sum(HC_5e6)), by = trait][order(-HC)], file.path(O, "summary_by_trait.tsv"), sep = "\t")
fwrite(m[, .(signals = uniqueN(paste(tissue, promoter, type, k)), tests = .N, HC = sum(HC)), by = .(stat, rule)], file.path(O, "summary_by_statistic.tsv"), sep = "\t")
## gene ranking: distinct HC DISEASE colocalizations (tissue x promoter x trait), the metric of scripts/04
source(file.path(NDIR, "figure6/scripts/trait_maps.R"))
DISCATS <- c("Cardiovascular","Lipids","Glycemic & metabolic","Immune & inflammatory","Neuropsychiatric","Endocrine")
sym <- fread("/home/e1327989/project/GTEx/QTL_figure/figure/data/ensg_symbol.tsv"); setnames(sym, c("gid", "symbol"))
rank_genes <- function(d, pp, gcol) {
  d <- d[!is.na(get(pp)) & get(pp) >= 0.8]; d[, cat := trait_cat[trait_map[trait]]]; d <- d[cat %in% DISCATS]
  d[, gid := sub("\\..*", "", get(gcol))]
  r <- d[, .(n_tissue_promoter_trait = uniqueN(paste(tissue, promoter, trait)), n_tissues = uniqueN(tissue),
             n_traits = uniqueN(trait_map[trait])), by = gid]
  r <- merge(r, sym, by = "gid", all.x = TRUE)[order(-n_tissue_promoter_trait)]
  r[, rank := frank(-n_tissue_promoter_trait, ties.method = "min")]; r[] }
stopifnot(all(unique(m$trait) %in% names(trait_map)))
rk <- list(final_1e5 = rank_genes(m, "PP.H4", "gene"), final_5e6 = rank_genes(m, "PP.H4_p12_5e6", "gene"))
for (k in names(rk)) fwrite(rk[[k]], file.path(O, sprintf("gene_rank_disease_HC_%s.tsv", k)), sep = "\t")
mr <- rbindlist(lapply(names(rk), function(k) { x <- rk[[k]]; y <- x[gid == "ENSG00000158186"]
  data.table(set = k, n_genes_with_HC = nrow(x), MRAS_rank = if (nrow(y)) y$rank else NA_integer_,
             MRAS_HC = if (nrow(y)) y$n_tissue_promoter_trait else 0L, MRAS_tissues = if (nrow(y)) y$n_tissues else 0L,
             MRAS_traits = if (nrow(y)) y$n_traits else 0L) }))
fwrite(mr, file.path(O, "MRAS_genome_wide_rank.tsv"), sep = "\t")
options(width = 200); print(tot); print(fread(file.path(O, "summary_by_type.tsv"))); print(mr); print(head(rk$final_1e5, 12))
