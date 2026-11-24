#!/usr/bin/env Rscript
# r5e80 MRAS candidate set (S9a rule): variants with r2 >= 0.8 (over each tissue's donors, ld_v3) to the lead of any
# significant r5e80 MRAS puQTL signal; compared with the 161 candidates of the earlier version and the 170 AlphaGenome-scored variants.
# Read-only; writes module checks/ only. No AlphaGenome call.
suppressPackageStartupMessages(library(data.table)); options(width = 250)
G  <- "/projects/group/bliu/e1327989/GTEx"; NF <- file.path(G, "nofloor_full/GTEx")
ND <- file.path(NF, "new_GWAS/coloc_final_nogate_v3_r5e80"); QR <- file.path(NF, "QTL_figure/rerun_r5e80/results")
V2 <- "/home/e1327989/project/GTEx/QTL_figure/rerun_v2/results"; PAPER <- "/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/results_v3c"
FZ <- file.path(G, "manuscript/BOHR_v5/FROZEN_coauthor_review_20261117/tables/supplementary")
CK <- "/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/r5e80_20261119/supp_S7_S9/checks"
RS <- "/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/delivery_2026-11-05_v3m_C22_C29/text/rsid_map_dbSNP156.tsv"
MRAS <- paste0("PromoterID_", 45111:45114); CLASS <- c(mu = "μ", p0 = "p0", p1 = "p1")
lab <- function(d) d[, sprintf("%s:%s:%s:k%d", tissue, sub("PromoterID_", "P", promoter), CLASS[type], k)]

sigs <- function(asf) { A <- readRDS(asf); setDT(A); A[promoter %in% MRAS & sig == TRUE, .(tissue, promoter, type, k, signal_snp)] }
R <- sigs(file.path(ND, "results_v3c/analysis_set.rds")); P <- sigs(file.path(PAPER, "analysis_set.rds"))
cat(sprintf("signals: r5e80 %d, paper %d\n", nrow(R), nrow(P)))
px_r <- fread(cmd = paste("zcat", file.path(QR, "ld_v3/signal_proxies.tsv.gz")))[promoter %in% MRAS & r2 >= 0.8]
px_p <- fread(cmd = paste("zcat", file.path(V2, "ld_v3/signal_proxies.tsv.gz")))[promoter %in% MRAS & r2 >= 0.8]
K <- c("tissue", "promoter", "type", "k", "signal_snp")
pr <- merge(px_r, R, by = K); pp <- merge(px_p, P, by = K)
stopifnot(uniqueN(pr[, ..K]) == nrow(R), uniqueN(pp[, ..K]) == nrow(P))

## same signal (tissue, promoter, class, k, lead) in both versions: are the proxy sets the same? (LD is genotype-only)
same <- merge(R, P, by = K)
cat("signals identical in both versions (same tissue, promoter, class, k and lead):", nrow(same), "\n")
cmp <- merge(pr[same, on = K, .(n_r = .N), by = K], pp[same, on = K, .(n_p = .N), by = K], by = K)
dif <- merge(pr[same, on = K][, .(tissue, promoter, type, k, signal_snp, proxy_snp, r2_r = r2)],
             pp[same, on = K][, .(tissue, promoter, type, k, signal_snp, proxy_snp, r2_p = r2)], by = c(K, "proxy_snp"), all = TRUE)
cat("  proxy rows r5e80 / paper / in both / max |r2 diff|:", sum(!is.na(dif$r2_r)), "/", sum(!is.na(dif$r2_p)), "/",
    dif[!is.na(r2_r) & !is.na(r2_p), .N], "/", signif(dif[!is.na(r2_r) & !is.na(r2_p), max(abs(r2_r - r2_p))], 3), "\n")
## same lead in the same tissue, any promoter: same proxies?
ld_r <- unique(px_r[, .(tissue, signal_snp, proxy_snp, r2)]); ld_p <- unique(px_p[, .(tissue, signal_snp, proxy_snp, r2)])
tl <- merge(unique(ld_r[, .(tissue, signal_snp)]), unique(ld_p[, .(tissue, signal_snp)]))
z <- merge(ld_r[tl, on = .(tissue, signal_snp)], ld_p[tl, on = .(tissue, signal_snp)], by = c("tissue", "signal_snp", "proxy_snp"), all = TRUE)
cat("  tissue x lead pairs in both:", nrow(tl), "; proxy rows only r5e80", z[is.na(r2.y), .N], ", only paper", z[is.na(r2.x), .N], "\n")

## candidate sets
cs <- function(p) p[, .(n_signals_proxied = .N, max_r2_to_signal_lead = max(r2), n_signals_led = sum(proxy_snp == signal_snp),
                        signals_proxied = paste(sort(lab(.SD)), collapse = "; ")), by = .(variant = proxy_snp), .SDcols = c("tissue", "promoter", "type", "k")]
cr <- cs(pr); cp <- cs(pp)
P9a <- fread(file.path(FZ, "Supplementary_Table_S9a_mras_candidate_variants.tsv"))
P9d <- fread(file.path(FZ, "Supplementary_Table_S9d_mras_alphagenome_peak_scores.tsv"))
stopifnot(setequal(cp$variant, P9a$variant_id))
cat("paper candidate set rebuilt from rerun_v2 ld_v3 = frozen S9a:", setequal(cp$variant, P9a$variant_id), "(", nrow(cp), ")\n")
scored <- unique(P9d$variant_id); rsm <- fread(RS)
all <- merge(cr[, .(variant, r5_n = n_signals_proxied, r5_led = n_signals_led, r5_maxr2 = max_r2_to_signal_lead, r5_signals = signals_proxied)],
             P9a[, .(variant = variant_id, paper_n = n_signals_proxied, paper_led = n_signals_led, paper_rank = rank_of_161, max_abs_log2FC)], by = "variant", all = TRUE)
all[, status := fifelse(!is.na(r5_n) & !is.na(paper_n), "both", fifelse(is.na(paper_n), "new in r5e80", "dropped in r5e80"))]
all[, scored_by_AlphaGenome := variant %in% scored]
all <- merge(all, rsm[, .(variant = variant_id, rsID)], by = "variant", all.x = TRUE)
sp <- tstrsplit(all$variant, "_", fixed = TRUE); all[, pos := as.integer(sp[[2]])]
setorder(all, status, paper_rank, pos)
fwrite(all, file.path(CK, "candidates_r5e80_vs_paper.tsv"), sep = "\t", na = "NA", quote = FALSE)
print(all[, .N, by = .(status, scored_by_AlphaGenome)])
cat("\nnew variants:\n"); print(all[status == "new in r5e80", .(variant, rsID, pos, r5_n, r5_led, r5_maxr2, r5_signals)])
nv <- all[status == "new in r5e80", variant]
for (v in nv) {
  s <- pr[proxy_snp == v, ..K]; cat("signals with", v, ":\n"); print(s)
  for (i in seq_len(nrow(s))) { x <- pr[s[i], on = K]; cat("  proxies of this signal:", nrow(x), "; in the r5e80 set elsewhere:", sum(x$proxy_snp %in% cr[n_signals_proxied > 1 | variant != v, variant]), "\n") }
  w <- merge(P, s[, .(tissue, promoter, type)], by = c("tissue", "promoter", "type"))
  cat("  paper signal(s) of the same unit:\n"); print(w)
}
cat("\nscored variants that are r5e80 candidates:", sum(all$status != "dropped in r5e80" & all$scored_by_AlphaGenome), "of", nrow(cr), "\n")
cat("r5e80 leads among candidates:", cr[n_signals_led > 0, .N], "; paper:", cp[n_signals_led > 0, .N], "\n")
cat("rank among the r5e80 candidates (paper |log2FC|, scored only): rs11713141 =",
    all[status == "both"][order(-max_abs_log2FC)][, which(variant == "chr3_138348784_C_T_b38")], "of", all[status == "both", .N], "\n")
cat("dropped variants among the paper's top 20 by |log2FC|:", paste(all[status == "dropped in r5e80" & paper_rank <= 20, sprintf("%s (rank %d)", rsID, paper_rank)], collapse = ", "), "\n")

## signals: r5e80 vs the earlier version ("paper"), by unit
cat("\nsignals by promoter x class, r5e80 vs paper:\n")
print(merge(R[, .(r5e80 = .N), by = .(promoter, type)], P[, .(paper = .N), by = .(promoter, type)], by = c("promoter", "type"), all = TRUE))
cat("r5e80 signals whose lead is a paper lead (any signal):", R[signal_snp %in% P$signal_snp, .N], "of", nrow(R), "\n")
fwrite(R[, .(tissue, promoter, class = CLASS[type], k, signal_snp, lead_in_paper_set = signal_snp %in% P$signal_snp,
             same_signal_in_paper = paste(tissue, promoter, type, k, signal_snp) %in% P[, paste(tissue, promoter, type, k, signal_snp)])],
       file.path(CK, "mras_signals_r5e80.tsv"), sep = "\t", quote = FALSE)
cat("done\n")
