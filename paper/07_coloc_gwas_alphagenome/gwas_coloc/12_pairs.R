#!/usr/bin/env Rscript
# LD-window GWAS coloc (author 2026-10-12), step 12: candidate puQTL signal x GWAS locus pairs.
# - puQTL windows: step 12a (Figure 4's shared table, tissue donors) + step 12b (the 122 non-significant MRAS leads).
# - GWAS loci: step 11 (NSCC): recursive-clumping leads at P < 1e-5 (ALL donors); window = [min(lead, partners) - 50 kb,
#   max(lead, partners) + 50 kb], partners = panel variants with r2 >= 0.1 (A'); a lead with no partner gets +/- 50 kb.
# - A pair is a candidate when the two windows overlap by position (intersection [max start, min end] non-empty).
#   Validity (>= 5 SNPs with statistics in both arms in the intersection) and the >= 10 shared SNPs in the union are
#   applied in step 14. The main threshold (lead P) is applied in step 15, so every tier comes from one run.
# Output: pairs/<tissue>.rds, gwas_loci.tsv.gz (all leads + windows), union_intervals/<trait>.rds (for step 13).
source("/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/scripts/coloc_final_lib.R")
D <- LDU; PAD <- 50000L
dir.create(file.path(D, "pairs"), showWarnings = FALSE); dir.create(file.path(D, "union_intervals"), showWarnings = FALSE)
## puQTL windows
W <- fread(file.path(D, "puqtl_windows.tsv"))
X <- fread(file.path(LDU0, "nscc_out/puqtl_extra_windows.tsv"))
if (UNITSET != "v3") for (xf in Sys.glob(file.path(LDU, "extra_windows*.tsv"))) X <- unique(rbind(X, fread(xf)), by = c("tissue", "snp"))   # v3m: windows of leads new to v3m
X[, `:=`(p_start = pmin(pos, min_partner_pos) - PAD, p_end = pmax(pos, max_partner_pos) + PAD)]
W <- merge(W, X[, .(tissue, signal_snp = snp, x_partners = n_partners, x_start = p_start, x_end = p_end)], by = c("tissue", "signal_snp"), all.x = TRUE)
W[is.na(p_start), `:=`(p_partners = x_partners, p_start = x_start, p_end = x_end)]; W[, c("x_partners", "x_start", "x_end") := NULL]
stopifnot(!anyNA(W$p_start), !anyNA(W$p_end), nrow(W) == nrow(readRDS(file.path(RES_DIR, "analysis_set.rds"))))
W[, chrn := sub("^chr", "", chr)]
## GWAS loci
TRB <- sub("\\.txt$", "", TRAITS)
G <- rbindlist(lapply(TRB, function(tb) {
  f <- file.path(LDU0, "nscc_out/clump", tb); stopifnot(file.exists(file.path(f, "DONE")))
  l <- fread(file.path(f, "leads.tsv"), colClasses = list(character = c("chr", "lead", "gwas_id")))
  if (!nrow(l)) return(NULL)
  q <- fread(file.path(f, "lead_partners.tsv"), colClasses = list(character = "lead"))
  stopifnot(!anyNA(q$pos), setequal(l$lead, q$lead))
  l <- merge(l, q, by = "lead")
  l[, pos := as.integer(sub("^chr[0-9]+_([0-9]+)_.*$", "\\1", lead))]
  l[n_partners == 0, `:=`(min_partner_pos = pos, max_partner_pos = pos)]
  l[, `:=`(g_start = pmin(pos, min_partner_pos) - PAD, g_end = pmax(pos, max_partner_pos) + PAD, trait = paste0(tb, ".txt"))]
  l }))
fwrite(G, file.path(D, "gwas_loci.tsv.gz"), sep = "\t", quote = FALSE)
cat("GWAS loci (P < 1e-5):", nrow(G), "; P < 5e-8:", G[lead_p < 5e-8, .N], "; P < 1e-7:", G[lead_p < 1e-7, .N], "; P < 1e-6:", G[lead_p < 1e-6, .N],
    "; A' leads (not in panel):", G[in_panel == FALSE, .N], "\n")
cat("GWAS window width (kb), P < 5e-8:\n"); print(G[lead_p < 5e-8, quantile((g_end - g_start) / 1e3, c(0, .1, .5, .9, .99, 1))])
print(G[, .(loci_5e8 = sum(lead_p < 5e-8), loci_1e7 = sum(lead_p < 1e-7), loci_1e6 = sum(lead_p < 1e-6), loci_1e5 = .N), by = trait][order(-loci_5e8)])
## candidate pairs: positional overlap of the two windows
setkey(G, chr, g_start, g_end)
Wk <- W[, .(tissue, promoter, type, k, sig, chr = chrn, signal_snp, signal_pos, p_start, p_end)]
P <- foverlaps(Wk, G[, .(chr, g_start, g_end, trait, lead, gwas_id, lead_p, lead_absz, rank, in_panel, n_partners)],
               by.x = c("chr", "p_start", "p_end"), by.y = c("chr", "g_start", "g_end"), type = "any", nomatch = 0L)
P[, `:=`(int_start = pmax(p_start, g_start), int_end = pmin(p_end, g_end), uni_start = pmin(p_start, g_start), uni_end = pmax(p_end, g_end))]
stopifnot(P[, all(int_start <= int_end)])
setnames(P, c("lead", "gwas_id", "rank", "in_panel", "n_partners"), c("gwas_locus_lead", "gwas_locus_id", "gwas_locus_rank", "gwas_lead_in_panel", "gwas_lead_partners"))
cat("\ncandidate pairs:", nrow(P), "; lead P < 5e-8:", P[lead_p < 5e-8, .N], "; signals with >= 1 pair at 5e-8:",
    uniqueN(P[lead_p < 5e-8 & sig == TRUE, paste(tissue, promoter, type, k)]), "of", W[sig == TRUE, .N], "\n")
cat("pairs per signal x trait (5e-8):\n"); print(P[lead_p < 5e-8, .N, by = .(tissue, promoter, type, k, trait)][, quantile(N, c(.5, .9, .99, 1))])
for (ti in unique(W$tissue)) saveRDS(P[tissue == ti], file.path(D, "pairs", paste0(ti, ".rds")))
for (tr in TRAITS) { u <- unique(P[trait == tr, .(chr, s = uni_start, e = uni_end)]); saveRDS(u, file.path(D, "union_intervals", sub("\\.txt$", ".rds", tr))) }
cat("done\n")
