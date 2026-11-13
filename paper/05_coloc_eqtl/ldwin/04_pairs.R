## Figure 4 panel e, LD-window coloc, step 4: the final signal x eQTL pairs and their regions (author 2026-10-12).
## Every analysis signal of coloc_eqtl_v3 (14,783; analysis_set.rds, read-only) is paired with the eQTL lead of its gene in
## that tissue (the GTEx v11 eGenes-file top variant; eqtl_leads.tsv). Windows (lab rule): r2 >= 0.1 partners (PLINK 1.9,
## tissue donors, +-1 Mb), farthest partner +/- 50 kb; a lead without a partner gets +/- 50 kb.
##   puQTL lead and eQTL leads in the MAF >= 0.05 panel: lead_windows.tsv (step 2);
##   eQTL leads absent from the panel (1,485, all MAF < 0.05): option A' (author 2026-10-12 "OK, use A'") -- the lead is
##   typed in the unfiltered v9 VCF and its partners are restricted to the panel (rarelead_windows.tsv, sAp/eAp; step 3).
## Intersection = [max start, min end] (empty when start > end); union = [min start, max end]. The validity (>= 5 SNPs with
## statistics in both datasets inside the intersection) and the >= 10 shared SNPs in the union are applied in step 6.
## Output: results/coloc_eqtl_ldwin/pairs.tsv.gz and eqtl_windows.tsv (per tissue x gene: the span of its pairs' unions).
suppressPackageStartupMessages(library(data.table))
V2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"
O  <- file.path(V2, "results/coloc_eqtl_ldwin")
A  <- readRDS(file.path(V2, "results/coloc_eqtl_v3/analysis_set.rds"))
E  <- fread(file.path(O, "eqtl_leads.tsv"))
L  <- fread(file.path(O, "lead_windows.tsv"))
R  <- fread(file.path(O, "rarelead_windows.tsv"))
stopifnot(nrow(R) == L[absent == TRUE, .N], R[, all(typed)])
## one window per tissue x lead SNP
W <- rbind(L[absent == FALSE, .(tissue, snp, n_partners, w_start = win_start, w_end = win_end, rare = FALSE)],
           R[, .(tissue, snp, n_partners = nAp, w_start = sAp, w_end = eAp, rare = TRUE)])
stopifnot(!anyDuplicated(W[, .(tissue, snp)]), !anyNA(W$w_start), !anyNA(W$w_end))
P <- merge(A[, .(tissue, promoter, type, k, gene_base, chr, signal_snp, signal_pos, egene)], E[, .(tissue, gene_base, eqtl_top = eqtl_lead)],
           by = c("tissue", "gene_base"), all.x = TRUE)
stopifnot(nrow(P) == nrow(A), !anyNA(P$eqtl_top))
P <- merge(P, W[, .(tissue, signal_snp = snp, p_partners = n_partners, p_start = w_start, p_end = w_end)], by = c("tissue", "signal_snp"), all.x = TRUE)
P <- merge(P, W[, .(tissue, eqtl_top = snp, e_rare = rare, e_partners = n_partners, e_start = w_start, e_end = w_end)], by = c("tissue", "eqtl_top"), all.x = TRUE)
stopifnot(nrow(P) == nrow(A), !anyNA(P$p_start), !anyNA(P$e_start))
P[, `:=`(int_start = pmax(p_start, e_start), int_end = pmin(p_end, e_end), uni_start = pmin(p_start, e_start), uni_end = pmax(p_end, e_end))]
P[, `:=`(eqtl_top_pos = as.integer(sub("^chr[0-9]+_([0-9]+)_.*$", "\\1", eqtl_top)), overlap_bp = pmax(0L, int_end - int_start + 1L))]
P[, lead_dist := abs(eqtl_top_pos - signal_pos)]
setcolorder(P, c("tissue", "promoter", "type", "k", "gene_base", "chr", "signal_snp", "signal_pos", "egene", "eqtl_top", "eqtl_top_pos", "e_rare",
                 "p_partners", "p_start", "p_end", "e_partners", "e_start", "e_end", "int_start", "int_end", "overlap_bp", "uni_start", "uni_end", "lead_dist"))
setorder(P, tissue, promoter, type, k)
fwrite(P, file.path(O, "pairs.tsv.gz"), sep = "\t")
## eQTL extraction windows: per tissue x gene, the span of its pairs' unions (the coloc region of each pair lies inside)
X <- P[, .(start = min(uni_start), end = max(uni_end)), by = .(tissue, gene_base, chr)]
fwrite(X, file.path(O, "eqtl_windows.tsv"), sep = "\t")
cat("pairs:", nrow(P), "; rare eQTL lead (A'):", P[e_rare == TRUE, .N], "; positional overlap:", P[overlap_bp > 0, .N], "\n")
q <- function(x) round(quantile(x / 1e3, c(.1, .5, .9, .99, 1)))
cat("union width (kb):\n"); print(q(P$uni_end - P$uni_start))
cat("eQTL extraction windows:", nrow(X), "; width (kb):\n"); print(q(X$end - X$start))
cat("v3 extraction windows (signal +-200 kb), total Mb:", round(fread(file.path(V2, "results/coloc_eqtl_v3/eqtl_windows.tsv"))[, sum(end - start)] / 1e6),
    "; now:", round(X[, sum(end - start)] / 1e6), "\n")
