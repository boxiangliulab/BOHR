## Mirror rule v3 coloc, step 2: inputs that do not exist yet for the new rows (step 1).
##   (a) LD windows (NSCC PLINK, rule of coloc_eqtl_ldwin/01_ldwin_nscc.pbs) for: new-row signal SNPs absent from
##       lead_windows.tsv, and the GTEx eGenes-file top variants of genes that have no top-variant window yet;
##   (b) recursive clumping (NSCC, sens_clump/s2_clump_nscc.R) for eGenes with no clumped leads yet (their GTEx significant
##       pairs, as sens_clump/s1_sig_pairs.py);
##   (c) eQTL extraction need: tissue x gene whose new coloc regions extend beyond the already extracted span.
## Output: results/coloc_eqtl_v3m/new_inputs/ (snps/<tissue>.<chr>.txt, egene_top.tsv, clump_genes.tsv, extract_genes.tsv)
suppressPackageStartupMessages(library(data.table))
R2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"; O <- file.path(R2, "results/coloc_eqtl_v3m"); LW <- file.path(R2, "results/coloc_eqtl_ldwin")
NI <- file.path(O, "new_inputs"); dir.create(file.path(NI, "snps"), recursive = TRUE, showWarnings = FALSE)
EQ <- "/home/e1327989/project/GTEx/Rawdata/gtex_v11/GTEx_Analysis_v11_eQTL"
A <- readRDS(file.path(O, "analysis_set.rds")); N <- A[carried == FALSE]
M <- fread(file.path(O, "missing.tsv"))
L <- fread(file.path(LW, "lead_windows.tsv")); RW <- fread(file.path(LW, "rarelead_windows.tsv"))
have_win <- unique(rbind(L[absent == FALSE, .(tissue, snp)], RW[, .(tissue, snp)]))
## GTEx top variant of every new-row gene
EL <- fread(file.path(LW, "eqtl_leads.tsv"))
G <- unique(N[, .(tissue, gene_base, chr)])
G0 <- merge(G, EL[, .(tissue, gene_base, eqtl_lead)], by = c("tissue", "gene_base"), all.x = TRUE)
need_top <- G0[is.na(eqtl_lead)]
E <- rbindlist(lapply(unique(need_top$tissue), function(t) {
  e <- fread(cmd = paste("zcat", file.path(EQ, paste0(t, ".v11.eGenes.txt.gz"))), select = c("gene_id", "variant_id", "pval_nominal", "qval"))
  e[, .(tissue = t, gene_base = sub("\\..*$", "", gene_id), eqtl_lead = variant_id, eqtl_lead_p = pval_nominal, eqtl_qval = qval)] }))
NT <- merge(need_top[, .(tissue, gene_base, chr)], E, by = c("tissue", "gene_base"), all.x = TRUE)
cat("new-row genes:", nrow(G), "; without an eqtl_leads row:", nrow(need_top), "; of those without an eGenes-file row:", NT[is.na(eqtl_lead), .N], "\n")
fwrite(NT, file.path(NI, "egene_top.tsv"), sep = "\t")
TOP <- rbind(G0[!is.na(eqtl_lead), .(tissue, gene_base, eqtl_lead)], NT[!is.na(eqtl_lead), .(tissue, gene_base, eqtl_lead)])
## (a) SNPs needing an LD window
need_snp <- unique(rbind(M[, .(tissue, snp = signal_snp, role = "puQTL")],
                         TOP[!paste(tissue, eqtl_lead) %in% have_win[, paste(tissue, snp)], .(tissue, snp = eqtl_lead, role = "eQTL")]))
need_snp[, chr := sub("_.*$", "", snp)]
cat("SNPs needing an LD window:", nrow(need_snp), "(puQTL", need_snp[role == "puQTL", .N], "; eQTL top", need_snp[role == "eQTL", .N], ")\n")
fwrite(need_snp, file.path(NI, "ld_snps.tsv"), sep = "\t")
for (k in split(need_snp, by = c("tissue", "chr"))) writeLines(unique(k$snp), file.path(NI, "snps", paste0(k$tissue[1], ".", k$chr[1], ".txt")))
## (b) eGenes without clumped leads
CL <- unique(fread(cmd = paste("zcat", file.path(LW, "sens_clump/pairs.tsv.gz")))[, .(tissue, gene_base)])
CG <- unique(N[egene == TRUE & !paste(tissue, gene_base) %in% CL[, paste(tissue, gene_base)], .(tissue, gene_base, chr)])
cat("eGenes needing clumping:", nrow(CG), "\n"); fwrite(CG, file.path(NI, "clump_genes.tsv"), sep = "\t")
## (c) eQTL extraction coverage for the rows whose windows exist (pairs of step 1); new genes always need extraction
EW <- rbind(fread(file.path(LW, "eqtl_windows.tsv"))[, src := "top"], fread(file.path(LW, "sens_clump/eqtl_windows.tsv"))[, src := "clump"])
PT <- fread(cmd = paste("zcat", file.path(O, "pairs_top_new.tsv.gz"))); PC <- fread(cmd = paste("zcat", file.path(O, "pairs_clump_new.tsv.gz")))
sp <- rbind(PT[, .(tissue, gene_base, chr, start = uni_start, end = uni_end, src = "top")], PC[, .(tissue, gene_base, chr, start = uni_start, end = uni_end, src = "clump")])
sp <- sp[, .(start = min(start), end = max(end)), by = .(tissue, gene_base, chr, src)]
cv <- merge(sp, EW, by = c("tissue", "gene_base", "chr", "src"), all.x = TRUE, suffixes = c("", "_have"))
cv[, covered := !is.na(start_have) & start >= start_have & end <= end_have]
cat("eQTL coverage of new pair regions: top", cv[src == "top", sum(covered)], "/", cv[src == "top", .N], "; clump", cv[src == "clump", sum(covered)], "/", cv[src == "clump", .N], "\n")
X <- rbind(cv[covered == FALSE, .(tissue, gene_base, chr, src)], unique(N[paste(tissue, signal_snp) %in% M[, paste(tissue, signal_snp)], .(tissue, gene_base, chr)])[, src := "new gene/SNP"])
fwrite(unique(X), file.path(NI, "extract_genes.tsv"), sep = "\t"); cat("tissue x gene needing (re)extraction:", uniqueN(X[, .(tissue, gene_base)]), "\n")
