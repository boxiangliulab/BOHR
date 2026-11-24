## Figure 4 panel e, LD-window coloc, step 0: the lead SNPs whose LD windows are needed.
## Lab rule: each lead's LD partners at r2 >= 0.1 (PLINK 1.9, GTEx v9), window = farthest partner
## +/- 50 kb; the two windows of a pair are valid when their intersection holds >= 5 SNPs with statistics in both
## datasets; coloc input = the union of the two windows (>= 10 shared SNPs kept). Leads:
##   puQTL: every analysis signal of coloc_eqtl_v3 (14,783; the LOO/marginal scan of each signal is unchanged);
##   eQTL : the GTEx v11 eGenes-file top variant of the unit's gene in that tissue.
## LD donors follow the ld_v3 standard: the tissue's freeze donors. Output: results/coloc_eqtl_ldwin/leads.tsv and
## per tissue x chr SNP lists for the NSCC PLINK run (ldwin_inputs/snps/<tissue>.<chrN>.txt).
suppressPackageStartupMessages(library(data.table))
V2  <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"
OUT <- file.path(V2, "results/coloc_eqtl_ldwin"); dir.create(file.path(OUT, "ldwin_inputs/snps"), recursive = TRUE, showWarnings = FALSE)
EQ  <- "/home/e1327989/project/GTEx/Rawdata/gtex_v11/GTEx_Analysis_v11_eQTL"
A <- readRDS(file.path(V2, "results/coloc_eqtl_v3/analysis_set.rds"))
E <- rbindlist(lapply(unique(A$tissue), function(t) {
  e <- fread(cmd = paste("zcat", file.path(EQ, paste0(t, ".v11.eGenes.txt.gz"))), select = c("gene_id", "variant_id", "pval_nominal", "qval"))
  e[, .(tissue = t, gene_base = sub("\\..*$", "", gene_id), eqtl_lead = variant_id, eqtl_lead_p = pval_nominal, eqtl_qval = qval)] }))
G <- merge(unique(A[, .(tissue, gene_base, chr)]), E, by = c("tissue", "gene_base"), all.x = TRUE)
cat("gene x tissue:", nrow(G), "; without an eGenes-file row:", G[is.na(eqtl_lead), .N], "; lead on another chr:",
    G[!is.na(eqtl_lead) & sub("_.*$", "", eqtl_lead) != chr, .N], "\n")
L <- unique(rbind(A[, .(tissue, chr, snp = signal_snp, role = "puQTL")],
                  G[!is.na(eqtl_lead), .(tissue, chr = sub("_.*$", "", eqtl_lead), snp = eqtl_lead, role = "eQTL")]))
L <- L[, .(role = paste(sort(unique(role)), collapse = "+")), by = .(tissue, chr, snp)]
fwrite(L, file.path(OUT, "leads.tsv"), sep = "\t"); fwrite(G, file.path(OUT, "eqtl_leads.tsv"), sep = "\t")
for (k in split(L, by = c("tissue", "chr"))) writeLines(unique(k$snp), file.path(OUT, "ldwin_inputs/snps", paste0(k$tissue[1], ".", k$chr[1], ".txt")))
cat("leads:", nrow(L), "(tissue x SNP); by role:\n"); print(L[, .N, by = role])
cat("SNP-list files:", length(list.files(file.path(OUT, "ldwin_inputs/snps"))), "\n")
