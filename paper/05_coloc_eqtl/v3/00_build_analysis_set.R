## Figure 4 panel e: puQTL x eQTL coloc on units v3, with Figure 6's coloc machinery (the same coloc as the GWAS analysis).
## Step 0: the analysis set = every significant primary unit x its backward-kept signals, built with Figure 6's rules
## (new_GWAS/coloc_final_nogate_v3/scripts/01_build_analysis_set.R, significant units only):
##   n_kept >= 2 -> LOO-conditional scan of each kept signal; n_kept == 1 -> marginal scan centred on that signal.
## Adds gene_base, the eGene flag of the unit's gene in that tissue (GTEx v11 eGenes, qval <= 0.05), the eQTL N (the
## tissue's freeze donors, results/ld_v3/donors) and writes the +-200 kb windows per (tissue, gene) for the eQTL extraction.
source("/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/scripts/coloc_final_lib.R")
OUT <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/results/coloc_eqtl_v3"
EQ  <- "/home/e1327989/project/GTEx/Rawdata/gtex_v11/GTEx_Analysis_v11_eQTL"
U <- read_units_final()
S <- fread(cmd = paste("zcat", file.path(V2, "conditional_final_v3/freeze_real_signals.tsv.gz")))[, .(tissue, promoter, type, k, snp, backward_keep)]
A <- merge(U[sig == TRUE], S, by = c("tissue", "promoter", "type"))
stopifnot(uniqueN(A[, .(tissue, promoter, type)]) == U[sig == TRUE, .N])
A[, n_signals := .N, by = .(tissue, promoter, type)]
A[, n_kept := sum(backward_keep), by = .(tissue, promoter, type)]
A <- A[backward_keep == TRUE]
A[, stat := fifelse(n_kept >= 2, "loo", "marginal")]
A[, rule := fcase(n_signals == 1, "single", n_kept == 1, "kept1_marginal", n_kept < n_signals, "loo_conditions_on_failed", default = "loo")]
A[, `:=`(signal_snp = snp, signal_pos = pos_of(snp), chr = chr_of(snp), gene_base = sub("\\..*$", "", gene))]
## eGenes (GTEx v11, qval <= 0.05) and eQTL N
EG <- rbindlist(lapply(unique(A$tissue), function(t) { e <- fread(cmd = paste("zcat", file.path(EQ, paste0(t, ".v11.eGenes.txt.gz"))), select = c("gene_id", "qval"))
  e[, .(tissue = t, gene_base = sub("\\..*$", "", gene_id), eqtl_qval = qval)] }))
A <- merge(A, EG, by = c("tissue", "gene_base"), all.x = TRUE)
A[, `:=`(egene = !is.na(eqtl_qval) & eqtl_qval <= 0.05, eqtl_tested = !is.na(eqtl_qval))]
NE <- data.table(tissue = unique(A$tissue))[, n_eqtl := vapply(tissue, function(t) length(readLines(file.path(V2, "ld_v3/donors", paste0(t, ".txt")))), 1L)]
A <- merge(A, NE, by = "tissue")
setorder(A, tissue, promoter, type, k)
saveRDS(A, file.path(OUT, "analysis_set.rds"))
## windows for the eQTL extraction: one row per (tissue, gene, chr) with the union of the signal windows
W <- A[, .(start = min(signal_pos) - WIN, end = max(signal_pos) + WIN), by = .(tissue, gene_base, chr)]
fwrite(W, file.path(OUT, "eqtl_windows.tsv"), sep = "\t")
cat("analysis signals:", nrow(A), "; units:", uniqueN(A[, .(tissue, promoter, type)]), "; genes x tissues:", uniqueN(A[, .(tissue, gene_base)]), "\n")
print(A[, .N, by = .(stat, rule)][order(-N)])
cat("units on eGenes:", uniqueN(A[egene == TRUE, .(tissue, promoter, type)]), "; on genes tested but not eGenes:", uniqueN(A[eqtl_tested & !egene, .(tissue, promoter, type)]),
    "; on genes not in the eQTL tables:", uniqueN(A[eqtl_tested == FALSE, .(tissue, promoter, type)]), "\n")
## Figure 6's own analysis set must agree on the significant rows
F6 <- readRDS(file.path(NDIR, "results/analysis_set.rds"))[sig == TRUE]
cat("rows vs Figure 6 analysis set:", nrow(A), "/", nrow(F6), "; same signal rows:", nrow(fintersect(A[, .(tissue, promoter, type, k, signal_snp, stat)], F6[, .(tissue, promoter, type, k, signal_snp, stat)])), "\n")
