#!/usr/bin/env Rscript
# LD-window GWAS coloc, step 12a: LD window of every puQTL analysis signal (tissue donors), taken from the shared table
# that Figure 4 computed for the same 14,783 signals (QTL_figure/rerun_v2/results/coloc_eqtl_ldwin/lead_windows.tsv,
# any role; PLINK 1.9 --r2 >= 0.1, panel, tissue donors, +-1 Mb; window = farthest partner +/- 50 kb).
# Checks that Figure 4's analysis set equals ours, and lists the leads with no window (the non-significant MRAS units
# of the MRAS panel) for an NSCC run (step 12b).
# Output: puqtl_windows.tsv (one row per analysis-set row), puqtl_leads_missing.tsv (tissue, snp).
source("/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/scripts/coloc_final_lib.R")
D  <- LDU
F4 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/results"
A  <- readRDS(file.path(RES_DIR, "analysis_set.rds"))
A4 <- readRDS(file.path(F4, paste0("coloc_eqtl_", UNITSET), "analysis_set.rds"))
k1 <- A[sig == TRUE, paste(tissue, promoter, type, k, signal_snp)]; k4 <- A4[, paste(tissue, promoter, type, k, signal_snp)]
cat("our significant analysis signals:", length(k1), "; Figure 4's:", length(k4), "; identical key sets:", setequal(k1, k4), "\n")
stopifnot(setequal(k1, k4))
L <- fread(file.path(F4, "coloc_eqtl_ldwin/lead_windows.tsv"))[absent == FALSE,   # roles puQTL and eQTL+puQTL (a window depends on tissue x SNP only)
          .(tissue, snp, p_partners = n_partners, p_start = win_start, p_end = win_end)]
stopifnot(!anyDuplicated(L[, .(tissue, snp)]))
W <- merge(A[, .(tissue, promoter, type, k, sig, rule, chr, signal_snp, signal_pos)], L, by.x = c("tissue", "signal_snp"),
           by.y = c("tissue", "snp"), all.x = TRUE)
stopifnot(nrow(W) == nrow(A))
miss <- unique(W[is.na(p_start), .(tissue, snp = signal_snp, chr)])
cat("analysis rows:", nrow(W), "; with a window:", W[!is.na(p_start), .N], "; missing:", W[is.na(p_start), .N],
    "(significant:", W[is.na(p_start) & sig == TRUE, .N], ") ; distinct missing tissue x lead:", nrow(miss), "\n")
print(W[is.na(p_start), .N, by = .(sig, rule)])
fwrite(W, file.path(D, "puqtl_windows.tsv"), sep = "\t", quote = FALSE)
fwrite(miss, file.path(D, "puqtl_leads_missing.tsv"), sep = "\t", quote = FALSE)
cat("window width (kb), significant signals:\n"); print(W[sig == TRUE & !is.na(p_start), quantile((p_end - p_start) / 1e3, c(0, .1, .5, .9, .99, 1))])
