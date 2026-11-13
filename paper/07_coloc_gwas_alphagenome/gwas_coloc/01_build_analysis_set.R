#!/usr/bin/env Rscript
# One row per coloc "analysis signal": every significant primary unit x its independent signals, plus every
# non-significant (primary, estimable) MRAS unit (marginal, for the MRAS tissue panel). Also writes the union of
# signal windows used to pre-cut the GWAS files.
# Backward-check rule (25 signals fail in v3; treat as not independent):
#   drop failed signals; a unit left with ONE kept signal uses the MARGINAL scan centred on that signal
#   ("kept1_marginal"); a unit left with >= 2 uses the LOO scans of the kept signals, flagged
#   "loo_conditions_on_failed" (those LOO scans still condition on the failed signal).
source("/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/scripts/coloc_final_lib.R")
U <- read_units_final()
S <- fread(cmd = paste("zcat", file.path(V2, COND_ROOT, "freeze_real_signals.tsv.gz")))[, .(tissue, promoter, type, k, snp, backward_keep)]
A <- merge(U[sig == TRUE], S, by = c("tissue", "promoter", "type"))
stopifnot(uniqueN(A[, .(tissue, promoter, type)]) == U[sig == TRUE, .N])
A[, n_signals := .N, by = .(tissue, promoter, type)]
A[, n_kept := sum(backward_keep), by = .(tissue, promoter, type)]
A <- A[backward_keep == TRUE]
A[, stat := fifelse(n_kept >= 2, "loo", "marginal")]
A[, rule := fcase(n_signals == 1, "single", n_kept == 1, "kept1_marginal", n_kept < n_signals, "loo_conditions_on_failed", default = "loo")]
M <- U[promoter %in% MRAS & primary == TRUE & sig == FALSE]
M[, `:=`(k = 1L, snp = leader, backward_keep = TRUE, n_signals = NA_integer_, n_kept = NA_integer_, stat = "marginal", rule = "mras_nonsig")]
A <- rbind(A, M, fill = TRUE)
A[, `:=`(signal_snp = snp, signal_pos = pos_of(snp), chr = chr_of(snp))]
saveRDS(A, file.path(RES_DIR, "analysis_set.rds"))
iv <- unique(A[, .(chr = sub("^chr", "", chr), s = signal_pos - WIN, e = signal_pos + WIN)])
saveRDS(iv, file.path(RES_DIR, "analysis_windows.rds"))
cat("analysis signals:", nrow(A), "\n"); print(A[, .N, by = .(sig, stat, rule)][order(-N)])
cat("significant units:", uniqueN(A[sig == TRUE, .(tissue, promoter, type)]), " signals:", A[sig == TRUE, .N], "\n")
print(A[sig == TRUE, .N, by = type])
cat("MRAS rows:\n"); print(A[promoter %in% MRAS, .N, by = .(sig)])
