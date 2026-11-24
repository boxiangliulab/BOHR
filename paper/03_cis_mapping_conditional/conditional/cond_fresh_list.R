## Stage S3 of the final run (authored; a fresh conditional analysis in place of the earlier incremental runs):
## every significant unit of units v3m (q <= 0.05) is run from its leader with run_conditional_v3m.R ("new"), and the
## stubs that let merge_conditional_v3m.R and merge_conditional_v3c.R run without an earlier conditional analysis.
## Usage (cwd = the run root): Rscript cond_fresh_list.R [tissue ...]   (default: every tissue with units v3m)
## Writes: work/cond_v3m/{units_to_run.tsv, jobs.tsv}; results/conditional_final_v3/ (stub: one placeholder unit row per
##   tissue, header-only signal and round tables, empty loo/<tissue>.rds); work/cond_v3c/units_to_run.tsv (header only).
suppressMessages(library(data.table))
a <- commandArgs(TRUE)
fs <- Sys.glob("results/final_nogate_v3m/freeze/*.real.units.rds")
U <- rbindlist(lapply(fs, function(f) as.data.table(readRDS(f))[, .(tissue, promoter, type, p, q)]))
TIS <- if (length(a)) a else sort(unique(U$tissue)); U <- U[tissue %in% TIS]
thr <- U[, .(thr_new = suppressWarnings(max(p[q <= 0.05], na.rm = TRUE))), by = tissue]
S <- merge(U[q <= 0.05], thr, by = "tissue")
chr <- rbindlist(lapply(TIS, function(t) { td <- readRDS(file.path("work/final/freeze", paste0(t, ".rds")))
  data.table(tissue = t, promoter = rownames(td$info), chr = as.character(td$info$seqnames)) }))
S <- merge(S, chr, by = c("tissue", "promoter"), all.x = TRUE); stopifnot(!anyNA(S$chr))
dir.create("work/cond_v3m", recursive = TRUE, showWarnings = FALSE); dir.create("work/cond_v3c", recursive = TRUE, showWarnings = FALSE)
fwrite(S[, .(tissue, chr, promoter, type, reason = "new", p, q, thr_new, thr_old = NA_real_)][order(tissue, chr, promoter, type)], "work/cond_v3m/units_to_run.tsv", sep = "\t")
J <- S[, .(n_units = .N), by = .(tissue, chr)][order(-n_units)]
fwrite(J, "work/cond_v3m/jobs.tsv", sep = "\t")
fwrite(data.table(tissue = character(), chr = character(), promoter = character(), type = character(), reason = character()), "work/cond_v3c/units_to_run.tsv", sep = "\t")
## stub of the earlier conditional analysis (no unit is carried: every significant unit is in units_to_run.tsv; the merge
## keeps only rows whose unit is significant in units v3m, so the placeholders never reach the output)
O <- "results/conditional_final_v3"; dir.create(file.path(O, "loo"), recursive = TRUE, showWarnings = FALSE)
fwrite(data.table(tissue = TIS, promoter = "(none)", type = "(none)", stop = NA_character_, n_signals = NA_integer_, chr = NA_character_,
                  n_signals_backward = NA_integer_), file.path(O, "freeze_real_units.tsv.gz"), sep = "\t")
## placeholder rows (promoter "(none)") carry the column types: header-only tables would be read back as logical columns,
## and rbindlist() would then coerce the tissue names of the new rows to NA
fwrite(data.table(tissue = TIS, promoter = "(none)", type = "(none)", k = 1L, snp = "(none)", loo_est = NA_real_, loo_acat = NA_real_,
                  loo_top_snp = NA_character_, backward_keep = NA), file.path(O, "freeze_real_signals.tsv.gz"), sep = "\t")
fwrite(data.table(tissue = TIS, promoter = "(none)", type = "(none)", round = 1L, acat_p = NA_real_), file.path(O, "freeze_real_rounds.tsv.gz"), sep = "\t")
for (t in TIS) saveRDS(data.table(), file.path(O, "loo", paste0(t, ".rds")))
writeLines("STUB written by cond_fresh_list.R for a fresh conditional run (no carried units).", file.path(O, "README_STUB.txt"))
cat("tissues", length(TIS), "; significant units", nrow(S), "(", paste(S[, .N, by = type][, paste(type, N)], collapse = ", "), "); (tissue, chr) jobs", nrow(J), "; max units per job", max(J$n_units), "\n")
