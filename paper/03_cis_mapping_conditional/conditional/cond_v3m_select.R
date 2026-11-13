## Mirror rule v3 (2026-11-04): which significant units need a (re)run of the forward conditional analysis.
## Current result = results/conditional_final_v3 (units v3, mirror v2). New primary set = results/final_nogate_v3m.
##   new        : significant under v3m with no current result (mostly p0(partner) twins of v2 p1(canonical))
##   continue   : current forward run stopped on "window ACAT above threshold" but that ACAT <= the new (looser) threshold
##   truncate   : a signal was added at a round whose ACAT lies in (thr_v3m, thr_v3]  (stricter threshold)
## All three are re-run from the leader (forward is deterministic, so continued units reproduce their prefix exactly).
## Output: work/cond_v3m/units_to_run.tsv (tissue, chr, promoter, type, reason) ; work/cond_v3m/jobs.tsv (tissue, chr)
suppressPackageStartupMessages(library(data.table))
V2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"; CF <- file.path(V2, "results/conditional_final_v3")
U3 <- fread(cmd = paste("zcat", file.path(CF, "freeze_real_units.tsv.gz"))); R3 <- fread(cmd = paste("zcat", file.path(CF, "freeze_real_rounds.tsv.gz")))
N <- rbindlist(lapply(Sys.glob(file.path(V2, "results/final_nogate_v3m/freeze/*.real.units.rds")), function(f) as.data.table(readRDS(f))[, .(tissue, promoter, type, p, q, q_mirror_v2)]))
thr <- N[, .(thr_new = max(p[q <= 0.05], na.rm = TRUE), thr_old = max(p[q_mirror_v2 <= 0.05], na.rm = TRUE)), by = tissue]
S <- N[q <= 0.05]; S <- merge(S, thr, by = "tissue")
cur <- U3[, .(tissue, promoter, type, stop, n_signals, chr)]
X <- merge(S, cur, by = c("tissue", "promoter", "type"), all.x = TRUE)
last <- R3[, .SD[which.max(round)], by = .(tissue, promoter, type)][, .(tissue, promoter, type, last_acat = acat_p)]
X <- merge(X, last, by = c("tissue", "promoter", "type"), all.x = TRUE)
## rounds where a signal was added = all rounds except a final one that stopped on the threshold
r2 <- merge(R3, U3[, .(tissue, promoter, type, stop, n_signals)], by = c("tissue", "promoter", "type"))
add <- r2[round < n_signals, .(max_added_acat = max(acat_p)), by = .(tissue, promoter, type)]
X <- merge(X, add, by = c("tissue", "promoter", "type"), all.x = TRUE)
X[, reason := fcase(is.na(stop), "new",
                    stop == "window ACAT above threshold" & is.finite(last_acat) & last_acat <= thr_new, "continue",
                    is.finite(max_added_acat) & max_added_acat > thr_new, "truncate", default = "keep")]
print(X[, .N, by = .(reason, type)][order(reason, type)])
RUN <- X[reason != "keep"]
## chr for new units: from any tissue's promoter info (promoter coordinates are tissue-independent)
allinfo <- unique(rbindlist(lapply(Sys.glob(file.path(V2, "work/final/freeze/*.rds")), function(f) { td <- readRDS(f); data.table(promoter = rownames(td$info), chr2 = as.character(td$info$seqnames)) })))
RUN <- merge(RUN, allinfo, by = "promoter", all.x = TRUE); RUN[is.na(chr), chr := chr2]; stopifnot(!anyNA(RUN$chr))
fwrite(RUN[, .(tissue, chr, promoter, type, reason, p, q, thr_new, thr_old)], file.path(V2, "work/cond_v3m/units_to_run.tsv"), sep = "\t")
J <- RUN[, .(n_units = .N), by = .(tissue, chr)][order(-n_units)]
fwrite(J, file.path(V2, "work/cond_v3m/jobs.tsv"), sep = "\t")
cat("units to run:", nrow(RUN), "; (tissue, chr) jobs:", nrow(J), "; max units per job:", max(J$n_units), "\n")
cat("units dropped from the primary significant set (current result discarded at merge):", cur[!paste(tissue, promoter, type) %in% S[, paste(tissue, promoter, type)], .N], "\n")
