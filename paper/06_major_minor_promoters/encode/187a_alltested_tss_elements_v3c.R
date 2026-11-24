#!/usr/bin/env Rscript
## 187a_alltested_tss_elements_v3c.R -- TSS elements of every promoter tested in units v3c. Copy of
## Figure5_rerun_v3/scripts/25a_alltested_tss_elements.R with ONE change of input: the tested units come from units v3c
## (results/final_nogate_v3c/freeze, mirror rule v3) instead of units v3 (final_nogate_v3/freeze, mirror rule v2), so every promoter
## tested in v3c (the set 122a uses) gets its TSS-element lookup. Rule, cCRE file and window are unchanged (TSS +-200 bp, the
## overlapping cCRE closest to the TSS). Outputs go to a NEW folder; the old data/alltested_tss_elements.tsv is not touched.
## Output: <B1>/data/alltested_tss_elements_v3c.tsv, <B1>/data/zscore_accessions_alltested_v3c.txt
## (both linked into the Fig. 5 stage folder, in/fig5_v3data, under the names 122a reads).
suppressPackageStartupMessages(library(data.table)); setDTthreads(8)
B1 <- "/home/e1327989/project/GTEx/TATA_Enhancer/Figure5_rerun_v4/b1_tss_elements_v3c_20261113"
FZ <- "/home/e1327989/project/GTEx/QTL_figure/rerun_v2/results/final_nogate_v3c/freeze"
A <- rbindlist(lapply(list.files(FZ, pattern = "real.units.rds$", full.names = TRUE), function(f)
  as.data.table(readRDS(f))[, .(tissue, promoter, gene, type, q, estimable, dedup_drop)]))
A <- A[estimable == TRUE & dedup_drop == FALSE]
P <- A[, .(gene0 = sub("\\..*$", "", gene[1]), tested_types = paste(sort(unique(type)), collapse = ","),
           mu = any(type == "mu" & !is.na(q) & q <= .05), p0 = any(type == "p0" & !is.na(q) & q <= .05),
           p1 = any(type == "p1" & !is.na(q) & q <= .05)), by = promoter]
PI <- fread("/home/e1327989/project/GTEx/v11/Promoter_info_filter_5.txt"); setnames(PI, 1, "promoter")
n_tested <- nrow(P)
P <- merge(P, PI[, .(promoter, chr = seqnames, tss = start, strand)], by = "promoter")
stopifnot(nrow(P) == n_tested)   # every tested promoter has a TSS in Promoter_info
C <- fread("/home/e1327989/project/GTEx/TATA_Enhancer/ccre_pELS_dELS/GRCh38-cCREs.bed", header = FALSE, col.names = c("chr", "s", "e", "dhs", "acc", "type"))
C[, grp := fcase(grepl("dELS", type), "dELS", grepl("pELS", type), "pELS", grepl("PLS", type), "PLS", default = type)]; setkey(C, chr, s, e)
Q <- P[, .(promoter, chr, s = tss - 201L, e = tss + 200L, tss)]; setkey(Q, chr, s, e)
O <- foverlaps(Q, C, nomatch = 0L)[, dist := abs((s + 1 + e) / 2 - tss)][order(promoter, dist)][, .SD[1], by = promoter]
P <- merge(P, O[, .(promoter, acc, grp)], by = "promoter", all.x = TRUE); P[is.na(grp), grp := "none"]
fwrite(P, file.path(B1, "data/alltested_tss_elements_v3c.tsv"), sep = "\t")
writeLines(unique(P[!is.na(acc)]$acc), file.path(B1, "data/zscore_accessions_alltested_v3c.txt"))
cat(sprintf("tested promoters (v3c): %d; with a mu / p0 / p1 QTL: %d / %d / %d; TSS element found: %d\n",
            nrow(P), sum(P$mu), sum(P$p0), sum(P$p1), P[grp != "none", .N]))

## ---- check against the old (v3) table: shared promoters must get the same element ----
OLD <- fread("/home/e1327989/f5_stage_20261104_v4/in/fig5_v3data/alltested_tss_elements.tsv")
cm <- intersect(OLD$promoter, P$promoter)
x <- merge(OLD[promoter %in% cm, .(promoter, acc_old = acc, grp_old = grp, tss_old = tss)], P[promoter %in% cm, .(promoter, acc_new = acc, grp_new = grp, tss_new = tss)], by = "promoter")
cat(sprintf("old table: %d promoters; new: %d; shared: %d; only old: %d; only new: %d\n", nrow(OLD), nrow(P), length(cm),
            length(setdiff(OLD$promoter, P$promoter)), length(setdiff(P$promoter, OLD$promoter))))
cat(sprintf("shared promoters with identical TSS / element / element type: %d / %d / %d of %d\n", x[tss_old == tss_new, .N],
            x[(is.na(acc_old) & is.na(acc_new)) | (!is.na(acc_old) & !is.na(acc_new) & acc_old == acc_new), .N], x[grp_old == grp_new, .N], nrow(x)))
NEWP <- P[!promoter %in% OLD$promoter]
cat(sprintf("promoters new in v3c: %d; with a TSS element: %d; element type:\n", nrow(NEWP), NEWP[grp != "none", .N])); print(NEWP[, .N, by = grp][order(-N)])
cat("tested types of the new promoters:\n"); print(NEWP[, .N, by = tested_types])
## units (tissue, promoter) tested in v3c whose promoter was missing from the old table
U <- unique(A[, .(tissue, promoter, type, sig = !is.na(q) & q <= .05)])
U[, in_old := promoter %in% OLD$promoter]
cat(sprintf("tested (tissue, promoter) units: %d in %d promoters; without an old lookup: %d units in %d promoters\n",
            uniqueN(U[, .(tissue, promoter)]), uniqueN(U$promoter), uniqueN(U[in_old == FALSE, .(tissue, promoter)]), uniqueN(U[in_old == FALSE]$promoter)))
cat("significant units by class, without an old lookup / all:\n"); print(U[sig == TRUE, .(no_old_lookup = sum(!in_old), all = .N), by = type])
fwrite(NEWP, file.path(B1, "data/promoters_new_in_v3c.tsv"), sep = "\t")
