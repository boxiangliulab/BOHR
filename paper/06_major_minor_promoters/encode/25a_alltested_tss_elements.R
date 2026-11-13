#!/usr/bin/env Rscript
# 25a_alltested_tss_elements.R -- TSS elements for EVERY tested promoter (reference group for the enhancer-like-promoter
# panel), plus the promoter's QTL classes.
# Tested promoters: Figure 4's v3 freeze units with estimable == TRUE and dedup_drop == FALSE (the set entering BH).
# A promoter "has a <type>-QTL" when any of its units of that type has q <= 0.05 (= the v3 delivery).
# TSS element: the cCRE overlapping TSS +-200 bp closest to the TSS (as 22a).
# Output: data/alltested_tss_elements.tsv, data/zscore_accessions_alltested.txt (for 22b, suffix _alltested).
suppressPackageStartupMessages(library(data.table)); setDTthreads(8)
W <- "/home/e1327989/project/GTEx/TATA_Enhancer/Figure5_rerun_v3"
FZ <- "/home/e1327989/project/GTEx/QTL_figure/rerun_v2/results/final_nogate_v3/freeze"
A <- rbindlist(lapply(list.files(FZ, pattern = "real.units.rds$", full.names = TRUE), function(f)
  as.data.table(readRDS(f))[, .(tissue, promoter, gene, type, q, estimable, dedup_drop)]))
A <- A[estimable == TRUE & dedup_drop == FALSE]
P <- A[, .(gene0 = sub("\\..*$", "", gene[1]), tested_types = paste(sort(unique(type)), collapse = ","),
           mu = any(type == "mu" & !is.na(q) & q <= .05), p0 = any(type == "p0" & !is.na(q) & q <= .05),
           p1 = any(type == "p1" & !is.na(q) & q <= .05)), by = promoter]
PI <- fread("/home/e1327989/project/GTEx/v11/Promoter_info_filter_5.txt"); setnames(PI, 1, "promoter")
P <- merge(P, PI[, .(promoter, chr = seqnames, tss = start, strand)], by = "promoter")
C <- fread("/home/e1327989/project/GTEx/TATA_Enhancer/ccre_pELS_dELS/GRCh38-cCREs.bed", header = FALSE, col.names = c("chr", "s", "e", "dhs", "acc", "type"))
C[, grp := fcase(grepl("dELS", type), "dELS", grepl("pELS", type), "pELS", grepl("PLS", type), "PLS", default = type)]; setkey(C, chr, s, e)
Q <- P[, .(promoter, chr, s = tss - 201L, e = tss + 200L, tss)]; setkey(Q, chr, s, e)
O <- foverlaps(Q, C, nomatch = 0L)[, dist := abs((s + 1 + e) / 2 - tss)][order(promoter, dist)][, .SD[1], by = promoter]
P <- merge(P, O[, .(promoter, acc, grp)], by = "promoter", all.x = TRUE); P[is.na(grp), grp := "none"]
fwrite(P, file.path(W, "data/alltested_tss_elements.tsv"), sep = "\t")
writeLines(unique(P[!is.na(acc)]$acc), file.path(W, "data/zscore_accessions_alltested.txt"))
cat(sprintf("tested promoters: %d; with a mu / p0 / p1 QTL: %d / %d / %d; TSS element found: %d\n",
            nrow(P), sum(P$mu), sum(P$p0), sum(P$p1), P[grp != "none", .N]))
