#!/usr/bin/env Rscript
# LD-window GWAS coloc, step 13 (array over traits): cut one GWAS file to the union regions of its candidate pairs (step 12).
# Same reader as the main run (coloc_final_lib.R::read_gwas: numeric P, SE > 0, first row per chr:pos).
# Usage: Rscript 13_subset_gwas.R <trait_index>    (idempotent)
source("/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/scripts/coloc_final_lib.R")
D <- LDU; dir.create(file.path(D, "gwas_subset"), showWarnings = FALSE)
tr <- TRAITS[as.integer(commandArgs(TRUE)[1])]
of <- file.path(D, "gwas_subset", sub("\\.txt$", ".rds", tr)); if (file.exists(of)) quit(save = "no")
dir.create(dirname(of), showWarnings = FALSE)
iv <- readRDS(file.path(D, "union_intervals", sub("\\.txt$", ".rds", tr)))
if (!nrow(iv)) { saveRDS(data.table(gwas_id = character(), BETA = numeric(), SE = numeric(), Pvalue = numeric()), of); cat(tr, ": no pairs\n"); quit(save = "no") }
setorder(iv, chr, s)                                                   # merge overlapping unions into disjoint blocks
iv[, grp := cumsum(c(1L, as.integer(s[-1] > cummax(e)[-.N]))), by = chr]
iv <- iv[, .(s = min(s), e = max(e)), by = .(chr, grp)][, grp := NULL]
setkey(iv, chr, s, e)
gw <- read_gwas(tr)
gw[, `:=`(chr = sub(":.*", "", gwas_id), s = as.integer(sub(".*:", "", gwas_id)))]; gw[, e := s]; setkey(gw, chr, s, e)
hit <- unique(foverlaps(gw, iv, type = "within", nomatch = 0L, which = TRUE)$xid)
saveRDS(gw[hit, .(gwas_id, BETA, SE, Pvalue)], of); cat(sprintf("%s: %d -> %d SNPs\n", tr, nrow(gw), length(hit)))
