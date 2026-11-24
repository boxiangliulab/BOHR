#!/usr/bin/env Rscript
## 187c_check_zscores_v3c.R -- check: the new z-score extracts (187b) must reproduce the earlier extracts (22b2) exactly for every
## accession in both, and the biosample column files must be identical.
suppressPackageStartupMessages(library(data.table)); setDTthreads(4)
B1 <- "/home/e1327989/project/GTEx/TATA_Enhancer/Figure5_rerun_v4/b1_tss_elements_v3c_20261113/data"
OD <- "/home/e1327989/f5_stage_20261104_v4/in/fig5_v3data"
for (m in c("h3k4me3", "h3k27ac", "dnase")) {
  o <- fread(cmd = sprintf("zcat %s/zscores_%s_alltested.tsv.gz", OD, m)); n <- fread(cmd = sprintf("zcat %s/zscores_%s_alltested_v3c.tsv.gz", B1, m))
  oc <- fread(file.path(OD, sprintf("zscores_%s_alltested_columns.tsv", m))); nc <- fread(file.path(B1, sprintf("zscores_%s_alltested_v3c_columns.tsv", m)))
  cm <- intersect(o$acc, n$acc); setkey(o, acc); setkey(n, acc)
  same <- isTRUE(all.equal(as.matrix(o[cm, -1]), as.matrix(n[cm, -1]), check.attributes = FALSE, tolerance = 0))
  cat(sprintf("%s: old %d acc, new %d acc, shared %d, old-only %d; values identical on shared: %s; columns file identical: %s\n", m, nrow(o), nrow(n),
              length(cm), length(setdiff(o$acc, n$acc)), same, isTRUE(all.equal(oc, nc))))
  stopifnot(same, isTRUE(all.equal(oc, nc)))
}
