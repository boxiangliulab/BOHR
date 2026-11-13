#!/usr/bin/env Rscript
### REALJUNC copy of /data/projects/11003054/YuhaoDong/GTEx/simulation/sim_proactiv.R (2026-10-06).
### Author: the real junction files keep junctions with >= 5 UNIQUE reads (no multi-mapping clause).
### Only changes vs the original:
###   input  = star_output/*SJ.out_realfilter.tab  (awk '$7 >= 5' on the unfiltered *_SJ.out.tab; made by pbs_realjunc_M2.sh)
###   output = <sim_dir>/realjunc/{proactiv_result,recal_att_f,promoter_info}.RData  (the original files are untouched)
###   + realjunc/count_check.txt: minimum nonzero promoter count and share of 1-4-read counts, new vs original result.
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) stop("Usage: Rscript sim_proactiv_realjunc.R <sim_dir>")
sim_dir <- args[1]; setwd(sim_dir); dir.create("realjunc", showWarnings = FALSE)
cat("=== proActiv (real junction rule: unique reads >= 5) ===\nWorking directory:", sim_dir, "\n")
library(proActiv); library(stringr); library(SummarizedExperiment)
promoterAnnotation.gencode.v47 <- preparePromoterAnnotation(
    file = "/home/users/nus/e1327989/YuhaoDong/reference/gencode.v47.annotation.gtf.gz", species = 'Homo_sapiens')
files <- Sys.glob('star_output/*SJ.out_realfilter.tab')
cat("Found", length(files), "real-rule SJ files\n"); if(length(files) == 0) stop("No real-rule SJ files")
condition <- rep("test", length(files))
result <- proActiv(files = files, promoterAnnotation = promoterAnnotation.gencode.v47, condition = condition)
result <- result[complete.cases(assays(result)$promoterCounts),]
save(result, file='realjunc/proactiv_result.RData'); cat("proActiv result saved:", nrow(result), "promoters\n")
library(dplyr); library(tidyr); library(tibble); library(data.table); library(Matrix)
rdata <- rowData(result); row.names(rdata) <- paste0('PromoterID_', row.names(rdata))
abs_av <- as.data.frame(assays(result)[[3]]); row.names(abs_av) <- row.names(rdata)
rel_av <- as.data.frame(assays(result)[[4]]); row.names(rel_av) <- row.names(rdata)
trel <- t(rel_av); trel_att_raw <- trel
irdata <- rdata[rdata$internalPromoter == 'FALSE',]; in_list <- row.names(irdata); in_info <- as.data.frame(irdata)
in_trel <- trel_att_raw[, colnames(trel_att_raw) %in% in_list]
abs_att <- abs_av[row.names(abs_av) %in% colnames(in_trel), ]; abs_info <- cbind(in_info[2], abs_att)
calculate_promoter_ratios_fast <- function(abs_info) {
  dt <- as.data.table(as.data.frame(abs_info)); promoterID <- rownames(abs_info); dt[, promoterID := promoterID]
  setcolorder(dt, c("promoterID", "geneId")); mat <- as.matrix(dt[, -c(1:2)]); row.names(mat) <- dt$promoterID
  genes <- as.factor(dt$geneId); n <- nrow(mat); p <- length(levels(genes))
  group_mat <- sparseMatrix(i = seq_len(n), j = as.integer(genes), x = rep(1, n), dims = c(n, p))
  gene_totals <- as.matrix(t(group_mat) %*% mat); ratios <- mat / gene_totals[genes, ]
  result <- data.table(promoterID = promoterID, geneId = dt$geneId, as.data.table(ratios)); setkey(result, promoterID); result }
result_ratios <- calculate_promoter_ratios_fast(abs_info)
recal <- as.data.frame(result_ratios); row.names(recal) <- recal$promoterID; recal <- recal[-c(1:2)]
recal_att <- as.data.frame(t(recal))
recal_att_f <- recal_att[, !(colSums(!is.na(recal_att)) == 0 | (colSums(recal_att == 0, na.rm = TRUE) == colSums(!is.na(recal_att))))]
save(recal_att_f, file='realjunc/recal_att_f.RData')
info <- rdata[row.names(rdata) %in% colnames(recal_att_f), 1:8]; save(info, file='realjunc/promoter_info.RData')
## count check: new (real rule) vs the original result in sim_dir
chk <- function(r, lab) { pc <- as.matrix(assays(r)$promoterCounts); nz <- pc[is.finite(pc) & pc > 0]
  sprintf("%s: min nonzero promoter count %g | share of nonzero counts in [1, 5) %.6f | nonzero counts %d", lab, min(nz), mean(nz < 5), length(nz)) }
old <- new.env(); load("proactiv_result.RData", envir = old)
writeLines(c(chk(result, "real rule (unique >= 5)"), chk(old$result, "original filter ($7 > 5 || $8 > 5)")), "realjunc/count_check.txt")
cat(readLines("realjunc/count_check.txt"), sep = "\n")
cat("Promoters retained:", ncol(recal_att_f), " Samples:", nrow(recal_att_f), "\nREALJUNC_DONE\n")
