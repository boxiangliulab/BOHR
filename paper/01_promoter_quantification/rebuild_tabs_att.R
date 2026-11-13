#!/usr/bin/env Rscript
# Rebuild the CORRECT filter_5_tabs_att.RData (absolute promoter activity).
# The original file's sample rows were misaligned (cbind-by-position bug).
# Correct absolute source: tabs_att_filter_In.rds $abs (name-indexed, never went through the bad cbind).
# Target layout must match the intended tabs_att: rows = samples (trel_att order),
#   cols = sample/SMTS/SMTSD + the same promoter set & order as trel_att.
suppressMessages({library(Matrix)})
setwd("/home/e1327989/project/GTEx/v11")

trel_att <- get(load("filter_5_trel_att.RData"))          # samples x (3 meta + promoters), reference layout
meta  <- trel_att[, 1:3]                                   # sample, SMTS, SMTSD (correct, internally consistent)
proms <- colnames(trel_att)[-(1:3)]                        # promoter columns, in canonical order
ab    <- readRDS("tabs_att_filter_In.rds")$abs             # samples x all-non-internal promoters (CORRECT abs)

## checks before building
stopifnot(all(rownames(trel_att) %in% rownames(ab)))       # every sample present in abs source
stopifnot(all(proms %in% colnames(ab)))                    # every promoter present in abs source

## name-indexed extraction -> guarantees correct sample & promoter alignment
absmat <- as.matrix(ab[rownames(trel_att), proms])         # rows ordered like trel_att, cols like proms
tabs_att <- cbind(meta, as.data.frame(absmat))             # now meta and data share the SAME row order

## ---- verification: rebuilt absolute must round-trip to abs_all (cor == 1) ----
set.seed(1); chk_s <- sample(rownames(tabs_att), 2000); chk_p <- sample(proms, 50)
rebuilt <- as.matrix(tabs_att[chk_s, chk_p])
source0 <- as.matrix(ab[chk_s, chk_p])
maxdiff <- max(abs(rebuilt - source0))
cat(sprintf("verification max|diff| rebuilt vs abs_all on 2000x50 spot grid: %.3g\n", maxdiff))
stopifnot(maxdiff < 1e-9)

## cross-check against relative (recal): abs_i/sum(MRAS 4) must equal recal for MRAS
tryCatch({
  P4 <- paste0("PromoterID_", 45111:45114)
  recal <- get(load("filter_5_recal_att_filter.RData"))    # recalculated relative (downstream file)
  common <- head(intersect(rownames(tabs_att), rownames(recal)), 5000)
  A4 <- as.matrix(tabs_att[common, P4]); tot <- rowSums(A4)
  ok <- is.finite(tot) & tot > 0 & rowSums(is.na(as.matrix(recal[common, P4]))) == 0
  cat(sprintf("MRAS recal cross-check max|diff|: %.3g (should be ~1e-15)\n",
              max(abs(A4[ok, ] / tot[ok] - as.matrix(recal[common, P4])[ok, ]))))
}, error = function(e) cat("recal cross-check skipped:", conditionMessage(e), "\n"))

## back up the corrupted file, then overwrite
if (!file.exists("filter_5_tabs_att.BAD.RData"))
  file.copy("filter_5_tabs_att.RData", "filter_5_tabs_att.BAD.RData")
cat("backed up old (corrupted) file -> filter_5_tabs_att.BAD.RData\n")

save(tabs_att, file = "filter_5_tabs_att.RData")
cat(sprintf("WROTE corrected filter_5_tabs_att.RData : %d samples x %d cols (%d promoters)\n",
            nrow(tabs_att), ncol(tabs_att), ncol(tabs_att) - 3))
