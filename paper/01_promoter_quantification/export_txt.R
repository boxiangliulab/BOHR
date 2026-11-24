## Authored, not an original script (ran in the final run). The text exports of a 0103 output folder that later steps read:
##   filter_5_recal_att_filter.txt   relative activity, rows = samples (row names), columns sample (donor), SMTS, SMTSD, promoters
##   Promoter_info_filter_5.txt      promoter table (txId list joined by "|")
##   filter_5_recal_att_filter_SMTSD_mean.txt   PEER input: mean over a donor's samples within each tissue (na.rm), NA -> 0,
##                                   columns sample, SMTSD, SMTS, promoters (recipe of 0104_peer.R lines 48-57)
## Usage: Rscript export_txt.R <0103 output dir> <out dir for the two exports> <out dir for the PEER input>
suppressMessages(library(data.table))
a <- commandArgs(TRUE); IN <- a[1]; OUT <- a[2]; OUTP <- a[3]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE); dir.create(OUTP, recursive = TRUE, showWarnings = FALSE)
e <- new.env(); load(file.path(IN, "filter_5_recal_att_filter.RData"), envir = e); X <- e$recal_att_f; rm(e)
e <- new.env(); load(file.path(IN, "Promoter_info_filter_5.RData"), envir = e); info <- e$info; rm(e)
pr <- grep("^PromoterID_", names(X), value = TRUE)
cat("samples", nrow(X), " promoters", length(pr), " info rows", nrow(info), "\n")
stopifnot(setequal(pr, rownames(info)))
fwrite(X[, c("sample", "SMTS", "SMTSD", pr)], file.path(OUT, "filter_5_recal_att_filter.txt"), sep = "\t", row.names = TRUE)
fwrite(info, file.path(OUT, "Promoter_info_filter_5.txt"), sep = "\t", row.names = TRUE, sep2 = c("", "|", ""))
## PEER input
M <- as.matrix(X[, pr]); g <- paste(X$sample, X$SMTSD, sep = "\r")
ok <- !is.na(M); M[!ok] <- 0
S <- rowsum(M, g, reorder = TRUE); N <- rowsum(ok * 1, g, reorder = TRUE)
R <- S / N; R[N == 0] <- 0
key <- do.call(rbind, strsplit(rownames(R), "\r", fixed = TRUE))
smts <- X$SMTS[match(rownames(R), g)]
D <- data.table(sample = key[, 1], SMTSD = key[, 2], SMTS = smts)
D <- cbind(D, as.data.table(R)); setorder(D, sample, SMTSD)
fwrite(D, file.path(OUTP, "filter_5_recal_att_filter_SMTSD_mean.txt"), sep = "\t")
cat("PEER input rows", nrow(D), "; done", format(Sys.time()), "\n")
