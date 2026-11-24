## =============================================================================
## prep.R -- cache the donor x promoter activity matrix (freeze definition) + donor metadata for each tissue
##   samples: SMAFRZE == "RNASEQ", keeping for each donor the row with the highest SMRIN (ties: the first row in the file),
##         and present in v10_cov_pca_full.txt
##   promoters: promoters of genes with GTEx v11 mean TPM >= 1 in the tissue (TPM_v11/TPM1_genes.txt);
##     (the previous round trait_regression_20261001 used v10 TPM, whose gene versions are older than the promoter annotation, missing about 18% of promoters)
##           CML has no TPM file; as in v13 all promoters are used
## Output: work/<tissue>.rds = list(tissue, U, row_id, act, info); work/meta.rds; work/tissues.tsv
## Usage: Rscript prep.R   (reads the 843 MB activity table; run it under PBS)
## =============================================================================
source("/projects/group/bliu/e1327989/GTEx/tissue_regression/trait_regression_20261007/code/data.R")
setDTthreads(8)

TIS  <- readLines(PATHS$tissues); TIS <- TIS[nzchar(TIS) & !grepl("^#", TIS)]
meta <- read_meta(); saveRDS(meta, file.path(PATHS$work, "meta.rds"))
info <- fread(PATHS$info, data.table = FALSE); rownames(info) <- info[[1]]
att  <- fread(PATHS$attr, select = c("SAMPID", "SMAFRZE", "SMRIN"), quote = "")
T11  <- fread(PATHS$tpm11)                                   # GTEx v11 TPM >= 1 gene table

cat("reading activity table...\n")
act <- fread(PATHS$act, data.table = FALSE); names(act)[1] <- "row_id"
act$frz <- att$SMAFRZE[match(act$row_id, att$SAMPID)]
act$rin <- att$SMRIN[match(act$row_id, att$SAMPID)]
pcols <- grep("^PromoterID_", names(act), value = TRUE)
cat("activity rows", nrow(act), " promoters", length(pcols), "\n")

rows <- list()
for (t in TIS) {
  sub0 <- act[act$SMTSD == t & !is.na(act$SMTSD), c("row_id", "sample", "frz", "rin"), drop = FALSE]
  s <- sub0[sub0$frz %in% "RNASEQ", , drop = FALSE]
  s <- s[order(s$sample, -s$rin, seq_len(nrow(s))), , drop = FALSE]
  s <- s[!duplicated(s$sample) & s$sample %in% rownames(meta), , drop = FALSE]
  g  <- T11[SMTSD == t, Name]                       # GTEx v11 genes with mean TPM >= 1 (versions agree with the promoter annotation; matched exactly by versioned ID)
  pr <- if (length(g)) intersect(pcols, rownames(info)[info$geneId %in% g]) else pcols
  A <- as.matrix(act[match(s$row_id, act$row_id), pr, drop = FALSE]); rownames(A) <- s$sample
  td <- list(tissue = t, U = s$sample, row_id = s$row_id, act = A,
             info = info[pr, c("seqnames", "start", "geneId", "promoterPosition"), drop = FALSE])
  saveRDS(td, file.path(PATHS$work, paste0(t, ".rds")))
  rows[[t]] <- data.table(tissue = t, rows_all = nrow(sub0), rows_excluded = sum(!sub0$frz %in% "RNASEQ"),
                          donors = nrow(s), promoters = length(pr), tpm_filter = if (length(g)) "GTEx_v11_TPM1_genes" else "none")
  cat(t, nrow(s), "donors", length(pr), "promoters\n")
}
TT <- rbindlist(rows); fwrite(TT, file.path(PATHS$work, "tissues.tsv"), sep = "\t")
print(TT)
