## Authored for the repository (not an original script). Per-tissue cis-mapping promoter sets by the paper's rule. On the
## paper's data it reproduces all 131,361 tissue x promoter inputs of cis mapping, 50 of 50 tissues exactly.
## Rule = the old GAMLSS scan's selection, which prep_all.R took over through the published result files:
##   samples : every sample of the tissue (not one per donor) whose donor has genotype PCs (v10_cov_pca.txt, PC1-5) and,
##             when k > 0, PEER factors (peer/v11/v11_filter_5_peer_factors_<T>.txt)
##   filter  : gamlss_QTL_optimized.R lines 115-137: SD >= 0.005, NA share < 0.4, gene mean TPM >= 1 in the tissue
##             (TPM_v11/TPM1_genes.txt; no TPM filter when the tissue has no TPM rows)
##   skip    : orig_QTL_regression (gamlss_functions_optimized.R lines 548-557) on the samples with a value:
##             n >= 20 and n >= 5 + k, share of 0 < 0.95, share of 1 < 0.95, SD > 0.05
##   chr1-22 only.
## Usage: Rscript promoter_sets_rule.R <gtex_root> <activity .txt or .RData> <peer_dir> <out.tsv>
##   activity: the 0103 relative-activity matrix (rows = samples; columns sample (donor), SMTSD, PromoterID_*)
## Tissues: env TISSUES = a file with one tissue per line, or a folder of <tissue>.rds input objects (default: the
##   paper's 50 cis-mapping tissues, the lab copy of work/final/freeze; ../../07_coloc_gwas_alphagenome/config/tissues_50.txt
##   lists the same 50).
suppressMessages(library(data.table))
a <- commandArgs(TRUE); G <- a[1]; ACT <- a[2]; PEER <- a[3]; OUT <- a[4]
best <- fread(file.path(G, "QTL/tissue_list_best_peer_table.txt"), header = FALSE, fill = TRUE, col.names = c("tissue", "N", "peer"))
best[, peer := as.integer(gsub("[^0-9]", "", peer))]
TF <- Sys.getenv("TISSUES", "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/work/final/freeze")
TIS <- if (dir.exists(TF)) sub("\\.rds$", "", list.files(TF, "\\.rds$")) else readLines(TF)   # the 50 tissues
mc <- read.table(file.path(G, "v10_cov_pca.txt"), sep = "\t", header = TRUE, row.names = 1)[20:24]
if (grepl("\\.RData$", ACT)) { e <- new.env(); load(ACT, envir = e); A <- as.data.frame(get("recal_att_f", envir = e)) } else {
  A <- fread(ACT, data.table = FALSE); rownames(A) <- A[[1]] }
info <- local({ e <- new.env(); f <- file.path(dirname(ACT), "Promoter_info_filter_5.RData"); load(f, envir = e); x <- as.data.frame(get("info", envir = e))
  if (!is.null(x$promoterId)) rownames(x) <- paste0("PromoterID_", x$promoterId); x })
tpm <- fread(file.path(G, "TPM_v11/TPM1_genes.txt"))
pc <- grep("^PromoterID_", names(A), value = TRUE); pc <- pc[as.character(info[pc, "seqnames"]) %in% paste0("chr", 1:22)]
res <- list()
for (t in TIS) { k <- best[tissue == t, peer]; k <- if (length(k) && !is.na(k)) k else 0L
  donors <- rownames(mc); if (k > 0) donors <- intersect(donors, rownames(read.table(file.path(PEER, paste0("v11_filter_5_peer_factors_", t, ".txt")), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)))
  d <- A[A$SMTSD == t & A$sample %in% donors, pc, drop = FALSE]
  sdv <- vapply(d, function(x) sd(x[!is.na(x)]), 0); nas <- colMeans(is.na(d))
  keep <- pc[!is.na(sdv) & sdv >= 0.005 & nas < 0.4]
  g <- tpm[SMTSD == t, Name]; if (length(g)) keep <- keep[info[keep, "geneId"] %in% g]
  sk <- vapply(keep, function(p) { x <- d[[p]]; x <- x[!is.na(x)]; n <- length(x)
    n < 20 || n < 5 + k || mean(x == 0) >= 0.95 || mean(x == 1) >= 0.95 || sd(x) <= 0.05 }, TRUE)
  res[[t]] <- data.table(tissue = t, promoter = keep[!sk]); cat(t, length(keep[!sk]), "\n") }
R <- rbindlist(res); fwrite(R, OUT, sep = "\t"); cat("tissue x promoter:", nrow(R), "; promoters:", uniqueN(R$promoter), "\n")
