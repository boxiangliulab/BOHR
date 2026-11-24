## depth_extract.R -- raw junction counts behind the relative activities (gene read depth of the final run).
## Without a read floor a zero is 0 reads, so a gene eQTL can move zeros without a change in usage.
## Per tissue of the v4 set: the raw proActiv promoter counts (assay 1 of the per-tissue proActiv object, the input of
## 0103_proActiv_data_process_tissue.R) of the tissue's v4-set promoters, and per gene the total raw count over ALL of its
## non-internal promoters (the denominator of the relative activity in 0103), for the scan samples (freeze row_id).
## Checks per tissue: activity 0 <=> promoter count 0, and activity NA <=> gene count 0 (0103's ratio of activities).
## Output (private: sample-level counts, no genotypes): $N/depth/<tissue>.rds = list(tissue, row_id, R = samples x
## promoters, D = samples x genes, gene_of = named promoter -> gene), and $N/depth/extract_checks.tsv.
## Usage: Rscript depth_extract.R <ncores> [tissue ...]   (NSCC compute node)
suppressPackageStartupMessages({ library(SummarizedExperiment); library(data.table); library(parallel) })
N <- "/scratch/users/nus/e1327989/nofloor_full"; NC <- as.integer(commandArgs(TRUE)[1])
PA <- "/scratch/users/nus/e1327989/junction_floor/nofloor/proActiv"; OUT <- file.path(N, "depth"); dir.create(OUT, showWarnings = FALSE)
e <- new.env(); load("/scratch/users/nus/e1327989/junction_floor/nofloor/Promoter_info_all.RData", envir = e)
IA <- as.data.table(get(ls(e)[1], envir = e)); IA[, pid := paste0("PromoterID_", promoterId)]
NI <- IA[internalPromoter %in% FALSE, .(pid, geneId)]
fz <- list.files(file.path(N, "rerun_v4/work/final/freeze"), "[.]rds$", full.names = TRUE)
if (length(commandArgs(TRUE)) > 1) fz <- fz[sub("[.]rds$", "", basename(fz)) %in% commandArgs(TRUE)[-1]]   # a subset, for a test
pf <- list.files(PA, "[.]rds$", full.names = TRUE)
key <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
pkey <- key(sub("_n[0-9]+[.]rds$", "", sub("^proActiv_", "", basename(pf))))
clean <- function(x) { x <- gsub("\\.SJ\\.out.*$", "", x); gsub("\\.", "-", x) }
one <- function(f) {
  x <- readRDS(f); t <- x$tissue; j <- which(pkey == key(t))
  if (length(j) != 1) return(data.table(tissue = t, status = paste("proActiv files matched:", length(j))))
  se <- readRDS(pf[j]); C <- assays(se)[[1]]; C <- as.matrix(C); storage.mode(C) <- "double"
  rownames(C) <- paste0("PromoterID_", rowData(se)$promoterId); colnames(C) <- clean(colnames(C)); rm(se); gc()
  rid <- x$row_id; miss <- sum(!rid %in% colnames(C)); if (miss) return(data.table(tissue = t, status = paste("samples missing:", miss)))
  P <- colnames(x$act); gene_of <- setNames(IA$geneId[match(P, IA$pid)], P)
  R <- t(C[P, rid, drop = FALSE])
  ng <- NI[geneId %in% unique(gene_of)]; G <- sort(unique(ng$geneId))
  M <- Matrix::sparseMatrix(i = match(ng$pid, rownames(C)), j = match(ng$geneId, G), x = 1, dims = c(nrow(C), length(G)))
  D <- as.matrix(Matrix::crossprod(M, C[, rid, drop = FALSE])); D <- t(D); dimnames(D) <- list(rid, G); rm(C); gc()
  whole <- all(abs(R - round(R)) < 1e-9) && all(abs(D - round(D)) < 1e-9); if (whole) { storage.mode(R) <- "integer"; storage.mode(D) <- "integer" }
  y <- x$act[, P, drop = FALSE]; dg <- D[, gene_of[P], drop = FALSE]
  chk <- data.table(tissue = t, status = "ok", samples = length(rid), promoters = length(P), genes = length(G), integer_counts = whole,
    zero_value_with_reads = sum(!is.na(y) & y == 0 & R > 0), positive_value_without_reads = sum(!is.na(y) & y > 0 & R == 0),
    na_value_with_gene_reads = sum(is.na(y) & dg > 0), value_without_gene_reads = sum(!is.na(y) & dg == 0),
    cells_with_value = sum(!is.na(y)))
  saveRDS(list(tissue = t, row_id = rid, R = R, D = D, gene_of = gene_of), file.path(OUT, paste0(t, ".rds")))
  chk }
res <- rbindlist(mclapply(fz, function(f) tryCatch(one(f), error = function(e) data.table(tissue = basename(f), status = paste("ERROR", conditionMessage(e)))),
                          mc.cores = NC, mc.preschedule = FALSE), fill = TRUE)
fwrite(res, file.path(OUT, "extract_checks.tsv"), sep = "\t"); print(res, nrows = 60)
cat(if (all(res$status == "ok")) "DEPTH EXTRACT DONE" else "DEPTH EXTRACT INCOMPLETE", "\n")
