## Full rerun, step 0: cached inputs for all 50 tissues under two sample definitions.
##   freeze : GTEx v11 analysis freeze (SMAFRZE == "RNASEQ"); one row per donor (highest SMRIN, then first)
##   pilot  : first activity row per donor (the pilot's definition; used for like-for-like comparisons
##            with the old results and the QTLtools benchmark)
## Both: donors with complete PC1-5 and PEER 1..k_now, present in the VCF. Promoter set per tissue =
## promoters of the published run (QTL/new_v11_LRT names, chr1-22), as in the pilot.
## K = C3 = max(0, min(k_now, floor(n/30) - 5)), n = donors of that definition.
## Output: work/final/<def>/<T>.rds, work/final/tissues.tsv, perm_maps/final/<def>/<T>.perm0.perm_map.tsv
suppressPackageStartupMessages({ library(data.table); library(parallel) })
source("/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/code/bohr2_data.R")
setDTthreads(8)
OUTW <- file.path(V2, "work", "final"); for (d in c("freeze", "pilot")) dir.create(file.path(OUTW, d), recursive = TRUE, showWarnings = FALSE)
OUTM <- file.path(V2, "perm_maps", "final"); for (d in c("freeze", "pilot")) dir.create(file.path(OUTM, d), recursive = TRUE, showWarnings = FALSE)
QLCP <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/QTLtools_linear_compare/perm"
TIS <- sort(sub("\\.perm_map\\.tsv$", "", basename(Sys.glob(file.path(QLCP, "*.perm_map.tsv")))))
best <- fread(file.path(GTEX, "QTL/tissue_list_best_peer_table.txt"), header = FALSE, fill = TRUE, col.names = c("tissue", "N", "peer"))
best[, peer := as.integer(gsub("[^0-9]", "", peer))]
stopifnot(length(TIS) == 50, all(TIS %in% best$tissue))

## promoters of the published run, per tissue (chr1-22)
pub <- mclapply(TIS, function(t) { k <- best[tissue == t, peer]
  unique(unlist(lapply(1:22, function(c) { f <- sprintf("%s/QTL/new_v11_LRT/%s/chr%d/sensitivity_test_10k.%s.chr%d.factor_peer_%d_pca_5.results.rds", GTEX, t, c, t, c, k)
    if (!file.exists(f)) return(character(0)); p <- sub(paste0("^", t, "_"), "", names(readRDS(f))); p[grepl("^PromoterID_", p)] }))) }, mc.cores = 8)
names(pub) <- TIS
cat("published promoters per tissue: median", median(lengths(pub)), " range", range(lengths(pub)), "\n")

pcs <- as.matrix(read.table(PATHS$pcs, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)[, paste0("PC", 1:N_PC)])
info <- fread(PATHS$info, data.table = FALSE); rownames(info) <- info[[1]]
hdr <- names(fread(PATHS$act, nrows = 0)); allp <- intersect(hdr, unique(unlist(pub)))
act <- fread(PATHS$act, select = c(hdr[1], "sample", "SMTSD", allp), data.table = FALSE); names(act)[1] <- "row_id"
att <- fread("/home/e1327989/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt", select = c("SAMPID", "SMAFRZE", "SMRIN"), quote = "")
act$frz <- att$SMAFRZE[match(act$row_id, att$SAMPID)]; act$rin <- att$SMRIN[match(act$row_id, att$SAMPID)]
vcf_s <- samples(scanVcfHeader(file.path(PATHS$vcf, "chr1.vcf.gz")))
cat("activity rows", nrow(act), " promoter columns", length(allp), "\n")

rows <- list()
for (i in seq_along(TIS)) { t <- TIS[i]; k <- best[tissue == t, peer]
  peer <- as.matrix(read.table(file.path(PATHS$peer, paste0("v11_filter_5_peer_factors_", t, ".txt")), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
  sub0 <- act[act$SMTSD == t, , drop = FALSE]
  for (def in c("freeze", "pilot")) {
    if (def == "freeze") { s <- sub0[sub0$frz %in% "RNASEQ", , drop = FALSE]; s <- s[order(s$sample, -s$rin), , drop = FALSE] }
    else s <- sub0
    s <- s[!duplicated(s$sample), , drop = FALSE]
    ok <- s$sample %in% rownames(pcs) & s$sample %in% rownames(peer) & s$sample %in% vcf_s
    U <- s$sample[ok]
    U <- U[complete.cases(pcs[U, , drop = FALSE]) & complete.cases(peer[U, seq_len(k), drop = FALSE])]
    pr <- intersect(pub[[t]], names(s))
    A <- as.matrix(s[match(U, s$sample), pr, drop = FALSE]); rownames(A) <- U
    n <- length(U); K3 <- max(0L, min(k, floor(n / 30) - N_PC))
    td <- list(tissue = t, k_now = k, K = K3, def = def, U = U, row_id = s$row_id[match(U, s$sample)], act = A,
               pcs = pcs[U, , drop = FALSE], peer = peer[U, , drop = FALSE],
               info = info[colnames(A), c("seqnames", "start", "geneId"), drop = FALSE])
    saveRDS(td, file.path(OUTW, def, paste0(t, ".rds")))
    ## permutation map (perm0): pilot definition reuses the published Fig 2e map; freeze draws a new one
    mf <- file.path(OUTM, def, paste0(t, ".perm0.perm_map.tsv"))
    if (!file.exists(mf)) {
      if (def == "pilot") { m <- fread(file.path(QLCP, paste0(t, ".perm_map.tsv")))
        if (!setequal(m$slot, U)) stop("pilot slot set differs from the rep0 map for ", t); fwrite(m, mf, sep = "\t") }
      else { set.seed(20261001L + i); fwrite(data.table(slot = U, donor = sample(U)), mf, sep = "\t") } }
    rows[[length(rows) + 1]] <- data.table(tissue = t, def = def, n = n, k_now = k, K_C3 = K3, promoters = ncol(A),
      excluded_rows_dropped = if (def == "freeze") sum(!sub0$frz %in% "RNASEQ") else 0L, dup_rows_dropped = nrow(if (def == "freeze") sub0[sub0$frz %in% "RNASEQ", ] else sub0) - nrow(s))
  }
  cat(t, "done\n")
}
TT <- rbindlist(rows); fwrite(TT, file.path(OUTW, "tissues.tsv"), sep = "\t")
print(dcast(TT, tissue + k_now ~ def, value.var = c("n", "K_C3")))
