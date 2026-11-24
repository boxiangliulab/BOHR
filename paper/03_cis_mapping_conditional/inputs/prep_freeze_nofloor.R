## Cis-mapping inputs (freeze definition) of the final run, on all junctions. Copy of prep_all.R (this folder)
## with two changes, both needed because no earlier run exists for these data:
##  1. per-tissue promoters = promoter_sets_rule.R output (the rule that reproduces all 131,361 inputs of the earlier run),
##     instead of the result files of the old GAMLSS scan (prep_all.R lines 21-25);
##  2. freeze definition only (the pilot definition fed the pilot comparison and the QTLtools pilot maps only).
## Everything else (donor choice, covariates, K = C3, object layout) is prep_all.R's code. Permutation maps: the earlier
## run's perm0-2 maps are copied (donor sets are checked to be identical).
## Env: BOHR2_GTEX = GTEx tree with the final run's v11/ and peer/v11/ files, BOHR2_ROOT = run root.
## Usage: Rscript prep_freeze_nofloor.R <promoter_sets.tsv> <stored run root for perm maps>
suppressPackageStartupMessages({ library(data.table); library(parallel) })
CODE <- Sys.getenv("BOHR2_CODE"); source(file.path(CODE, "bohr2_data.R"))
a <- commandArgs(TRUE); SETS <- fread(a[1]); STORED <- a[2]
setDTthreads(8)
OUTW <- file.path(V2, "work", "final", "freeze"); dir.create(OUTW, recursive = TRUE, showWarnings = FALSE)
OUTM <- file.path(V2, "perm_maps", "final", "freeze"); dir.create(OUTM, recursive = TRUE, showWarnings = FALSE)
TIS <- sort(unique(SETS$tissue))
best <- fread(file.path(GTEX, "QTL/tissue_list_best_peer_table.txt"), header = FALSE, fill = TRUE, col.names = c("tissue", "N", "peer"))
best[, peer := as.integer(gsub("[^0-9]", "", peer))]
## NF_PARTIAL=1: a subset of tissues (the promoter-set file holds only those); NF_APPEND=1 also keeps the other tissues'
## rows of an existing tissues.tsv (second batch)
PARTIAL <- nzchar(Sys.getenv("NF_PARTIAL")); APPEND <- nzchar(Sys.getenv("NF_APPEND"))
stopifnot(length(TIS) == 50 || PARTIAL, all(TIS %in% best$tissue))
pub <- split(SETS$promoter, SETS$tissue)
cat("promoters per tissue: median", median(lengths(pub)), " range", range(lengths(pub)), " total", nrow(SETS), "\n")

pcs <- as.matrix(read.table(PATHS$pcs, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)[, paste0("PC", 1:N_PC)])
info <- fread(PATHS$info, data.table = FALSE); rownames(info) <- info[[1]]
hdr <- names(fread(PATHS$act, nrows = 0)); allp <- intersect(hdr, unique(unlist(pub)))
stopifnot(length(allp) == length(unique(unlist(pub))))
act <- fread(PATHS$act, select = c(hdr[1], "sample", "SMTSD", allp), data.table = FALSE); names(act)[1] <- "row_id"
att <- fread("/home/e1327989/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt", select = c("SAMPID", "SMAFRZE", "SMRIN"), quote = "")
act$frz <- att$SMAFRZE[match(act$row_id, att$SAMPID)]; act$rin <- att$SMRIN[match(act$row_id, att$SAMPID)]
vcf_s <- samples(scanVcfHeader(file.path(PATHS$vcf, "chr1.vcf.gz")))
cat("activity rows", nrow(act), " promoter columns", length(allp), "\n")

rows <- list()
for (i in seq_along(TIS)) { t <- TIS[i]; k <- best[tissue == t, peer]; def <- "freeze"
  peer <- as.matrix(read.table(file.path(PATHS$peer, paste0("v11_filter_5_peer_factors_", t, ".txt")), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
  sub0 <- act[act$SMTSD == t, , drop = FALSE]
  s <- sub0[sub0$frz %in% "RNASEQ", , drop = FALSE]; s <- s[order(s$sample, -s$rin), , drop = FALSE]
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
  saveRDS(td, file.path(OUTW, paste0(t, ".rds")))
  ## the earlier run's maps, when the donor set is the same
  for (arm in c("perm0", "perm1", "perm2")) { src <- file.path(STORED, "perm_maps/final/freeze", paste0(t, ".", arm, ".perm_map.tsv"))
    m <- fread(src); if (!setequal(m$slot, U)) stop("donor set differs from the paper's map: ", t, " ", arm)
    dst <- file.path(OUTM, basename(src)); if (!file.exists(dst)) file.copy(src, dst) }
  rows[[length(rows) + 1]] <- data.table(tissue = t, def = def, n = n, k_now = k, K_C3 = K3, promoters = ncol(A))
  cat(t, "done: n", n, " K", K3, " promoters", ncol(A), "\n")
}
TT <- rbindlist(rows); ftt <- file.path(V2, "work", "final", "tissues.tsv")
if (APPEND && file.exists(ftt)) { TT0 <- fread(ftt); TT <- rbind(TT0[!tissue %in% TT$tissue], TT); setorder(TT, tissue) }
fwrite(TT, ftt, sep = "\t")
OLD <- fread(file.path(STORED, "work/final/tissues.tsv"))[def == "freeze"]
cmp <- merge(TT, OLD[, .(tissue, n_paper = n, K_paper = K_C3, promoters_paper = promoters)], by = "tissue")
cat("donors identical to the paper in", cmp[n == n_paper, .N], "of", nrow(cmp), "tissues; K identical in", cmp[K_C3 == K_paper, .N], "\n")
print(cmp[, .(tissue, n, n_paper, K_C3, K_paper, promoters, promoters_paper)])
