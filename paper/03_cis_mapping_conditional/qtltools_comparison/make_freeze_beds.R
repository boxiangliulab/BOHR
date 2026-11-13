## All-tissue BOHR vs QTLtools comparison (author approved 2026-10-01, deferred; started 2026-10-02).
## QTLtools inputs are built from BOHR's own freeze tissue objects (rerun_v2/work/final/freeze/<T>.rds), so the donors,
## the author's filtered promoters, the covariates and the perm0 permutation are exactly BOHR's:
##   real : phenotype = td$act[U, ]; covariates = genotype PC1-5 + PEER 1..K
##   perm0: slot s carries the phenotype AND PEER row of donor map[s]; genotype PCs stay with the slot (= run_final2.R)
## Phenotype rows: the tissue's promoters on chr1-22; NA mean-imputed and constant rows dropped (QTLtools needs complete
## rows; same as build_pheno_cov.py). BED coordinates as build_pheno_cov.py: start = TSS - 1, end = TSS, strand from
## Promoter_info_filter_5.txt. Two BEDs: raw (run with and without --normal) and rint (z = qnorm((rank - 0.5) / n), ties
## averaged; the GTEx order, QTLtools then regresses the covariates).
## Output: QTLtools_linear_compare/freeze_all/inputs/<T>.<arm>.{raw,rint}.bed.gz(.tbi) and <T>.<arm>.cov.txt
## Usage: ARMS="perm1 perm2" Rscript make_freeze_beds.R [tissue ...]   (default ARMS = "real perm0"; permK uses perm_maps/final/freeze/<T>.permK.perm_map.tsv)
suppressPackageStartupMessages(library(data.table))
V2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"; QLC <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/QTLtools_linear_compare"
OUT <- file.path(QLC, "freeze_all/inputs"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
HTS <- "/home/e1327989/.local/share/mamba/envs/htsenv/bin"
INFO <- fread("/home/e1327989/project/GTEx/v11/Promoter_info_filter_5.txt"); setnames(INFO, 1, "pid")
args <- commandArgs(TRUE); TIS <- if (length(args)) args else sub("[.]rds$", "", list.files(file.path(V2, "work/final/freeze"), pattern = "[.]rds$"))
rint <- function(x) qnorm((rank(x, ties.method = "average") - 0.5) / length(x))
write_bed <- function(P, coord, f) {                     # P: promoters x samples
  b <- data.table(`#Chr` = coord$seqnames, start = coord$start - 1L, end = coord$start, pid = rownames(P), gid = rownames(P), strand = coord$strand)
  b <- cbind(b, as.data.table(round(P, 6))); setorder(b, `#Chr`, start)
  fwrite(b, f, sep = "\t", quote = FALSE)
  stopifnot(system2(file.path(HTS, "bgzip"), c("-f", f)) == 0, system2(file.path(HTS, "tabix"), c("-f", "-p", "bed", paste0(f, ".gz"))) == 0) }
for (t in TIS) {
  td <- readRDS(file.path(V2, "work/final/freeze", paste0(t, ".rds"))); U <- td$U; K <- td$K
  for (arm in strsplit(Sys.getenv("ARMS", "real perm0"), " ")[[1]]) {
    if (arm == "real") { A <- td$act[U, , drop = FALSE]; Pe <- td$peer[U, , drop = FALSE] } else {
      m <- fread(file.path(V2, "perm_maps/final/freeze", paste0(t, ".", arm, ".perm_map.tsv")))
      stopifnot(setequal(m$slot, U), setequal(m$donor, U), !anyDuplicated(m$donor))
      d <- setNames(m$donor, m$slot)[U]; A <- td$act[d, , drop = FALSE]; Pe <- td$peer[d, , drop = FALSE]; rownames(A) <- rownames(Pe) <- U }
    A <- A[, td$info[colnames(A), "seqnames"] %in% paste0("chr", 1:22), drop = FALSE]
    na_share <- mean(is.na(A))
    A <- apply(A, 2, function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x })
    A <- A[, apply(A, 2, function(x) all(is.finite(x)) && length(unique(x)) > 1), drop = FALSE]
    coord <- data.table(seqnames = td$info[colnames(A), "seqnames"], start = as.integer(td$info[colnames(A), "start"]),
                        strand = INFO$strand[match(colnames(A), INFO$pid)])
    stopifnot(!anyNA(coord$strand))
    write_bed(t(A), coord, file.path(OUT, paste(t, arm, "raw.bed", sep = ".")))
    R <- apply(A, 2, rint); dimnames(R) <- dimnames(A)
    write_bed(t(R), coord, file.path(OUT, paste(t, arm, "rint.bed", sep = ".")))
    C <- cbind(td$pcs[U, , drop = FALSE], Pe[, seq_len(K), drop = FALSE]); colnames(C) <- c(paste0("PC", 1:ncol(td$pcs)), if (K > 0) paste0("peer_", seq_len(K)))
    stopifnot(!anyNA(C))
    cv <- data.table(id = colnames(C), t(C)); setnames(cv, c("id", U))
    fwrite(cv, file.path(OUT, paste(t, arm, "cov.txt", sep = ".")), sep = "\t", quote = FALSE)
    cat(t, arm, ": samples", length(U), "; promoters", ncol(A), "; covariates", ncol(C), "; NA share imputed", round(na_share, 3), "\n") } }
