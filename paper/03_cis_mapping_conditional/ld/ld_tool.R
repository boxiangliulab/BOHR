## LD standard shared by Figures 4, 5 and 6 = PLINK 1.9 (the PLINK approach is the unified standard).
##   genotypes : the model's VCF (GTEx v9 WGS, MAF >= 0.05); GT -> ALT-allele dosage 0/1/2; missing calls stay MISSING
##   per pair  : donors missing either SNP are dropped from that pair (PLINK 1.9 --r / --r2). Verified identical to PLINK
##               (work/plink_v3: 477,984 shared proxy pairs, signs 100% equal with --keep-allele-order)
##   donors    : the tissue's freeze donors (donors/<tissue>.txt = the model's td$U); "ALL" = the union over the 50
##               tissues (943 donors), only for statements not tied to one tissue
##   r         : Pearson correlation of ALT dosages (signed; PLINK --r with --keep-allele-order); r2 = r^2
##   proxy     : r2 >= 0.8; independent: r2 < 0.2
##   EUR       : optional sensitivity for GWAS-facing statements only (RACE == 3; Figure 3), donors/<set>.EUR.txt
## Reference PLINK run (NSCC, PLINK v1.90b7.1): bfiles from the VCF with --double-id --keep-allele-order --vcf-half-call m;
##   then --keep <donors> --ld-snp-list <signal SNPs> --r --keep-allele-order --ld-window 99999 --ld-window-kb 1000
##   --ld-window-r2 0.8 (work/plink_v3). The earlier files in Rawdata/v9_vcf_MAF_0.05/LD used the default --ld-window 10.
## NOT LD: ld_dosage_modal() returns the model's genotypes (missing -> most common genotype), for phenotype summaries only.
## Library: source("code/ld_tool.R"); ld_pairs(data.table(tissue, snp_a, snp_b), eur = FALSE)
## CLI    : Rscript code/ld_tool.R <in.tsv: tissue, snp_a, snp_b> <out.tsv> [EUR]   (run on PBS, not on the login node)
##          Output adds r, r2, n_pair (donors called at both SNPs), n_donors, in_vcf_a, in_vcf_b.
suppressPackageStartupMessages(library(data.table))
if (!exists("read_geno")) source("/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/code/bohr2_data.R")
LD_DIR <- file.path(V2, "results/ld_v3")
ld_donors <- function(tissue, eur = FALSE) readLines(file.path(LD_DIR, "donors", paste0(tissue, if (eur) ".EUR" else "", ".txt")))
ld_pos <- function(id) as.integer(vapply(strsplit(id, "_", fixed = TRUE), `[`, "", 2))
ld_chr <- function(id) vapply(strsplit(id, "_", fixed = TRUE), `[`, "", 1)
## raw dosages (donors x SNPs, missing calls NA) for a region; record selection = readVcf overlap, as read_geno
ld_raw_region <- function(chr, from, to, donors) {
  prm <- ScanVcfParam(which = GRanges(chr, IRanges(max(1L, from), to)), samples = donors, fixed = NA, info = NA, geno = "GT")
  v <- readVcf(TabixFile(file.path(PATHS$vcf, paste0(chr, ".vcf.gz"))), "hg38", prm)
  if (nrow(v) == 0) return(NULL)
  t(gt_dosage(geno(v)$GT[, donors, drop = FALSE])) }
## raw dosages for a set of SNP IDs (one tabix query per SNP); IDs absent from the VCF are dropped
ld_raw_dosage <- function(ids, donors) {
  ids <- unique(ids); out <- list()
  for (id in ids) { D <- tryCatch(ld_raw_region(ld_chr(id), ld_pos(id), ld_pos(id), donors), error = function(e) NULL)
    if (!is.null(D) && id %in% colnames(D)) out[[id]] <- D[, id] }
  if (!length(out)) return(matrix(numeric(0), length(donors), 0, dimnames = list(donors, NULL)))
  do.call(cbind, out) }
## PLINK's per-pair rule: donors missing either SNP are dropped (x a vector; y a vector or a donors x SNPs matrix)
ld_r <- function(x, y) suppressWarnings(drop(cor(x, y, use = "pairwise.complete.obs")))
## the model's genotypes (missing -> most common genotype): for phenotype-by-genotype summaries, NOT for LD
ld_dosage_modal <- function(ids, donors) {
  ids <- unique(ids); out <- list()
  for (id in ids) { g <- tryCatch(read_geno(ld_chr(id), ld_pos(id), ld_pos(id), donors), error = function(e) NULL)
    if (!is.null(g) && id %in% colnames(g$G)) out[[id]] <- g$G[, id] }
  if (!length(out)) return(matrix(numeric(0), length(donors), 0, dimnames = list(donors, NULL)))
  do.call(cbind, out) }
ld_dosage <- ld_dosage_modal                      # old name (fig5_joint_moments.R); modal-filled, NOT for LD
ld_pairs <- function(P, eur = FALSE) {
  P <- as.data.table(P)[, row_id := .I]; res <- list()
  for (t in unique(P$tissue)) { x <- P[tissue == t]; d <- ld_donors(t, eur); D <- ld_raw_dosage(c(x$snp_a, x$snp_b), d)
    x[, `:=`(in_vcf_a = snp_a %in% colnames(D), in_vcf_b = snp_b %in% colnames(D), n_donors = length(d))]
    ok2 <- function(a, b) a %in% colnames(D) && b %in% colnames(D)
    x[, r := as.numeric(mapply(function(a, b) if (ok2(a, b)) ld_r(D[, a], D[, b]) else NA_real_, snp_a, snp_b))]
    x[, n_pair := as.integer(mapply(function(a, b) if (ok2(a, b)) sum(!is.na(D[, a]) & !is.na(D[, b])) else NA_integer_, snp_a, snp_b))]
    res[[t]] <- x }
  R <- rbindlist(res)[, r2 := r^2]; setorder(R, row_id); R[, row_id := NULL][] }
if (sys.nframe() == 0L) {
  a <- commandArgs(TRUE); if (!length(a) %in% 2:3) stop("usage: Rscript code/ld_tool.R <in.tsv> <out.tsv> [EUR]")
  fwrite(ld_pairs(fread(a[1]), eur = length(a) == 3 && a[3] == "EUR"), a[2], sep = "\t") }
