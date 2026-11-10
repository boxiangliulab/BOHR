## Inputs as R objects or file paths, sample tables and tissue loops, and the TSV writers (author decisions of
## 2026-11-10: no reader/writer machinery beyond this; data.table is used when it is installed).

has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)
is_path <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && is.null(dim(x))
is_vcf_path <- function(x) is_path(x) && grepl("\\.vcf(\\.b?gz)?$", x, ignore.case = TRUE)

## a tab-separated file with a header line; the first column is read as character (IDs), or every column when
## all_character = TRUE (tables). fread() needs R.utils for compressed files, so those go to read.delim().
read_tsv <- function(path, what, all_character = FALSE) {
  if (!file.exists(path)) stop("`", what, "`: file not found: ", path, call. = FALSE)
  if (dir.exists(path)) stop("`", what, "`: ", path, " is a folder, not a file.", call. = FALSE)
  compressed <- grepl("\\.(gz|bz2|xz)$", path, ignore.case = TRUE)
  if (has_pkg("data.table") && (!compressed || (grepl("\\.(gz|bz2)$", path, ignore.case = TRUE) && has_pkg("R.utils")))) {
    d <- data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE, na.strings = c("NA", ""),
                           colClasses = if (all_character) "character" else list(character = 1L), showProgress = FALSE)
  } else {
    hdr <- strsplit(readLines(path, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1]]
    d <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", ""), quote = "", comment.char = "",
                           colClasses = if (all_character) "character" else c("character", rep(NA, length(hdr) - 1L)))
  }
  if (!ncol(d)) stop("`", what, "`: no columns in ", path, ".", call. = FALSE)
  d
}
## a numeric matrix with row IDs in the first column (or the object itself)
read_matrix <- function(x, what) {
  if (!is_path(x)) return(x)
  d <- read_tsv(x, what)
  m <- as.matrix(d[, -1, drop = FALSE])
  if (!is.numeric(m) && !(is.logical(m) && all(is.na(m))))
    stop("`", what, "`: every column after the first (IDs) must be numeric in ", x, ".", call. = FALSE)
  storage.mode(m) <- if (is.integer(m)) "integer" else "double"
  rownames(m) <- as.character(d[[1]])
  m
}
## a table (or the object itself); from a file every column is character, converted by the checks
read_table <- function(x, what) if (is_path(x)) read_tsv(x, what, all_character = TRUE) else x
## covariates of a single-tissue call: an object or a file
read_covariates <- function(x) {
  if (is_path(x) && dir.exists(x)) stop("`covariates` is a folder; give `tissue` to read one file per tissue, or give one file.", call. = FALSE)
  read_matrix(x, "covariates")
}
## covariates of tissue tt in a tissue loop: a list named by tissue, a folder with <tissue>.tsv files, or one table with
## row names (sample IDs) for every sample (cis_one() takes the tissue's rows)
tissue_covariates <- function(covariates, tt) {
  if (is.null(covariates)) return(NULL)
  if (is.list(covariates) && !is.data.frame(covariates)) {
    if (is.null(names(covariates)) || !tt %in% names(covariates))
      stop("`covariates` has no element for tissue ", tt, "; give a list named by tissue.", call. = FALSE)
    return(read_matrix(covariates[[tt]], paste0("covariates$", tt)))
  }
  if (is_path(covariates)) {
    f <- file.path(covariates, paste0(tt, c(".tsv", ".txt", ".tsv.gz", ".txt.gz")))
    f <- f[file.exists(f)]
    if (!length(f)) stop("No covariate file for tissue ", tt, " in ", covariates, " (expected ", tt, ".tsv).", call. = FALSE)
    return(read_matrix(f[1], paste0("covariates (", tt, ")")))
  }
  if (is.null(rownames(covariates)))
    stop("With `tissue`, `covariates` must be a list named by tissue, a folder with one file per tissue, or a table with row names (sample IDs).",
         call. = FALSE)
  covariates
}

## tissue labels and the genotype column (individual) of each sample, from vectors or a sample table
align_samples <- function(v, samples, what) {
  if (is.list(v) || !is.null(dim(v))) stop("`", what, "` must be a vector, one value per sample.", call. = FALSE)
  if (!is.null(names(v)) && all(samples %in% names(v))) v <- v[samples]
  else if (length(v) != length(samples)) stop("`", what, "` must have one value per column of `phenotype` (or be named by sample ID).", call. = FALSE)
  v <- as.character(unname(v))
  if (anyNA(v) || any(!nzchar(v))) stop("`", what, "` has missing values.", call. = FALSE)
  v
}
sample_info <- function(tissue, individual, samples) {
  lev <- NULL
  if (is_path(tissue) || is.data.frame(tissue)) {
    tab <- read_table(tissue, "tissue")
    if (!all(c("sample", "tissue") %in% names(tab)))
      stop("A sample table given as `tissue` needs the columns sample and tissue (and optionally individual).", call. = FALSE)
    i <- match(samples, as.character(tab$sample))
    if (anyNA(i)) stop("The sample table has no row for ", sum(is.na(i)), " sample(s) of `phenotype`, e.g. ", samples[is.na(i)][1], ".", call. = FALSE)
    if (is.factor(tab$tissue)) lev <- levels(droplevels(tab$tissue[i]))
    tissue <- stats::setNames(as.character(tab$tissue[i]), samples)
    if (is.null(individual) && "individual" %in% names(tab)) individual <- stats::setNames(as.character(tab$individual[i]), samples)
  } else if (is.factor(tissue)) lev <- levels(droplevels(tissue))
  if (!is.null(tissue)) {
    tissue <- align_samples(tissue, samples, "tissue")
    lev <- if (is.null(lev)) sort(unique(tissue), method = "radix") else intersect(lev, tissue)
  }
  individual <- if (is.null(individual)) samples else align_samples(individual, samples, "individual")
  list(tissue = tissue, individual = individual, levels = lev)
}

## bind per-tissue bohr_cis results with a first column `tissue`
bind_tissues <- function(res, levels) {
  tag <- function(d, tt) if (is.null(d)) NULL else data.frame(tissue = rep(tt, nrow(d)), d, check.names = FALSE, stringsAsFactors = FALSE)
  tab <- function(k) { b <- bind_df(Map(function(r, tt) tag(r[[k]], tt), res, levels)); if (!is.null(b)) rownames(b) <- NULL; b }
  s <- res[[1]]$settings
  s$tissues <- levels
  s$chromosomes <- sort(unique(unlist(lapply(res, function(r) r$settings$chromosomes))))
  s$n_samples <- stats::setNames(vapply(res, function(r) r$settings$n_samples, 0L), levels)
  s$n_covariates <- stats::setNames(vapply(res, function(r) r$settings$n_covariates, 0L), levels)
  structure(list(snp = tab("snp"), units = tab("units"), promoters = tab("promoters"), mirror = tab("mirror"),
                 threshold = stats::setNames(vapply(res, function(r) r$threshold, 0), levels), settings = s,
                 seconds = sum(vapply(res, function(r) r$seconds, 0))), class = "bohr_cis")
}

## ---- writers ------------------------------------------------------------------------------------------
write_tsv <- function(d, path) {
  if (has_pkg("data.table")) data.table::fwrite(d, path, sep = "\t", quote = FALSE, na = "NA", showProgress = FALSE)
  else utils::write.table(d, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}
write_cis_tables <- function(out, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out_dir)) stop("Cannot create `out_dir`: ", out_dir, call. = FALSE)
  for (k in c("snp", "units", "promoters", "mirror")) if (!is.null(out[[k]]) && nrow(out[[k]])) write_tsv(out[[k]], file.path(out_dir, paste0(k, ".tsv")))
  invisible(out_dir)
}

#' Read allele dosages from a VCF file
#'
#' Reads the genotype calls (`GT`) of a VCF file with the Bioconductor package
#' VariantAnnotation and turns them into allele dosages: the number of
#' alternative alleles, 0, 1 or 2, as in the BOHR paper's pipeline. The calls
#' `0/0`, `0/1`, `1/0`, `1/1` and their phased forms (`0|1` and so on) are
#' coded; every other call (missing `./.`, multi-allelic, haploid) becomes
#' `NA`, which [bohr_cis()] fills with the variant's most common dosage among a
#' tissue's samples. [bohr_cis()] calls this function when `genotype` is a
#' path ending in `.vcf` or `.vcf.gz`.
#'
#' The whole file is read; to save memory, give one VCF per chromosome or
#' restrict `samples`. Variant IDs are the VCF's `ID` column, or
#' `chrom:pos_ref/alt` when it is missing (VariantAnnotation's rule).
#'
#' @param file Path of a VCF file (plain or bgzip-compressed).
#' @param samples Optional character vector of the samples (VCF columns) to
#'   read; by default all.
#' @return A list with `genotype`, an integer matrix of dosages (variants in
#'   rows, samples in columns), and `variants`, a data frame with `id`, `chr`,
#'   `pos` and `end` (`pos` + length of the reference allele - 1), the
#'   formats [bohr_cis()] takes.
#' @seealso [bohr_cis()], [simulate_bohr()] (which writes a small VCF).
#' @examples
#' if (requireNamespace("VariantAnnotation", quietly = TRUE)) {
#'   d <- file.path(tempdir(), "bohr_vcf_example")
#'   sim <- simulate_bohr(n_individuals = 30, n_tissues = 1, n_genes = 2,
#'                        variant_spacing = 2e5, out_dir = d)
#'   g <- read_vcf_dosage(file.path(d, "genotypes.vcf"))
#'   g$genotype[1:3, 1:5]
#'   head(g$variants)
#'   unlink(d, recursive = TRUE)
#' }
#' @export
read_vcf_dosage <- function(file, samples = NULL) {
  if (!is_path(file)) stop("`file` must be the path of one VCF file.", call. = FALSE)
  if (!file.exists(file)) stop("VCF file not found: ", file, call. = FALSE)
  if (!has_pkg("VariantAnnotation") || !has_pkg("SummarizedExperiment"))
    stop("Reading a VCF file needs the Bioconductor packages VariantAnnotation and SummarizedExperiment: install them with ",
         "BiocManager::install(\"VariantAnnotation\"), or pass a dosage matrix (variants x samples) as `genotype` and the ",
         "variant positions as `variants`.", call. = FALSE)
  if (!is.null(samples)) {
    samples <- unique(as.character(samples))
    have <- VariantAnnotation::samples(VariantAnnotation::scanVcfHeader(file))
    miss <- setdiff(samples, have)
    if (length(miss)) stop("The VCF file has no column for ", length(miss), " sample(s), e.g. ", miss[1], ".", call. = FALSE)
  }
  param <- VariantAnnotation::ScanVcfParam(info = NA, geno = "GT", samples = if (is.null(samples)) character() else samples)
  vcf <- VariantAnnotation::readVcf(file, param = param)
  gt <- VariantAnnotation::geno(vcf)$GT
  rr <- as.data.frame(SummarizedExperiment::rowRanges(vcf))
  ids <- rownames(gt)
  if (is.null(ids)) ids <- paste0(rr$seqnames, ":", rr$start)
  if (anyDuplicated(ids)) {
    dup <- duplicated(ids) | duplicated(ids, fromLast = TRUE)
    warning(sum(dup), " records share a variant ID; they are renamed ID_1, ID_2, ... in file order.", call. = FALSE)
    ids[dup] <- stats::ave(ids[dup], ids[dup], FUN = function(v) paste0(v, "_", seq_along(v)))
  }
  D <- gt_dosage(gt)
  dimnames(D) <- list(ids, colnames(gt))
  list(genotype = D, variants = data.frame(id = ids, chr = as.character(rr$seqnames), pos = as.numeric(rr$start),
                                           end = as.numeric(rr$end), stringsAsFactors = FALSE))
}
