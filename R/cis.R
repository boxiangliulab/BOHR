#' cis mapping of promoter-usage QTLs
#'
#' Tests every variant within a window around each promoter's transcription
#' start site for association with each component of the promoter's relative
#' activity (mu, p0 and p1; see [bohr_test()]), combines the variant P values
#' of each promoter and component with ACAT, and controls the false-discovery
#' rate with the Benjamini-Hochberg procedure across all tests of a tissue.
#' With `tissue`, it loops over tissues and binds the results.
#'
#' The steps follow the analysis engine of the BOHR paper; each tissue is
#' mapped on its own samples.
#' \enumerate{
#'   \item Covariates are standardised over the tissue's samples and an
#'     intercept is added. Missing genotypes are filled, per variant, with the
#'     most common dosage over the tissue's samples (ties go to the smaller
#'     dosage). Variant filtering (for example by allele frequency) is left to
#'     the user; the paper tested variants with a minor-allele frequency of at
#'     least 0.05 over all 944 GTEx donors. A warning is given when there are
#'     more covariates than one per 30 samples (\eqn{\lfloor n/30 \rfloor});
#'     nothing is dropped.
#'   \item Each promoter is prepared and fitted without variants as in
#'     [bohr_fit()] (skips, endpoint rule, family, orthonormal design).
#'   \item Variants overlapping \[`pos` - `window`, `pos` + `window` - 1\] on
#'     the promoter's chromosome are tested when at least
#'     `min_genotype_count` of the promoter's samples carry each dosage (0, 1
#'     and 2). Each variant is tested with [bohr_test()] as an allele dosage:
#'     mu with the F-corrected likelihood-ratio test, p0 and p1 with the
#'     Firth-penalised likelihood-ratio test.
#'   \item A unit is a promoter and one of its components: mu for every fitted
#'     promoter, p0 when zeros remain after the endpoint rule (families BEINF
#'     and BEZI), p1 when ones remain (BEINF and BEOI). The unit P value
#'     combines the P values of all tested variants in the window with
#'     [acat()]. The leader is the variant with the smallest P value; ties are
#'     broken by variant ID (C-locale order), as in the paper's pipeline.
#'   \item Estimability: a mu unit needs at least `min_n_mid` intermediate
#'     samples; a p0 or p1 unit needs a boundary share (samples at that
#'     boundary over all samples used) within `boundary_share`. These rules
#'     were set for GTEx sample sizes.
#'   \item Exact mirror pairs ([mirror_pairs()], judged in each tissue; needs
#'     `promoters$gene`): for a `boundary_exact` pair the p1 units of both
#'     members are dropped (p1 of one promoter is p0 of the other, which is
#'     kept: each boundary event is tested as p0 at the promoter that goes
#'     silent); for an `interior_exact` pair the mu unit of the partner
#'     (non-canonical member) is dropped.
#'   \item q values: Benjamini-Hochberg over the units of the tissue that are
#'     estimable and not dropped, all components together. Other units have
#'     `q = NA`.
#' }
#' To match the paper's within-tissue false-discovery control, run all
#' chromosomes of a tissue in one call, or run chromosomes separately (for
#' example as cluster jobs) and combine the results with [bohr_combine()].
#' Results do not depend on `ncores` or `chunk_size`.
#'
#' **Inputs as objects or files.** Each of `phenotype`, `genotype`,
#' `variants`, `promoters`, `covariates` and `tissue` is either an R object or
#' the path of a tab-separated file with a header line, read with
#' [data.table::fread()] when data.table is installed and with
#' [utils::read.delim()] otherwise (`NA` or an empty field is missing):
#' \describe{
#'   \item{matrices (`phenotype`, `genotype`, `covariates`)}{The first column
#'     holds the row IDs (promoter, variant or sample IDs); the header holds the
#'     column IDs. Covariate files have samples in rows.}
#'   \item{tables (`variants`, `promoters`, sample table)}{One row per item with
#'     the named columns; other columns are ignored.}
#'   \item{VCF}{A `genotype` path ending in `.vcf` or `.vcf.gz` is read with
#'     [read_vcf_dosage()] (needs the Bioconductor package VariantAnnotation);
#'     the variant positions then come from the VCF and `variants` must be
#'     `NULL`.}
#' }
#' [simulate_bohr()] writes a complete example input set in these formats.
#'
#' **Tissues.** With `tissue`, the samples (columns of `phenotype`) are
#' split by tissue and each tissue is mapped as above, on its own samples,
#' covariates and genotype imputation; q values and `threshold` stay within
#' each tissue. The tables are bound with a first column `tissue`. The
#' canonical member of a mirror pair is chosen by the mean relative activity
#' over all supplied samples pooled across tissues, the rule of the paper (one
#' choice for every tissue); a single-tissue call uses that tissue's mean.
#' Each sample's genotypes are looked up in the column `individual` of
#' `genotype` (by default the sample ID itself), so samples of one donor in
#' several tissues share the donor's genotypes.
#'
#' With `all_stats = TRUE`, the `snp` table also holds every statistic of
#' the analysis engine under its original name: unpenalised likelihood-ratio
#' and Wald log P values (`mu_lrt_logp`, `nu_lrt_logp`, `tau_lrt_logp`,
#' `*_wald_logp`), the joint test, estimates of the unpenalised fits (`nu_est`,
#' `tau_est`), iteration counts and log-likelihoods. Here nu and tau are the
#' p0 and p1 parameters of the model.
#'
#' @param phenotype Relative activities in \[0, 1\], promoters in rows and
#'   samples in columns, with row names (promoter IDs) and unique column names
#'   (sample IDs); a numeric matrix or the path of a TSV file. `NA` marks a
#'   missing value. Should hold all tested promoters, so that two-promoter
#'   genes are recognised.
#' @param genotype Dosages (0, 1, 2 or `NA`), variants in rows (row names =
#'   variant IDs) and individuals in columns: a numeric matrix, the path of a
#'   TSV file, or the path of a VCF file (`.vcf` or `.vcf.gz`; see
#'   [read_vcf_dosage()]). Its columns must include the genotype column of
#'   every sample (see `individual`).
#' @param variants Positions of the rows of `genotype`: a data frame or the
#'   path of a TSV file with columns `chr` and `pos` (start position), and
#'   optionally `id` (matched to the row names of `genotype`; otherwise rows
#'   are taken in the same order) and `end` (end position, default `pos`; for
#'   an indel, `pos` + length of the reference allele - 1). `NULL` when
#'   `genotype` is a VCF file.
#' @param promoters One row per promoter, a data frame or the path of a TSV
#'   file: `id`, `chr`, `pos` (transcription start site) and optionally `gene`
#'   (needed for mirror pairs). Rows for promoters not in `phenotype` (for
#'   example internal or untested promoters of a TSS table) are ignored.
#' @param covariates `NULL`, or covariates with one row per sample and no
#'   intercept column: a numeric matrix or data frame (row names = sample IDs,
#'   or rows in the order of the columns of `phenotype`), or the path of a TSV
#'   file (first column = sample IDs). With `tissue`: a list named by tissue
#'   (each element a matrix, data frame or file path), the path of a folder
#'   with one file per tissue named `<tissue>.tsv` (or `.txt`, optionally
#'   gzipped), or one table whose row names cover the samples of every tissue.
#' @param tissue `NULL` (one tissue), or the tissue of each sample: a vector
#'   with one label per column of `phenotype` (or named by sample ID), or a
#'   sample table (a data frame or the path of a TSV file) with columns
#'   `sample` and `tissue`, and optionally `individual`. See Details.
#' @param individual `NULL`, or the genotype column of each sample: a vector
#'   with one value per column of `phenotype` (or named by sample ID). By
#'   default the `individual` column of a sample table given as `tissue`, or
#'   else the sample ID itself.
#' @param window Half-width of the cis window in base pairs: variants
#'   overlapping \[`pos` - `window`, `pos` + `window` - 1\] are tested. Default
#'   5e5 (+-500 kb; Methods Section 4).
#' @param chromosomes Optional vector of chromosome names; only promoters on
#'   these chromosomes are mapped.
#' @param min_genotype_count A variant is tested only when at least this many
#'   of the promoter's samples carry each dosage (0, 1 and 2). Default 5 (SI
#'   Section 2).
#' @param min_n_mid A mu unit is estimable when the promoter has at least this
#'   many intermediate samples. Default 80 (Methods Section 2).
#' @param boundary_share A p0 (p1) unit is estimable when the share of the
#'   promoter's samples at 0 (at 1) lies in this range, ends included. Default
#'   `c(0.10, 0.90)` (Methods Section 2).
#' @param fdr q-value threshold of the `significant` column and of
#'   `threshold`. Default 0.05 (Methods Section 2).
#' @param acat_cap Cap each variant P value at \eqn{1 - 1/M}, for the \eqn{M}
#'   tested variants of the window, before ACAT (see [acat()]). Default `TRUE`,
#'   as in the paper's analysis code (not stated in the Methods or SI).
#' @param mirror Apply the exact-mirror-pair rule. Default `TRUE` (Methods
#'   Section 2; SI Section 2).
#' @param mirror_tol Tolerance of \eqn{y_1 + y_2 = 1} for `interior_exact`
#'   pairs. Default 1e-6 (SI Section 2).
#' @param major_score Optional numeric vector named by promoter ID that picks
#'   the canonical member of a mirror pair (see [mirror_pairs()]). Defaults to
#'   the mean relative activity in `phenotype`: over all samples pooled across
#'   tissues when `tissue` is given (the paper's rule: all samples of all 50
#'   tissues, SI Section 2), otherwise over the tissue's samples.
#' @param control Fitting rules, optimiser limits and the choice of the p0/p1
#'   test, from [bohr_control()]; the defaults are the values of the BOHR
#'   paper.
#' @param ncores Number of cores (forked processes; 1 on Windows). When BLAS
#'   is multithreaded (e.g. OpenBLAS), setting the environment variable
#'   `OPENBLAS_NUM_THREADS=1` before starting R avoids oversubscription.
#' @param chunk_size Number of variants per parallel task.
#' @param all_stats Logical: keep every statistic of the analysis engine in the
#'   `snp` table (see Details).
#' @param out_dir Optional folder (created if needed). When given, the tables
#'   `snp`, `units`, `promoters` and `mirror` are also written to it as
#'   `snp.tsv`, `units.tsv`, `promoters.tsv` and `mirror.tsv` (tab-separated,
#'   `NA` for missing values); a table without rows is not written. The tables
#'   are returned either way.
#' @return An object of class `"bohr_cis"`, a list with
#'   \describe{
#'     \item{`snp`}{One row per tested promoter-variant pair: `promoter`,
#'       `variant`, `chr`, `pos`, `alt_af` (frequency of the counted allele
#'       among the promoter's samples), `maf`, `min_geno` (smallest dosage
#'       class), and the columns of [bohr_test()].}
#'     \item{`units`}{One row per unit with a finite P value: `promoter`,
#'       `gene`, `component`, `family`, `n`, `n_mid`, `n0`, `n1`,
#'       `boundary_share`, `n_variants`, `p` (ACAT), `leader`, `leader_pos`,
#'       `leader_log_p`, `leader_est`, `leader_effect`, the other components'
#'       log P values and estimates at the leader, `estimable`, the mirror
#'       columns (`n_prom_gene`, `mirror_pair_id`, `boundary_exact`,
#'       `interior_exact`, `mirror_role`, `mirror_partner`), `dedup_drop`,
#'       `duplicate_of`, `q` and `significant`.}
#'     \item{`promoters`}{One row per promoter: family or skip reason, sample
#'       counts, design size, null-fit convergence and log-likelihoods, and
#'       variant counts (`n_variants_window`, `n_variants_tested`).}
#'     \item{`mirror`}{The pairs found by [mirror_pairs()].}
#'     \item{`threshold`}{Largest unit P value with `q <= fdr` (`NA` if
#'       none); with `tissue`, one value per tissue.}
#'     \item{`settings`, `seconds`}{Settings used and elapsed time.}
#'   }
#'   With `tissue`, every table starts with a column `tissue`.
#' @seealso [bohr_test()], [bohr_combine()], [mirror_pairs()], [acat()],
#'   [read_vcf_dosage()], [simulate_bohr()].
#' @examples
#' set.seed(4)
#' n <- 120
#' samples <- paste0("S", 1:n)
#' geno <- matrix(rbinom(30 * n, 2, 0.35), 30, n,
#'                dimnames = list(paste0("v", 1:30), samples))
#' vars <- data.frame(id = rownames(geno), chr = "chr1", pos = seq(1000, 30000, 1000))
#' y1 <- simulate_response(n, x = geno[5, ], mu = 0.4, p0 = 0.3, p1 = 0,
#'                         effect_p0 = 0.15)$y
#' y2 <- simulate_response(n, x = geno[20, ], mu = 0.6, p0 = 0, p1 = 0.2,
#'                         effect_mu = 0.08)$y
#' ph <- rbind(P1 = y1, P2 = y2); colnames(ph) <- samples
#' proms <- data.frame(id = c("P1", "P2"), chr = "chr1", pos = c(5000, 20000),
#'                     gene = c("G1", "G2"))
#' res <- bohr_cis(ph, geno, vars, proms, window = 3000)
#' res
#' res$units[, c("promoter", "component", "p", "q", "leader")]
#'
#' ## files and a tissue loop: a small synthetic input set
#' d <- file.path(tempdir(), "bohr_example")
#' sim <- simulate_bohr(n_individuals = 90, n_tissues = 2, n_genes = 3,
#'                      variant_spacing = 1e5, out_dir = d)
#' res2 <- bohr_cis(file.path(d, "relative_activity.tsv"),
#'                  genotype = file.path(d, "genotypes.tsv"),
#'                  variants = file.path(d, "variants.tsv"),
#'                  promoters = file.path(d, "tss.tsv"),
#'                  covariates = file.path(d, "covariates"),
#'                  tissue = file.path(d, "samples.tsv"), window = 2e5)
#' res2
#' unlink(d, recursive = TRUE)
#' @export
bohr_cis <- function(phenotype, genotype, variants = NULL, promoters, covariates = NULL, tissue = NULL, individual = NULL,
                     window = 5e5, chromosomes = NULL, min_genotype_count = 5L,
                     min_n_mid = 80L, boundary_share = c(0.10, 0.90), fdr = 0.05,
                     acat_cap = TRUE, mirror = TRUE, mirror_tol = 1e-6, major_score = NULL,
                     control = bohr_control(), ncores = 1L, chunk_size = 100L, all_stats = FALSE, out_dir = NULL) {
  t0 <- proc.time()[["elapsed"]]
  check_control(control); check_flag(acat_cap, "acat_cap")
  if (!is.numeric(mirror_tol) || length(mirror_tol) != 1L || !(mirror_tol >= 0)) stop("`mirror_tol` must be a single number >= 0.", call. = FALSE)
  phenotype <- check_phenotype(read_matrix(phenotype, "phenotype"))
  samples <- colnames(phenotype)
  if (!is.numeric(window) || length(window) != 1L || !is.finite(window) || window < 1) stop("`window` must be a single number >= 1.", call. = FALSE)
  min_genotype_count <- check_count(min_genotype_count, "min_genotype_count")
  min_n_mid <- check_count(min_n_mid, "min_n_mid")
  ncores <- check_count(ncores, "ncores", min = 1); chunk_size <- check_count(chunk_size, "chunk_size", min = 1)
  check_flag(mirror, "mirror"); check_flag(all_stats, "all_stats")
  if (!is.numeric(boundary_share) || length(boundary_share) != 2L || anyNA(boundary_share) || boundary_share[1] > boundary_share[2])
    stop("`boundary_share` must be two increasing numbers.", call. = FALSE)
  if (!is.numeric(fdr) || length(fdr) != 1L || !(fdr > 0 && fdr < 1)) stop("`fdr` must lie strictly between 0 and 1.", call. = FALSE)
  if (!is.null(major_score) && (!is.numeric(major_score) || is.null(names(major_score))))
    stop("`major_score` must be a numeric vector named by promoter ID.", call. = FALSE)
  if (!is.null(out_dir) && (!is.character(out_dir) || length(out_dir) != 1L || is.na(out_dir) || !nzchar(out_dir)))
    stop("`out_dir` must be the path of one folder.", call. = FALSE)
  st <- sample_info(tissue, individual, samples)
  prom <- check_promoters(read_table(promoters, "promoters"), rownames(phenotype))
  sel <- if (is.null(chromosomes)) rep(TRUE, nrow(prom)) else prom$chr %in% as.character(chromosomes)
  if (!any(sel)) stop("No promoter lies on the requested chromosomes; promoter chromosomes are: ",
                      paste(utils::head(unique(prom$chr), 10), collapse = ", "), call. = FALSE)
  if (is_vcf_path(genotype)) {
    if (!is.null(variants)) stop("`variants` must be NULL when `genotype` is a VCF file: the positions are read from the VCF.", call. = FALSE)
    vc <- read_vcf_dosage(genotype, samples = unique(st$individual)); genotype <- vc$genotype; variants <- vc$variants
  } else {
    genotype <- read_matrix(genotype, "genotype")
    if (is.null(variants)) stop("`variants` is needed (a data frame or TSV file with columns id, chr and pos) unless `genotype` is a VCF file.", call. = FALSE)
    variants <- read_table(variants, "variants")
  }
  opt <- list(window = window, chromosomes = chromosomes, sel = sel, min_genotype_count = min_genotype_count, min_n_mid = min_n_mid,
              boundary_share = boundary_share, fdr = fdr, acat_cap = acat_cap, mirror = mirror, mirror_tol = mirror_tol,
              control = control, ncores = ncores, chunk_size = chunk_size, all_stats = all_stats)
  if (is.null(st$tissue)) {
    G_all <- check_genotype(genotype, st$individual, samples, individual_given = !identical(st$individual, samples))
    V <- check_variants(variants, rownames(G_all))
    out <- cis_one(phenotype, G_all, V, prom, read_covariates(covariates), major_score, opt, t0)
  } else {
    G_all <- check_genotype(genotype, unique(st$individual), individual_given = !identical(st$individual, samples))
    V <- check_variants(variants, rownames(G_all))
    if (is_path(covariates) && !dir.exists(covariates)) covariates <- read_matrix(covariates, "covariates")
    ## the canonical member of mirror pairs: mean over all supplied samples pooled across tissues (paper's rule)
    if (is.null(major_score)) major_score <- rowMeans(phenotype, na.rm = TRUE)
    res <- lapply(st$levels, function(tt) {
      w <- st$tissue == tt; s <- samples[w]
      G <- G_all[, st$individual[w], drop = FALSE]; colnames(G) <- s
      cis_one(phenotype[, s, drop = FALSE], G, V, prom, tissue_covariates(covariates, tt), major_score, opt,
              proc.time()[["elapsed"]], label = tt)
    })
    out <- bind_tissues(res, st$levels)
    out$seconds <- proc.time()[["elapsed"]] - t0
  }
  if (!is.null(out_dir)) write_cis_tables(out, out_dir)
  out
}

## cis mapping of one tissue (run_final2.R, run_muF.R, discovery_nogate.R, units_v3_estimable.R, units_v3m_mirror.R).
## phenotype: promoters x samples of the tissue; G_all: checked dosages (integer), variants x the same samples;
## V: checked variants (rows of G_all); prom/opt$sel: checked promoters (all rows of phenotype) and the selection;
## covariates: NULL, matrix or data frame (rows: sample IDs or the order of the samples). label: tissue for messages.
cis_one <- function(phenotype, G_all, V, prom, covariates, major_score, opt, t0, label = NULL) {
  window <- opt$window; chromosomes <- opt$chromosomes; sel <- opt$sel; min_genotype_count <- opt$min_genotype_count
  min_n_mid <- opt$min_n_mid; boundary_share <- opt$boundary_share; fdr <- opt$fdr; acat_cap <- opt$acat_cap; mirror <- opt$mirror
  mirror_tol <- opt$mirror_tol; control <- opt$control; ncores <- opt$ncores; chunk_size <- opt$chunk_size; all_stats <- opt$all_stats
  samples <- colnames(phenotype)
  Z <- as_covariates(covariates, length(samples))
  if (!is.null(covariates)) {
    if (!is.null(rownames(Z))) {
      miss <- setdiff(samples, rownames(Z))
      if (length(miss)) stop("`covariates` has no row for ", length(miss), " sample(s)", if (!is.null(label)) paste0(" of tissue ", label),
                             ", e.g. ", miss[1], ".", call. = FALSE)
      Z <- Z[samples, , drop = FALSE]
    } else if (nrow(Z) != length(samples)) stop("`covariates` must have one row per sample of `phenotype`.", call. = FALSE)
  }
  rownames(Z) <- samples
  C_all <- make_design(Z)
  ## covariate count (author 2026-11-10): warn, never drop; the paper kept at most one covariate per 30 donors where
  ## possible (SI Section 4: K = max(0, min(K0, floor(n/30) - 5)) PEER factors plus 5 genotype PCs)
  if (ncol(Z) > floor(length(samples) / 30))
    warning(sprintf(paste0("%s%d covariates for %d samples: more than floor(n/30) = %d. Consider fewer covariates (for example fewer hidden ",
                           "factors). The BOHR paper kept at most one covariate per 30 donors where possible. Nothing was dropped."),
                    if (is.null(label)) "" else paste0("Tissue ", label, ": "), ncol(Z), length(samples), floor(length(samples) / 30)), call. = FALSE)

  ## mirror pairs over the whole matrix (promoter counts per gene), flags for the selected chromosomes
  has_gene <- "gene" %in% names(prom) && any(!is.na(prom$gene))
  n_prom_gene <- if (has_gene) { tb <- table(prom$gene[!is.na(prom$gene)]); as.integer(tb[prom$gene]) } else rep(NA_integer_, nrow(prom))
  pairs <- NULL
  if (mirror && has_gene) {
    if (is.null(major_score)) major_score <- rowMeans(phenotype, na.rm = TRUE)
    pairs <- pair_table(phenotype, prom$gene, major_score, tol = mirror_tol, genes = unique(prom$gene[sel]))
  }

  ## variants: chromosome index; fill missing dosages of the variants that can enter a window
  todo <- which(sel); todo <- todo[order(match(prom$chr[todo], unique(prom$chr[todo])), prom$pos[todo], method = "radix")]
  vchr <- split(seq_len(nrow(V)), V$chr)
  win <- lapply(todo, function(i) { vi <- vchr[[prom$chr[i]]]; if (is.null(vi)) integer(0) else vi[cis_index(V$start[vi], V$end[vi], prom$pos[i], window)] })
  used <- sort(unique(unlist(win)))
  if (length(used)) { Gi <- impute_modal(G_all[used, , drop = FALSE]); G_all <- NULL } else Gi <- NULL
  pos_in <- match(seq_len(nrow(V)), used)

  ## per-promoter preparation (run_final2.R:37-68)
  prep_one <- function(k) {
    i <- todo[k]; p <- prom$id[i]
    base <- list(promoter = p, chr = prom$chr[i], pos = prom$pos[i], gene = prom$gene[i], family = NA_character_, skip = NA_character_,
                 n = NA_integer_, n0 = NA_integer_, n1 = NA_integer_, n_mid = NA_integer_, n_cov = NA_integer_, n_cov_dropped = NA_integer_,
                 rank_mid = NA_integer_, null_beta_conv = NA, null_cat_conv = NA, null_ll_beta = NA_real_, null_ll_cat = NA_real_,
                 n_variants_window = NA_integer_, n_variants_tested = 0L, n_errors = 0L, n_firth_errors = 0L)
    fit <- fit_promoter(phenotype[p, ], C_all, control)
    for (k2 in c("family", "skip", "n", "n0", "n1", "n_mid", "n_cov", "n_cov_dropped", "rank_mid")) base[[k2]] <- fit[[k2]]
    if (!is.na(fit$skip)) return(list(prom = base))
    nul <- fit$null
    base$null_beta_conv <- if (!is.null(nul$beta)) nul$beta$converged else NA
    base$null_cat_conv  <- if (!is.null(nul$cat)) nul$cat$converged else NA
    base$null_ll_beta   <- if (!is.null(nul$beta)) nul$beta$ll else NA_real_
    base$null_ll_cat    <- if (!is.null(nul$cat)) nul$cat$ll else NA_real_
    idx <- win[[k]]
    if (!length(idx)) { base$skip <- "no_snp"; return(list(prom = base)) }
    base$n_variants_window <- length(idx)
    G <- Gi[pos_in[idx], fit$samples, drop = FALSE]
    cnt <- apply(G, 1, function(x) tabulate(x + 1L, nbins = 3L)); dim(cnt) <- c(3L, length(idx))
    keep <- which(apply(cnt, 2, min) >= min_genotype_count)
    if (!length(keep)) { base$skip <- "no_snp"; return(list(prom = base)) }
    list(prom = base, fit = fit, G = G[keep, , drop = FALSE], cnt = cnt[, keep, drop = FALSE], vid = idx[keep])
  }
  prep <- par_apply(seq_along(todo), prep_one, ncores, preschedule = TRUE)

  ## variant tests in chunks, heaviest first (run_final2.R:75-95)
  tasks <- list()
  for (i in seq_along(prep)) if (!is.null(prep[[i]]$G)) {
    m <- nrow(prep[[i]]$G); fam <- prep[[i]]$prom$family
    w <- c(BEINF = 3, BEZI = 1, BEOI = 1, BE = 0.7)[[fam]] * prep[[i]]$prom$n * ncol(prep[[i]]$fit$design)^1.5 * (if (control$firth && fam != "BE") 4 else 1)
    for (ch in split(seq_len(m), ceiling(seq_len(m) / chunk_size))) tasks[[length(tasks) + 1L]] <- list(i = i, j = ch, w = w * length(ch))
  }
  if (length(tasks)) tasks <- tasks[order(-vapply(tasks, `[[`, 0, "w"))]
  run_chunk <- function(t) {
    P <- prep[[t$i]]
    rows <- lapply(t$j, function(j) {
      e <- engine_one(P$fit, unname(P$G[j, ]), control)
      v <- P$vid[j]; cn <- P$cnt[, j]; af <- (cn[2] + 2 * cn[3]) / (2 * length(P$fit$y))
      c(list(promoter = P$prom$promoter, variant = V$id[v], chr = V$chr[v], pos = V$start[v], alt_af = af, maf = min(af, 1 - af),
             min_geno = min(cn)), engine_to_row(e, all_stats, control$firth))
    })
    list(i = t$i, df = rows_to_df(rows))
  }
  res <- par_apply(tasks, run_chunk, ncores, preschedule = FALSE)
  for (r in res) { b <- prep[[r$i]]$prom
    b$n_variants_tested <- b$n_variants_tested + nrow(r$df); b$n_errors <- b$n_errors + sum(!is.na(r$df$error))
    b$n_firth_errors <- b$n_firth_errors + sum(r$df$firth_error); prep[[r$i]]$prom <- b }
  snp <- bind_df(lapply(res, `[[`, "df"))
  prom_out <- rows_to_df(lapply(prep, `[[`, "prom"))
  if (!is.null(snp)) { o <- order(match(snp$promoter, prom_out$promoter), snp$pos, snp$variant, method = "radix"); snp <- snp[o, , drop = FALSE]; rownames(snp) <- NULL }

  units <- make_units(snp, prom_out, n_prom_gene[match(prom_out$promoter, prom$id)], mirror_flags(pairs),
                      min_n_mid, boundary_share, acat_cap)
  units <- units_bh(units, fdr)
  out <- list(snp = snp, units = units, promoters = prom_out, mirror = pairs,
              threshold = unit_threshold(units, fdr),
              settings = list(window = window, chromosomes = if (is.null(chromosomes)) sort(unique(prom$chr[sel])) else as.character(chromosomes),
                              min_genotype_count = min_genotype_count, min_n_mid = min_n_mid, boundary_share = boundary_share,
                              fdr = fdr, acat_cap = acat_cap, mirror = mirror, mirror_tol = mirror_tol, control = control,
                              n_samples = length(samples), n_covariates = ncol(Z)),
              seconds = proc.time()[["elapsed"]] - t0)
  structure(out, class = "bohr_cis")
}

## ---- units: ACAT per promoter x component (discovery_nogate.R:43-49) -----------------------------
UNIT_COMP <- list(mu = c(col = "log_p_mu", est = "est_mu", eff = "effect_mu"),
                  p0 = c(col = "log_p_p0", est = "est_p0", eff = "effect_p0"),
                  p1 = c(col = "log_p_p1", est = "est_p1", eff = "effect_p1"))
make_units <- function(snp, prom, n_prom_gene, flags, min_n_mid, boundary_share, acat_cap = TRUE) {
  fitted <- prom[is.na(prom$skip) & !is.na(prom$family), , drop = FALSE]
  if (is.null(snp) || !nrow(snp) || !nrow(fitted)) return(empty_units())
  rows_of <- split(seq_len(nrow(snp)), snp$promoter)
  has <- list(mu = fitted$promoter, p0 = fitted$promoter[fitted$family %in% c("BEZI", "BEINF")],
              p1 = fitted$promoter[fitted$family %in% c("BEOI", "BEINF")])
  U <- list()
  for (ty in names(UNIT_COMP)) for (p in has[[ty]]) {
    r <- rows_of[[p]]; if (is.null(r)) next
    r <- r[order(snp$variant[r], method = "radix")]                        # = row order of the original (merge on promoter, snp)
    lp <- snp[[UNIT_COMP[[ty]][["col"]]]][r]
    j <- which.min(replace(lp, !is.finite(lp), Inf)); jr <- r[j]
    U[[length(U) + 1L]] <- list(promoter = p, component = ty, n_variants = sum(is.finite(lp)), p = acat_log(lp, cap = acat_cap),
      leader = snp$variant[jr], leader_pos = snp$pos[jr], leader_log_p = lp[j], leader_est = snp[[UNIT_COMP[[ty]][["est"]]]][jr],
      leader_effect = snp[[UNIT_COMP[[ty]][["eff"]]]][jr],
      log_p_mu_at_leader = snp$log_p_mu[jr], log_p_p0_at_leader = snp$log_p_p0[jr], log_p_p1_at_leader = snp$log_p_p1[jr],
      est_mu_at_leader = snp$est_mu[jr], est_p0_at_leader = snp$est_p0[jr], est_p1_at_leader = snp$est_p1[jr])
  }
  U <- rows_to_df(U)
  if (is.null(U)) return(empty_units())
  U <- U[is.finite(U$p), , drop = FALSE]
  if (!nrow(U)) return(empty_units())
  i <- match(U$promoter, fitted$promoter)
  U <- data.frame(promoter = U$promoter, gene = fitted$gene[i], component = U$component, family = fitted$family[i],
                  n = fitted$n[i], n_mid = fitted$n_mid[i], n0 = fitted$n0[i], n1 = fitted$n1[i],
                  boundary_share = ifelse(U$component == "p0", fitted$n0[i] / fitted$n[i], ifelse(U$component == "p1", fitted$n1[i] / fitted$n[i], NA_real_)),
                  U[, setdiff(names(U), c("promoter", "component")), drop = FALSE], stringsAsFactors = FALSE)
  ## estimability (units_v3_estimable.R:15-17)
  U$estimable <- !(U$component == "mu" & U$n_mid < min_n_mid) &
    (U$component == "mu" | (U$boundary_share >= boundary_share[1] & U$boundary_share <= boundary_share[2]))
  ## mirror flags and the v3m de-duplication (discovery_nogate.R:50-57, units_v3m_mirror.R:35-37)
  U$n_prom_gene <- n_prom_gene[match(U$promoter, prom$promoter)]
  fi <- if (is.null(flags)) rep(NA_integer_, nrow(U)) else match(U$promoter, flags$promoter)
  get_f <- function(k, na) { v <- if (is.null(flags)) rep(na, nrow(U)) else flags[[k]][fi]; v[is.na(fi)] <- na; v }
  U$mirror_pair_id <- get_f("mirror_pair_id", NA_character_)
  U$boundary_exact <- get_f("boundary_exact", FALSE); U$interior_exact <- get_f("interior_exact", FALSE)
  U$mirror_role <- get_f("mirror_role", NA_character_); U$mirror_partner <- get_f("mirror_partner", NA_character_)
  drop_mu <- get_f("drop_mu", FALSE)
  U$dedup_drop <- ifelse(U$component == "mu", drop_mu, U$boundary_exact & U$component == "p1")
  lab <- c(mu = "mu", p0 = "p1", p1 = "p0")
  U$duplicate_of <- ifelse(U$dedup_drop, paste0(U$mirror_partner, ":", lab[U$component]), NA_character_)
  rownames(U) <- NULL
  U
}
empty_units <- function() NULL

## BH within the tissue on the primary units (units_v3m_mirror.R:44)
units_bh <- function(U, fdr) {
  if (is.null(U) || !nrow(U)) return(U)
  prim <- !U$dedup_drop & U$estimable
  U$q <- NA_real_
  U$q[prim] <- stats::p.adjust(U$p[prim], "BH")
  U$significant <- !is.na(U$q) & U$q <= fdr
  U
}
unit_threshold <- function(U, fdr) if (is.null(U) || !any(U$significant)) NA_real_ else max(U$p[U$significant])
## BH within each tissue of a tissue-tagged unit table (one call of units_bh() per tissue; rows keep their order)
units_bh_tissue <- function(U, fdr) {
  if (is.null(U) || !nrow(U) || is.null(U$tissue)) return(units_bh(U, fdr))
  U$q <- NA_real_; U$significant <- FALSE
  for (tt in unique(U$tissue)) { i <- which(U$tissue == tt); b <- units_bh(U[i, , drop = FALSE], fdr); U$q[i] <- b$q; U$significant[i] <- b$significant }
  U
}

#' Combine cis-mapping results run in parts
#'
#' Combines [bohr_cis()] results of the same samples run separately, for
#' example one chromosome per cluster job, and recomputes the q values over
#' all units together, so that the false-discovery control is within the
#' whole tissue. Results of a tissue loop (`bohr_cis(tissue = )`) are combined
#' the same way, with the q values recomputed within each tissue.
#'
#' @param ... `bohr_cis` objects, or one list of them. They must use the same
#'   samples, tissues and settings and cover different promoters (for example
#'   different chromosomes).
#' @param fdr False-discovery rate for the `significant` column; defaults to
#'   that of the inputs.
#' @return A `bohr_cis` object.
#' @seealso [bohr_cis()].
#' @examples
#' set.seed(5)
#' n <- 100; s <- paste0("S", 1:n)
#' geno <- matrix(rbinom(20 * n, 2, 0.4), 20, n, dimnames = list(paste0("v", 1:20), s))
#' vars <- data.frame(id = rownames(geno), chr = rep(c("chr1", "chr2"), each = 10),
#'                    pos = rep(seq(1000, 10000, 1000), 2))
#' ph <- rbind(P1 = simulate_response(n, p1 = 0)$y, P2 = simulate_response(n, p0 = 0)$y)
#' colnames(ph) <- s
#' proms <- data.frame(id = c("P1", "P2"), chr = c("chr1", "chr2"), pos = 5000)
#' r1 <- bohr_cis(ph, geno, vars, proms, window = 2000, chromosomes = "chr1")
#' r2 <- bohr_cis(ph, geno, vars, proms, window = 2000, chromosomes = "chr2")
#' bohr_combine(r1, r2)$units[, c("promoter", "component", "p", "q")]
#' @export
bohr_combine <- function(..., fdr = NULL) {
  L <- list(...)
  if (length(L) == 1L && is.list(L[[1]]) && !inherits(L[[1]], "bohr_cis")) L <- L[[1]]
  if (!length(L) || !all(vapply(L, inherits, TRUE, what = "bohr_cis"))) stop("Supply bohr_cis objects.", call. = FALSE)
  key <- function(s) s[c("window", "min_genotype_count", "min_n_mid", "boundary_share", "acat_cap", "mirror", "mirror_tol", "control",
                         "n_samples", "n_covariates", "tissues")]
  if (length(unique(lapply(L, function(x) key(x$settings)))) > 1L) stop("The results use different settings or samples.", call. = FALSE)
  P <- bind_df(lapply(L, `[[`, "promoters"))
  if (anyDuplicated(if (is.null(P$tissue)) P$promoter else paste(P$tissue, P$promoter, sep = "\r")))
    stop("The results overlap (a promoter appears more than once).", call. = FALSE)
  if (is.null(fdr)) fdr <- L[[1]]$settings$fdr
  U <- bind_df(lapply(L, `[[`, "units"))
  if (!is.null(U)) U$q <- U$significant <- NULL
  tis <- L[[1]]$settings$tissues
  U <- if (is.null(tis)) units_bh(U, fdr) else units_bh_tissue(U, fdr)
  M <- bind_df(lapply(L, `[[`, "mirror"))
  s <- L[[1]]$settings; s$fdr <- fdr; s$chromosomes <- sort(unique(unlist(lapply(L, function(x) x$settings$chromosomes))))
  thr <- if (is.null(tis)) unit_threshold(U, fdr) else
    stats::setNames(vapply(tis, function(tt) unit_threshold(if (is.null(U)) NULL else U[U$tissue == tt, , drop = FALSE], fdr), 0), tis)
  structure(list(snp = bind_df(lapply(L, `[[`, "snp")), units = U, promoters = P, mirror = M,
                 threshold = thr, settings = s, seconds = sum(vapply(L, `[[`, 0, "seconds"))),
            class = "bohr_cis")
}

#' @export
print.bohr_cis <- function(x, ...) {
  P <- x$promoters; U <- x$units
  if (!is.null(x$settings$tissues)) {
    tis <- x$settings$tissues
    cat("BOHR cis mapping:", length(tis), "tissues,", length(unique(P$promoter)), "promoters on", length(x$settings$chromosomes), "chromosome(s)\n")
    cnt <- function(v) vapply(tis, function(tt) sum(v[P$tissue == tt]), 0L)
    tab <- data.frame(tissue = tis, samples = unname(x$settings$n_samples), fitted = cnt(is.na(P$skip)), skipped = cnt(!is.na(P$skip)),
                      variant_tests = vapply(tis, function(tt) if (is.null(x$snp)) 0L else sum(x$snp$tissue == tt), 0L),
                      units = vapply(tis, function(tt) if (is.null(U)) 0L else sum(U$tissue == tt), 0L),
                      tested = vapply(tis, function(tt) if (is.null(U)) 0L else sum(U$tissue == tt & !U$dedup_drop & U$estimable), 0L),
                      significant = vapply(tis, function(tt) if (is.null(U)) 0L else sum(U$tissue == tt & U$significant), 0L),
                      row.names = NULL, stringsAsFactors = FALSE)
    print(tab, row.names = FALSE)
    cat(sprintf("q values within each tissue; significant at q <= %g. Tested units: estimable and not mirror duplicates.\n", x$settings$fdr))
    return(invisible(x))
  }
  cat("BOHR cis mapping:", x$settings$n_samples, "samples,", nrow(P), "promoters on", length(x$settings$chromosomes), "chromosome(s)\n")
  cat("  fitted promoters:", sum(is.na(P$skip)), "; skipped:", sum(!is.na(P$skip)), "\n")
  cat("  variant tests:", if (is.null(x$snp)) 0L else nrow(x$snp), "\n")
  if (!is.null(U) && nrow(U)) {
    prim <- !U$dedup_drop & U$estimable
    cat("  units:", nrow(U), "; tested (estimable, not mirror duplicates):", sum(prim), "\n")
    sig <- table(factor(U$component[U$significant], levels = c("mu", "p0", "p1")))
    cat(sprintf("  significant at q <= %g: %d (mu %d, p0 %d, p1 %d)\n", x$settings$fdr, sum(U$significant), sig[["mu"]], sig[["p0"]], sig[["p1"]]))
  } else cat("  units: 0\n")
  invisible(x)
}

## ---- parallel map with one serial retry of failed tasks ------------------------------------------
par_apply <- function(X, FUN, ncores, preschedule = TRUE) {
  if (!length(X)) return(list())
  if (ncores > 1L && .Platform$OS.type != "windows") {
    res <- parallel::mclapply(X, FUN, mc.cores = ncores, mc.preschedule = preschedule)
    bad <- which(vapply(res, function(r) is.null(r) || inherits(r, "try-error"), TRUE))
    for (b in bad) res[b] <- list(FUN(X[[b]]))                               # a worker died or failed: retry once here
    res
  } else lapply(X, FUN)
}

## ---- input checks ------------------------------------------------------------------------------------
check_phenotype <- function(phenotype, sample_names = TRUE) {
  if (is.data.frame(phenotype)) phenotype <- as.matrix(phenotype)
  if (!is.matrix(phenotype) || !is.numeric(phenotype)) stop("`phenotype` must be a numeric matrix (promoters x samples).", call. = FALSE)
  if (is.null(rownames(phenotype)) || anyDuplicated(rownames(phenotype)) || anyNA(rownames(phenotype)))
    stop("`phenotype` needs unique row names (promoter IDs).", call. = FALSE)
  if (sample_names && (is.null(colnames(phenotype)) || anyDuplicated(colnames(phenotype)) || anyNA(colnames(phenotype))))
    stop("`phenotype` needs unique column names (sample IDs).", call. = FALSE)
  check_activity(phenotype, "phenotype")
  phenotype
}
## cols: the genotype columns to take (one per sample, or unique individuals); names: the column names of the result
check_genotype <- function(genotype, cols, names = cols, individual_given = FALSE) {
  if (is.data.frame(genotype)) genotype <- as.matrix(genotype)
  if (!is.matrix(genotype) || !is.numeric(genotype)) stop("`genotype` must be a numeric matrix (variants x samples).", call. = FALSE)
  if (is.null(rownames(genotype)) || anyDuplicated(rownames(genotype)) || anyNA(rownames(genotype)))
    stop("`genotype` needs unique row names (variant IDs).", call. = FALSE)
  miss <- setdiff(cols, colnames(genotype))
  if (length(miss)) stop("`genotype` has no column for ", length(miss), if (individual_given) " individual(s) of `individual`" else
                         " sample(s) of `phenotype`", ", e.g. ", miss[1], ".", call. = FALSE)
  G <- genotype[, cols, drop = FALSE]
  if (any(!is.na(G) & !(G %in% c(0, 1, 2)))) stop("`genotype` must hold dosages 0, 1, 2 or NA.", call. = FALSE)
  storage.mode(G) <- "integer"
  colnames(G) <- names
  G
}
check_variants <- function(variants, ids) {
  if (!is.data.frame(variants) || !all(c("chr", "pos") %in% names(variants))) stop("`variants` must be a data frame with columns chr and pos.", call. = FALSE)
  if ("id" %in% names(variants)) {
    i <- match(ids, as.character(variants$id))
    if (anyNA(i)) stop("`variants$id` has no entry for ", sum(is.na(i)), " row(s) of `genotype`, e.g. ", ids[is.na(i)][1], ".", call. = FALSE)
    variants <- variants[i, , drop = FALSE]
  } else if (nrow(variants) != length(ids)) stop("`variants` must have one row per row of `genotype` (or an id column).", call. = FALSE)
  start <- as.numeric(variants$pos); end <- if ("end" %in% names(variants)) as.numeric(variants$end) else start
  if (anyNA(start) || anyNA(end) || any(end < start)) stop("`variants$pos` (and `end`) must be non-missing with end >= pos.", call. = FALSE)
  data.frame(id = ids, chr = as.character(variants$chr), start = start, end = end, stringsAsFactors = FALSE)
}
check_promoters <- function(promoters, ids) {
  if (!is.data.frame(promoters) || !all(c("id", "chr", "pos") %in% names(promoters))) stop("`promoters` must be a data frame with columns id, chr and pos.", call. = FALSE)
  i <- match(ids, as.character(promoters$id))
  if (anyNA(i)) stop("`promoters` has no row for ", sum(is.na(i)), " promoter(s) of `phenotype`, e.g. ", ids[is.na(i)][1], ".", call. = FALSE)
  pr <- promoters[i, , drop = FALSE]
  pos <- as.numeric(pr$pos)
  if (anyNA(pos)) stop("`promoters$pos` must be non-missing.", call. = FALSE)
  data.frame(id = ids, chr = as.character(pr$chr), pos = pos,
             gene = if ("gene" %in% names(pr)) as.character(pr$gene) else NA_character_, stringsAsFactors = FALSE)
}
