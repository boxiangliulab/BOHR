## proActiv adapter (layer 3): from_proactiv() reads proActiv results, relative_activity() computes the relative activity
## over a gene's external promoters and applies the paper's promoter filter.
##
## Follows code/promoter_quantification/0103_proActiv_data_process_tissue.R of the paper (md5 6fe1bc29...; SI Section 1,
## code map rows S1.3-S1.4): proActiv classes (column "<condition>.class"; complete in every condition and not "Inactive" in
## at least one) select the active promoters (:96-98); promoters of genes with one annotated promoter go (:105-111); in each
## tissue the share of samples with a relative activity above 0, the share of missing values and the mean relative activity
## are computed on proActiv's relative activity (:118-160) and a promoter needs, in at least one tissue (not necessarily the
## same one), a share above 0 greater than 0.25 (:163), a missing share below 0.40 (:168) and a mean greater than 0.1
## (:173); genes left with one promoter go (:180-189). The relative activity used downstream is recomputed as the absolute
## activity over the sum of the absolute activities of the gene's external promoters, all of them, before any filter
## (:214-268), and restricted to the kept promoters (:303). proActiv's class rule (categorizePromoters(), proActiv 1.16.0):
## a promoter is "Inactive" in a condition when its mean absolute activity there is below 0.25; internal promoters get NA.

PROACTIV_ASSAYS <- c(absolute = "absolutePromoterActivity", relative = "relativePromoterActivity")

## the object both functions share: matrices for the external promoters, the table and classes for every promoter
new_proactiv <- function(absolute, relative, promoters, class = NULL, condition = NULL) {
  structure(list(absolute = absolute, relative = relative, promoters = promoters, class = class, condition = condition),
            class = "bohr_proactiv")
}

#' @export
print.bohr_proactiv <- function(x, ...) {
  p <- x$promoters
  cat("proActiv results for BOHR:", nrow(p), "promoters (", sum(!p$internal), "external ) of", length(unique(p$gene)),
      "genes in", ncol(x$absolute), "samples")
  if (!is.null(x$class)) cat(";", ncol(x$class), "conditions with promoter classes")
  cat("\n")
  invisible(x)
}

as_dense <- function(m, what) {
  if (is.data.frame(m)) m <- as.matrix(m[, vapply(m, is.numeric, TRUE), drop = FALSE])
  if (!is.matrix(m)) m <- as.matrix(m)                                   # Matrix classes and DelayedArray
  if (!is.numeric(m) && !(is.logical(m) && all(is.na(m)))) stop("`", what, "` must be numeric.", call. = FALSE)
  storage.mode(m) <- "double"
  m
}

#' Read proActiv results
#'
#' Collects what [relative_activity()] needs from the output of
#' `proActiv::proActiv()`: the absolute and relative promoter activities of
#' the external promoters, the promoter annotation (gene and internal flag of
#' every promoter) and proActiv's promoter classes.
#'
#' proActiv is often run once per condition (in the BOHR paper, once per
#' GTEx tissue). Give the results as a list (or as paths of `.rds` files);
#' they must hold the same promoters in the same order, and their samples are
#' bound together. Each result's class columns (`"<condition>.class"` in the
#' row data, written by proActiv when it is given a `condition`) are kept,
#' one per condition.
#'
#' Internal promoters are listed in the promoter table but their activities
#' are not kept: proActiv cannot quantify them, they are not in the
#' denominator of [relative_activity()] and they are never tested.
#'
#' @param se A `SummarizedExperiment` returned by `proActiv::proActiv()`, a
#'   list of them, or the path(s) of `.rds` files holding them. Needs the
#'   package SummarizedExperiment.
#' @param id_prefix Prefix of the promoter IDs, which are built from
#'   proActiv's `promoterId`. Default `"PromoterID_"`, the IDs of the paper.
#' @return An object of class `"bohr_proactiv"`, a list with `absolute` and
#'   `relative` (proActiv's absolute and relative activity of the external
#'   promoters; rows = promoter IDs, columns = samples), `promoters` (one row
#'   per promoter: `id`, `gene`, `chr`, `pos`, `strand`, `internal`), `class`
#'   (proActiv's class of every promoter in every condition: `"Major"`,
#'   `"Minor"`, `"Inactive"` or `NA`; `NULL` when there are none) and
#'   `condition` (each sample's condition, when proActiv recorded one).
#' @seealso [relative_activity()].
#' @examples
#' \dontrun{
#' pa <- from_proactiv(list.files("proActiv", pattern = "\\.rds$", full.names = TRUE))
#' ra <- relative_activity(pa, group = sample_tissue)
#' }
#' @export
from_proactiv <- function(se, id_prefix = "PromoterID_") {
  if (!has_pkg("SummarizedExperiment")) stop("from_proactiv() needs the package SummarizedExperiment.", call. = FALSE)
  if (!is.character(id_prefix) || length(id_prefix) != 1L || is.na(id_prefix)) stop("`id_prefix` must be a single string.", call. = FALSE)
  if (is.character(se)) {
    miss <- se[!file.exists(se)]
    if (length(miss)) stop("`se`: file not found: ", miss[1], call. = FALSE)
    se <- lapply(se, readRDS)
  }
  if (!is.list(se) || methods::is(se, "SummarizedExperiment")) se <- list(se)
  if (!length(se)) stop("`se` is empty.", call. = FALSE)
  ok <- vapply(se, function(s) methods::is(s, "SummarizedExperiment"), TRUE)
  if (!all(ok)) stop("`se` must hold SummarizedExperiment objects returned by proActiv().", call. = FALSE)
  rd <- lapply(se, function(s) as.data.frame(SummarizedExperiment::rowData(s)))
  need <- c("promoterId", "geneId", "internalPromoter")
  for (k in seq_along(rd)) {
    if (!all(need %in% names(rd[[k]]))) stop("Result ", k, " of `se` has no row data columns ", paste(need, collapse = ", "),
                                              "; is it a proActiv result?", call. = FALSE)
    if (!identical(as.character(rd[[k]]$promoterId), as.character(rd[[1]]$promoterId)))
      stop("Result ", k, " of `se` does not hold the same promoters in the same order as the first.", call. = FALSE)
    an <- SummarizedExperiment::assayNames(se[[k]])
    if (!all(PROACTIV_ASSAYS %in% an)) stop("Result ", k, " of `se` has no assays ", paste(PROACTIV_ASSAYS, collapse = " and "), ".", call. = FALSE)
  }
  r1 <- rd[[1]]
  ids <- paste0(id_prefix, r1$promoterId)
  if (anyDuplicated(ids)) stop("proActiv's promoterId is not unique.", call. = FALSE)
  internal <- as.logical(r1$internalPromoter)
  if (anyNA(internal)) internal[is.na(internal)] <- TRUE                 # proActiv gives these no class either
  col_or_na <- function(k) if (k %in% names(r1)) r1[[k]] else rep(NA, nrow(r1))
  promoters <- data.frame(id = ids, gene = as.character(r1$geneId), chr = as.character(col_or_na("seqnames")),
                          pos = suppressWarnings(as.numeric(col_or_na("start"))), strand = as.character(col_or_na("strand")),
                          internal = internal, stringsAsFactors = FALSE)
  ext <- which(!internal)
  get_assay <- function(k, a) {
    m <- as_dense(SummarizedExperiment::assay(se[[k]], a), a)[ext, , drop = FALSE]
    rownames(m) <- ids[ext]
    m
  }
  absolute <- do.call(cbind, lapply(seq_along(se), get_assay, a = PROACTIV_ASSAYS[["absolute"]]))
  relative <- do.call(cbind, lapply(seq_along(se), get_assay, a = PROACTIV_ASSAYS[["relative"]]))
  if (is.null(colnames(absolute)) || anyDuplicated(colnames(absolute)))
    stop("The samples of `se` need unique names (column names).", call. = FALSE)
  cls <- NULL
  for (k in seq_along(rd)) {
    cc <- grep("\\.class$", names(rd[[k]]), value = TRUE)
    for (c1 in cc) {
      nm <- sub("\\.class$", "", c1)
      if (!is.null(cls) && nm %in% colnames(cls)) next
      cls <- cbind(cls, matrix(as.character(rd[[k]][[c1]]), ncol = 1L, dimnames = list(NULL, nm)))
    }
  }
  if (!is.null(cls)) rownames(cls) <- ids
  cond <- unlist(lapply(se, function(s) {
    cd <- SummarizedExperiment::colData(s)
    if ("condition" %in% names(cd)) stats::setNames(as.character(cd$condition), colnames(s)) else NULL
  }))
  if (!is.null(cond) && length(cond) != ncol(absolute)) cond <- NULL
  new_proactiv(absolute, relative, promoters, cls, cond)
}

## the gene, internal flag and ID of each row of a plain matrix: a vector per row, or a TSS table (id, gene, internal)
ra_annotation <- function(x, gene, internal) {
  ids <- rownames(x)
  if (is.null(ids) || anyDuplicated(ids)) stop("`x` needs unique row names (promoter IDs).", call. = FALSE)
  if (is.data.frame(gene) || is_path(gene)) {
    tab <- read_table(gene, "gene")
    if (!all(c("id", "gene") %in% names(tab))) stop("A table given as `gene` needs the columns id and gene (and optionally internal).", call. = FALSE)
    i <- match(ids, as.character(tab$id))
    if (anyNA(i)) stop("The table `gene` has no row for ", sum(is.na(i)), " promoter(s) of `x`, e.g. ", ids[is.na(i)][1], ".", call. = FALSE)
    if (is.null(internal) && "internal" %in% names(tab)) internal <- as.logical(toupper(as.character(tab$internal[i])))
    gene <- as.character(tab$gene[i])
  } else {
    if (is.null(gene) || length(gene) != length(ids)) stop("`gene` must give the gene of each row of `x` (or be a table with columns id and gene).", call. = FALSE)
    if (!is.null(names(gene)) && all(ids %in% names(gene))) gene <- gene[ids]
    gene <- as.character(unname(gene))
  }
  if (anyNA(gene) || any(!nzchar(gene))) stop("`gene` has missing values.", call. = FALSE)
  if (is.null(internal)) internal <- rep(FALSE, length(ids))
  if (!is.logical(internal) || length(internal) != length(ids) || anyNA(internal))
    stop("`internal` must be TRUE or FALSE for each row of `x`.", call. = FALSE)
  data.frame(id = ids, gene = gene, internal = unname(internal), stringsAsFactors = FALSE)
}

## per row: is the criterion met in at least one group? (stat(): rows x samples of one group -> one value per row)
any_group <- function(R, groups, stat, keep) {
  hit <- rep(FALSE, nrow(R)); best <- rep(NA_real_, nrow(R))
  for (g in groups) {
    v <- stat(R[, g, drop = FALSE])
    hit <- hit | (!is.na(v) & keep(v))
    best <- pmax(best, v, na.rm = TRUE)
  }
  list(hit = hit, best = best)
}

#' Relative promoter activity and the promoter filter of the BOHR paper
#'
#' Computes the relative activity of each external promoter, its absolute
#' activity over the sum of the absolute activities of all external promoters
#' of its gene, and, with `filter = TRUE`, keeps the promoters that pass the
#' paper's filter (SI Section 1).
#'
#' **Relative activity.** The denominator is the sum over all external
#' promoters of the gene, fixed before any promoter is filtered, so a
#' promoter dropped by the filter still counts in its gene's denominator.
#' Internal promoters are left out: proActiv does not quantify them. Where
#' the gene's sum is 0 the relative activity is `NA`.
#'
#' **Filter.** In this order (each step reported in `steps`):
#' \enumerate{
#'   \item `require_active`: a promoter must be active in at least one
#'     condition. For a `from_proactiv()` result these are proActiv's classes:
#'     not `"Inactive"` in at least one condition, and classified in every
#'     condition (as in the paper). For a matrix, active means a mean
#'     absolute activity of at least `min_abs` in at least one group of
#'     `group`, proActiv's rule.
#'   \item Promoters of genes with one annotated promoter (internal promoters
#'     included) are dropped.
#'   \item In at least one group (not necessarily the same one for each
#'     criterion), a promoter must have a share of samples with a relative
#'     activity above 0 greater than `min_nonzero`, a share of missing values
#'     less than `max_missing`, and a mean relative activity (over the samples
#'     with a value) greater than `min_mean`. The shares are over all samples
#'     of the group.
#'   \item Genes left with fewer than `min_promoters` promoters are dropped.
#' }
#' As in the paper, the criteria of step 3 are evaluated on proActiv's own
#' relative activity when `x` comes from [from_proactiv()] (its denominator
#' also holds the internal promoters), and on the relative activity computed
#' here when `x` is a matrix. The values returned are always the ones
#' computed here. Steps 2 and 4 together keep genes with at least two kept
#' promoters; step 2 only changes the counts reported in between.
#'
#' The defaults are the paper's values. With them, on the paper's proActiv
#' results (GTEx v8 RNA-seq, 22,865 samples, 57 tissues as `group`), the
#' 92,557 external promoters of 52,037 genes come down to 10,826 promoters in
#' 4,744 genes.
#'
#' @param x A [from_proactiv()] result, or a matrix (or the path of a TSV
#'   file, first column = promoter IDs) of absolute promoter activities,
#'   promoters in rows (row names = promoter IDs) and samples in columns.
#' @param gene For a matrix: the gene of each row, or a TSS table (a data
#'   frame or the path of a TSV file) with columns `id`, `gene` and
#'   optionally `internal`, such as `simulate_bohr()$tss`. Ignored for a
#'   `from_proactiv()` result.
#' @param internal For a matrix: `TRUE` for each internal promoter (default:
#'   none, or the table's column `internal`).
#' @param group `NULL` (all samples form one group), or the group (tissue)
#'   of each sample: a vector with one label per column (or named by sample
#'   ID), or a sample table with columns `sample` and `tissue`, as in
#'   [bohr_cis()].
#' @param filter `FALSE` computes the relative activity only.
#' @param min_nonzero,max_missing,min_mean The thresholds of step 3:
#'   defaults 0.25, 0.40 and 0.1.
#' @param require_active Apply step 1. Default `TRUE`.
#' @param min_abs For a matrix, the mean absolute activity needed to count
#'   as active in a group. Default 0.25 (proActiv).
#' @param min_promoters Genes with fewer kept promoters are dropped. Default 2.
#' @param out_dir Optional folder (created if needed); when given, the
#'   results are also written to it: `relative_activity.tsv` (promoter IDs in
#'   the first column), `promoters.tsv` and `steps.tsv`.
#' @return A list of class `"bohr_relative"`:
#'   \describe{
#'     \item{`relative`}{The relative activity of the kept promoters (all
#'       external promoters with `filter = FALSE`), promoters x samples, in
#'       the order of `x`: the `phenotype` of [bohr_cis()].}
#'     \item{`promoters`}{One row per promoter of the annotation: `id`,
#'       `gene`, `internal`, `kept`, `reason` (the first step that dropped it:
#'       `"internal"`, `"inactive"`, `"single_promoter_gene"`, `"nonzero"`,
#'       `"missing"`, `"mean"` or `"few_promoters"`), and the best value over
#'       the groups of each criterion: `max_nonzero`, `min_missing`,
#'       `max_mean`.}
#'     \item{`steps`}{The promoters and genes left after each step
#'       (`annotated`, `external`, then `active`, `multi_promoter_gene`,
#'       `nonzero`, `missing`, `mean` and `min_promoters`), and the genes
#'       with at least `min_promoters` of them (`genes_multi`).}
#'     \item{`settings`}{The arguments.}
#'   }
#' @seealso [from_proactiv()], [bohr_cis()].
#' @examples
#' set.seed(1)
#' sim <- simulate_bohr(n_individuals = 60, n_tissues = 2, n_genes = 4)
#' ra <- relative_activity(sim$absolute, gene = sim$tss, group = sim$samples)
#' ra$steps
#' head(ra$promoters)
#' ## without the filter, the relative activity of every external promoter
#' all_ext <- relative_activity(sim$absolute, gene = sim$tss, filter = FALSE)
#' dim(all_ext$relative)
#' @export
relative_activity <- function(x, gene = NULL, internal = NULL, group = NULL, filter = TRUE,
                              min_nonzero = 0.25, max_missing = 0.40, min_mean = 0.1,
                              require_active = TRUE, min_abs = 0.25, min_promoters = 2L, out_dir = NULL) {
  filter <- check_flag(filter, "filter"); require_active <- check_flag(require_active, "require_active")
  min_promoters <- check_count(min_promoters, "min_promoters", min = 1)
  for (a in c("min_nonzero", "max_missing", "min_mean", "min_abs")) {
    v <- get(a)
    if (!is.numeric(v) || length(v) != 1L || is.na(v)) stop("`", a, "` must be a single number.", call. = FALSE)
  }
  pa <- inherits(x, "bohr_proactiv")
  if (pa) {
    ann <- x$promoters[, c("id", "gene", "internal")]
    A <- x$absolute; Rc <- x$relative
    if (!identical(rownames(A), ann$id[!ann$internal]) || !identical(dimnames(A), dimnames(Rc)))
      stop("`x`: the matrices do not match each other or the external promoters of its table.", call. = FALSE)
  } else {
    A <- as_dense(read_matrix(x, "x"), "x")
    ann <- ra_annotation(A, gene, internal)
    A <- A[!ann$internal, , drop = FALSE]
    Rc <- NULL
  }
  samples <- colnames(A)
  if (is.null(samples) || anyDuplicated(samples)) stop("`x` needs unique column names (sample IDs).", call. = FALSE)
  ext <- ann[!ann$internal, ]
  ## relative activity over all external promoters of the gene
  gf <- factor(ext$gene, levels = unique(ext$gene))
  tot <- rowsum(A, gf, reorder = FALSE)
  R <- A / tot[as.integer(gf), , drop = FALSE]
  R[is.nan(R)] <- NA_real_
  dimnames(R) <- dimnames(A)

  count <- function(step, ids) {
    g <- table(ann$gene[match(ids, ann$id)])
    data.frame(step = step, promoters = length(ids), genes = length(g), genes_multi = sum(g >= min_promoters), stringsAsFactors = FALSE)
  }
  steps <- list(count("annotated", ann$id), count("external", ext$id))
  reason <- ifelse(ann$internal, "internal", NA_character_)
  out_tab <- data.frame(ann, kept = !ann$internal, reason = reason, max_nonzero = NA_real_, min_missing = NA_real_,
                        max_mean = NA_real_, stringsAsFactors = FALSE)
  rownames(out_tab) <- NULL
  si <- sample_info(group, NULL, samples)
  if (!filter) {
    res <- structure(list(relative = R, promoters = out_tab, steps = do.call(rbind, steps), settings = list(filter = FALSE)),
                     class = "bohr_relative")
    return(ra_write(res, out_dir))
  }
  tissue <- if (is.null(si$tissue)) rep("all", length(samples)) else si$tissue
  groups <- split(seq_along(samples), factor(tissue, levels = unique(tissue)))
  ei <- match(ext$id, out_tab$id)
  alive <- rep(TRUE, nrow(ext))
  drop_step <- function(fail, reason, step) {
    out_tab$reason[ei[alive & fail]] <<- reason
    alive <<- alive & !fail
    steps[[length(steps) + 1L]] <<- count(step, ext$id[alive])
  }
  ## 1. active
  if (require_active) {
    if (pa && !is.null(x$class)) {
      cl <- x$class[ext$id, , drop = FALSE]
      act <- stats::complete.cases(cl) & rowSums(cl != "Inactive", na.rm = TRUE) > 0
    } else {
      act <- any_group(A, groups, function(m) rowMeans(m, na.rm = TRUE), function(v) v >= min_abs)$hit
    }
    drop_step(!act, "inactive", "active")
  }
  ## 2. genes with one annotated promoter
  n_ann <- table(ann$gene)
  drop_step(ext$gene %in% names(n_ann)[n_ann == 1L], "single_promoter_gene", "multi_promoter_gene")
  ## 3. the criteria, on proActiv's relative activity when there is one
  C <- if (pa) Rc else R
  nz <- any_group(C, groups, function(m) rowSums(m > 0, na.rm = TRUE) / ncol(m), function(v) v > min_nonzero)
  ms <- any_group(C, groups, function(m) -rowSums(is.na(m)) / ncol(m), function(v) -v < max_missing)
  mn <- any_group(C, groups, function(m) { s <- rowMeans(m, na.rm = TRUE); s[is.nan(s)] <- NA; s }, function(v) v > min_mean)
  out_tab$max_nonzero[ei] <- nz$best; out_tab$min_missing[ei] <- -ms$best; out_tab$max_mean[ei] <- mn$best
  drop_step(!nz$hit, "nonzero", "nonzero")
  drop_step(!ms$hit, "missing", "missing")
  drop_step(!mn$hit, "mean", "mean")
  ## 4. genes left with too few promoters
  n_left <- table(ext$gene[alive])
  drop_step(!(ext$gene %in% names(n_left)[n_left >= min_promoters]), "few_promoters", "min_promoters")
  out_tab$kept[ei] <- alive
  res <- structure(list(relative = R[alive, , drop = FALSE], promoters = out_tab, steps = do.call(rbind, steps),
                        settings = list(filter = TRUE, min_nonzero = min_nonzero, max_missing = max_missing, min_mean = min_mean,
                                        require_active = require_active, min_abs = min_abs, min_promoters = min_promoters,
                                        criteria_on = if (pa) "proActiv relative activity" else "relative activity (external promoters)",
                                        active_from = if (!require_active) NA_character_ else if (pa && !is.null(x$class)) "proActiv classes" else "mean absolute activity",
                                        groups = names(groups))),
                   class = "bohr_relative")
  ra_write(res, out_dir)
}

ra_write <- function(res, out_dir) {
  if (is.null(out_dir)) return(res)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write_tsv(data.frame(id = rownames(res$relative), res$relative, check.names = FALSE), file.path(out_dir, "relative_activity.tsv"))
  write_tsv(res$promoters, file.path(out_dir, "promoters.tsv"))
  write_tsv(res$steps, file.path(out_dir, "steps.tsv"))
  res
}

#' @export
print.bohr_relative <- function(x, ...) {
  cat("Relative promoter activity:", nrow(x$relative), "promoters of", length(unique(x$promoters$gene[x$promoters$kept])),
      "genes in", ncol(x$relative), "samples\n")
  print(x$steps, row.names = FALSE)
  invisible(x)
}
