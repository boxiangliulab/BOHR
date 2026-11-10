## Numeric part of a promoter ID for the tie-break of mirror_pairs.R:25 (num(): "PromoterID_5712" -> 5712).
## General IDs: the trailing digits when every ID has them, otherwise the ID itself (C-locale order).
id_order <- function(ids) {
  if (all(grepl("[0-9]$", ids))) order(as.numeric(sub("^.*[^0-9]", "", ids)), ids, method = "radix") else order(ids, method = "radix")
}

#' Find exact mirror pairs of features
#'
#' When a group has exactly two features whose proportions are shares of the
#' same total (for example the two promoters of a gene, with relative
#' activities that sum to 1 in every sample), their association tests describe
#' the same events: p0 of one feature is p1 of the other, and the two mu tests
#' are the same test with opposite effect directions. This function finds such
#' pairs so that [bohr_cis()] can count each event once.
#'
#' Following the BOHR paper (rules of `mirror_pairs.R` in the analysis code;
#' SI Section 2), for each group with exactly two features among the rows of
#' `phenotype`, using the samples where both are observed:
#' \describe{
#'   \item{`boundary_exact`}{One feature is at 0 exactly when the other is at
#'     1, and vice versa, in every sample. Then p0 of one feature and p1 of
#'     the other test the same event, even when the two values do not sum to
#'     1 (for example when other, untested promoters of the gene stay in the
#'     denominator).}
#'   \item{`interior_exact`}{`boundary_exact`, the same samples are
#'     intermediate for both, and the two values sum to 1 within `tol` in
#'     each of them. Then the two mu tests are the same test.}
#' }
#' The canonical member of a pair is the feature with the higher
#' `major_score` (by default its mean in `phenotype`); ties go to the feature
#' with the smaller identifier number. The paper used the mean relative
#' activity over all analysed samples of all 50 tissues. In [bohr_cis()], the
#' canonical member decides only which of the two mu tests of an
#' `interior_exact` pair is kept.
#'
#' @param phenotype Numeric matrix of proportions, features in rows (row names
#'   = feature IDs) and samples in columns (same samples in every row; missing
#'   values allowed). Should hold all features of each group, so that
#'   two-feature groups are recognised.
#' @param group Group of each feature (for promoters, the gene): a vector with
#'   one value per row of `phenotype`, or a vector named by feature ID.
#' @param major_score Optional numeric vector named by feature ID used to
#'   choose the canonical member (higher wins). Defaults to the row means of
#'   `phenotype`.
#' @param tol Tolerance of \eqn{y_1 + y_2 = 1} on the intermediate samples of an
#'   `interior_exact` pair. Default 1e-6 (SI Section 2).
#' @return A data frame with one row per group that has exactly two features:
#'   `group`, `feature_1` and `feature_2` (ordered by identifier number),
#'   `n_samples` (samples where both are observed), `boundary_exact`,
#'   `interior_exact`, `exact_mirror` (both), `max_dev` (largest
#'   \eqn{|y_1 + y_2 - 1|}), `n_both0` (samples where both are 0), `mean_1`,
#'   `mean_2`, `score_1`, `score_2`, `canonical` and `partner`.
#' @seealso [bohr_cis()].
#' @examples
#' set.seed(3)
#' a <- simulate_response(60, p0 = 0.2, p1 = 0.2)$y
#' b <- runif(60)
#' ph <- rbind(P1 = a, P2 = 1 - a, P3 = b, P4 = 1 - b, P5 = runif(60))
#' colnames(ph) <- paste0("S", 1:60)
#' mirror_pairs(ph, group = c("G1", "G1", "G2", "G2", "G3"))
#' @export
mirror_pairs <- function(phenotype, group, major_score = NULL, tol = 1e-6) {
  phenotype <- check_phenotype(phenotype, sample_names = FALSE)
  group <- align_gene(group, rownames(phenotype))
  if (!is.numeric(tol) || length(tol) != 1L || !(tol >= 0)) stop("`tol` must be a single number >= 0.", call. = FALSE)
  if (is.null(major_score)) major_score <- rowMeans(phenotype, na.rm = TRUE)
  else {
    if (!is.numeric(major_score) || is.null(names(major_score))) stop("`major_score` must be a numeric vector named by feature ID.", call. = FALSE)
  }
  pair_table(phenotype, group, major_score, tol)
}

## mirror_pairs.R:26-42 for one tissue matrix (promoters in rows here, in columns there).
pair_table <- function(phenotype, gene, major_score, tol = 1e-6, genes = NULL) {
  ok_g <- !is.na(gene)
  np <- table(gene[ok_g])
  two <- names(np)[np == 2]
  if (!is.null(genes)) two <- intersect(two, genes)
  cols <- c("group", "feature_1", "feature_2", "n_samples", "boundary_exact", "interior_exact", "exact_mirror",
            "max_dev", "n_both0", "mean_1", "mean_2", "score_1", "score_2", "canonical", "partner")
  if (!length(two)) {
    out <- data.frame(group = character(0), feature_1 = character(0), feature_2 = character(0), n_samples = integer(0),
                      boundary_exact = logical(0), interior_exact = logical(0), exact_mirror = logical(0), max_dev = numeric(0),
                      n_both0 = integer(0), mean_1 = numeric(0), mean_2 = numeric(0), score_1 = numeric(0), score_2 = numeric(0),
                      canonical = character(0), partner = character(0), stringsAsFactors = FALSE)
    return(out)
  }
  ids <- rownames(phenotype)
  rows <- lapply(two, function(x) {
    pr <- ids[ok_g & gene == x]; pr <- pr[id_order(pr)]; a <- phenotype[pr[1], ]; b <- phenotype[pr[2], ]
    ok <- !is.na(a) & !is.na(b); d <- abs(a[ok] + b[ok] - 1); aa <- a[ok]; bb <- b[ok]
    za <- aa == 0; oa <- aa == 1; zb <- bb == 0; ob <- bb == 1; ia <- !za & !oa; ib <- !zb & !ob
    bex <- any(ok) && all(za == ob) && all(oa == zb)
    iex <- bex && all(ia == ib) && (!any(ia) || max(abs(aa[ia] + bb[ia] - 1)) <= tol)
    s1 <- unname(major_score[pr[1]]); s2 <- unname(major_score[pr[2]])
    if (is.na(s1) || is.na(s2)) stop("`major_score` is missing for feature ", if (is.na(s1)) pr[1] else pr[2], ".", call. = FALSE)
    canon <- if (s1 >= s2) pr[1] else pr[2]
    list(group = x, feature_1 = pr[1], feature_2 = pr[2], n_samples = sum(ok), boundary_exact = bex, interior_exact = iex,
         exact_mirror = iex && bex, max_dev = if (any(ok)) max(d) else NA_real_, n_both0 = sum(a[ok] == 0 & b[ok] == 0),
         mean_1 = mean(a, na.rm = TRUE), mean_2 = mean(b, na.rm = TRUE), score_1 = s1, score_2 = s2,
         canonical = canon, partner = if (canon == pr[1]) pr[2] else pr[1])
  })
  rows_to_df(rows, cols)
}

## per-promoter flags (mirror_pairs.R:44-50) and the v3m drop rule (units_v3m_mirror.R:23-26):
## drop_mu = interior_exact & partner; drop_p1 = boundary_exact (both members).
mirror_flags <- function(pairs) {
  if (is.null(pairs) || !nrow(pairs)) return(NULL)
  f <- function(role, me, other) data.frame(promoter = pairs[[me]], mirror_pair_id = pairs$group, boundary_exact = pairs$boundary_exact,
                                            interior_exact = pairs$interior_exact, mirror_role = role, mirror_partner = pairs[[other]],
                                            stringsAsFactors = FALSE)
  P <- rbind(f("canonical", "canonical", "partner"), f("partner", "partner", "canonical"))
  P$drop_mu <- P$interior_exact & P$mirror_role == "partner"
  P$drop_p1 <- P$boundary_exact
  P
}

## gene vector aligned with the promoter IDs
align_gene <- function(gene, ids) {
  if (!is.null(names(gene)) && all(ids %in% names(gene))) gene <- gene[ids]
  if (length(gene) != length(ids)) stop("`group` must have one value per feature (or be named by feature ID).", call. = FALSE)
  as.character(unname(gene))
}
