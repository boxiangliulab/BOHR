## Representative mean of each feature in each tissue (the input of the paper's tissue-specificity analysis).
##
## Follows raw_fit() of the paper's analysis code (code_orig/gamlss_functions_optimized.R, md5 d81d14fd..., lines 324-444;
## code map rows M6.1.1 and S6.1.3-S6.1.4): fewer than 10 values -> no mean (:330-334); all values equal -> sample mean (:335-346);
## a share of 0s or of 1s of at least 0.90, or SD <= 0.05 -> sample mean (:348-368); fit error -> sample mean (:371-397);
## otherwise the overall mean mu (1 - p0 - p1) + p1 of an intercept-only fit (:399-444), whose family follows the boundaries
## present (select_model(), :34-48). The paper fitted with gamlss (BEINF, BEZI, BEOI, BE; n.cyc = 200); here the fit is the
## BOHR engine's promoter_null() with an intercept-only design; its beta part is fitted from 2 values strictly between 0 and 1
## (promoter_null() itself needs 4 = columns + 3, a rule for regressions). With 0 or 1 such values the fitted overall mean is
## the sample mean (p0 and p1 are the shares of 0s and 1s, and mu is the one value), which is used. A fit that fails or does
## not converge also falls back to the sample mean (gamlss reported non-converged fits as "Success" and used them).

## one feature in one tissue: y = its values (NA allowed)
tm_one <- function(y, min_samples, max_boundary, min_sd, maxit) {
  y <- y[!is.na(y)]; n <- length(y); mid <- y > 0 & y < 1
  r <- list(mean = NA_real_, method = NA_character_, reason = NA_character_, family = NA_character_, n = n, n0 = sum(y == 0),
            n1 = sum(y == 1), n_mid = sum(mid), mu = NA_real_, phi = NA_real_, p0 = NA_real_, p1 = NA_real_, converged = NA)
  if (n < min_samples) { r$reason <- "too_few"; return(r) }
  r$p0 <- r$n0 / n; r$p1 <- r$n1 / n
  if (r$n_mid) r$mu <- mean(y[mid])
  sample_mean <- function(reason) { r$mean <- mean(y); r$method <- "sample_mean"; r$reason <- reason; r }
  if (all(y == y[1])) return(sample_mean("all_equal"))
  if (r$p0 >= max_boundary || r$p1 >= max_boundary) return(sample_mean("boundary"))
  if (stats::sd(y) <= min_sd) return(sample_mean("low_sd"))
  if (r$n_mid < 2) return(sample_mean("few_mid"))                        # = the fitted overall mean
  fam <- if (r$n0 > 0 && r$n1 > 0) "BEINF" else if (r$n0 > 0) "BEZI" else if (r$n1 > 0) "BEOI" else "BE"
  r$family <- fam
  C <- matrix(1, n, 1L, dimnames = list(NULL, "(Intercept)"))
  nul <- tryCatch({ f <- promoter_null(y, C, fam, maxit = maxit)
    if (is.null(f$beta)) f$beta <- beta_fit(y[mid], C[mid, , drop = FALSE], maxit = maxit)     # 2 or 3 values in between
    f }, error = function(e) NULL)
  ok <- !is.null(nul) && isTRUE(nul$beta$converged) && (fam == "BE" || isTRUE(nul$cat$converged))
  if (ok) {
    e <- nul$cat$coef
    p <- switch(fam, BEINF = { m <- max(0, e[1], e[2]); d <- exp(-m) + exp(e[1] - m) + exp(e[2] - m); c(exp(e[1] - m), exp(e[2] - m)) / d },
                BEZI = c(stats::plogis(e[1]), 0), BEOI = c(0, stats::plogis(e[1])), BE = c(0, 0))
    mu <- stats::plogis(nul$beta$coef[1]); m <- mu * (1 - p[1] - p[2]) + p[2]
    ok <- all(is.finite(c(p, mu, m)))
  }
  if (!ok) { r <- sample_mean("fit_failed"); r$converged <- FALSE; return(r) }
  r$mu <- unname(mu); r$phi <- unname(exp(nul$beta$coef[2])); r$p0 <- unname(p[1]); r$p1 <- unname(p[2])
  r$mean <- unname(m); r$method <- "fit"; r$converged <- TRUE
  r
}

#' Representative mean of each feature in each tissue
#'
#' Summarises each feature (for example a promoter's relative activity) in
#' each tissue by one mean, the input of tissue-specificity methods such as
#' ROKU (package TCC). The mean is the overall mean
#' \eqn{\mu (1 - p_0 - p_1) + p_1} of an intercept-only zero-and-one-inflated
#' beta fit, or the sample mean where such a fit is not useful, as in the
#' tissue-specificity analysis of the BOHR paper (SI Section 6).
#'
#' For each feature and tissue, over the samples with a value:
#' \enumerate{
#'   \item fewer than `min_samples` values: no mean (`NA`, reason
#'     `"too_few"`);
#'   \item all values equal: the sample mean (`"all_equal"`);
#'   \item a share of 0s or of 1s of at least `max_boundary`
#'     (`"boundary"`), or a standard deviation of `min_sd` or less
#'     (`"low_sd"`): the sample mean;
#'   \item otherwise the model is fitted with the BOHR engine without
#'     covariates and without the fitting rules of [bohr_control()] (of which
#'     only `maxit` is used): the family follows the boundaries present (BEINF
#'     with 0s and 1s, BEZI with 0s only, BEOI with 1s only, BE with neither);
#'     mu is the mean of the beta part, fitted to the values strictly between
#'     0 and 1, and p0 and p1 are the probabilities of 0 and 1. The mean is
#'     \eqn{\mu (1 - p_0 - p_1) + p_1} (`method = "fit"`). With fewer than 2
#'     values strictly between 0 and 1 the fitted mean is the sample mean,
#'     which is used (`"few_mid"`); when the fit fails or does not converge,
#'     the sample mean is used too (`"fit_failed"`).
#' }
#' Without covariates, p0 and p1 are the shares of 0s and 1s, so the fitted
#' mean differs from the sample mean only through mu: the maximum-likelihood
#' mean of a beta distribution instead of the arithmetic mean of the values
#' between 0 and 1.
#'
#' **The paper used gamlss.** The tissue means of the BOHR paper were fitted
#' with the gamlss package (families BEINF, BEZI, BEOI and BE, intercept only,
#' at most 200 cycles), with the same rules and the same mean. The package
#' fits the same model with its own engine, so the numbers differ slightly:
#' gamlss stops its iterations at a looser tolerance, which leaves its p0 off
#' the share of 0s by up to about 0.002, and it also fitted pairs with 0 or 1
#' values between 0 and 1 (where the exact fitted mean is the sample mean,
#' used here) and kept fits that had not converged. On the paper's data the
#' two agree on which pairs have too few values and on every sample-mean rule, and
#' their means agree to about 0.001 wherever both report a converged fit.
#'
#' **Tissue specificity.** The paper applied ROKU to the feature x tissue
#' matrix of these means (`TCC::ROKU(x, upper.limit = 0.25, sort = TRUE)`)
#' and kept a call when the tissue's mean differed by at least 0.10 from the
#' median of the tissues without a call; see the vignette.
#'
#' @param phenotype Proportions in \[0, 1\] (for example relative
#'   activities), features in rows (row names = feature IDs) and samples in
#'   columns (unique column names); a numeric matrix or the path of a TSV file
#'   (first column = feature IDs). `NA` marks a missing value.
#' @param tissue `NULL` (all samples form one group, labelled `"all"`), or the
#'   tissue of each sample: a vector with one label per column of `phenotype`
#'   (or named by sample ID), or a sample table (a data frame or the path of a
#'   TSV file) with columns `sample` and `tissue`, as in [bohr_cis()].
#' @param min_samples A feature gets no mean in a tissue with fewer values
#'   than this. Default 10 (SI Section 6).
#' @param max_boundary The sample mean is used when the share of 0s, or of 1s,
#'   is at least this. Default 0.90 (SI Section 6).
#' @param min_sd The sample mean is used when the standard deviation is this
#'   or less. Default 0.05 (SI Section 6).
#' @param control [bohr_control()]; only `maxit` (Newton steps of each fit) is
#'   used.
#' @param ncores Number of cores (forked processes; 1 on Windows).
#' @param out_dir Optional folder (created if needed); when given, the table
#'   is also written to it as `tissue_means.tsv`.
#' @return A data frame (a tidy table) with one row per feature and tissue:
#'   `feature`, `tissue`, `mean` (the representative mean; `NA` when there
#'   are too few values), `method` (`"fit"` or `"sample_mean"`), `reason` (why
#'   the sample mean was used, or `"too_few"`), `family`, `n`, `n0`, `n1`,
#'   `n_mid` (values at 0, at 1 and in between), `mu`, `phi` (the fitted beta
#'   mean and precision; for the sample mean, `mu` is the mean of the values
#'   in between), `p0`, `p1` and `converged`. Rows are ordered by tissue, then
#'   by feature as in `phenotype`.
#' @seealso [bohr_cis()], [bohr_trait()]; the vignette shows the call of
#'   ROKU on the result.
#' @examples
#' sim <- simulate_bohr(n_individuals = 60, n_tissues = 2, n_genes = 3)
#' tm <- tissue_means(sim$relative, tissue = sim$samples)
#' head(tm)
#' table(tm$tissue, tm$method)
#' ## features x tissues, the input of TCC::ROKU()
#' M <- tapply(tm$mean, list(tm$feature, tm$tissue), c)
#' head(M)
#' @export
tissue_means <- function(phenotype, tissue = NULL, min_samples = 10L, max_boundary = 0.90, min_sd = 0.05,
                         control = bohr_control(), ncores = 1L, out_dir = NULL) {
  check_control(control)
  min_samples <- check_count(min_samples, "min_samples", min = 1)
  if (!is.numeric(max_boundary) || length(max_boundary) != 1L || is.na(max_boundary) || max_boundary <= 0 || max_boundary > 1)
    stop("`max_boundary` must be a single number in (0, 1].", call. = FALSE)
  if (!is.numeric(min_sd) || length(min_sd) != 1L || is.na(min_sd) || min_sd < 0) stop("`min_sd` must be a single number >= 0.", call. = FALSE)
  ncores <- check_count(ncores, "ncores", min = 1)
  if (!is.null(out_dir) && (!is.character(out_dir) || length(out_dir) != 1L || is.na(out_dir) || !nzchar(out_dir)))
    stop("`out_dir` must be the path of one folder.", call. = FALSE)
  phenotype <- check_phenotype(read_matrix(phenotype, "phenotype"))
  samples <- colnames(phenotype); feats <- rownames(phenotype)
  st <- sample_info(tissue, NULL, samples)
  lev <- if (is.null(st$tissue)) "all" else st$levels
  cols <- split(seq_along(samples), factor(if (is.null(st$tissue)) rep("all", length(samples)) else st$tissue, levels = lev))
  chunks <- split(seq_along(feats), ceiling(seq_along(feats) / max(1L, ceiling(length(feats) / (4L * ncores)))))
  res <- par_apply(chunks, function(ix) lapply(ix, function(i) lapply(cols, function(j)
    tm_one(phenotype[i, j], min_samples, max_boundary, min_sd, control$maxit))), ncores, preschedule = FALSE)
  res <- unlist(res, recursive = FALSE)                                   # features, each a list over tissues
  rows <- unlist(lapply(seq_along(lev), function(t) lapply(seq_along(feats), function(i) res[[i]][[t]])), recursive = FALSE)
  out <- rows_to_df(rows)
  out <- data.frame(feature = rep(feats, length(lev)), tissue = rep(lev, each = length(feats)), out, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(out_dir)) stop("Cannot create `out_dir`: ", out_dir, call. = FALSE)
    write_tsv(out, file.path(out_dir, "tissue_means.tsv"))
  }
  out
}
