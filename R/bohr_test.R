## Coding of the predictor for the reported effects. The engine reports each effect as the prediction at
## x' = 2 minus x' = 0 (bohr2_engine.R:165-166, 181, 208-213). A contrast c(from, to) on the scale of x is
## turned into that form by x' = 2 (x - from) / (to - from). The defaults reproduce the paper:
##   dosage      c(0, 2):            x' = x (no arithmetic; identical to the cis-mapping engine)
##   binary      c(lo, hi):          x' = 2 (x - lo) / (hi - lo); for 0/1 coding x' = 2 x (= trait_contrast(),
##                                   run_trait data.R:100-102)
##   continuous  c(m, m + s):        x' = 2 (x - m) / s          (= trait_contrast(); m, s = mean, sd of x)
## Returns the coded predictor, the contrast, the factor that converts coefficients on x' to per unit of x
## (per_unit), and the offset of the coding (x' = offset + per_unit * x), used to map coefficients back.
## type = "auto" (coordinator decision 2026-11-10): exactly two distinct values -> binary, anything else ->
## continuous. A dosage is never detected; the user sets type = "dosage" (bohr_cis() always codes dosages).
detect_type <- function(x) {
  v <- unique(x[!is.na(x)])
  if (length(v) == 2L) "binary" else "continuous"
}
code_predictor <- function(x, type, contrast = NULL) {
  if (is.null(contrast)) {
    if (type == "dosage") return(list(x = x, contrast = c(0, 2), per_unit = 1, offset = 0))
    if (type == "binary") {
      v <- sort(unique(x[!is.na(x)]))
      if (length(v) > 2L) stop("`type = \"binary\"` needs a predictor with exactly two distinct values; this one has ",
                               length(v), ".", call. = FALSE)
      if (length(v) < 2L) return(NULL)
      contrast <- v
    } else {
      m <- mean(x, na.rm = TRUE); s <- stats::sd(x, na.rm = TRUE)
      if (!is.finite(s) || s <= 0) return(NULL)
      return(list(x = 2 * (x - m) / s, contrast = c(m, m + s), per_unit = 2 / s, offset = -2 * m / s))
    }
  }
  if (identical(as.numeric(contrast), c(0, 2))) return(list(x = x, contrast = c(0, 2), per_unit = 1, offset = 0))
  list(x = 2 * (x - contrast[1]) / (contrast[2] - contrast[1]), contrast = as.numeric(contrast),
       per_unit = 2 / (contrast[2] - contrast[1]), offset = -2 * contrast[1] / (contrast[2] - contrast[1]))
}

#' Test a predictor for association with a proportion response
#'
#' Tests whether a numeric predictor `x` is associated with each component of
#' a response `y` in \[0, 1\] that has exact 0s and 1s: the probability of 0
#' (p0), the probability of 1 (p1) and the mean of the values in between (mu).
#' The predictor can be anything numeric, for example an allele dosage, a
#' binary exposure or a continuous measurement.
#'
#' For each predictor, the model of [bohr_fit()] is refitted with the
#' predictor added, starting from the fit without it, and each component is
#' tested with one degree of freedom:
#' \describe{
#'   \item{mu}{Likelihood-ratio statistic \eqn{W} of the beta part, with a
#'     small-sample correction: \eqn{F = (n_{mid} - p)(e^{W/n_{mid}} - 1)} is
#'     compared with an \eqn{F(1, n_{mid} - p)} distribution, where
#'     \eqn{n_{mid}} is the number of intermediate observations and \eqn{p}
#'     the rank of the design on them plus 1 (the predictor). `NA` when the
#'     beta part cannot be fitted or when \eqn{n_{mid} \le p + 1}.}
#'   \item{p0, p1}{Firth-penalised likelihood-ratio test of the predictor's
#'     coefficient for \eqn{\nu} (p0) or \eqn{\tau} (p1) against
#'     \eqn{\chi^2_1}. Full and restricted penalised fits are restarted from
#'     each other, for at most `control$max_rounds` rounds, until neither
#'     improves; a penalised fit that has not converged after
#'     `control$maxit` Newton steps is continued with BFGS, and the BFGS
#'     result is kept only if it does not lower the penalised
#'     log-likelihood.}
#' }
#' **Effects** are on the scale of the component: for each observation the
#' value is predicted at `x = to` and at `x = from`, with that observation's
#' covariates, and the difference is averaged over observations. For mu it is
#' the change in the mean among intermediate values; for p0 and p1 it is the
#' change in probability, from the unpenalised fit. The default contrast
#' depends on `type` (as in the BOHR paper): `"dosage"`, 2 versus 0 (two
#' copies of an allele versus none); `"binary"`, the larger of the two values
#' versus the smaller (1 versus 0 for a 0/1 predictor); `"continuous"`, one
#' standard deviation above the mean versus the mean (mean and standard
#' deviation over the observations with a non-missing `x`). A `contrast`
#' overrides the default of any type.
#'
#' **Predictor type.** `type = "auto"` (the default) chooses only between two
#' types: a predictor with exactly two distinct values is `"binary"`, any other
#' is `"continuous"`. An allele dosage is never detected automatically, because
#' a dosage vector without one of the three values would look like another
#' type: set `type = "dosage"` for allele dosages. [bohr_cis()] always codes
#' its variants as dosages. Internally the predictor is coded as
#' \eqn{x' = 2(x - from)/(to - from)}, so that the contrast is \eqn{x'} from 0
#' to 2; for a dosage (contrast 0 to 2) \eqn{x' = x}, the coding of the paper's
#' cis mapping.
#'
#' **Scales.** Coefficients and standard errors (`est_*`, `se_*`) are reported
#' per unit of `x`, as in [stats::glm()], whatever the type and contrast.
#' Effects (`effect_*`) are for the contrast. P values do not depend on the
#' type or the contrast.
#'
#' @param y Numeric vector, the response in \[0, 1\]; `NA` allowed.
#' @param x Numeric predictor: a vector with one value per element of `y`, or a
#'   matrix with one row per element of `y` and one column per predictor
#'   (each tested separately). For a vector `x`, observations where `x` is `NA`
#'   are left out (before the fit without `x`). A matrix `x` must have no
#'   missing values on the observations used.
#' @param covariates `NULL`, or a numeric matrix or data frame with one row per
#'   element of `y` and no intercept column.
#' @param fit Optional fit from [bohr_fit()] on the same `y` and `covariates`,
#'   to reuse it across many predictors.
#' @param type How `x` is treated for the effect contrast: `"auto"` (default;
#'   `"binary"` for exactly two distinct values, otherwise `"continuous"`),
#'   `"dosage"`, `"binary"` or `"continuous"` (see Details). Set `"dosage"` for
#'   allele dosages; it is never chosen automatically.
#' @param contrast Optional length-2 numeric vector `c(from, to)` on the scale of
#'   `x`; the effects compare `x = to` with `x = from`. Overrides the default of
#'   `type`.
#' @param control Fitting rules, optimiser limits and the choice of the p0/p1
#'   test, from [bohr_control()]; the defaults are the values of the BOHR
#'   paper. When `fit` is supplied, its fitting rules were already applied and
#'   `control` sets only the tests.
#' @param all_stats Logical. Also return every statistic of the analysis
#'   engine under its original name (see Details of [bohr_cis()]); these are on
#'   the scale of the coded predictor \eqn{x'}.
#' @return A data frame with one row per predictor:
#'   \describe{
#'     \item{`predictor`, `type`, `contrast_from`, `contrast_to`}{Name of the
#'       predictor (column name of `x`, or `"x"`), how it was treated, and the
#'       contrast of the effects.}
#'     \item{`family`, `n`, `n_mid`, `n0`, `n1`}{From the fit without `x`.}
#'     \item{`log_p_mu`, `log_p_p0`, `log_p_p1`}{Natural logarithm of the
#'       P values: mu F-corrected likelihood-ratio test, p0 and p1
#'       Firth-penalised likelihood-ratio tests (unpenalised with
#'       `bohr_control(firth = FALSE)`). `NA` when a component does not exist
#'       in the family or cannot be tested.}
#'     \item{`p_mu`, `p_p0`, `p_p1`}{The same P values, `exp(log_p_*)`.}
#'     \item{`effect_mu`, `effect_p0`, `effect_p1`}{Effects of the contrast
#'       (see Details). For `"BEZI"`, `effect_p1` is 0; for `"BEOI"`,
#'       `effect_p0` is 0.}
#'     \item{`est_mu`, `se_mu`}{Coefficient of `x`, per unit of `x`, in the
#'       logit of mu, and its standard error.}
#'     \item{`est_p0`, `se_p0`, `est_p1`, `se_p1`}{Firth-penalised coefficients
#'       of `x`, per unit of `x`, for \eqn{\log\nu} (p0) and \eqn{\log\tau} (p1),
#'       the log odds of 0 and of 1 against an intermediate value (for
#'       `"BEZI"` and `"BEOI"`, the logit of being at the boundary), and their
#'       standard errors.}
#'     \item{`stat_mu`, `stat_p0`, `stat_p1`}{Likelihood-ratio statistic
#'       \eqn{W} of mu, and the penalised likelihood-ratio statistics.}
#'     \item{`conv_mu`, `conv_p0`, `conv_p1`}{Convergence of the beta fit and of
#'       the penalised fits.}
#'     \item{`sep_mu`, `sep_p0`, `sep_p1`}{Separation flags of the unpenalised
#'       fits: not converged, |estimate| > 15 or standard error > 50 (on the
#'       scale of \eqn{x'}).}
#'     \item{`complete_sep`}{`TRUE` when the covariates alone separate the
#'       categories completely (the unpenalised boundary tests are then not
#'       defined; the Firth tests still are).}
#'     \item{`error`, `firth_error`}{Error message of the unpenalised fits
#'       (`NA` if none) and whether the penalised fits failed.}
#'   }
#' @seealso [bohr_fit()], [bohr_control()], [simulate_bohr()]; [bohr_cis()]
#'   for genome-wide cis mapping.
#' @examples
#' set.seed(2)
#' d <- simulate_response(200, mu = 0.4, p0 = 0.25, p1 = 0.1,
#'                        effect_mu = 0.05, effect_p0 = -0.08)
#' ## an allele dosage: set the type (effects for 2 copies versus 0)
#' bohr_test(d$y, d$x, type = "dosage")[, c("type", "family", "p_mu", "p_p0",
#'                                          "p_p1", "effect_mu", "effect_p0")]
#' ## a continuous predictor: effects for +1 standard deviation
#' z <- rnorm(200)
#' bohr_test(d$y, z)[, c("type", "contrast_from", "contrast_to", "p_mu", "effect_mu")]
#' @export
bohr_test <- function(y, x, covariates = NULL, fit = NULL, type = c("auto", "dosage", "binary", "continuous"),
                      contrast = NULL, control = bohr_control(), all_stats = FALSE) {
  check_activity(y); check_control(control)
  check_flag(all_stats, "all_stats")
  type <- match.arg(type)
  if (!is.null(contrast) && (!is.numeric(contrast) || length(contrast) != 2L || anyNA(contrast) || !all(is.finite(contrast)) ||
                             contrast[1] == contrast[2])) stop("`contrast` must be two different finite numbers c(from, to).", call. = FALSE)
  y <- as.vector(y)
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  X <- if (is.null(dim(x))) matrix(x, ncol = 1L, dimnames = list(NULL, "x")) else x
  if (nrow(X) != length(y)) stop("`x` must have one value (row) per element of `y`.", call. = FALSE)
  if (is.null(colnames(X))) colnames(X) <- if (ncol(X) == 1L) "x" else paste0("x", seq_len(ncol(X)))
  if (any(is.infinite(X))) stop("`x` must not contain infinite values.", call. = FALSE)
  if (is.null(fit)) {
    if (is.null(dim(x)) && anyNA(X)) y[is.na(X[, 1])] <- NA
    fit <- bohr_fit(y, covariates, control)
  } else {
    if (!inherits(fit, "bohr_fit")) stop("`fit` must come from bohr_fit().", call. = FALSE)
    if (is.na(fit$skip) && (max(fit$samples) > length(y) || !isTRUE(all(fit$y == y[fit$samples]))))
      stop("`fit` does not match `y`; refit with bohr_fit(y, covariates).", call. = FALSE)
  }
  rows <- lapply(seq_len(ncol(X)), function(j) test_column(fit, X[, j], colnames(X)[j], type, contrast, control, all_stats)$row)
  out <- rows_to_df(rows)
  rownames(out) <- NULL
  out
}

## One predictor (column xj of length(y)) against a bohr_fit: the row of bohr_test(); with keep_fits = TRUE also
## the engine's full fits and the predictor coding (used by bohr()).
test_column <- function(fit, xj, name, type, contrast, control, all_stats = FALSE, keep_fits = FALSE) {
  tj <- if (type == "auto") detect_type(xj) else type
  cp <- code_predictor(xj, tj, contrast)
  e <- stats::setNames(rep(list(NA), length(ENGINE_COLS)), ENGINE_COLS); e$firth_error <- FALSE; scale <- 1
  if (!is.na(fit$skip)) e$error <- paste("skipped:", fit$skip)
  else if (is.null(cp)) e$error <- "x is constant"
  else {
    g <- cp$x[fit$samples]
    if (anyNA(g)) stop("`x` has missing values on observations used by the fit (column ", name, ").", call. = FALSE)
    if (length(unique(g)) < 2L) e$error <- "x is constant on the observations used"
    else { e <- engine_one(fit, unname(g), control, keep_fits = keep_fits); scale <- cp$per_unit }
  }
  ctr <- if (is.null(cp)) c(NA_real_, NA_real_) else cp$contrast
  row <- c(list(predictor = name, type = tj, contrast_from = ctr[1], contrast_to = ctr[2]),
           list(family = fit$family, n = fit$n, n_mid = fit$n_mid, n0 = fit$n0, n1 = fit$n1),
           engine_to_row(e, all_stats, control$firth, per_unit = scale))
  list(row = row, fits = if (keep_fits) e$fits else NULL, coding = cp, type = tj)
}
