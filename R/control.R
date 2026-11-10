#' Settings of the BOHR model fits
#'
#' Collects the fitting rules and optimiser limits of the BOHR model of a
#' proportion response with exact 0s and 1s. Every
#' default is the value used in the BOHR paper; the section of the paper's
#' Supplementary Methods (SI) or Methods that states it is given for each.
#' [bohr_fit()], [bohr_test()] and [bohr_cis()] take the result as `control`.
#'
#' @param min_samples A response is skipped when fewer observations than this
#'   have a non-missing value, or fewer observations than covariates. Default
#'   20 (SI Section 2).
#' @param max_boundary A response is skipped when at least this share of its
#'   observations lies at 0, or at 1. Default 0.95 (SI Section 2).
#' @param min_sd A response is skipped when its standard deviation is this
#'   value or less. Default 0.05 (SI Section 2).
#' @param endpoint Before fitting, the observations at a boundary (0 or 1) are
#'   removed when that boundary holds this share of the observations or fewer;
#'   the family then follows from the boundaries that remain. Default 0.05 (SI
#'   Section 2).
#' @param maxit Maximum number of Newton (Fisher-scoring) iterations of each
#'   fit, unpenalised and penalised. Default 200 (SI Section 2; Methods Section
#'   2).
#' @param bfgs_maxit Maximum number of BFGS iterations that continue a
#'   penalised fit not converged after `maxit` Newton steps; the BFGS result is
#'   kept only if it does not lower the penalised log-likelihood. Default 1000
#'   (Methods Section 2; SI Section 2 describes the step).
#' @param max_rounds Maximum number of rounds in which the penalised fits with
#'   and without the predictor are restarted from each other's estimates,
#'   until neither improves. Default 6 (SI Section 2).
#' @param firth Test p0 and p1 with Firth-penalised likelihood-ratio tests.
#'   Default `TRUE` (SI Section 2). With `FALSE`, the unpenalised
#'   likelihood-ratio tests (and estimates) are reported for p0 and p1
#'   instead; these are not finite under separation.
#' @param firth_null_start Also start each restricted penalised fit (the
#'   tested coefficient held at 0) from the fit without the predictor, when
#'   that fit converged, and keep the start that reaches the highest
#'   penalised log-likelihood. The penalised likelihood can have more than one
#'   local maximum. Near separation, the other two starts (the full fit with
#'   the tested coefficient set to 0, and the default start) can stop at a
#'   lower maximum of the restricted fit, which makes the likelihood-ratio
#'   statistic too large. `TRUE` (the default) finds the higher maximum.
#'   `FALSE` reproduces the analysis engine of the BOHR paper exactly. The
#'   paper's published conditional results include one refit that reached the
#'   higher maximum, Brain_Putamen PromoterID_47842 p0 (one signal); with
#'   `TRUE` the package reproduces it.
#' @return A list of class `"bohr_control"`.
#' @seealso [bohr_fit()], [bohr_test()], [bohr_cis()].
#' @examples
#' bohr_control()
#' bohr_control(maxit = 500)
#' bohr_control(firth_null_start = FALSE)   # the paper's engine exactly
#' @export
bohr_control <- function(min_samples = 20L, max_boundary = 0.95, min_sd = 0.05, endpoint = 0.05,
                         maxit = 200L, bfgs_maxit = 1000L, max_rounds = 6L, firth = TRUE, firth_null_start = TRUE) {
  share <- function(v, what) if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < 0 || v > 1) stop("`", what, "` must be a single number in [0, 1].", call. = FALSE)
  share(max_boundary, "max_boundary"); share(endpoint, "endpoint")
  if (!is.numeric(min_sd) || length(min_sd) != 1L || is.na(min_sd) || min_sd < 0) stop("`min_sd` must be a single number >= 0.", call. = FALSE)
  structure(list(min_samples = check_count(min_samples, "min_samples"), max_boundary = max_boundary, min_sd = min_sd, endpoint = endpoint,
                 maxit = check_count(maxit, "maxit", min = 1), bfgs_maxit = check_count(bfgs_maxit, "bfgs_maxit"),
                 max_rounds = check_count(max_rounds, "max_rounds", min = 1), firth = check_flag(firth, "firth"),
                 firth_null_start = check_flag(firth_null_start, "firth_null_start")),
            class = "bohr_control")
}

#' @export
print.bohr_control <- function(x, ...) {
  cat("BOHR settings\n")
  cat(sprintf("  skip: fewer than %d samples, >= %g at one boundary, or SD <= %g\n", x$min_samples, x$max_boundary, x$min_sd))
  cat(sprintf("  endpoint rule: remove a boundary holding <= %g of samples\n", x$endpoint))
  cat(sprintf("  optimiser: %d Newton steps; BFGS up to %d; at most %d restart rounds\n", x$maxit, x$bfgs_maxit, x$max_rounds))
  cat(sprintf("  p0/p1 tests: %s\n", if (!x$firth) "unpenalised likelihood ratio" else
    paste0("Firth-penalised likelihood ratio", if (isTRUE(x$firth_null_start)) "; restricted fits also started from the null fit"
           else " (the paper's engine: no null start)")))
  invisible(x)
}

check_control <- function(control) {
  if (!inherits(control, "bohr_control")) stop("`control` must come from bohr_control().", call. = FALSE)
  control
}
