## Internal: the per-feature preparation of run_final2.R:44-58 (prep_promoter) and run_muF.R:27-30;
## control = bohr_control() (its defaults are the original constants).
## yv: activities over the analysis samples, named by (unique) sample IDs; C_all: the standardised design
## (cov_matrix) over the same samples, rownames = names(yv).
fit_promoter <- function(yv, C_all, control) {
  out <- list(family = NA_character_, skip = NA_character_, n = NA_integer_, n0 = NA_integer_, n1 = NA_integer_,
              n_mid = NA_integer_, n_cov = NA_integer_, n_cov_dropped = NA_integer_, rank_mid = NA_integer_,
              p_mean = NA_integer_, samples = integer(0), y = numeric(0), design = NULL, null = NULL)
  z <- prepare_y(yv, ncol(C_all) - 1L, min_n = control$min_samples, max_boundary = control$max_boundary,
                 min_sd = control$min_sd, endpoint = control$endpoint)
  if (!is.na(z$skip)) { out$skip <- z$skip; return(out) }
  if (!any(z$y > 0 & z$y < 1)) { out$skip <- "no_interior"; return(out) }
  s <- names(z$y); y <- unname(z$y); Cs <- orth_design(C_all[s, , drop = FALSE])
  n <- length(y); n0 <- sum(y == 0); n1 <- sum(y == 1); mid <- y > 0 & y < 1
  out$family <- z$family; out$n <- n; out$n0 <- n0; out$n1 <- n1; out$n_mid <- sum(mid)
  out$n_cov <- ncol(Cs) - 1L; out$n_cov_dropped <- length(attr(Cs, "dropped"))
  out$rank_mid <- rank_mid(Cs, mid); out$p_mean <- out$rank_mid + 1L
  out$samples <- match(s, names(yv)); out$y <- y; out$design <- Cs
  out$null <- promoter_null(y, Cs, z$family, maxit = control$maxit)
  out
}

#' Fit the BOHR model of a proportion response without the tested predictor
#'
#' Fits the zero-and-one-inflated beta model of one response `y` in
#' \[0, 1\] with exact 0s and 1s (for example the relative usage of one
#' promoter or isoform across samples), with covariates only. Every
#' association test of the response starts from this fit (see
#' [bohr_test()]).
#'
#' The fitting rules are those of the BOHR paper (Supplementary Methods,
#' Section 2); their thresholds are set with [bohr_control()].
#' \enumerate{
#'   \item Missing values of `y` are removed. The response is skipped
#'     (`skip = "n"`) when fewer than `min_samples` (20) observations remain,
#'     or fewer observations than covariates. It is also skipped
#'     (`skip = "prop_or_sd"`) when at least `max_boundary` (95\%) of the
#'     observations lie at 0 or at 1, or when the standard deviation is
#'     `min_sd` (0.05) or less.
#'   \item Boundary rule: the observations at a boundary are removed when that
#'     boundary holds `endpoint` (5\%) of the observations or fewer.
#'   \item The family follows from the boundaries that remain: `"BEINF"` (0s
#'     and 1s; components p0, p1 and mu), `"BEZI"` (0s only; p0 and mu),
#'     `"BEOI"` (1s only; p1 and mu) or `"BE"` (neither; mu only). A response
#'     without values strictly between 0 and 1 is skipped
#'     (`skip = "no_interior"`).
#'   \item The covariates are standardised over all observations (mean 0,
#'     standard deviation 1) and an intercept is added. The design of the
#'     remaining observations is then replaced by an orthonormal basis of its
#'     column space; columns collinear at tolerance 1e-7 are dropped. This
#'     changes neither the maximised likelihood nor any test statistic.
#'   \item The two parts of the likelihood share no parameters and are fitted
#'     separately: a beta regression of the intermediate values (logit-linked
#'     mean, constant precision), fitted when there are more intermediate
#'     observations than design columns plus 2; and the categorical part, a
#'     multinomial logit of 0 / intermediate / 1 with the intermediate state as
#'     reference (`"BEINF"`) or a logistic regression of being at the
#'     boundary (`"BEZI"`, `"BEOI"`). The maximiser is a monotone Newton /
#'     Fisher-scoring algorithm with step halving and Levenberg-Marquardt
#'     damping, at most `maxit` (200) iterations.
#' }
#'
#' @param y Numeric vector, the response in \[0, 1\], one value per
#'   observation. `NA` marks a missing value.
#' @param covariates `NULL`, or a numeric matrix or data frame with one row per
#'   element of `y` and no intercept column. Must not contain missing values.
#' @param control Fitting rules and optimiser limits from [bohr_control()];
#'   the defaults are the values of the BOHR paper.
#' @return An object of class `"bohr_fit"`, a list with elements
#'   \describe{
#'     \item{`family`}{`"BEINF"`, `"BEZI"`, `"BEOI"` or `"BE"`; `NA` when skipped.}
#'     \item{`skip`}{`NA` when fitted, otherwise the reason (`"n"`,
#'       `"prop_or_sd"` or `"no_interior"`).}
#'     \item{`n`, `n0`, `n1`, `n_mid`}{Observations used after the boundary
#'       rule; those at 0, at 1 and in between.}
#'     \item{`n_cov`, `n_cov_dropped`}{Covariate columns kept in the design and
#'       dropped as collinear.}
#'     \item{`rank_mid`, `p_mean`}{Rank of the design on the intermediate
#'       observations, and `rank_mid + 1`, the number of mean parameters used by
#'       the F correction of the mu test.}
#'     \item{`samples`}{Indices (into `y`) of the observations used.}
#'     \item{`y`, `design`}{The response values and the orthonormal design
#'       used.}
#'     \item{`control`}{The settings used.}
#'     \item{`null`}{The fitted parts: `null$beta` (beta regression; `NULL`
#'       when not fitted) and `null$cat` (categorical part; `NULL` for `"BE"`),
#'       each with `ll` (log-likelihood), `coef`, `vcov`, `converged`, `iter`
#'       and `maxgrad`.}
#'     \item{`input`}{Counts before the fitting rules (observations given,
#'       missing, at 0 and at 1, standard deviation, covariate names and
#'       scaling), from which `print()` reports every rule that changed the
#'       data or the model.}
#'   }
#' @seealso [bohr_test()] for association tests, [bohr()] for the formula
#'   interface, [bohr_control()] for the settings.
#' @examples
#' set.seed(1)
#' d <- simulate_response(200, mu = 0.4, p0 = 0.2, p1 = 0.1)
#' fit <- bohr_fit(d$y)
#' fit
#' fit$null$beta$coef
#' @export
bohr_fit <- function(y, covariates = NULL, control = bohr_control()) {
  check_activity(y); check_control(control)
  y <- as.vector(y)
  Z <- as_covariates(covariates, length(y))
  if (nrow(Z) != length(y)) stop("`covariates` must have one row per element of `y`.", call. = FALSE)
  nm <- paste0("s", seq_along(y))
  rownames(Z) <- nm; yv <- stats::setNames(y, nm)
  fit <- fit_promoter(yv, make_design(Z), control)
  fit$control <- control
  fit$input <- input_summary(y, Z)
  structure(fit, class = "bohr_fit")
}

## counts before the fitting rules (for the rule messages) and the covariate scaling of cov_matrix()
## (for mapping coefficients back to the covariates' own scale in bohr())
input_summary <- function(y, Z) {
  ok <- !is.na(y)
  sc <- if (ncol(Z)) scale(Z) else NULL
  list(n_input = length(y), n_missing = sum(!ok), n_at0 = sum(y[ok] == 0), n_at1 = sum(y[ok] == 1),
       sd = if (sum(ok) > 1L) stats::sd(y[ok]) else NA_real_, covariates = colnames(Z),
       center = if (is.null(sc)) numeric(0) else attr(sc, "scaled:center"),
       scale = if (is.null(sc)) numeric(0) else attr(sc, "scaled:scale"))
}

#' @export
print.bohr_fit <- function(x, ...) {
  if (!is.na(x$skip)) {
    cat("BOHR null fit: skipped (", x$skip, ")\n", sep = "")
  } else {
    cat("BOHR null fit, family ", x$family, "\n", sep = "")
    cat(sprintf("  samples %d: %d at 0, %d at 1, %d in between\n", x$n, x$n0, x$n1, x$n_mid))
    cat(sprintf("  covariates in the design: %d (%d dropped as collinear)\n", x$n_cov, x$n_cov_dropped))
    b <- x$null$beta; k <- x$null$cat
    if (!is.null(b)) cat(sprintf("  beta part: log-likelihood %.4f, converged %s\n", b$ll, b$converged))
    else cat("  beta part: not fitted (too few intermediate samples)\n")
    if (!is.null(k)) cat(sprintf("  categorical part: log-likelihood %.4f, converged %s\n", k$ll, k$converged))
  }
  print_rules(fit_rules(x))
  invisible(x)
}
