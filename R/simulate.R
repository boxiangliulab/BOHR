#' Simulate a proportion response with exact 0s and 1s
#'
#' Draws a predictor and a response from a zero-and-one-inflated beta
#' distribution whose three components depend on the predictor. Uses base R
#' only.
#'
#' For an observation with predictor value \eqn{x}, the components are
#' \eqn{\mu(x) = \mu + x\,\delta_\mu}, \eqn{p_0(x) = p_0 + x\,\delta_0} and
#' \eqn{p_1(x) = p_1 + x\,\delta_1}, with effects \eqn{\delta} per unit of
#' \eqn{x} on the probability scale (as in the simulations of the BOHR paper,
#' where \eqn{x} was an allele dosage). The response is 0 with probability
#' \eqn{p_0(x)}, 1 with probability \eqn{p_1(x)}, and otherwise a beta value
#' with mean \eqn{\mu(x)} and dispersion \eqn{\sigma} (variance
#' \eqn{\sigma^2\mu(1-\mu)}, precision \eqn{\phi = (1 - \sigma^2)/\sigma^2};
#' the BEINF parameterisation of gamlss), kept within `range`. By default
#' \eqn{x} is drawn as the dosage (0, 1 or 2 copies) of an allele with
#' frequency `freq` under Hardy-Weinberg equilibrium.
#'
#' The effect that [bohr_test()] estimates for a contrast from \eqn{a} to
#' \eqn{b} is \eqn{(b - a)\,\delta}: \eqn{2\delta} for a dosage (2 versus 0).
#'
#' Random numbers come from the current R session; call [set.seed()] first for
#' reproducible draws.
#'
#' @param n Number of observations.
#' @param x Optional numeric vector of `n` predictor values. By default a
#'   dosage drawn with allele frequency `freq`.
#' @param mu,sigma Mean (0 < `mu` < 1) and dispersion (0 < `sigma` < 1) of the
#'   intermediate values at \eqn{x = 0}.
#' @param p0,p1 Probabilities of 0 and of 1 at \eqn{x = 0}.
#' @param effect_mu,effect_p0,effect_p1 Change in mu, p0 and p1 per unit of
#'   \eqn{x}. The components must stay valid at every value of \eqn{x}:
#'   0 < mu < 1, p0 >= 0, p1 >= 0 and p0 + p1 < 1.
#' @param freq Allele frequency (0 < `freq` < 1) of the default dosage
#'   predictor; ignored when `x` is given.
#' @param range Length-2 vector; intermediate values are kept within it (the
#'   paper used 0.005 and 0.995).
#' @return A data frame with columns `x` and `y`. Attribute `"effects"` holds
#'   the true effects per unit of \eqn{x}: `effect_mu`, `effect_p0` and
#'   `effect_p1`.
#' @seealso [bohr_test()], [bohr_fit()].
#' @examples
#' set.seed(1)
#' d <- simulate_bohr(300, mu = 0.5, p0 = 0.2, p1 = 0.1, effect_p0 = 0.1)
#' table(x = d$x, at_zero = d$y == 0)
#' attr(d, "effects")
#' ## a continuous predictor
#' d2 <- simulate_bohr(300, x = runif(300), p0 = 0.1, effect_mu = 0.2)
#' summary(d2$y)
#' @export
simulate_bohr <- function(n, x = NULL, mu = 0.5, sigma = 0.3, p0 = 0.2, p1 = 0.1,
                          effect_mu = 0, effect_p0 = 0, effect_p1 = 0,
                          freq = 0.3, range = c(0.005, 0.995)) {
  n <- check_count(n, "n", min = 1)
  num1 <- function(v, what) if (!is.numeric(v) || length(v) != 1L || !is.finite(v)) stop("`", what, "` must be a single finite number.", call. = FALSE)
  for (a in c("mu", "sigma", "p0", "p1", "effect_mu", "effect_p0", "effect_p1")) num1(get(a), a)
  if (is.null(x)) {
    num1(freq, "freq")
    if (freq <= 0 || freq >= 1) stop("`freq` must lie strictly between 0 and 1.", call. = FALSE)
    x <- stats::rbinom(n, 2L, freq)
  } else if (!is.numeric(x) || length(x) != n || anyNA(x) || any(!is.finite(x))) stop("`x` must be a vector of n finite numbers.", call. = FALSE)
  if (sigma <= 0 || sigma >= 1) stop("`sigma` must lie strictly between 0 and 1.", call. = FALSE)
  if (!is.numeric(range) || length(range) != 2L || anyNA(range) || range[1] < 0 || range[2] > 1 || range[1] >= range[2])
    stop("`range` must be two increasing values in [0, 1].", call. = FALSE)
  m <- mu + effect_mu * x; q0 <- p0 + effect_p0 * x; q1 <- p1 + effect_p1 * x
  if (any(m <= 0 | m >= 1)) stop("mu + effect_mu * x must lie strictly between 0 and 1 for every x.", call. = FALSE)
  if (any(q0 < 0 | q1 < 0 | q0 + q1 >= 1)) stop("p0 and p1 (with their effects) must be >= 0 and sum to less than 1 for every x.", call. = FALSE)
  phi <- (1 - sigma^2) / sigma^2
  y <- stats::rbeta(n, m * phi, (1 - m) * phi)
  y <- pmin(pmax(y, range[1]), range[2])
  u <- stats::runif(n)
  y[u < q0] <- 0
  y[u >= q0 & u < q0 + q1] <- 1
  structure(data.frame(x = x, y = y), effects = c(effect_mu = effect_mu, effect_p0 = effect_p0, effect_p1 = effect_p1))
}
