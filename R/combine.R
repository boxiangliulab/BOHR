## ---- ACAT -----------------------------------------------------------------------------------------
## acat_log: verbatim from the analysis engine of the BOHR paper,
##   Figure4_puQTL_landscape_eQTL_coloc/rerun_v2/code/discovery_nogate.R:18-21 (= final_discovery.R:16-19).
## Input: natural-log P values; non-finite entries (NA, -Inf) are ignored. Package change (author decision
## 2026-11-10, every setting an argument): argument `cap` (default TRUE = the original) switches the 1 - 1/M cap.
acat_log <- function(logp, cap = TRUE) { p <- exp(logp[is.finite(logp)]); if (!length(p)) return(NA_real_); if (any(p == 0)) return(0)   # = final_discovery.R
  if (cap && length(p) > 1) p <- pmin(p, 1 - 1 / length(p))
  sm <- p < 1e-15; t <- numeric(length(p)); t[!sm] <- tan((0.5 - p[!sm]) * pi); t[sm] <- 1 / (p[sm] * pi); st <- mean(t)
  if (st > 1e15) return(1 / (st * pi)); pcauchy(st, lower.tail = FALSE) }

#' Aggregated Cauchy association test (ACAT)
#'
#' Combines P values with the Cauchy combination test of Liu and Xie (2020),
#' with equal weights. It stays valid under arbitrary dependence between the
#' tests, for example when many correlated predictors are tested against one
#' response (in [bohr_cis()], the variants of one cis window).
#'
#' Following the BOHR analysis engine, each P value is capped at
#' \eqn{1 - 1/M} when there are \eqn{M > 1} of them (`cap = TRUE`), P values below
#' \eqn{10^{-15}} use the approximation \eqn{\tan((0.5 - p)\pi) \approx
#' 1/(p\pi)}, and a combined statistic above \eqn{10^{15}} gives
#' \eqn{P = 1/(T\pi)}.
#'
#' @param p Numeric vector of P values in \[0, 1\], or of natural-log P values
#'   when `log = TRUE`. Missing values are ignored; with `log = TRUE`, so are
#'   `-Inf` values.
#' @param log Logical: are `p` natural-log P values? The cis mapping of
#'   [bohr_cis()] combines log P values.
#' @param cap Cap each P value at \eqn{1 - 1/M} (\eqn{M} = number of P values)
#'   before combining. Default `TRUE`, as in the paper's analysis code (not
#'   stated in the Methods or SI).
#' @return The combined P value (not logged); `NA` when no P value is usable,
#'   and 0 when any P value is 0.
#' @references Liu Y, Xie J (2020). Cauchy combination test: a powerful test
#'   with analytic p-value calculation under arbitrary dependency structures.
#'   *Journal of the American Statistical Association* 115, 393-402.
#'   \doi{10.1080/01621459.2018.1554485}
#' @seealso [simes()], [bb_select()].
#' @examples
#' acat(c(0.01, 0.2, 0.5, 0.9))
#' acat(log(c(0.01, 0.2, 0.5, 0.9)), log = TRUE)
#' @export
acat <- function(p, log = FALSE, cap = TRUE) {
  check_flag(log, "log"); check_flag(cap, "cap")
  if (!is.numeric(p) && !(is.logical(p) && all(is.na(p)))) stop("`p` must be numeric.", call. = FALSE)
  p <- as.numeric(p)
  if (log) {
    if (any(p > 0, na.rm = TRUE)) stop("log P values must be <= 0.", call. = FALSE)
    return(acat_log(p, cap = cap))
  }
  if (any(p < 0 | p > 1, na.rm = TRUE)) stop("`p` must lie in [0, 1].", call. = FALSE)
  p <- p[is.finite(p)]
  if (!length(p)) return(NA_real_)
  if (any(p == 0)) return(0)
  acat_log(log(p), cap = cap)
}

## ---- Simes ----------------------------------------------------------------------------------------
## simes_p: verbatim from Figure4_puQTL_landscape_eQTL_coloc/rerun_v2/code/fig2e_deliver_freeze_v3m.R:18
## (= fig2e_deliver_freeze_v3m_pairs.R:24). sort() drops missing values.
simes_p <- function(p) { p <- sort(p); k <- length(p); min(1, min(k * p / seq_len(k))) }

#' Simes combination of P values
#'
#' Combines a family of P values, for example those of the components of one
#' response, into one P value with the Simes method:
#' \eqn{\min(1, \min_i k\,p_{(i)}/i)} over the \eqn{k} sorted P values. The
#' BOHR paper used it only to compare BOHR with linear QTL mapping (one P value
#' per feature).
#'
#' @param p Numeric vector of P values in \[0, 1\]. Missing values are ignored.
#' @return The Simes P value; `NA` when `p` has no non-missing value.
#' @references Simes RJ (1986). An improved Bonferroni procedure for multiple
#'   tests of significance. *Biometrika* 73, 751-754.
#'   \doi{10.1093/biomet/73.3.751}
#' @seealso [acat()], [bb_select()].
#' @examples
#' simes(c(0.01, 0.04, 0.3))
#' @export
simes <- function(p) {
  if (!is.numeric(p) && !(is.logical(p) && all(is.na(p)))) stop("`p` must be numeric.", call. = FALSE)
  p <- as.numeric(p)
  if (any(p < 0 | p > 1, na.rm = TRUE)) stop("`p` must lie in [0, 1].", call. = FALSE)
  if (all(is.na(p))) return(NA_real_)
  simes_p(p)
}

## ---- Benjamini-Bogomolov -------------------------------------------------------------------------
#' Select significant components with the Benjamini-Bogomolov procedure
#'
#' Two-stage selection of significant groups and of the significant members of
#' each selected group (Benjamini and Bogomolov 2014), for example responses
#' (groups) and their component tests mu, p0 and p1 (members). The BOHR paper
#' used it to compare BOHR with linear QTL mapping, with features (or exact
#' mirror pairs) as groups.
#'
#' Stage 1 combines the P values of each group with [simes()] and applies the
#' Benjamini-Hochberg procedure across groups; the \eqn{R} groups with
#' \eqn{q \le} `alpha`, of \eqn{m} groups, are selected. Stage 2 applies the
#' Benjamini-Hochberg procedure within each selected group and rejects the
#' members whose adjusted P value is at most `alpha` \eqn{\times R/m}.
#'
#' @param p Numeric vector of P values in \[0, 1\] (members of all groups).
#' @param group Vector of the same length as `p` giving each P value's group.
#' @param alpha Target level of both stages.
#' @return A data frame with one row per element of `p`: `group`, `p`,
#'   `p_group` (Simes P value of the group), `q_group` (its
#'   Benjamini-Hochberg q value across groups), `selected` (group selected at
#'   stage 1) and `rejected` (member significant at stage 2). Attributes `R`,
#'   `m` and `level` (the stage-2 level `alpha * R / m`).
#' @references Benjamini Y, Bogomolov M (2014). Selective inference on
#'   multiple families of hypotheses. *Journal of the Royal Statistical
#'   Society Series B* 76, 297-318. \doi{10.1111/rssb.12028}
#' @seealso [simes()], [acat()].
#' @examples
#' p <- c(1e-6, 0.03, 0.4, 0.2, 0.5, 0.01, 0.002)
#' g <- c("A", "A", "A", "B", "B", "C", "C")
#' bb_select(p, g)
#' @export
bb_select <- function(p, group, alpha = 0.05) {
  if (!is.numeric(p) || anyNA(p) || any(p < 0 | p > 1)) stop("`p` must be P values in [0, 1] without missing values.", call. = FALSE)
  if (length(group) != length(p) || anyNA(group)) stop("`group` must have one non-missing value per P value.", call. = FALSE)
  if (!is.numeric(alpha) || length(alpha) != 1L || !(alpha > 0 && alpha < 1)) stop("`alpha` must lie strictly between 0 and 1.", call. = FALSE)
  ## fig2e_deliver_freeze_v3m.R:29-34: Simes per group, BH over groups, BH within each selected group at
  ## level alpha * R / m.
  gr <- unique(group)
  pg <- vapply(gr, function(k) simes_p(p[group == k]), 0)
  qg <- stats::p.adjust(pg, "BH")
  sel <- qg <= alpha
  m <- length(gr); R <- sum(sel); lvl <- alpha * R / m
  rej <- rep(FALSE, length(p))
  for (k in gr[sel]) { i <- which(group == k); rej[i] <- stats::p.adjust(p[i], "BH") <= lvl }
  j <- match(group, gr)
  structure(data.frame(group = group, p = p, p_group = unname(pg[j]), q_group = unname(qg[j]), selected = unname(sel[j]),
                       rejected = rej, stringsAsFactors = FALSE),
            R = R, m = m, level = lvl)
}
