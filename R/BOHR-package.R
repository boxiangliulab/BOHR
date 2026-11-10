#' BOHR: beta one-zero hybrid regression for proportions with exact 0s and 1s
#'
#' Many proportions are measured with exact 0s and 1s: the share of a gene's
#' activity that comes from one promoter, the usage of an isoform or exon,
#' the fraction of a cell population, a methylation level. BOHR models such a
#' response with a zero-and-one-inflated beta distribution and tests a
#' predictor's effect on its three components separately: p0, the
#' probability of 0; p1, the probability of 1; and mu, the mean of the values
#' in between. Alternative promoter usage is its first application.
#'
#' General regression:
#' \itemize{
#'   \item [bohr_fit()]: fit of one response without the tested predictor;
#'   \item [bohr_test()]: association of one response with any numeric
#'     predictor (allele dosage, binary or continuous), with effects for a
#'     chosen contrast;
#'   \item [bohr_control()]: fitting rules and optimiser limits (defaults are
#'     the values of the BOHR paper);
#'   \item [simulate_bohr()]: simulated predictor and response.
#' }
#'
#' @references
#' Ospina R, Ferrari SLP (2010). Inflated beta distributions. *Statistical
#' Papers* 51, 111-126. \doi{10.1007/s00362-008-0125-4}
#'
#' Ferrari SLP, Cribari-Neto F (2004). Beta regression for modelling rates and
#' proportions. *Journal of Applied Statistics* 31, 799-815.
#' \doi{10.1080/0266476042000214501}
#'
#' Firth D (1993). Bias reduction of maximum likelihood estimates.
#' *Biometrika* 80, 27-38. \doi{10.1093/biomet/80.1.27}
#' @keywords internal
#' @importFrom stats optim p.adjust pcauchy pchisq pf plogis pnorm qlogis rbeta rbinom runif sd setNames var
"_PACKAGE"
