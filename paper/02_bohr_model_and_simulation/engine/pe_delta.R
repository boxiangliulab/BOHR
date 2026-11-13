## Probability-scale effects with delta-method SEs, from the SAME fits as the tests:
##   mu    : the beta regression of snp_tests (ML, logit mean, interior donors; same start values)
##   p0/p1 : the full Firth fit that firth_snp_tests selects (multi-start + alternating rounds). firth_snp_tests_ff is
##           firth_snp_tests with that fit returned as r$ff (body copied programmatically, so the search is identical).
## pe = mean over the model's donors of P(dosage 2) - P(dosage 0), covariates at their observed values (the snp_tests
## definition: mu over all donors' covariate profiles; p0/p1 from the logit (BEZI/BEOI) or multinomial (BEINF) fit).
## SE = sqrt(grad' V grad), V = inverse expected information of the same fit (the vcov behind est / se).
## Sourced after bohr2_engine.R and bohr2_addon3.R.
firth_snp_tests_ff <- firth_snp_tests
local({ b <- as.list(body(firth_snp_tests_ff)); n <- length(b); stopifnot(identical(b[[n]], quote(r)))
  body(firth_snp_tests_ff) <<- as.call(c(b[-n], quote(if (exists("ff", inherits = FALSE)) r$ff <- ff), quote(r))) })
dlogis1 <- function(x) { q <- plogis(x); q * (1 - q) }
quad_se <- function(gr, V) if (is.null(V) || any(!is.finite(V))) NA_real_ else sqrt(max(drop(crossprod(gr, V %*% gr)), 0))
## mu: refit exactly as snp_tests; returns est, se (link) and pe, pe_se
pe_mu_fit <- function(y, C, g, nul) {
  out <- list(est = NA_real_, se = NA_real_, pe = NA_real_, pe_se = NA_real_, conv = NA)
  mid <- y > 0 & y < 1; pC <- ncol(C)
  if (is.null(nul$beta) || length(unique(g[mid])) < 2) return(out)
  st <- c(nul$beta$coef[1:pC], 0, nul$beta$coef[pC + 1])
  f <- beta_fit(y[mid], cbind(C[mid, , drop = FALSE], g = g[mid]), start = st)
  b <- f$coef[1:pC]; gm <- f$coef[pC + 1]; eta <- drop(C %*% b)
  gr <- c(colMeans((dlogis1(eta + 2 * gm) - dlogis1(eta)) * C), mean(2 * dlogis1(eta + 2 * gm)), 0)   # (b, gamma, log phi)
  list(est = unname(gm), se = if (!is.null(f$vcov)) sqrt(f$vcov[pC + 1, pC + 1]) else NA_real_,
       pe = mean(plogis(eta + 2 * gm) - plogis(eta)), pe_se = quad_se(gr, f$vcov), conv = f$converged)
}
## p0/p1 from the Firth full fit ff (family BEZI / BEOI: one logit; BEINF: multinomial, baseline = interior)
pe_firth <- function(ff, C, fam) {
  out <- list(p0_pe = NA_real_, p0_pe_se = NA_real_, p1_pe = NA_real_, p1_pe_se = NA_real_)
  if (is.null(ff)) return(out); th <- ff$coef; V <- ff$vcov; pC <- ncol(C)
  if (fam %in% c("BEZI", "BEOI")) {
    eta <- drop(C %*% th[1:pC]); a <- th[pC + 1]
    pe <- mean(plogis(eta + 2 * a) - plogis(eta)); gr <- c(colMeans((dlogis1(eta + 2 * a) - dlogis1(eta)) * C), mean(2 * dlogis1(eta + 2 * a)))
    if (fam == "BEZI") { out$p0_pe <- pe; out$p0_pe_se <- quad_se(gr, V) } else { out$p1_pe <- pe; out$p1_pe_se <- quad_se(gr, V) }
  } else if (fam == "BEINF") {
    e0 <- drop(C %*% th[1:pC]); a0 <- th[pC + 1]; e1 <- drop(C %*% th[pC + 1 + 1:pC]); a1 <- th[2 * (pC + 1)]
    pr <- function(d) { l0 <- e0 + a0 * d; l1 <- e1 + a1 * d; m <- pmax(0, l0, l1); den <- exp(-m) + exp(l0 - m) + exp(l1 - m)
      list(q0 = exp(l0 - m) / den, q1 = exp(l1 - m) / den) }
    P0 <- pr(0); P2 <- pr(2)
    v00 <- function(P) P$q0 * (1 - P$q0); v11 <- function(P) P$q1 * (1 - P$q1); v01 <- function(P) -P$q0 * P$q1
    g0 <- c(colMeans((v00(P2) - v00(P0)) * C), mean(2 * v00(P2)), colMeans((v01(P2) - v01(P0)) * C), mean(2 * v01(P2)))
    g1 <- c(colMeans((v01(P2) - v01(P0)) * C), mean(2 * v01(P2)), colMeans((v11(P2) - v11(P0)) * C), mean(2 * v11(P2)))
    out$p0_pe <- mean(P2$q0 - P0$q0); out$p0_pe_se <- quad_se(g0, V); out$p1_pe <- mean(P2$q1 - P0$q1); out$p1_pe_se <- quad_se(g1, V)
  }
  out
}
## interior-only moments by genotype: donors with 0 < y < 1 of the RAW usage
interior_moments <- function(y, g) { ok <- !is.na(y) & !is.na(g) & y > 0 & y < 1; y <- y[ok]; g <- g[ok]
  f <- function(d, fun) { v <- y[g == d]; if (length(v)) fun(v) else NA_real_ }                  # as genotype_pheno's gstat
  c(n0_mid = sum(g == 0), n1_mid = sum(g == 1), n2_mid = sum(g == 2), m0_mid = f(0, mean), m1_mid = f(1, mean), m2_mid = f(2, mean),
    sd0_mid = f(0, sd), sd1_mid = f(1, sd), sd2_mid = f(2, sd)) }
