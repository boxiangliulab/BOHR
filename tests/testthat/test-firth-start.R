## The null start of the restricted Firth fits (bohr_control(firth_null_start = TRUE), stage 3b-1).

## A synthetic BEINF case with two maxima of the restricted penalised likelihood of the p0 test, built like the conditional
## test of Brain_Putamen PromoterID_47842 p0 in the paper: the tested variant g equals a variant L among the covariates
## except in one sample, p0 falls and p1 rises with the dosage, and no sample with dosage 0 is at 1 (the p1 equation is
## nearly separated). Found by a seed search; in about 1.5% of such data sets the paper's engine stops at the lower maximum.
two_maxima_case <- function() with_seed(655, {
  n <- 150; k <- 8
  Z <- matrix(stats::rnorm(n * k), n, k, dimnames = list(NULL, paste0("c", 1:k)))
  L <- stats::rbinom(n, 2, 0.4)
  g <- L; m <- sample(1:3, 1); i <- sample.int(n, m); g[i] <- ifelse(g[i] == 0, 1, ifelse(g[i] == 2, 1, sample(c(0, 2), m, TRUE)))
  p0 <- c(0.97, 0.58, 0.06)[L + 1]; p1 <- c(0, 0.07, 0.16)[L + 1]; pm <- 1 - p0 - p1
  l0 <- log(p0 / pm) + 0.7 * Z[, 1]; l1 <- ifelse(p1 > 0, log(pmax(p1, 1e-300) / pm) + 0.7 * Z[, 2], -Inf)
  den <- 1 + exp(l0) + exp(l1); q0 <- exp(l0) / den; q1 <- exp(l1) / den
  mu <- stats::plogis(stats::qlogis(0.4) + 0.3 * Z[, 3])
  y <- stats::rbeta(n, mu * 6, (1 - mu) * 6); u <- stats::runif(n)
  y[u < q0] <- 0; y[u >= q0 & u < q0 + q1] <- 1
  list(y = y, g = g, X = cbind(Z, L = L))
})

test_that("firth_null_start = TRUE returns the higher of two restricted maxima", {
  d <- two_maxima_case()
  fit <- bohr_fit(d$y, d$X)
  expect_identical(fit$family, "BEINF")
  g <- d$g[fit$samples]
  rT <- BOHR:::firth_snp_tests(fit$y, fit$design, g, fit$null, keep_fits = TRUE, null_start = TRUE)
  rF <- BOHR:::firth_snp_tests(fit$y, fit$design, g, fit$null, keep_fits = TRUE, null_start = FALSE)
  ## restricted maximum = full penalised log-likelihood - statistic / 2
  llr <- function(r) r$fits$cat$llpen - r$nu_firth_stat / 2
  ## the restricted fit started from the null fit (the tested p0 coefficient held at 0) reaches the maximum TRUE returns
  pC <- ncol(fit$design); Cg <- cbind(fit$design, g = g)
  sn <- c(BOHR:::na0(fit$null$cat$coef[1:pC]), 0, BOHR:::na0(fit$null$cat$coef[pC + 1:pC]), 0)
  direct <- BOHR:::firth_mnl3(as.numeric(fit$y == 0), as.numeric(fit$y == 1), Cg, Cg, start = sn, fixed = pC + 1L)
  expect_equal(llr(rT), direct$llpen, tolerance = 1e-8)
  expect_gte(llr(rT), llr(rF) - 1e-9)
  expect_equal(rT$fits$cat$llpen, rF$fits$cat$llpen, tolerance = 1e-8)        # the full fit is the same
  ## the public functions: the default (TRUE) gives the small statistic of the higher maximum, |z| about 1
  a <- bohr_test(d$y, d$g, d$X, type = "dosage", all_stats = TRUE)
  b <- bohr_test(d$y, d$g, d$X, type = "dosage", all_stats = TRUE, control = engine_control())
  expect_identical(a$stat_p0, rT$nu_firth_stat); expect_identical(b$stat_p0, rF$nu_firth_stat)
  expect_lt(a$stat_p0, 5); expect_gt(a$p_p0, 0.01)
  expect_lt(abs(a$est_p0 / a$se_p0), 2)
  expect_equal(a$stat_p1, b$stat_p1, tolerance = 1e-8)                         # the p1 test has one maximum here
  ## on this platform the paper's engine stops at the lower maximum (log P about -18 with |z| about 1, as in the paper's
  ## Brain_Putamen refit); other platforms may reach the higher one, so this is checked locally only
  if (identical(Sys.getenv("NOT_CRAN"), "true")) {
    expect_gt(b$stat_p0, 25); expect_gt(llr(rT) - llr(rF), 10)
  }
})

test_that("with one maximum, the null start changes nothing beyond optimiser precision", {
  cases <- parity_cases()
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    a <- bohr_test(cs$y, cs$g, cs$Z, type = "dosage", all_stats = TRUE)
    b <- bohr_test(cs$y, cs$g, cs$Z, type = "dosage", all_stats = TRUE, control = engine_control())
    for (k in c(BOHR:::FIRTH_COLS, "log_p_p0", "log_p_p1", "est_p0", "se_p1")) {
      if (is.logical(b[[k]])) expect_identical(a[[k]], b[[k]], info = paste(nm, k))
      else expect_equal(a[[k]], b[[k]], tolerance = 1e-6, info = paste(nm, k))
    }
    for (k in c("log_p_mu", "est_mu", "effect_p0", "effect_p1")) expect_identical(a[[k]], b[[k]], info = paste(nm, k))
  }
  ## with firth = FALSE (unpenalised p0/p1 tests) the setting has no effect
  cs <- cases$BEINF_cov
  expect_identical(bohr_test(cs$y, cs$g, cs$Z, control = bohr_control(firth = FALSE)),
                   bohr_test(cs$y, cs$g, cs$Z, control = bohr_control(firth = FALSE, firth_null_start = FALSE)))
})
