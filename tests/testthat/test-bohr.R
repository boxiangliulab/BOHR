## The formula interface bohr(): agreement with bohr_test(), coefficients on the user's scale, the methods, and
## one test per fitting-rule message.

dat <- function(seed = 61, n = 220, p0 = 0.25, p1 = 0.12, ...) {
  with_seed(seed, {
    d <- data.frame(x = stats::rbinom(n, 2, 0.35), age = stats::rnorm(n, 50, 10), bmi = stats::rnorm(n, 25, 4),
                    sex = factor(sample(c("female", "male"), n, replace = TRUE)))
    d$y <- simulate_response(n, x = d$x, p0 = p0, p1 = p1, ...)$y
    d
  })
}
rule_codes <- function(m) c(m$rules$rule, m$notes$rule)

test_that("bohr() equals bohr_test() on the same data", {
  d <- dat(effect_p0 = 0.06, effect_mu = 0.04)
  Z <- stats::model.matrix(~ age + bmi + sex, d)[, -1]
  for (ty in c("dosage", "continuous", "auto")) {
    m <- bohr(y ~ x + age + bmi + sex, data = d, test = "x", type = ty)
    b <- bohr_test(d$y, cbind(x = d$x), covariates = Z, type = ty, all_stats = TRUE)
    expect_identical(m$result, b, info = ty)
  }
  ## a contrast and bohr_control() settings are passed on
  m <- bohr(y ~ x + age, data = d, type = "dosage", contrast = c(0, 1), control = bohr_control(firth = FALSE))
  b <- bohr_test(d$y, cbind(x = d$x), covariates = d["age"], type = "dosage", contrast = c(0, 1),
                 control = bohr_control(firth = FALSE), all_stats = TRUE)
  expect_identical(m$result, b)
  ## a two-level factor as the tested term is the 0/1 indicator of its second level
  m <- bohr(y ~ sex + age, data = d, test = "sex")
  b <- bohr_test(d$y, cbind(sexmale = as.numeric(d$sex == "male")), covariates = d["age"], all_stats = TRUE)
  expect_identical(m$result, b)
  expect_identical(m$test$name, "sexmale"); expect_identical(m$result$type, "binary")
  expect_identical(c(m$result$contrast_from, m$result$contrast_to), c(0, 1))
  ## a logical tested term
  d$flag <- d$sex == "male"
  expect_identical(bohr(y ~ flag + age, data = d, test = "flag")$result[, -1], m$result[, -1])
  ## a function of a variable as the tested term
  m <- bohr(y ~ log(bmi) + age, data = d, test = "log(bmi)")
  b <- bohr_test(d$y, cbind(`log(bmi)` = log(d$bmi)), covariates = d["age"], all_stats = TRUE)
  expect_identical(m$result, b)
  ## the first term is tested by default; variables can come from the formula's environment
  yy <- d$y; xx <- d$x
  expect_identical(bohr(yy ~ xx, type = "dosage")$result[, -1], bohr_test(yy, xx, type = "dosage", all_stats = TRUE)[, -1])
})

test_that("bohr() with missing values equals bohr_test() on the complete cases", {
  d <- dat(seed = 62, effect_mu = 0.05)
  d$y[c(3, 10)] <- NA; d$x[c(5, 11)] <- NA; d$age[c(7, 12, 13)] <- NA
  m <- bohr(y ~ x + age, data = d, type = "dosage")
  cc <- stats::complete.cases(d[c("y", "x", "age")])
  b <- bohr_test(d$y[cc], cbind(x = d$x[cc]), covariates = d["age"][cc, , drop = FALSE], type = "dosage", all_stats = TRUE)
  expect_identical(m$result, b)
  expect_identical(nobs(m), m$fit$n)
  ## bohr_test() on all rows standardises the covariates over the rows with a response too: same model
  d2 <- d[!is.na(d$age), ]
  b2 <- bohr_test(d2$y, d2$x, covariates = d2["age"], type = "dosage")
  for (k in c("log_p_mu", "log_p_p0", "log_p_p1", "est_mu", "est_p0", "effect_mu", "effect_p1"))
    expect_equal(m$result[[k]], b2[[k]], tolerance = 1e-7, info = k)
})

test_that("coefficients are on the covariates' scale and per unit of the tested term", {
  ## BEZI and BEOI with unpenalised p0/p1: the categorical part is a logistic regression
  for (cs in list(list(p0 = 0.3, p1 = 0, f = function(y) y == 0, k = "p0"), list(p0 = 0, p1 = 0.3, f = function(y) y == 1, k = "p1"))) {
    d <- dat(seed = 63, p0 = cs$p0, p1 = cs$p1, effect_mu = 0.04)
    m <- bohr(y ~ x + age + bmi + sex, data = d, type = "continuous", control = bohr_control(firth = FALSE))
    g <- stats::glm(cs$f(y) ~ x + age + bmi + sex, family = stats::binomial, data = d)
    expect_equal(unname(coef(m, cs$k)[names(coef(g))]), unname(coef(g)), tolerance = 1e-7)
    expect_equal(unname(sqrt(diag(vcov(m, cs$k)))[names(coef(g))]), unname(sqrt(diag(stats::vcov(g)))), tolerance = 1e-5)
    expect_equal(m$result$ll_cat_full, as.numeric(stats::logLik(g)), tolerance = 1e-9)
    expect_null(m$coefficients[[setdiff(c("p0", "p1"), cs$k)]])
  }
  ## BEINF with Firth (default): the tested term's coefficient and SE are bohr_test()'s est and se, per unit of x
  d <- dat(seed = 64, effect_p0 = 0.06)
  for (ty in c("dosage", "binary", "continuous")) {
    xx <- if (ty == "binary") as.numeric(d$x > 0) else d$x
    d$xx <- xx
    m <- bohr(y ~ xx + age + sex, data = d, type = ty)
    for (k in c("mu", "p0", "p1")) {
      expect_equal(coef(m, k)[["xx"]], m$result[[paste0("est_", k)]], tolerance = 1e-12, info = paste(ty, k))
      expect_equal(sqrt(vcov(m, k)["xx", "xx"]), m$result[[paste0("se_", k)]], tolerance = 1e-10, info = paste(ty, k))
    }
  }
  ## the linear predictors on the user's scale equal the engine's on the observations used
  m <- bohr(y ~ x + age + bmi, data = d, type = "dosage")
  f <- m$fit; s <- f$samples
  fits <- BOHR:::engine_one(f, d$x[s], bohr_control(), keep_fits = TRUE)$fits
  X <- cbind(1, d$x, d$age, d$bmi)[s, ]; Cg <- cbind(f$design, d$x[s]); r1 <- ncol(Cg)
  expect_equal(drop(X %*% coef(m, "mu")), unname(drop(Cg %*% fits$beta$coef[seq_len(r1)])), tolerance = 1e-10)
  expect_equal(drop(X %*% coef(m, "p0")), unname(drop(Cg %*% fits$firth$coef[seq_len(r1)])), tolerance = 1e-10)
  expect_equal(drop(X %*% coef(m, "p1")), unname(drop(Cg %*% fits$firth$coef[r1 + seq_len(r1)])), tolerance = 1e-10)
  expect_equal(coef(m, "phi")[[1]], fits$beta$coef[[r1 + 1]])
  ## coef() without a component: every component, named component:term
  expect_identical(names(coef(m)), c(paste0("mu:", c("(Intercept)", "x", "age", "bmi")), "phi:(Intercept)",
                                     paste0("p0:", c("(Intercept)", "x", "age", "bmi")), paste0("p1:", c("(Intercept)", "x", "age", "bmi"))))
})

test_that("vcov(), logLik() and nobs() are consistent", {
  d <- dat(seed = 65)
  m <- bohr(y ~ x + age, data = d, type = "dosage")
  V <- vcov(m)
  expect_identical(rownames(V), names(coef(m)))
  expect_true(isSymmetric(unname(V)))
  expect_true(all(V[grep("^mu|^phi", rownames(V)), grep("^p0|^p1", rownames(V))] == 0))   # the two parts share no parameter
  expect_true(all(eigen(V, only.values = TRUE)$values > 0))
  ll <- logLik(m)
  expect_equal(as.numeric(ll), m$result$ll_beta_full + m$result$ll_cat_full)
  expect_identical(attr(ll, "df"), (3L + 1L) + 2L * 3L)                    # beta: 3 + phi; categorical: 2 equations of 3
  expect_identical(nobs(m), m$fit$n)
  expect_equal(stats::AIC(m), -2 * as.numeric(ll) + 2 * attr(ll, "df"))
  ## collinear covariates: NA coefficients and NA rows/columns
  d$age2 <- 3 * d$age - 2
  m2 <- bohr(y ~ x + age + age2, data = d, type = "dosage")
  expect_true(is.na(coef(m2, "mu")[["age2"]])); expect_true(all(is.na(vcov(m2, "p0")["age2", ])))
  expect_equal(coef(m2, "mu")[c("x", "age")], coef(m, "mu")[c("x", "age")], tolerance = 1e-8)
})

test_that("predict() gives mu, p0, p1 and the mean", {
  d <- dat(seed = 66, effect_p1 = 0.05)
  m <- bohr(y ~ x + sex, data = d, type = "dosage")
  P <- sapply(c("mu", "p0", "p1", "mean"), function(t) predict(m, type = t))
  expect_equal(P[, "mean"], P[, "p1"] + (1 - P[, "p0"] - P[, "p1"]) * P[, "mu"])
  expect_true(all(P[, c("mu", "p0", "p1")] > 0 & P[, c("mu", "p0", "p1")] < 1))
  ## newdata, by levels of the factor covariate
  nd <- data.frame(x = c(0, 2), sex = c("male", "female"))
  b <- coef(m, "mu")
  expect_equal(unname(predict(m, nd, type = "mu")), stats::plogis(c(b[["(Intercept)"]] + b[["sexmale"]], b[["(Intercept)"]] + 2 * b[["x"]])))
  expect_equal(unname(predict(m, d[1:5, ], type = "p0")), unname(P[1:5, "p0"]))
  ## BEZI: p1 is 0; the mean is (1 - p0) mu
  d0 <- dat(seed = 67, p0 = 0.3, p1 = 0)
  m0 <- bohr(y ~ x, data = d0, type = "dosage")
  expect_true(all(predict(m0, type = "p1") == 0))
  expect_equal(predict(m0), (1 - predict(m0, type = "p0")) * predict(m0, type = "mu"))
  ## na.exclude pads the predictions
  d$y[4] <- NA
  me <- bohr(y ~ x, data = d, type = "dosage", na.action = stats::na.exclude)
  pe <- predict(me); expect_length(pe, nrow(d)); expect_true(is.na(pe[4]))
  expect_error(predict(bohr(y ~ sex + x, data = d, test = "sex"), data.frame(sex = "other", x = 1)), "new level")
})

test_that("print() and summary() show the tables, the rules and the notes", {
  d <- dat(seed = 68, effect_mu = 0.05)
  m <- bohr(y ~ x + age + sex, data = d, type = "dosage")
  expect_output(print(m), "Tested term: x \\(dosage; effects compare x = 2 with x = 0\\)")
  expect_output(print(m), "Family: BEINF")
  s <- summary(m)
  expect_s3_class(s, "summary.bohr")
  expect_named(s$tables, c("mu", "p0", "p1"))
  expect_identical(rownames(s$tables$p0), c("(Intercept)", "x", "age", "sexmale"))
  expect_equal(s$tables$mu["x", "P value"], m$result$p_mu)
  expect_equal(s$tables$p1["x", "Effect"], m$result$effect_p1)
  expect_true(all(is.na(s$tables$mu[c("age", "sexmale"), "Statistic"])))
  expect_identical(unname(s$n), c(m$fit$n_mid, m$fit$n, m$fit$n))
  expect_output(print(s), "p1: probability of 1")
  expect_output(print(s), "Precision phi")
  expect_output(print(bohr(y ~ age + x, data = d, test = "age")), "continuous; effects compare age = .* \\(mean \\+ 1 SD\\)")
})

## ---- one test per rule message ------------------------------------------------------------------------
test_that("rule: observations left out by na.action", {
  d <- dat(seed = 69); d$y[1:3] <- NA; d$x[4] <- NA
  m <- bohr(y ~ x, data = d, type = "dosage")
  expect_true("na_action" %in% rule_codes(m))
  expect_output(print(m), "4 observations with missing values left out by na.action \\(omit\\)")
})

test_that("rule: missing responses left out (na.pass, and bohr_fit())", {
  d <- dat(seed = 70); d$y[1:2] <- NA
  m <- bohr(y ~ x, data = d, type = "dosage", na.action = stats::na.pass)
  expect_true("missing" %in% rule_codes(m)); expect_false("na_action" %in% rule_codes(m))
  expect_output(print(m), "2 observations with a missing response")
  expect_output(print(bohr_fit(d$y)), "2 observations with a missing response")
})

test_that("rule: the boundary rule removes 0s or 1s (with counts)", {
  d <- dat(seed = 71, p0 = 0.25, p1 = 0)
  d$y[1:4] <- 1                                                    # 4 of 220 at 1 (1.8%)
  m <- bohr(y ~ x, data = d, type = "dosage")
  expect_identical(m$fit$family, "BEZI")
  expect_true("boundary_1" %in% rule_codes(m))
  expect_output(print(m), "4 observations at 1 \\(1.8% of 220\\) removed by the boundary rule")
  expect_output(print(m), "p1 is not modelled")
  d2 <- dat(seed = 72, p0 = 0, p1 = 0.25); d2$y[1:3] <- 0
  m2 <- bohr(y ~ x, data = d2, type = "dosage")
  expect_true("boundary_0" %in% rule_codes(m2))
  expect_output(print(m2), "3 observations at 0 \\(1.4% of 220\\) removed by the boundary rule")
  expect_false(any(c("boundary_0", "boundary_1") %in% rule_codes(bohr(y ~ x, data = dat(seed = 73), type = "dosage"))))
})

test_that("rule: fitted mu only", {
  d <- dat(seed = 74, p0 = 0, p1 = 0); d$y[1:2] <- 0
  m <- bohr(y ~ x, data = d, type = "dosage")
  expect_identical(m$fit$family, "BE")
  expect_true(all(c("boundary_0", "mu_only") %in% rule_codes(m)))
  expect_output(print(m), "Fitted mu only: no observations at 0 or 1 remain")
  expect_named(m$coefficients, c("mu", "phi"))
})

test_that("rule: skipped for too few observations", {
  d <- dat(seed = 75, n = 15)
  m <- bohr(y ~ x, data = d, type = "dosage")
  expect_true("skip" %in% rule_codes(m))
  expect_output(print(m), "Skipped: 15 observations with a value, fewer than the 20 required")
  expect_true(is.na(m$result$p_mu)); expect_length(coef(m), 0L); expect_identical(nobs(m), 0L)
  expect_true(all(is.na(predict(m))))
})

test_that("rule: skipped for fewer observations than covariates", {
  d <- dat(seed = 76, n = 24)
  Z <- matrix(with_seed(1, stats::rnorm(24 * 26)), 24, 26)
  m <- bohr(y ~ x + Z, data = cbind(d, Z = I(Z)), type = "dosage")
  expect_output(print(m), "Skipped: 24 observations, fewer than the 26 covariate columns")
})

test_that("rule: skipped for the boundary share or the standard deviation", {
  d <- dat(seed = 77); d$y[1:212] <- 0                                 # 96.4% at 0
  expect_output(print(bohr(y ~ x, data = d, type = "dosage")), "Skipped: 96.4% of the observations are at 0 \\(rule: at least 95% at one boundary")
  d$y <- 0.5 + with_seed(2, stats::runif(nrow(d), -0.03, 0.03))
  expect_output(print(bohr(y ~ x, data = d, type = "dosage")), "Skipped: the standard deviation is 0.01")
})

test_that("rule: skipped without values between 0 and 1", {
  d <- dat(seed = 78); d$y <- with_seed(3, stats::rbinom(nrow(d), 1, 0.5))
  m <- bohr(y ~ x, data = d, type = "dosage")
  expect_output(print(m), "Skipped: no observations strictly between 0 and 1")
})

test_that("rule: collinear covariate columns dropped", {
  d <- dat(seed = 79); d$age2 <- 2 * d$age + 1
  m <- bohr(y ~ x + age + age2, data = d, type = "dosage")
  expect_true("collinear" %in% rule_codes(m))
  expect_output(print(m), "1 covariate column dropped as collinear with the others \\(tolerance 1e-7\\): age2")
})

test_that("rule: mu not fitted with too few values in between", {
  y <- c(rep(0, 28), rep(1, 28), 0.2, 0.5, 0.7)
  m <- bohr(y ~ x, data = data.frame(y = y, x = rep(0:2, length.out = 59)), type = "dosage")
  expect_true("mu_not_fitted" %in% rule_codes(m))
  expect_output(print(m), "mu not fitted: 3 observations between 0 and 1, but more than 3 \\(design columns \\+ 2\\) are needed")
  expect_null(m$coefficients$mu); expect_true(is.na(logLik(m)))
  expect_true(all(is.na(predict(m, type = "mu")))); expect_true(all(is.finite(predict(m, type = "p0"))))
})

test_that("rule: a constant tested term is not tested", {
  d <- dat(seed = 80); d$k <- 1
  m <- bohr(y ~ k, data = d)
  expect_true("x_constant" %in% rule_codes(m))
  expect_output(print(m), "The tested term is constant on the observations used")
})

test_that("rule: failures of the fits are reported", {
  R <- BOHR:::test_rules(list(error = "singular system", firth_error = FALSE), "BEINF")
  expect_identical(R$rule, "fit_error"); expect_match(R$message, "The fits with the tested term failed \\(singular system\\)")
  R <- BOHR:::test_rules(list(error = NA_character_, firth_error = TRUE, complete_sep = FALSE, sep_p0 = FALSE, sep_p1 = FALSE,
                              sep_mu = FALSE, conv_mu = TRUE, conv_p0 = FALSE, conv_p1 = FALSE), "BEINF")
  expect_true("firth_error" %in% R$rule); expect_match(R$message[R$rule == "firth_error"], "Firth-penalised fits failed")
})

test_that("note: complete separation by the covariates", {
  cs <- parity_cases()$complete_sep
  d <- data.frame(y = cs$y, g = cs$g, c1 = cs$Z[, 1], c2 = cs$Z[, 2])
  m <- bohr(y ~ g + c1 + c2, data = d, type = "dosage")
  expect_true("complete_separation" %in% m$notes$rule)
  expect_output(print(m), "The covariates alone separate the boundary categories completely")
})

test_that("note: separation of an unpenalised boundary fit", {
  cs <- parity_cases()$sep_beinf
  d <- data.frame(y = cs$y, g = cs$g, c1 = cs$Z[, 1], c2 = cs$Z[, 2])
  m <- bohr(y ~ g + c1 + c2, data = d, type = "dosage")
  expect_true("separation_p0" %in% m$notes$rule)
  expect_output(print(m), "p0: the unpenalised fit shows separation .* the reported test is Firth-penalised")
  m2 <- bohr(y ~ g + c1 + c2, data = d, type = "dosage", control = bohr_control(firth = FALSE))
  expect_output(print(m2), "its unpenalised test \\(firth = FALSE\\) is unreliable")
  R <- BOHR:::test_rules(list(error = NA_character_, firth_error = FALSE, complete_sep = FALSE, sep_p0 = FALSE, sep_p1 = TRUE,
                              sep_mu = FALSE, conv_mu = TRUE, conv_p0 = TRUE, conv_p1 = TRUE), "BEOI")
  expect_identical(R$rule, "separation_p1"); expect_match(R$message, "p1: the unpenalised fit shows separation")
})

test_that("note: extreme mu estimate and non-convergence", {
  base <- list(error = NA_character_, firth_error = FALSE, complete_sep = FALSE, sep_p0 = FALSE, sep_p1 = FALSE)
  R <- BOHR:::test_rules(c(base, sep_mu = TRUE, conv_mu = TRUE, conv_p0 = TRUE, conv_p1 = TRUE), "BEINF")
  expect_identical(R$rule, "separation_mu"); expect_match(R$message, "mu: extreme estimate")
  R <- BOHR:::test_rules(c(base, sep_mu = TRUE, conv_mu = FALSE, conv_p0 = FALSE, conv_p1 = TRUE), "BEZI")
  expect_identical(R$rule, c("not_converged_mu", "not_converged_p0"))
  expect_match(R$message[2], "p0: the penalised fits did not converge")
  R <- BOHR:::test_rules(c(base, sep_mu = FALSE, conv_mu = TRUE, conv_p0 = TRUE, conv_p1 = FALSE), "BEOI", firth = FALSE)
  expect_match(R$message, "p1: the unpenalised fits did not converge")
  ## a real fit stopped after one Newton step
  d <- dat(seed = 81)
  m <- bohr(y ~ x, data = d, type = "dosage", control = bohr_control(maxit = 1))
  expect_true("not_converged_mu" %in% m$notes$rule)
})

test_that("bohr() rejects terms it cannot test", {
  d <- dat(seed = 82); d$grp <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  expect_error(bohr(y ~ x, data = d, test = "z"), "not a term of the formula")
  expect_error(bohr(y ~ x * age, data = d, test = "x"), "also appears in an interaction")
  expect_error(bohr(y ~ x * age, data = d, test = "x:age"), "main effect")
  expect_error(bohr(y ~ grp + x, data = d, test = "grp"), "factor with two levels; `grp` has 3 levels")
  expect_error(bohr(y ~ poly(age, 2) + x, data = d, test = "poly(age, 2)"), "a matrix with 2 columns")
  expect_error(bohr(y ~ 0 + x, data = d), "intercept")
  expect_error(bohr(y ~ x + offset(age), data = d), "Offsets")
  expect_error(bohr(y ~ sex, data = d, type = "dosage"), "two-level factor is tested as binary")
  expect_error(bohr(y ~ x, data = d, contrast = 1), "contrast")
  expect_error(bohr(~ x, data = d), "two-sided")
  d$y[1] <- 1.5; expect_error(bohr(y ~ x, data = d), "\\[0, 1\\]")
  expect_error(bohr(y ~ x, data = transform(d, y = 0.5, x = rep(0:3, length.out = nrow(d))), type = "binary"),
               "exactly two distinct values")
})

test_that("keep_fits returns the engine's fits without changing any statistic", {
  for (nm in c("BEINF_cov", "BEZI_cov", "BEOI_nocov", "BE_cov", "sep_beinf", "few_mid")) {
    cs <- parity_cases()[[nm]]
    f <- bohr_fit(cs$y, cs$Z); g <- cs$g[f$samples]
    a <- BOHR:::snp_tests(f$y, f$design, g, f$null); b <- BOHR:::snp_tests(f$y, f$design, g, f$null, keep_fits = TRUE)
    expect_identical(b[names(a)], a, info = nm); expect_identical(setdiff(names(b), names(a)), "fits", info = nm)
    expect_identical(is.null(b$fits$beta), is.null(f$null$beta), info = nm)
    if (f$family != "BE") {
      a <- BOHR:::firth_snp_tests(f$y, f$design, g, f$null); b <- BOHR:::firth_snp_tests(f$y, f$design, g, f$null, keep_fits = TRUE)
      expect_identical(b[names(a)], a, info = nm); expect_false(is.null(b$fits$cat), info = nm)
    }
    e1 <- BOHR:::engine_one(f, g, bohr_control()); e2 <- BOHR:::engine_one(f, g, bohr_control(), keep_fits = TRUE)
    expect_identical(e2[names(e1)], e1, info = nm)
  }
})

test_that("type = 'auto' chooses binary (two values) or continuous, never dosage", {
  d <- dat(seed = 83)
  expect_identical(bohr_test(d$y, d$x)$type, "continuous")                # 0/1/2 is not detected as a dosage
  x02 <- 2 * as.numeric(d$x > 0); x01 <- as.numeric(d$x > 0); x12 <- x01 + 1
  r02 <- bohr_test(d$y, x02); r01 <- bohr_test(d$y, x01); r12 <- bohr_test(d$y, x12)
  expect_identical(c(r02$type, r01$type, r12$type), rep("binary", 3))
  expect_identical(c(r02$contrast_from, r02$contrast_to), c(0, 2))       # contrast = the two values
  expect_identical(c(r12$contrast_from, r12$contrast_to), c(1, 2))
  ## {0, 2} as binary is the dosage coding itself; {0, 1} and {1, 2} give the same model
  rd <- bohr_test(d$y, x02, type = "dosage")
  expect_identical(r02[setdiff(names(r02), "type")], rd[setdiff(names(rd), "type")])
  for (k in c("log_p_mu", "log_p_p0", "log_p_p1", "effect_mu", "effect_p0", "effect_p1"))
    expect_equal(r12[[k]], r01[[k]], tolerance = 1e-8, info = k)
  ## coefficients are per unit of x, whatever the type: x01 and x02 differ by a factor 2
  expect_equal(r01$est_mu, 2 * r02$est_mu, tolerance = 1e-8); expect_equal(r01$se_p0, 2 * r02$se_p0, tolerance = 1e-8)
  expect_equal(r01$effect_mu, r02$effect_mu, tolerance = 1e-8)
  ## dosage, continuous: same P values, coefficients per unit; effects follow the contrast
  rc <- bohr_test(d$y, d$x); rd <- bohr_test(d$y, d$x, type = "dosage")
  expect_equal(rc$log_p_p0, rd$log_p_p0, tolerance = 1e-7); expect_equal(rc$est_p0, rd$est_p0, tolerance = 1e-6)
  expect_identical(c(rd$contrast_from, rd$contrast_to), c(0, 2))
  expect_equal(c(rc$contrast_from, rc$contrast_to), c(mean(d$x), mean(d$x) + stats::sd(d$x)))
  expect_error(bohr_test(d$y, d$x, type = "binary"), "exactly two distinct values")
})
