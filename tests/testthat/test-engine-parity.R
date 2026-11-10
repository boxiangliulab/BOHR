`%||%` <- function(a, b) if (is.null(a)) b else a
## Engine parity: the package must reproduce the original analysis engine of the BOHR paper exactly
## (identical values, or all.equal at tolerance 1e-10) for every statistic of bohr_fit() and bohr_test().
## The original engine is sourced read-only from the paper's storage repository; skipped when absent.

fit_parts <- c("ll", "coef", "vcov", "converged", "iter", "maxgrad")

test_that("transcribed lines of the original scripts are unchanged", {
  skip_if_no_engine()
  expect_code_line("run_muF.R", "rk <- if (sum(mid) > 0) qr(Cs[mid, , drop = FALSE])$rank else NA_integer_")
  expect_code_line("run_muF.R", "pf((n_mid - p_mean) * expm1(pmax(W, 0) / n_mid), 1, n_mid - p_mean, lower.tail = FALSE, log.p = TRUE)")
  expect_code_line("run_final2.R", "z <- prepare_y(A[, p], ncol(C_all) - 1L)")
  expect_code_line("run_final2.R", "s <- names(z$y); y <- unname(z$y); Cs <- orth_design(C_all[s, , drop = FALSE])")
  expect_code_line("run_final2.R", "nul <- promoter_null(y, Cs, z$family)")
  expect_code_line("run_final2.R", "r <- tryCatch(snp_tests(P$y, P$Cs, P$G[, j], P$nul), error = function(e) NULL)")
  expect_code_line("run_final2.R", "ff <- tryCatch(firth_snp_tests(P$y, P$Cs, P$G[, j], P$nul), error = function(e) NULL)")
  expect_code_line("discovery_nogate.R", 'PCOL <- c(mu = "mu_F_logp", p0 = "nu_firth_logp", p1 = "tau_firth_logp")')
})

test_that("bohr_fit() reproduces the engine's null fits (all four families, separation, special cases)", {
  skip_if_no_engine()
  eng <- load_engine()
  cases <- parity_cases()
  fams <- character(0)
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    ref <- ref_promoter(eng, cs$y, cs$Z)
    fit <- bohr_fit(cs$y, cs$Z)
    expect_identical(fit$skip, ref$skip, info = nm)
    expect_identical(fit$family, ref$family, info = nm)
    fams <- c(fams, fit$family)
    expect_identical(fit$y, ref$y, info = nm)
    expect_identical(fit$design, ref$Cs, info = nm)
    expect_identical(fit$samples, match(ref$s, paste0("s", seq_along(cs$y))), info = nm)
    expect_identical(fit$rank_mid, ref$rank_mid, info = nm)
    expect_identical(fit$n_mid, ref$n_mid, info = nm)
    for (part in c("beta", "cat")) {
      a <- fit$null[[part]]; b <- ref$nul[[part]]
      expect_identical(is.null(a), is.null(b), info = paste(nm, part))
      if (!is.null(b)) for (k in fit_parts) expect_true(isTRUE(all.equal(a[[k]], b[[k]], tolerance = 1e-10)) || identical(a[[k]], b[[k]]),
                                                         info = paste(nm, part, k))
    }
    expect_identical(fit$null[c("family", "n", "n_mid", "n0", "n1")], ref$nul[c("family", "n", "n_mid", "n0", "n1")], info = nm)
  }
  expect_setequal(unique(fams), c("BEINF", "BEZI", "BEOI", "BE"))
})

test_that("bohr_test() reproduces every engine statistic (mu F test, Firth p0/p1, effects, flags)", {
  skip_if_no_engine()
  eng <- load_engine()
  cases <- parity_cases()
  snp_cols <- BOHR:::SNP_COLS; firth_cols <- BOHR:::FIRTH_COLS
  n_cmp <- 0L; n_exact <- 0L
  seen <- list(sep = FALSE, complete_sep = FALSE, beta_null = FALSE)
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    ref <- ref_promoter(eng, cs$y, cs$Z, cs$g)
    tst <- bohr_test(cs$y, cs$g, cs$Z, contrast = c(0, 2), all_stats = TRUE, control = engine_control())   # the engine's own coding of x
    expect_identical(nrow(tst), 1L)
    for (k in snp_cols) { n_cmp <- n_cmp + 1L; n_exact <- n_exact + exact_value(tst[[k]], ref$r[[k]])
      expect_true(same_value(tst[[k]], ref$r[[k]]), info = paste(nm, k)) }
    for (k in firth_cols) { b <- if (is.null(ref[["f"]])) NA else ref[["f"]][[k]]; n_cmp <- n_cmp + 1L; n_exact <- n_exact + exact_value(tst[[k]], b)
      expect_true(same_value(tst[[k]], b), info = paste(nm, k)) }
    n_cmp <- n_cmp + 1L; n_exact <- n_exact + exact_value(tst$mu_F_logp, ref$mu_F_logp)
    expect_true(same_value(tst$mu_F_logp, ref$mu_F_logp), info = paste(nm, "mu_F_logp"))
    ## the package's named columns are the engine's published statistics
    map <- BOHR:::CLEAR_MAP
    for (k in names(map)) {
      b <- if (map[[k]] == "mu_F_logp") ref$mu_F_logp else if (map[[k]] %in% firth_cols) (if (is.null(ref[["f"]])) NA else ref[["f"]][[map[[k]]]]) else ref$r[[map[[k]]]]
      expect_true(same_value(tst[[k]], b), info = paste(nm, k))
    }
    expect_true(same_value(tst$p_mu, exp(ref$mu_F_logp)), info = nm)
    if (isTRUE(ref$r$nu_sep) || isTRUE(ref$r$tau_sep)) seen$sep <- TRUE
    if (isTRUE(ref$r$cat_complete_sep)) seen$complete_sep <- TRUE
    if (is.null(ref$nul$beta)) seen$beta_null <- TRUE
  }
  expect_true(seen$sep)            # a separated case was covered
  expect_true(seen$complete_sep)   # covariates separating the categories completely
  expect_true(seen$beta_null)      # too few intermediate samples for the beta part
  expect_gt(n_exact / n_cmp, 0.99) # in one R session the values are bit-identical
})

test_that("bohr_control(firth = FALSE) reports the engine's unpenalised p0/p1 tests", {
  skip_if_no_engine()
  eng <- load_engine()
  cs <- parity_cases()$BEINF_cov
  ref <- ref_promoter(eng, cs$y, cs$Z, cs$g)
  tst <- bohr_test(cs$y, cs$g, cs$Z, control = bohr_control(firth = FALSE))
  expect_true(same_value(tst$log_p_p0, ref$r$nu_lrt_logp)); expect_true(same_value(tst$log_p_p1, ref$r$tau_lrt_logp))
  expect_true(same_value(tst$est_p0, ref$r$nu_est)); expect_true(same_value(tst$se_p1, ref$r$tau_se))
  expect_true(same_value(tst$log_p_mu, ref$mu_F_logp))
})

test_that("skipped promoters are skipped for the same reason as in the engine", {
  skip_if_no_engine()
  eng <- load_engine()
  for (nm in names(skip_cases())) {
    cs <- skip_cases()[[nm]]
    ref <- ref_promoter(eng, cs$y, cs$Z)
    fit <- bohr_fit(cs$y, cs$Z)
    expect_identical(fit$skip, ref$skip, info = nm)
    tst <- bohr_test(cs$y, cs$g, cs$Z)
    expect_true(is.na(tst$log_p_mu) && is.na(tst$log_p_p0) && is.na(tst$log_p_p1), info = nm)
  }
})

test_that("binary and continuous predictors follow the engine with the trait-analysis coding", {
  skip_if_no_engine()
  eng <- load_engine()
  ## run_trait.R codes a trait as trait_contrast(): 2 * g (binary) or 2 * (g - mean(g)) / sd(g) (continuous), over all
  ## samples with the trait, then calls snp_tests / firth_snp_tests; it reports coefficients x 2 (per contrast).
  cs <- parity_cases()$BEINF_cov
  n <- length(cs$y)
  b <- with_seed(31, stats::rbinom(n, 1, 0.4)); z <- with_seed(32, stats::rnorm(n, 50, 12))
  for (case in list(list(x = b, type = "binary", coded = 2 * b, per_unit = 2),
                    list(x = z, type = "continuous", coded = 2 * (z - mean(z)) / stats::sd(z), per_unit = 2 / stats::sd(z)))) {
    ref <- ref_promoter(eng, cs$y, cs$Z, case$coded)
    tst <- bohr_test(cs$y, case$x, cs$Z, all_stats = TRUE, control = engine_control())
    expect_identical(tst$type, case$type)
    for (k in c(BOHR:::SNP_COLS, BOHR:::FIRTH_COLS)) expect_true(same_value(tst[[k]], ref$r[[k]] %||% ref[["f"]][[k]]), info = paste(case$type, k))
    expect_true(same_value(tst$log_p_mu, ref$mu_F_logp)); expect_true(same_value(tst$log_p_p0, ref[["f"]]$nu_firth_logp))
    expect_true(same_value(tst$effect_mu, ref$r$mu_pe)); expect_true(same_value(tst$effect_p1, ref$r$p1_pe))
    expect_true(same_value(tst$est_mu, ref$r$mu_est * case$per_unit)); expect_true(same_value(tst$se_p0, ref[["f"]]$nu_firth_se * case$per_unit))
  }
  ## the dosage path is the engine itself (no recoding), and contrast c(0, 2) is the same
  d <- bohr_test(cs$y, cs$g, cs$Z, type = "dosage"); d2 <- bohr_test(cs$y, cs$g, cs$Z, contrast = c(0, 2))
  expect_identical(d, transform(d2, type = "dosage"))
  ## P values do not depend on the contrast; effects and coefficients follow it
  r1 <- bohr_test(cs$y, z, cs$Z); r2 <- bohr_test(cs$y, z, cs$Z, contrast = c(40, 60))
  expect_equal(r2$log_p_mu, r1$log_p_mu, tolerance = 1e-6); expect_equal(r2$est_mu, r1$est_mu, tolerance = 1e-6)
  expect_identical(c(r2$contrast_from, r2$contrast_to), c(40, 60))
})
