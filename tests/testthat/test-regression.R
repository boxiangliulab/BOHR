## Regression test on a frozen fixture (fixtures/engine_fixture.rds, built by data-raw/make_engine_fixture.R).
## The expected values were computed with the ORIGINAL analysis engine of the BOHR paper. Runs everywhere;
## the tolerance (1e-6, relative) allows for floating-point differences between platforms and BLAS libraries.

fx <- readRDS(test_path("fixtures", "engine_fixture.rds"))
close_vec <- function(a, b, tol = 1e-6) {
  a <- as.numeric(a); b <- as.numeric(b)
  (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & (a == b | abs(a - b) <= tol * pmax(1, abs(b))))
}

test_that("bohr_test() matches the frozen engine output (four families, separation)", {
  for (nm in names(fx$single)) {
    cs <- fx$single[[nm]]; e <- cs$engine
    r <- bohr_test(cs$input$y, cs$input$g, cs$input$Z, type = "dosage", control = engine_control())
    expect_identical(r$family, cs$family, info = nm)
    map <- BOHR:::CLEAR_MAP
    for (k in setdiff(names(map), c("conv_mu", "conv_p0", "conv_p1", "sep_mu", "sep_p0", "sep_p1", "complete_sep"))) {
      ev <- e[[map[[k]]]]; if (is.null(ev)) ev <- NA
      expect_true(close_vec(r[[k]], ev), info = paste(nm, k))
    }
    for (k in c("conv_mu", "conv_p0", "conv_p1", "complete_sep")) {
      ev <- e[[map[[k]]]]; if (is.null(ev)) ev <- NA
      expect_identical(as.logical(r[[k]]), as.logical(ev), info = paste(nm, k))
    }
  }
})

test_that("bohr_cis() matches the frozen engine output (variants, units, q)", {
  x <- fx$tissue; ref <- fx$tissue_engine
  expect_warning(res <- bohr_cis(x$phenotype, x$genotype, x$variants, x$promoters, covariates = x$covariates, all_stats = TRUE,
                                control = engine_control()),
                 "3 covariates for 80 samples")                                  # more than floor(80/30) = 2: warning only
  kr <- paste(ref$snp$promoter, ref$snp$snp); kp <- paste(res$snp$promoter, res$snp$variant)
  expect_setequal(kp, kr)
  S <- res$snp[match(kr, kp), ]
  for (k in c("mu_F_logp", "nu_firth_logp", "tau_firth_logp", "mu_stat", "nu_firth_stat", "tau_firth_stat", "mu_est", "mu_se",
              "nu_firth_est", "tau_firth_est", "mu_pe", "p0_pe", "p1_pe", "alt_af", "min_geno"))
    expect_true(all(close_vec(S[[k]], ref$snp[[k]])), info = k)
  ku <- paste(ref$units$promoter, ref$units$type); kq <- paste(res$units$promoter, res$units$component)
  expect_setequal(kq, ku)
  U <- res$units[match(ku, kq), ]
  expect_true(all(close_vec(U$p, ref$units$p)))
  expect_true(all(close_vec(U$q, ref$units$q)))
  expect_identical(U$leader, ref$units$leader)
  expect_identical(U$estimable, ref$units$estimable)
  expect_identical(U$dedup_drop, ref$units$dedup_drop)
})
