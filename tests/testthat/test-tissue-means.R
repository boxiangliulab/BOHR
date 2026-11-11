## tissue_means(): the paper's rules (raw_fit() of code_orig/gamlss_functions_optimized.R:324-444), the fitted mean, inputs.

## one tissue of 60 samples; each feature exercises one rule
tm_data <- function() with_seed(71, {
  n <- 60
  beinf <- stats::rbeta(n, 3, 4); beinf[1:12] <- 0; beinf[13:20] <- 1                  # 20% at 0, 13% at 1
  bezi <- stats::rbeta(n, 2, 5); bezi[1:15] <- 0
  beoi <- stats::rbeta(n, 5, 2); beoi[1:9] <- 1
  be <- stats::rbeta(n, 2, 2)
  few <- c(rep(0, 30), rep(1, 29), 0.4)                                                # one value in between
  bound <- c(rep(0, 55), stats::rbeta(5, 2, 2))                                        # 92% at 0
  flat <- 0.5 + stats::rnorm(n, 0, 0.01)                                               # SD 0.01
  same <- rep(0.25, n)
  sparse <- c(stats::rbeta(9, 2, 2), rep(NA, n - 9))                                   # 9 values
  rbind(BEINF = beinf, BEZI = bezi, BEOI = beoi, BE = be, FEW = few, BOUND = bound, FLAT = flat, SAME = same, SPARSE = sparse)
})
ph_named <- function(m) { colnames(m) <- paste0("S", seq_len(ncol(m))); m }

test_that("tissue_means() applies the paper's rules and gives the fitted overall mean", {
  ph <- ph_named(tm_data())
  tm <- tissue_means(ph)
  expect_identical(tm$feature, rownames(ph)); expect_true(all(tm$tissue == "all"))
  r <- stats::setNames(tm$reason, tm$feature); m <- stats::setNames(tm$method, tm$feature)
  expect_identical(unname(m[c("BEINF", "BEZI", "BEOI", "BE")]), rep("fit", 4))
  expect_identical(unname(tm$family[1:4]), c("BEINF", "BEZI", "BEOI", "BE"))
  expect_identical(unname(r[c("FEW", "BOUND", "FLAT", "SAME", "SPARSE")]), c("few_mid", "boundary", "low_sd", "all_equal", "too_few"))
  expect_true(is.na(tm$mean[tm$feature == "SPARSE"]))
  for (f in c("FEW", "BOUND", "FLAT", "SAME")) expect_equal(tm$mean[tm$feature == f], mean(ph[f, ]), tolerance = 1e-14, info = f)
  ## fitted: p0, p1 = shares of 0s and 1s; mu = the maximum-likelihood beta mean of the values in between (checked with optim)
  for (f in c("BEINF", "BEZI", "BEOI", "BE")) {
    y <- ph[f, ]; w <- tm$feature == f; mid <- y[y > 0 & y < 1]
    expect_equal(tm$p0[w], mean(y == 0), tolerance = 1e-8, info = f); expect_equal(tm$p1[w], mean(y == 1), tolerance = 1e-8, info = f)
    nll <- function(th) -sum(stats::dbeta(mid, stats::plogis(th[1]) * exp(th[2]), (1 - stats::plogis(th[1])) * exp(th[2]), log = TRUE))
    o <- stats::optim(c(0, 1), nll, method = "BFGS", control = list(reltol = 1e-14, maxit = 1000))
    expect_equal(tm$mu[w], stats::plogis(o$par[1]), tolerance = 1e-5, info = f)
    expect_equal(tm$mean[w], tm$mu[w] * (1 - tm$p0[w] - tm$p1[w]) + tm$p1[w], tolerance = 1e-14, info = f)
    expect_false(isTRUE(all.equal(tm$mean[w], mean(y))))                             # not the sample mean
    expect_true(tm$converged[w])
  }
  ## each rule is an argument
  t2 <- tissue_means(ph, min_samples = 5, max_boundary = 0.95, min_sd = 0)
  expect_identical(t2$method[t2$feature %in% c("SPARSE", "BOUND", "FLAT")], rep("fit", 3))
  expect_identical(t2$mean[1:4], tm$mean[1:4])
})

test_that("tissue_means() works per tissue, from files, in parallel, and writes out_dir", {
  ph <- ph_named(cbind(tm_data(), tm_data()[, 1:30]))
  lab <- rep(c("A", "B"), c(60, 30))
  tm <- tissue_means(ph, tissue = lab)
  expect_identical(unique(tm$tissue), c("A", "B")); expect_identical(nrow(tm), 2L * nrow(ph))
  a <- tissue_means(ph[, lab == "A"]); a$tissue <- "A"
  expect_identical(tm[tm$tissue == "A", ], a)
  st <- data.frame(sample = colnames(ph), tissue = factor(lab, levels = c("B", "A")))
  expect_identical(unique(tissue_means(ph, tissue = st)$tissue), c("B", "A"))      # factor levels give the order
  d <- withr::local_tempdir()
  f <- file.path(d, "ph.tsv"); utils::write.table(data.frame(id = rownames(ph), ph, check.names = FALSE), f, sep = "\t", quote = FALSE, row.names = FALSE)
  tf <- tissue_means(f, tissue = lab, out_dir = file.path(d, "out"))
  expect_equal(tf$mean, tm$mean, tolerance = 1e-12)
  back <- utils::read.delim(file.path(d, "out", "tissue_means.tsv"), stringsAsFactors = FALSE)
  expect_identical(back$feature, tm$feature); expect_equal(back$mean, tm$mean, tolerance = 1e-12)
  skip_on_cran(); skip_on_os("windows")
  expect_identical(tissue_means(ph, tissue = lab, ncores = 2), tm)
})

test_that("tissue_means() checks its inputs", {
  ph <- ph_named(tm_data())
  expect_error(tissue_means(ph, min_samples = 0), "min_samples")
  expect_error(tissue_means(ph, max_boundary = 0), "max_boundary")
  expect_error(tissue_means(ph, min_sd = -1), "min_sd")
  expect_error(tissue_means(ph, tissue = c("A", "B")), "tissue")
  expect_error(tissue_means(ph, control = list()), "bohr_control")
  expect_error(tissue_means(ph, out_dir = c("a", "b")), "out_dir")
})
