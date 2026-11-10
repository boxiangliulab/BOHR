## Behaviour of the small utilities and of simulate_bohr(), mirror_pairs() and bohr_cis() output.

test_that("acat() follows its definition", {
  expect_equal(acat(0.3), 0.3)
  expect_equal(acat(c(0.2, 0.2)), 0.2)                                     # equal P values combine to themselves
  expect_identical(acat(c(0, 0.5)), 0)
  expect_true(is.na(acat(c(NA, NA))))
  p <- c(0.01, 0.4, 0.9)
  t <- mean(tan((0.5 - pmin(p, 1 - 1 / 3)) * pi))
  expect_equal(acat(p), stats::pcauchy(t, lower.tail = FALSE))
  expect_equal(acat(p, cap = FALSE), stats::pcauchy(mean(tan((0.5 - p) * pi)), lower.tail = FALSE))
  expect_equal(acat(log(p), log = TRUE), acat(p))
  expect_true(acat(c(1e-20, 0.5)) < 1e-19)                                 # small-p branch
})

test_that("simes() and bb_select() follow their definitions", {
  expect_equal(simes(c(0.01, 0.04, 0.3)), min(1, min(3 * c(0.01, 0.04, 0.3) / 1:3)))
  expect_equal(simes(c(0.5, NA)), 0.5)
  expect_true(is.na(simes(NA_real_)))
  b <- bb_select(c(1e-6, 0.03, 0.4, 0.2, 0.5, 0.01, 0.002), c("A", "A", "A", "B", "B", "C", "C"))
  expect_identical(attr(b, "R"), 2L); expect_identical(attr(b, "m"), 3L)
  expect_equal(attr(b, "level"), 0.05 * 2 / 3)
  expect_identical(b$rejected, c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_identical(b$selected, c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE))
})

test_that("simulate_bohr() draws the intended distribution", {
  d <- with_seed(1, simulate_bohr(20000, freq = 0.3, mu = 0.4, sigma = 0.3, p0 = 0.2, p1 = 0.1, effect_p0 = 0.05))
  expect_named(d, c("x", "y"))
  expect_equal(mean(d$x) / 2, 0.3, tolerance = 0.02)
  expect_equal(mean(d$y[d$x == 0] == 0), 0.2, tolerance = 0.05)
  expect_equal(mean(d$y[d$x == 2] == 0), 0.3, tolerance = 0.1)
  expect_equal(mean(d$y == 1), 0.1, tolerance = 0.05)
  mid <- d$y > 0 & d$y < 1
  expect_true(all(d$y[mid] >= 0.005 & d$y[mid] <= 0.995))
  expect_equal(mean(d$y[mid & d$x == 0]), 0.4, tolerance = 0.03)
  expect_equal(attr(d, "effects")[["effect_p0"]], 0.05)
  d2 <- simulate_bohr(10, x = rep(0:1, 5))
  expect_identical(d2$x, rep(0:1, 5))
  d3 <- with_seed(2, simulate_bohr(5000, x = stats::runif(5000), p0 = 0.1, p1 = 0, effect_mu = 0.3))
  expect_gt(stats::cor(d3$x, d3$y), 0.2)
})

test_that("mirror_pairs() applies the exact and boundary-only rules", {
  a <- c(0, 1, 0.3, 0.6, 0, 1, 0.5)
  ph <- rbind(PromoterID_2 = a, PromoterID_1 = 1 - a,                      # exact pair (IDs ordered by number)
              PromoterID_3 = a, PromoterID_4 = ifelse(a %in% c(0, 1), 1 - a, 0.9 * (1 - a)),   # boundary-only
              PromoterID_5 = a, PromoterID_6 = c(0, 0, 0.5, 0.5, 0.5, 0, 0.5),               # not a mirror
              PromoterID_7 = a)                                            # single-promoter gene
  colnames(ph) <- paste0("S", 1:7)
  m <- mirror_pairs(ph, group = c("G1", "G1", "G2", "G2", "G3", "G3", "G4"))
  expect_identical(m$group, c("G1", "G2", "G3"))
  expect_identical(m$feature_1, c("PromoterID_1", "PromoterID_3", "PromoterID_5"))
  expect_identical(m$boundary_exact, c(TRUE, TRUE, FALSE))
  expect_identical(m$interior_exact, c(TRUE, FALSE, FALSE))
  expect_identical(m$canonical[1], if (mean(1 - a) >= mean(a)) "PromoterID_1" else "PromoterID_2")
  m2 <- mirror_pairs(ph, group = c("G1", "G1", "G2", "G2", "G3", "G3", "G4"), major_score = c(PromoterID_1 = 0, PromoterID_2 = 1,
                     PromoterID_3 = 1, PromoterID_4 = 1, PromoterID_5 = 1, PromoterID_6 = 0))
  expect_identical(m2$canonical, c("PromoterID_2", "PromoterID_3", "PromoterID_5"))   # ties go to the smaller ID
  ph[2, 3] <- ph[2, 3] + 1e-5
  expect_false(mirror_pairs(ph, group = c("G1", "G1", "G2", "G2", "G3", "G3", "G4"))$interior_exact[1])
  expect_true(mirror_pairs(ph, group = c("G1", "G1", "G2", "G2", "G3", "G3", "G4"), tol = 1e-4)$interior_exact[1])
})

test_that("bohr_cis() applies the estimability and mirror rules to its units", {
  x <- sim_tissue(seed = 14, n = 90)
  keep <- c("PromoterID_102", "PromoterID_103", "PromoterID_104", "PromoterID_105")
  r <- bohr_cis(x$phenotype[keep, ], x$genotype, x$variants, x$promoters, covariates = x$covariates)
  U <- r$units
  expect_true(all(U$dedup_drop[U$component == "p1" & U$boundary_exact]))
  expect_false(any(U$dedup_drop[U$component == "p0"]))
  expect_true(all(U$dedup_drop[U$component == "mu" & U$interior_exact & U$mirror_role == "partner"]))
  expect_false(any(U$dedup_drop[U$component == "mu" & !U$interior_exact]))
  expect_true(all(is.na(U$q[U$dedup_drop | !U$estimable])))
  expect_identical(U$estimable[U$component == "mu"], U$n_mid[U$component == "mu"] >= 80)
  b <- U$component != "mu"
  expect_identical(U$estimable[b], U$boundary_share[b] >= 0.1 & U$boundary_share[b] <= 0.9)
  r2 <- bohr_cis(x$phenotype[keep, ], x$genotype, x$variants, x$promoters, covariates = x$covariates, mirror = FALSE, min_n_mid = 10)
  expect_false(any(r2$units$dedup_drop))
  expect_true(all(r2$units$estimable[r2$units$component == "mu"]))
  expect_output(print(r), "BOHR cis mapping")
})

test_that("bohr_cis() warns about many covariates and keeps the genotype-count rule", {
  n <- 100; s <- paste0("S", 1:n)
  G <- with_seed(41, rbind(v1 = stats::rbinom(n, 2, 0.4), v2 = c(1, 1, rep(0, n - 2)), v3 = stats::rbinom(n, 2, 0.3)))
  colnames(G) <- s
  V <- data.frame(id = rownames(G), chr = "c1", pos = c(1000, 2000, 3000))
  y <- with_seed(42, simulate_bohr(n, x = G[1, ], p0 = 0.2, p1 = 0.1)$y)
  ph <- matrix(y, 1, n, dimnames = list("F1", s)); pr <- data.frame(id = "F1", chr = "c1", pos = 2000)
  r <- bohr_cis(ph, G, V, pr, window = 5000)
  expect_identical(r$promoters$n_variants_window, 3L)
  expect_false("v2" %in% r$snp$variant)                                       # fails the >= 5 per dosage rule
  Z <- matrix(with_seed(43, stats::rnorm(n * 4)), n, 4, dimnames = list(s, NULL))
  expect_warning(bohr_cis(ph, G, V, pr, covariates = Z, window = 5000), "4 covariates for 100 samples")
  expect_silent(bohr_cis(ph, G, V, pr, covariates = Z[, 1:3], window = 5000))
})
