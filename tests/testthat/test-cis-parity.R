## End-to-end parity of bohr_cis() with the original cis pipeline (run_final2.R, run_muF.R, discovery_nogate.R,
## units_v3_estimable.R, mirror_pairs.R, units_v3m_mirror.R) on a synthetic tissue. Skipped without the engine.

same_vec <- function(a, b, tol = 1e-10) vapply(seq_along(a), function(i) same_value(a[[i]], b[[i]], tol), TRUE)

test_that("the transcribed unit-step lines of the original scripts are unchanged", {
  skip_if_no_engine()
  for (f in names(ref_tissue_lines)) for (l in ref_tissue_lines[[f]]) expect_code_line(f, l)
})

test_that("bohr_cis() reproduces the original pipeline: promoters, every per-SNP statistic, units and q", {
  skip_if_no_engine()
  eng <- load_engine()
  x <- sim_tissue()
  ref <- ref_tissue(eng, x)
  res <- bohr_cis(x$phenotype, x$genotype, x$variants, x$promoters, covariates = x$covariates, all_stats = TRUE, control = engine_control())
  ## promoters: family, skip reason and sample counts
  P <- res$promoters[match(ref$prom$promoter, res$promoters$promoter), ]
  expect_identical(P$family, ref$prom$family)
  expect_identical(P$skip, ref$prom$skip)
  for (k in c("n", "n0", "n1")) expect_true(all(same_vec(P[[k]], ref$prom[[k]])), info = k)
  expect_true(all(c("BEINF", "BEZI", "BEOI", "BE") %in% P$family))
  expect_true(any(P$skip %in% "prop_or_sd"))
  ## per-SNP: same tested variants, every engine statistic identical
  kr <- paste(ref$snp$promoter, ref$snp$snp); kp <- paste(res$snp$promoter, res$snp$variant)
  expect_setequal(kp, kr)
  S <- res$snp[match(kr, kp), ]
  cols <- c(BOHR:::SNP_COLS, BOHR:::FIRTH_COLS, "mu_F_logp", "alt_af", "maf", "min_geno")
  for (k in cols) expect_true(all(same_vec(S[[k]], ref$snp[[k]])), info = k)
  expect_identical(S$pos, ref$snp$snp_pos)
  ## units: ACAT p, leader, estimability, mirror rule, q
  ku <- paste(ref$units$promoter, ref$units$type); kq <- paste(res$units$promoter, res$units$component)
  expect_setequal(kq, ku)
  U <- res$units[match(ku, kq), ]
  for (k in c("p", "leader_log_p", "leader_est", "leader_effect", "boundary_share", "q"))
    expect_true(all(same_vec(U[[k]], ref$units[[c(p = "p", leader_log_p = "leader_logp", leader_est = "leader_est",
                                                  leader_effect = "leader_pe", boundary_share = "boundary_share", q = "q")[[k]]]])), info = k)
  expect_identical(U$leader, ref$units$leader)
  expect_identical(U$n_variants, ref$units$n_snp)
  expect_identical(U$estimable, ref$units$estimable)
  expect_identical(U$dedup_drop, ref$units$dedup_drop)
  expect_identical(U$duplicate_of, ref$units$duplicate_of)
  expect_identical(U$boundary_exact, ref$units$boundary_exact)
  expect_identical(U$interior_exact, ref$units$interior_exact)
  ## the synthetic tissue exercises both mirror rules and both estimability rules
  expect_true(any(U$dedup_drop & U$component == "mu")); expect_true(any(U$dedup_drop & U$component == "p1"))
  expect_true(any(U$boundary_exact & !U$interior_exact))
  expect_true(any(!U$estimable & U$component == "mu")); expect_true(any(U$estimable & U$component == "mu"))
})

test_that("bohr_cis() results do not depend on ncores or chunk_size", {
  skip_on_cran(); skip_on_os("windows")
  x <- sim_tissue(seed = 12, n = 90)
  keep <- c("PromoterID_101", "PromoterID_102", "PromoterID_103", "PromoterID_106")
  ph <- x$phenotype[keep, ]
  r1 <- bohr_cis(ph, x$genotype, x$variants, x$promoters, covariates = x$covariates, ncores = 1, chunk_size = 100)
  r2 <- bohr_cis(ph, x$genotype, x$variants, x$promoters, covariates = x$covariates, ncores = 2, chunk_size = 7)
  expect_identical(r2$snp, r1$snp)
  expect_identical(r2$units, r1$units)
  expect_identical(r2$promoters, r1$promoters)
})

test_that("per-chromosome runs combined with bohr_combine() equal one run", {
  skip_on_cran()
  x <- sim_tissue(seed = 13, n = 90)
  keep <- c("PromoterID_101", "PromoterID_106", "PromoterID_111", "PromoterID_112")
  ph <- x$phenotype[keep, ]
  all <- bohr_cis(ph, x$genotype, x$variants, x$promoters, covariates = x$covariates)
  a <- bohr_cis(ph, x$genotype, x$variants, x$promoters, covariates = x$covariates, chromosomes = "chr1")
  b <- bohr_cis(ph, x$genotype, x$variants, x$promoters, covariates = x$covariates, chromosomes = "chr2")
  cb <- bohr_combine(a, b)
  o <- function(u) u[order(u$promoter, u$component), , drop = FALSE]
  expect_equal(o(cb$units), o(all$units), ignore_attr = TRUE)
  expect_identical(cb$threshold, all$threshold)
  expect_error(bohr_combine(a, a), "overlap")
})

test_that("acat() and simes() reproduce the engine's ACAT and Simes", {
  skip_if_no_engine()
  eng <- load_engine()
  with_seed(7, {
    for (m in c(1, 2, 5, 50, 400)) {
      lp <- log(stats::runif(m)); lp[sample.int(m, min(m, 2))] <- -c(40, 80)[seq_len(min(m, 2))]
      expect_identical(BOHR:::acat_log(lp), eng$acat_log(lp))
      expect_identical(acat(lp, log = TRUE), eng$acat_log(lp))
      p <- stats::runif(m)
      expect_identical(simes(p), eng$simes(p))
    }
  })
  expect_identical(BOHR:::acat_log(c(NA, -Inf, log(0.3))), eng$acat_log(c(NA, -Inf, log(0.3))))
  expect_identical(BOHR:::acat_log(c(-800, -1)), eng$acat_log(c(-800, -1)))          # exp underflow -> 0
  expect_identical(BOHR:::acat_log(-c(30, 35, 40)), eng$acat_log(-c(30, 35, 40)))    # small-p branch
})

test_that("single-tissue calls give the tables of the version before the tissue loop (ae071ae)", {
  fx <- readRDS(test_path("fixtures", "cis_ae071ae.rds"))
  cfg <- cis_v0_configs()
  for (nm in names(cfg)) {
    args <- cfg[[nm]]; if (is.null(args$control)) args$control <- engine_control()      # ae071ae: the paper's engine
    r <- suppressWarnings(do.call(bohr_cis, args))
    for (k in c("snp", "units", "promoters", "mirror", "threshold")) {
      if (identical(Sys.getenv("NOT_CRAN"), "true")) expect_identical(r[[k]], fx[[nm]][[k]], info = paste(nm, k))
      else expect_equal(r[[k]], fx[[nm]][[k]], tolerance = 1e-6, info = paste(nm, k))   # other platforms and BLAS libraries
    }
  }
})
