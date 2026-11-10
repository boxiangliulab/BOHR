## bohr_trait(): parity with the original trait regression, scales, rules, tissues, files and errors.

TRAITS4 <- c("SEX", "AGE", "BMI", "MHT2D")
run_trait <- function(x, ...) { inp <- trait_inputs(x); bohr_trait(x$phenotype, inp$traits, inp$covariates, promoters = x$promoters, ...) }

test_that("bohr_trait() equals the original trait regression (run_trait.R, data.R, summarize.R, mirror rule v3, option B)", {
  skip_if_no_trait_code()
  for (f in names(ref_trait_lines)) for (l in ref_trait_lines[[f]]) expect_trait_line(f, l)
  env <- load_trait_code()
  x <- sim_trait()
  res <- run_trait(x, absolute = x$absolute)
  R <- ref_trait(env, x, TRAITS4, A = x$absolute)
  P <- res$results; U <- R$units
  expect_identical(paste(P$trait, P$promoter, P$component), paste(U$trait, U$promoter, U$type))
  map <- c(p = "p", log_p = "logp", est = "est", se = "se", effect = "pe", q = "q_final", eligible = "final", significant = "sig_final",
           dedup_drop = "dedup_drop", bh_family = "bh_family", abs_relabel = "abs_relabel", n = "n", n0 = "n0", n1 = "n1", n_mid = "n_mid")
  for (k in names(map)) expect_identical(unname(P[[k]]), unname(U[[map[[k]]]]), info = k)
  ## the paper's rules are all exercised: mirror drops (exact and boundary-only pairs), relabel, ineligible units
  expect_true(any(P$dedup_drop & P$component == "p1")); expect_true(any(P$dedup_drop & P$component == "mu"))
  expect_true(any(!P$eligible & !P$dedup_drop))
  expect_true(all(is.na(P$mirror_role[!P$boundary_exact]))); expect_true(all(!is.na(P$mirror_role[P$boundary_exact])))
  expect_identical(sum(P$abs_relabel), 2L)
  L <- res$labels
  expect_identical(L$S, R$labels$S); expect_identical(L$z_S, unname(R$labels$z_S)); expect_identical(L$z_O, unname(R$labels$z_O))
  expect_identical(L$outcome, "O only")
  ## the relabelled event: p1 at the partner carries the q of the p0 test; p0 at the silent promoter is the duplicate
  o <- P[P$abs_relabel & P$component == "p1", ]; s <- P[P$abs_relabel & P$component == "p0", ]
  expect_identical(o$bh_family, "p0"); expect_true(o$significant); expect_false(s$eligible); expect_identical(s$duplicate_of, paste0(o$promoter, ":p1"))
  ## without absolute activity: the v3m labels, the same P values and q values
  plain <- run_trait(x)
  R0 <- ref_trait(env, x, TRAITS4)$units
  expect_identical(plain$results$q, R0$q_final); expect_false(any(plain$results$abs_relabel))
  expect_match(plain$settings$mirror_labels, "absolute activity not given"); expect_null(plain$labels)
})

test_that("bohr_trait() reports estimates per contrast; P values equal bohr_test()", {
  set.seed(5)
  n <- 160; s <- paste0("S", 1:n)
  age <- round(runif(n, 20, 70)); sex <- rbinom(n, 1, 0.5); pcs <- matrix(rnorm(2 * n), n, 2, dimnames = list(s, c("PC1", "PC2")))
  y <- rbeta(n, 3 + 0.04 * age, 4); y[runif(n) < plogis(-1.5 + 0.8 * sex)] <- 0
  ph <- rbind(P1 = y); colnames(ph) <- s
  res <- bohr_trait(ph, data.frame(AGE = age, SEX = sex, row.names = s), as.data.frame(pcs), min_n_mid = 30, min_donors = 20)
  r <- res$results
  b <- bohr_test(y, age, covariates = pcs, type = "continuous", all_stats = TRUE)
  ra <- r[r$trait == "AGE" & r$component == "mu", ]
  expect_identical(ra$p, b$p_mu)
  expect_equal(ra$est, b$est_mu * sd(age), tolerance = 1e-12); expect_equal(ra$se, b$se_mu * sd(age), tolerance = 1e-12)
  expect_identical(ra$est, 2 * b$mu_est)                                           # 2 x the coded estimate, as run_trait.R
  expect_identical(ra$effect, b$effect_mu)
  bs <- bohr_test(y, sex, covariates = pcs, type = "binary")
  rs <- r[r$trait == "SEX" & r$component == "p0", ]
  expect_identical(rs$p, bs$p_p0); expect_equal(rs$est, bs$est_p0, tolerance = 1e-12)  # 0 -> 1: per unit = per contrast
  expect_equal(res$traits$contrast_to - res$traits$contrast_from, c(sd(age), 1), tolerance = 1e-12)
})

test_that("bohr_trait() applies the trait, feature and eligibility rules", {
  x <- sim_trait(n = 120)
  inp <- trait_inputs(x)
  tr <- inp$traits
  tr$RARE <- 0; tr$RARE[1:4] <- 1                                            # 4 donors in the minority group
  tr$FEW <- NA_real_; tr$FEW[1:15] <- stats::rnorm(15)                       # fewer than 20 donors
  ph <- x$phenotype; t2 <- which(tr$MHT2D %in% 1); ph["PromoterID_310", t2[-(1:3)]] <- NA     # 3 exposed donors left for P310
  res <- bohr_trait(ph, tr, inp$covariates, promoters = x$promoters, min_donors = 150)
  T <- res$traits
  expect_identical(T$skip[T$trait == "RARE"], "trait_rare"); expect_identical(T$skip[T$trait == "FEW"], "trait_constant_or_n")
  expect_false(any(res$results$eligible))                                    # every tissue x trait has < 150 donors
  expect_true(all(is.na(res$results$q)))
  f <- res$features
  expect_identical(f$skip[f$trait == "MHT2D" & f$promoter == "PromoterID_310"], "trait_rare_in_feature")
  ## a two-level character trait is the indicator of its second level ("male"): SEX reversed, the same P values
  cs <- data.frame(SEX = tr$SEX, CAT = ifelse(tr$SEX == 1, "female", "male"), row.names = rownames(tr))
  rc <- bohr_trait(x$phenotype, cs, inp$covariates[, -1], min_donors = 0)$results
  a <- rc[rc$trait == "SEX", ]; b <- rc[rc$trait == "CAT", ]
  expect_equal(a$p, b$p, tolerance = 1e-8); expect_equal(a$est, -b$est, tolerance = 1e-5)
  ## the tested trait is not a covariate; Hardy-scale levels below min_level_count are merged
  expect_false(grepl("(^|;)SEX(;|$)", T$covariates[T$trait == "SEX"])); expect_true(grepl("(^|;)SEX(;|$)", T$covariates[T$trait == "AGE"]))
  r2 <- bohr_trait(x$phenotype, inp$traits[, "AGE", drop = FALSE], inp$covariates, min_level_count = 1)
  expect_gt(r2$traits$n_covariates, res$traits$n_covariates[res$traits$trait == "AGE"])
  ## a feature-level rare group skips the feature; trait_type overrides the automatic type
  r3 <- bohr_trait(x$phenotype, inp$traits[, "SEX", drop = FALSE], inp$covariates, trait_type = c(SEX = "continuous"))
  expect_identical(r3$traits$type, "continuous")
  ## exclude_separated and boundary_share change only eligibility
  r4 <- bohr_trait(x$phenotype, inp$traits[, "SEX", drop = FALSE], inp$covariates, promoters = x$promoters, boundary_share = c(0, 1), min_donors = 0)
  r5 <- bohr_trait(x$phenotype, inp$traits[, "SEX", drop = FALSE], inp$covariates, promoters = x$promoters, min_donors = 0)
  expect_identical(r4$results$p, r5$results$p); expect_gte(sum(r4$results$eligible), sum(r5$results$eligible))
})

test_that("bohr_trait() with tissues equals one call per tissue; traits by individual; out_dir", {
  skip_on_cran()
  x <- sim_trait(); inp <- trait_inputs(x)
  d <- colnames(x$phenotype)
  ## two tissues: the first and second half of the donors, samples renamed, traits and covariates given per individual
  smp <- data.frame(sample = paste0(ifelse(seq_along(d) <= 90, "A_", "B_"), d), tissue = ifelse(seq_along(d) <= 90, "A", "B"), individual = d,
                    stringsAsFactors = FALSE)
  ph <- x$phenotype; colnames(ph) <- smp$sample
  ph[c("PromoterID_305", "PromoterID_310"), smp$tissue == "B"] <- NA         # not measured in tissue B: no unit there
  td <- file.path(tempdir(), "bohr_trait_test"); on.exit(unlink(td, recursive = TRUE), add = TRUE)
  res <- bohr_trait(ph, inp$traits, inp$covariates, tissue = smp, promoters = x$promoters, min_donors = 20, out_dir = td)
  expect_false(any(res$results$tissue == "B" & res$results$promoter %in% c("PromoterID_305", "PromoterID_310")))
  ms <- rowMeans(ph, na.rm = TRUE)
  for (tt in c("A", "B")) {
    w <- smp$tissue == tt
    tw <- inp$traits[smp$individual[w], ]; cw <- inp$covariates[smp$individual[w], ]; rownames(tw) <- rownames(cw) <- smp$sample[w]
    one <- bohr_trait(ph[, w], tw, cw, promoters = x$promoters, min_donors = 20, major_score = ms)
    r <- res$results[res$results$tissue == tt, -1]; rownames(r) <- NULL
    expect_identical(r, one$results, info = tt)
  }
  expect_true(all(file.exists(file.path(td, c("results.tsv", "features.tsv", "traits.tsv")))))
  back <- utils::read.delim(file.path(td, "results.tsv"), stringsAsFactors = FALSE)
  expect_identical(nrow(back), nrow(res$results)); expect_equal(back$p, res$results$p, tolerance = 1e-12)
  expect_output(print(res), "within tissue x trait x component")
})

test_that("bohr_trait() matches the frozen output of the original trait regression", {
  fx <- readRDS(test_path("fixtures", "trait_fixture.rds"))
  x <- sim_trait()
  res <- run_trait(x, absolute = x$absolute)
  P <- res$results; U <- fx$units
  expect_identical(paste(P$trait, P$promoter, P$component), paste(U$trait, U$promoter, U$type))
  expect_equal(P$p, U$p, tolerance = 1e-6); expect_equal(P$est, U$est, tolerance = 1e-6); expect_equal(P$effect, U$pe, tolerance = 1e-6)
  expect_equal(P$q, U$q_final, tolerance = 1e-6)
  expect_identical(P$eligible, U$final); expect_identical(P$significant, U$sig_final); expect_identical(P$abs_relabel, U$abs_relabel)
})

test_that("bohr_trait() checks its inputs", {
  x <- sim_trait(n = 60); inp <- trait_inputs(x)
  expect_error(bohr_trait(x$phenotype, inp$traits[1:10, ], inp$covariates), "row for every sample")
  expect_error(bohr_trait(x$phenotype, inp$traits, inp$covariates, trait_type = c(SEX = "ordinal")), "binary")
  expect_error(bohr_trait(x$phenotype, inp$traits, inp$covariates, trait_type = c(AGE = "binary")), "declared binary")
  tr3 <- data.frame(G = sample(c("a", "b", "c"), 60, TRUE), row.names = rownames(inp$traits))
  expect_error(bohr_trait(x$phenotype, tr3, inp$covariates), "exactly two")
  expect_error(bohr_trait(x$phenotype, inp$traits, inp$covariates, promoters = data.frame(id = rownames(x$phenotype))), "id and gene")
  expect_error(bohr_trait(x$phenotype, inp$traits, inp$covariates, min_group = 0), "min_group")
  expect_error(bohr_trait(x$phenotype, inp$traits, inp$covariates, boundary_share = c(0.9, 0.1)), "boundary_share")
  expect_error(bohr_trait(x$phenotype, inp$traits, inp$covariates, fdr = 1), "fdr")
  expect_error(bohr_trait(x$phenotype, inp$traits, inp$covariates, abs_z = -1), "abs_z")
  expect_error(bohr_trait(x$phenotype, "no_such_file.tsv"), "not found")
  expect_error(bohr_trait(x$phenotype, data.frame(row.names = rownames(inp$traits))), "no columns")
})
