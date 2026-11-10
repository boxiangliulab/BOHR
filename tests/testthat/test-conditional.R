## bohr_conditional(): parity with the original conditional analysis, labels, invariances, tissues, files and errors.

cond_inputs <- function(x) list(x$phenotype, x$genotype, x$variants, x$promoters, covariates = x$covariates)
run_cis <- function(x, ...) bohr_cis(x$phenotype, x$genotype, x$variants, x$promoters, covariates = x$covariates, min_n_mid = 60, ...)
run_cond <- function(cis, x, ...) bohr_conditional(cis, x$phenotype, x$genotype, x$variants, x$promoters, covariates = x$covariates, ...)
## the reference's unit table from a bohr_cis() result (significant units, in the order of cis$units)
ref_units <- function(cis) {
  U <- cis$units[cis$units$significant, , drop = FALSE]
  data.frame(promoter = U$promoter, type = U$component, leader = U$leader, leader_logp = U$leader_log_p, p = U$p, stringsAsFactors = FALSE)
}

test_that("bohr_conditional() equals the original conditional analysis (run_conditional_v3c.R, merge_conditional_v3c.R)", {
  skip_if_no_engine()
  eng <- load_engine()
  for (f in names(ref_conditional_lines)) for (l in ref_conditional_lines[[f]]) expect_code_line(f, l)
  x <- sim_cond()
  cis <- run_cis(x)
  cond <- run_cond(cis, x)
  R <- ref_conditional(eng, x, ref_units(cis), cis$threshold)
  S <- cond$signals; RS <- R$signals
  expect_gt(nrow(S), nrow(R$units))                                   # some unit has two or more signals
  expect_true(any(S$component == "mu" & S$k == 2) && any(S$component == "p0" & S$k == 2))
  expect_identical(paste(S$promoter, S$component, S$k, S$variant), paste(RS$promoter, RS$type, RS$k, RS$snp))
  map <- c(selected_round = "selected_round", log_p = "logp_at_selection", log_p_chi = "logp_chi_at_selection", acat_p = "acat_at_selection",
           loo_est = "loo_est", loo_se = "loo_se", loo_log_p = "loo_logp", loo_log_p_chi = "loo_logp_chi", loo_acat_p = "loo_acat",
           loo_top_variant = "loo_top_snp", kept = "backward_keep")
  for (k in names(map)) expect_identical(unname(S[[k]]), unname(RS[[map[[k]]]]), info = k)
  ## units, rounds and the leave-one-out scans
  expect_identical(cond$units$n_signals, R$units$n_signals)
  expect_identical(cond$units$stop, R$units$stop)
  expect_identical(cond$units$n_kept, vapply(seq_len(nrow(R$units)), function(i)
    sum(RS$backward_keep[RS$promoter == R$units$promoter[i] & RS$type == R$units$type[i]]), 0L))
  cr <- cond$rounds; rr <- R$rounds                                   # both: round by round, units in the same order
  expect_identical(paste(cr$promoter, cr$component), paste(rr$promoter, rr$type))
  for (k in c("round", "n_tested", "n_failed", "acat_p", "acat_p_chi", "cols_dropped"))
    expect_identical(unname(cr[[k]]), unname(rr[[k]]), info = k)
  expect_identical(cr$best_log_p, rr$best_logp)
  expect_identical(cr$best_variant, rr$best_snp)
  expect_identical(cr$n_conditioned, rr$n_cond)
  expect_identical(paste(cond$loo$promoter, cond$loo$component, cond$loo$k, cond$loo$variant), paste(R$loo$promoter, R$loo$type, R$loo$k, R$loo$snp))
  for (k in c("est", "se")) expect_identical(cond$loo[[k]], R$loo[[k]], info = k)
  expect_identical(cond$loo$log_p, R$loo$logp); expect_identical(cond$loo$log_p_chi, R$loo$logp_chi)
  ## ties of the smallest conditional P go to the first variant in row order (a perfect-LD copy at an earlier position)
  expect_identical(S$variant[S$promoter == "PromoterID_201" & S$k == 2], "chr1_410000_v")
})

test_that("absolute-activity labels follow mirror_abs_labels.R and rename only O-only events", {
  skip_if_no_engine()
  eng <- load_engine()
  for (f in names(ref_abs_lines)) for (l in ref_abs_lines[[f]]) expect_code_line(f, l)
  x <- sim_cond(); cis <- run_cis(x)
  plain <- run_cond(cis, x)
  cond <- run_cond(cis, x, absolute = x$absolute)
  L <- cond$labels
  expect_identical(nrow(L), 1L)
  U <- cis$units[cis$units$significant & cis$units$component == "p0" & cis$units$boundary_exact, , drop = FALSE]
  kept <- plain$signals[plain$signals$component == "p0" & plain$signals$kept, c("promoter", "k", "variant")]
  names(kept)[3] <- "snp"
  R <- ref_abs_labels(eng, x, data.frame(S = U$promoter, O = U$mirror_partner, leader = U$leader, stringsAsFactors = FALSE), kept, x$absolute)
  expect_identical(L$z_S, R$z_S); expect_identical(L$z_O, R$z_O)
  expect_identical(L$outcome, R$outcome); expect_identical(L$label, R$label)
  expect_identical(L$outcome, "O only")
  ## the event is renamed p1 at the partner everywhere, with n0 and n1 swapped; statistics unchanged
  S <- L$S; O <- L$O
  expect_false(any(cond$units$promoter == S & cond$units$component == "p0"))
  u <- cond$units[cond$units$promoter == O & cond$units$component == "p1", ]
  u0 <- plain$units[plain$units$promoter == S & plain$units$component == "p0", ]
  expect_identical(c(u$n0, u$n1), c(u0$n1, u0$n0)); expect_identical(u$p, u0$p); expect_identical(u$relabelled_from, paste0(S, ":p0"))
  s1 <- cond$signals[cond$signals$promoter == O & cond$signals$component == "p1", ]
  s0 <- plain$signals[plain$signals$promoter == S & plain$signals$component == "p0", ]
  expect_identical(s1$log_p, s0$log_p); expect_identical(s1$variant, s0$variant)
  expect_identical(cond$signals[!(cond$signals$promoter == O), names(plain$signals)], plain$signals[!(plain$signals$promoter == S), ], ignore_attr = TRUE)
  ## when the silent promoter's absolute activity changes too, the label stays p0 at the silent promoter
  A2 <- x$absolute; A2["PromoterID_203", ] <- A2["PromoterID_203", ] - 0.5 * x$genotype["chr1_1310000_v", colnames(A2)] %in% 2
  c2 <- run_cond(cis, x, absolute = A2)
  expect_identical(c2$labels$outcome, "both"); expect_false(c2$labels$changed)
  expect_identical(c2$signals[, names(plain$signals)], plain$signals)
  expect_match(plain$settings$mirror_labels, "absolute activity not given")
  expect_null(plain$labels)
})

test_that("bohr_conditional() does not depend on ncores or chunk_size; chromosome subsets and thresholds", {
  skip_on_cran()
  x <- sim_cond(); cis <- run_cis(x)
  a <- run_cond(cis, x)
  b <- run_cond(cis, x, ncores = 2, chunk_size = 3)
  for (k in c("signals", "units", "rounds", "loo")) expect_identical(a[[k]], b[[k]], info = k)
  c2 <- run_cond(cis, x, chromosomes = "chr2")
  expect_true(all(c2$units$chr == "chr2"))
  expect_identical(c2$signals, {s <- a$signals[a$signals$promoter %in% c2$units$promoter, ]; rownames(s) <- NULL; s})
  ## a stricter threshold gives fewer signals; max_signals caps the forward step
  st <- run_cond(cis, x, threshold = cis$threshold / 1e3)
  expect_lte(nrow(st$signals), nrow(a$signals))
  expect_true(all(st$signals$acat_p[st$signals$k > 1] <= cis$threshold / 1e3))
  cap <- run_cond(cis, x, threshold = 1, max_signals = 2)
  expect_true(all(cap$units$n_signals <= 2)); expect_true(any(cap$units$stop == "cap"))
})

test_that("bohr_conditional() with tissues equals one call per tissue; out_dir writes the tables", {
  skip_on_cran()
  x1 <- sim_cond(); x2 <- sim_cond(seed = 32)
  ren <- function(m) { colnames(m) <- sub("^D", "E", colnames(m)); m }
  x2$phenotype <- ren(x2$phenotype); x2$genotype <- ren(x2$genotype); rownames(x2$covariates) <- sub("^D", "E", rownames(x2$covariates))
  all <- list(phenotype = cbind(x1$phenotype, x2$phenotype), genotype = cbind(x1$genotype, x2$genotype), variants = x1$variants,
              promoters = x1$promoters, covariates = rbind(x1$covariates, x2$covariates))
  tis <- rep(c("A", "B"), c(ncol(x1$phenotype), ncol(x2$phenotype)))
  cis <- bohr_cis(all$phenotype, all$genotype, all$variants, all$promoters, covariates = all$covariates, tissue = tis, min_n_mid = 60)
  d <- file.path(tempdir(), "bohr_cond_test"); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  cl <- bohr_conditional(cis, all$phenotype, all$genotype, all$variants, all$promoters, covariates = all$covariates, tissue = tis, out_dir = d)
  ms <- rowMeans(all$phenotype, na.rm = TRUE)
  for (tt in c("A", "B")) {
    x <- if (tt == "A") x1 else x2
    ci <- bohr_cis(x$phenotype, x$genotype, x$variants, x$promoters, covariates = x$covariates, min_n_mid = 60, major_score = ms)
    co <- run_cond(ci, x)
    s <- cl$signals[cl$signals$tissue == tt, -1]; rownames(s) <- NULL
    expect_identical(s, co$signals, info = tt)
    expect_identical(unname(cl$threshold[tt]), ci$threshold)
  }
  expect_true(all(file.exists(file.path(d, c("signals.tsv", "units.tsv", "rounds.tsv", "loo.tsv")))))
  back <- utils::read.delim(file.path(d, "signals.tsv"), stringsAsFactors = FALSE)
  expect_identical(back$variant, cl$signals$variant); expect_equal(back$log_p, cl$signals$log_p, tolerance = 1e-12)
  expect_output(print(cl), "kept")
})

test_that("bohr_conditional() matches the frozen output of the original pipeline", {
  fx <- readRDS(test_path("fixtures", "conditional_fixture.rds"))
  x <- sim_cond(); cis <- run_cis(x)
  cond <- run_cond(cis, x)
  S <- cond$signals; R <- fx$signals
  expect_identical(paste(S$promoter, S$component, S$k, S$variant), paste(R$promoter, R$type, R$k, R$snp))
  expect_equal(S$log_p, R$logp_at_selection, tolerance = 1e-6)
  expect_equal(S$acat_p, R$acat_at_selection, tolerance = 1e-6)
  expect_equal(S$loo_log_p, R$loo_logp, tolerance = 1e-6)
  expect_equal(S$loo_est, R$loo_est, tolerance = 1e-6)
  expect_equal(S$loo_acat_p, R$loo_acat, tolerance = 1e-6)
  expect_identical(S$kept, R$backward_keep)
  expect_identical(cond$units$n_signals, fx$units$n_signals)
  expect_equal(cond$loo$log_p, fx$loo$logp, tolerance = 1e-6)
})

test_that("bohr_conditional() checks its inputs", {
  set.seed(7)
  n <- 120; s <- paste0("S", 1:n)
  geno <- matrix(rbinom(20 * n, 2, 0.4), 20, n, dimnames = list(paste0("v", 1:20), s))
  vars <- data.frame(id = rownames(geno), chr = "chr1", pos = seq(1000, 20000, 1000))
  eta <- 0.5 * geno[5, ]; y <- rbeta(n, 8 * plogis(eta), 8 * (1 - plogis(eta)))
  ph <- rbind(P1 = y); colnames(ph) <- s
  proms <- data.frame(id = "P1", chr = "chr1", pos = 10000, gene = "G1")
  cis <- bohr_cis(ph, geno, vars, proms, window = 20000, min_n_mid = 50)
  expect_error(bohr_conditional(cis$units, ph, geno, vars, proms), "bohr_cis")
  expect_error(bohr_conditional(cis, ph, geno, vars, proms, threshold = c(0.1, 0.2)), "single number")
  expect_error(bohr_conditional(cis, ph, geno, vars, proms, threshold = -1), "threshold")
  expect_error(bohr_conditional(cis, ph, geno, vars, proms, max_signals = 0), "max_signals")
  expect_error(bohr_conditional(cis, ph, geno, vars, proms, tissue = rep("A", n)), "single-tissue")
  expect_error(bohr_conditional(cis, ph, geno, vars, proms, abs_z = 0), "abs_z")
  ph2 <- ph; ph2[1, 1:30] <- 0
  expect_error(bohr_conditional(cis, ph2, geno, vars, proms), "differ from the bohr_cis")
  expect_error(bohr_conditional(cis, ph, geno[-5, ], vars, proms), "variants")
  expect_error(bohr_conditional(cis, ph, geno, vars, proms, absolute = "no_such_file.tsv"), "not found")
  ok <- bohr_conditional(cis, ph, geno, vars, proms)
  expect_s3_class(ok, "bohr_conditional"); expect_true(all(ok$signals$kept))
})
