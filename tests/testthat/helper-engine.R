## Test helpers: the original analysis engine of the BOHR paper, loaded read-only into its own
## environment, and reference pipelines built from its functions. Tests that need the engine are
## skipped when it is not present (CRAN and other users); set BOHR_ENGINE_DIR to point at a copy of
## Figure4_puQTL_landscape_eQTL_coloc/rerun_v2/code to run them elsewhere.

engine_dir <- function() {
  d <- Sys.getenv("BOHR_ENGINE_DIR",
                  "/projects/group/bliu/e1327989/GTEx/TATA_Enhancer/_tssQTL_repo/Figure4_puQTL_landscape_eQTL_coloc/rerun_v2/code")
  f <- c("bohr2_engine.R", "bohr2_addon3.R", "bohr2_data.R", "run_muF.R", "run_final2.R", "discovery_nogate.R",
         "units_v3_estimable.R", "units_v3m_mirror.R", "mirror_pairs.R", "fig2e_deliver_freeze_v3m.R")
  if (nzchar(d) && all(file.exists(file.path(d, f)))) d else NULL
}
skip_if_no_engine <- function() {
  testthat::skip_on_cran()
  if (is.null(engine_dir())) testthat::skip("original BOHR engine not available")
}

## evaluate only the top-level assignments `name <- function(...)` (optionally restricted to `funs`) and the
## constants named in `consts` of a script; nothing else of the script runs (no library() calls, no I/O)
source_defs <- function(file, env, funs = NULL, consts = character(0)) {
  for (e in parse(file, keep.source = FALSE)) {
    if (!(is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]))) next
    nm <- as.character(e[[2]])
    is_fun <- is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))
    if ((is_fun && (is.null(funs) || nm %in% funs)) || nm %in% consts) eval(e, env)
  }
  invisible(env)
}
load_engine <- function() {
  d <- engine_dir(); env <- new.env(parent = globalenv())
  sys.source(file.path(d, "bohr2_engine.R"), envir = env)          # function definitions only
  sys.source(file.path(d, "bohr2_addon3.R"), envir = env)          # function definitions only
  source_defs(file.path(d, "bohr2_data.R"), env, consts = c("FLANK", "MIN_GENO", "N_PC"))
  source_defs(file.path(d, "discovery_nogate.R"), env, funs = "acat_log")
  source_defs(file.path(d, "fig2e_deliver_freeze_v3m.R"), env, funs = "simes")
  env
}
## the reference transcribes a few data.table lines of the original scripts; check they are still there
expect_code_line <- function(file, text) {
  testthat::expect_true(any(grepl(text, readLines(file.path(engine_dir(), file), warn = FALSE), fixed = TRUE)),
                        info = paste0(file, ": ", text))
}

## The settings of the paper's engine: the package defaults without the null start of the restricted Firth fits (stage 3b-1).
## Every comparison with the original engine, or with frozen output of the original pipeline, uses it.
engine_control <- function(...) bohr_control(firth_null_start = FALSE, ...)

## Checks that hold only on the machine the fixtures were made on (same CPU, maths and BLAS libraries): bit-for-bit
## equality with frozen output, and which local maximum the paper's engine reaches. They run when BOHR_EXACT_PARITY=true.
## testthat::test_local() and devtools::test() set NOT_CRAN=true on every machine, so NOT_CRAN cannot tell them apart.
exact_parity <- function() identical(Sys.getenv("BOHR_EXACT_PARITY"), "true")

## run expr with a fixed seed and restore the caller's RNG state
with_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, globalenv()))
  set.seed(seed); expr
}

## equal values: identical, or both missing, or all.equal at tolerance `tol`
same_value <- function(a, b, tol = 1e-10) {
  if (is.null(a) || !length(a)) a <- NA
  if (is.null(b) || !length(b)) b <- NA
  if (length(a) != 1L || length(b) != 1L) return(FALSE)
  if (is.na(a) && is.na(b)) return(TRUE)
  if (is.na(a) != is.na(b)) return(FALSE)
  if (is.character(a) || is.character(b)) return(identical(as.character(a), as.character(b)))
  a <- as.numeric(a); b <- as.numeric(b)
  if (a == b) return(TRUE)
  isTRUE(all.equal(a, b, tolerance = tol))
}
## identical in value and type up to NA type (logical NA vs NA_real_)
exact_value <- function(a, b) {
  if (is.null(a) || !length(a)) a <- NA
  if (is.null(b) || !length(b)) b <- NA
  if (is.na(a) && is.na(b)) return(TRUE)
  identical(unname(a), unname(b)) || (is.numeric(a) && is.numeric(b) && !is.na(a) && !is.na(b) && as.numeric(a) == as.numeric(b))
}

## ---- reference pipeline for one promoter: run_final2.R:44-58 + 84-92 and run_muF.R:27-34 ------------------
## y: activities over all samples; Z: samples x k covariates (with column names); g: predictor over all samples
ref_promoter <- function(eng, y, Z, g = NULL) {
  n <- length(y); nm <- paste0("s", seq_len(n)); rownames(Z) <- nm
  C_all <- eng$cov_matrix(list(pcs = Z, tissue = "test"), Z[, 0, drop = FALSE], 0L)   # cov_matrix(td, peer, K = 0)
  z <- eng$prepare_y(stats::setNames(y, nm), ncol(C_all) - 1L)
  if (!is.na(z$skip)) return(list(skip = z$skip))
  if (!any(z$y > 0 & z$y < 1)) return(list(skip = "no_interior"))
  s <- names(z$y); yy <- unname(z$y); Cs <- eng$orth_design(C_all[s, , drop = FALSE])
  nul <- eng$promoter_null(yy, Cs, z$family)
  mid <- yy > 0 & yy < 1
  rk <- if (sum(mid) > 0) qr(Cs[mid, , drop = FALSE])$rank else NA_integer_                 # run_muF.R:29
  out <- list(skip = NA_character_, family = z$family, s = s, y = yy, Cs = Cs, nul = nul, rank_mid = rk, n_mid = sum(mid))
  if (is.null(g)) return(out)
  gs <- g[match(s, nm)]
  out$r <- eng$snp_tests(yy, Cs, gs, nul)
  out$f <- if (z$family != "BE") eng$firth_snp_tests(yy, Cs, gs, nul) else NULL               # main run + Firth add-on
  W <- out$r$mu_stat; n_mid <- sum(mid); p_mean <- rk + 1L
  out$mu_F_logp <- if (is.finite(W) && n_mid > p_mean + 1L)                                   # run_muF.R:34
    stats::pf((n_mid - p_mean) * expm1(pmax(W, 0) / n_mid), 1, n_mid - p_mean, lower.tail = FALSE, log.p = TRUE) else NA_real_
  out
}

## ---- synthetic promoters covering every family and the special cases ----------------------------------
sim_case <- function(seed, n, k = 2, p0 = 0.2, p1 = 0.1, d0 = 0, d1 = 0, dmu = 0, maf = 0.35, sep = 0) {
  with_seed(seed, {
    g <- stats::rbinom(n, 2, maf); Z <- matrix(stats::rnorm(n * k), n, k, dimnames = list(NULL, if (k) paste0("c", seq_len(k))))
    eta <- if (k) 0.3 * Z[, 1] else 0
    mu <- stats::plogis(stats::qlogis(0.45) + dmu * g + eta)
    y <- stats::rbeta(n, mu * 8, (1 - mu) * 8); u <- stats::runif(n)
    q0 <- pmin(pmax(p0 + d0 * g, 0), 0.9); q1 <- pmin(pmax(p1 + d1 * g, 0), 0.9)
    if (sep == 1) q0[g == 2] <- 0
    y[u < q0] <- 0; y[u >= q0 & u < q0 + q1] <- 1
    list(y = y, g = g, Z = Z)
  })
}
parity_cases <- function() {
  cs <- list(
    BEINF_cov   = sim_case(101, 150, k = 3, p0 = 0.2, p1 = 0.15, d0 = 0.08, dmu = 0.3),
    BEZI_cov    = sim_case(102, 120, k = 2, p0 = 0.25, p1 = 0, d0 = -0.08, dmu = -0.2),
    BEOI_nocov  = sim_case(103, 100, k = 0, p0 = 0, p1 = 0.3, d1 = 0.1),
    BE_cov      = sim_case(104, 90, k = 2, p0 = 0, p1 = 0, dmu = 0.4),
    endpoint    = sim_case(105, 200, k = 2, p0 = 0.2, p1 = 0.02),
    sep_beinf   = sim_case(9, 60, k = 2, p0 = 0.3, p1 = 0.15, d0 = 0.2, d1 = -0.07, sep = 1),
    small_beinf = sim_case(24, 40, k = 2, p0 = 0.15, p1 = 0.1, d0 = 0.25, d1 = 0.1),
    few_mid     = sim_case(106, 60, k = 2, p0 = 0.45, p1 = 0.45, d0 = 0.05))
  ## a covariate that separates the zeros completely (categorical null at its supremum)
  cc <- sim_case(107, 100, k = 2, p0 = 0.3, p1 = 0)
  cc$Z[, 2] <- ifelse(cc$y == 0, 2, -2) + with_seed(1, stats::runif(100, -0.5, 0.5)); cs$complete_sep <- cc
  ## a continuous trait as predictor
  tr <- sim_case(108, 130, k = 2, p0 = 0.2, p1 = 0.2); tr$g <- 2 * as.numeric(scale(with_seed(2, stats::rnorm(130)) + 0.3 * tr$y)); cs$trait <- tr
  ## collinear covariates: the third column is dropped by orth_design
  cl <- sim_case(109, 110, k = 3, p0 = 0.2, p1 = 0.1); cl$Z[, 3] <- cl$Z[, 1] - 2 * cl$Z[, 2]; cs$collinear <- cl
  ## missing activities
  na <- sim_case(110, 140, k = 2, p0 = 0.2, p1 = 0.2, dmu = 0.2); na$y[c(3, 17, 40, 41, 99)] <- NA; cs$missing_y <- na
  cs
}
skip_cases <- function() {
  a <- sim_case(201, 15, k = 2)                                     # n < 20
  b <- sim_case(202, 100, k = 2); b$y[1:97] <- 0                      # >= 95% at 0
  d <- sim_case(203, 100, k = 2, p0 = 0, p1 = 0); d$y <- 0.5 + with_seed(3, stats::runif(100, -0.02, 0.02))   # sd <= 0.05
  e <- sim_case(204, 100, k = 2); e$y <- with_seed(4, stats::rbinom(100, 1, 0.5))                         # no interior
  f <- sim_case(205, 30, k = 40)                                      # fewer samples than covariates
  list(n = a, prop_or_sd = b, sd = d, no_interior = e, few_samples = f)
}
