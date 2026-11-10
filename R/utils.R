## Internal helpers: argument checks, the F correction of the mu test, and
## fast row binding. None of these is exported.

## ---- small-sample F correction of the mu likelihood-ratio test -------------------------------
## Transcribed from the analysis engine of the BOHR paper,
##   Figure4_puQTL_landscape_eQTL_coloc/rerun_v2/code/run_muF.R (md5 b1f21b17...), lines 28-34:
##   rk <- if (sum(mid) > 0) qr(Cs[mid, , drop = FALSE])$rank else NA_integer_;  p_mean = rk + 1L
##   mu_F_logp := fifelse(is.finite(W) & n_mid > p_mean + 1L,
##                        pf((n_mid - p_mean) * expm1(pmax(W, 0) / n_mid), 1, n_mid - p_mean,
##                           lower.tail = FALSE, log.p = TRUE), NA_real_)
## i.e. F = (n_mid - p) * (exp(W / n_mid) - 1) against F(1, n_mid - p), p = rank of the design on the
## interior samples + 1 (the predictor). Here pf() is evaluated only where the condition holds (the
## original evaluated it everywhere and discarded the rest); the values are identical.
mu_F_logp <- function(W, n_mid, p_mean) {
  out <- rep(NA_real_, length(W))
  ok <- is.finite(W) & !is.na(n_mid) & !is.na(p_mean) & n_mid > p_mean + 1L
  if (any(ok)) out[ok] <- pf((n_mid[ok] - p_mean[ok]) * expm1(pmax(W[ok], 0) / n_mid[ok]), 1, n_mid[ok] - p_mean[ok],
                             lower.tail = FALSE, log.p = TRUE)
  out
}
## interior rank of the per-promoter design (run_muF.R:29)
rank_mid <- function(Cs, mid) if (sum(mid) > 0) qr(Cs[mid, , drop = FALSE])$rank else NA_integer_

## ---- engine output columns ---------------------------------------------------------------------
## run_final2.R:29-35: the per-SNP columns kept from snp_tests() and firth_snp_tests().
SNP_COLS <- c("mu_stat", "nu_stat", "tau_stat", "cat_joint_stat", "joint_stat", "joint_df",
              "mu_lrt_logp", "nu_lrt_logp", "tau_lrt_logp", "joint_logp", "mu_wald_logp", "nu_wald_logp", "tau_wald_logp",
              "mu_est", "mu_se", "nu_est", "nu_se", "tau_est", "tau_se", "mu_pe", "p0_pe", "p1_pe",
              "mu_conv", "nu_conv", "tau_conv", "mu_sep", "nu_sep", "tau_sep", "mu_iter", "cat_iter",
              "ll_beta_full", "ll_cat_full", "cat_complete_sep")
FIRTH_COLS <- c("nu_firth_stat", "nu_firth_logp", "nu_firth_est", "nu_firth_se", "nu_firth_conv",
                "tau_firth_stat", "tau_firth_logp", "tau_firth_est", "tau_firth_se", "tau_firth_conv", "firth_rounds")
ENGINE_COLS <- c(SNP_COLS, FIRTH_COLS, "mu_F_logp")

## Package column name -> engine column name. The published tests are mu: F-corrected LRT;
## p0 and p1: Firth-penalised LRT (discovery_nogate.R:34-35).
CLEAR_MAP <- c(log_p_mu = "mu_F_logp", log_p_p0 = "nu_firth_logp", log_p_p1 = "tau_firth_logp",
               effect_mu = "mu_pe", effect_p0 = "p0_pe", effect_p1 = "p1_pe",
               est_mu = "mu_est", se_mu = "mu_se", est_p0 = "nu_firth_est", se_p0 = "nu_firth_se",
               est_p1 = "tau_firth_est", se_p1 = "tau_firth_se",
               stat_mu = "mu_stat", stat_p0 = "nu_firth_stat", stat_p1 = "tau_firth_stat",
               conv_mu = "mu_conv", conv_p0 = "nu_firth_conv", conv_p1 = "tau_firth_conv",
               sep_mu = "mu_sep", sep_p0 = "nu_sep", sep_p1 = "tau_sep",
               complete_sep = "cat_complete_sep")
## with control$firth = FALSE the p0/p1 columns come from the unpenalised fits of snp_tests()
CLEAR_MAP_LRT <- replace(CLEAR_MAP, c("log_p_p0", "log_p_p1", "est_p0", "se_p0", "est_p1", "se_p1", "stat_p0", "stat_p1", "conv_p0", "conv_p1"),
                         c("nu_lrt_logp", "tau_lrt_logp", "nu_est", "nu_se", "tau_est", "tau_se", "nu_stat", "tau_stat", "nu_conv", "tau_conv"))
## the per-predictor columns, in output order
TEST_COLS <- c("log_p_mu", "log_p_p0", "log_p_p1", "p_mu", "p_p0", "p_p1",
               "effect_mu", "effect_p0", "effect_p1", "est_mu", "se_mu", "est_p0", "se_p0", "est_p1", "se_p1",
               "stat_mu", "stat_p0", "stat_p1", "conv_mu", "conv_p0", "conv_p1", "sep_mu", "sep_p0", "sep_p1",
               "complete_sep", "error", "firth_error")

## ---- one predictor against a fitted null ------------------------------------------------------
## The per-SNP step of run_final2.R:84-92 (snp_tests, then firth_snp_tests when requested; a failed
## snp_tests gives no statistics and no Firth fit; a failed Firth fit gives NA Firth columns with
## nu/tau_firth_conv = FALSE), followed by the F correction of run_muF.R. Returns a named list of the
## engine columns plus `error` (character) and `firth_error` (logical). control: bohr_control().
## keep_fits = TRUE (bohr()) also returns `fits`: the full beta fit (`beta`), the full unpenalised categorical
## fit (`cat`) and the full Firth-penalised categorical fit (`firth`), each NULL when not fitted.
engine_one <- function(fit, g, control, keep_fits = FALSE) {
  out <- stats::setNames(rep(list(NA), length(ENGINE_COLS)), ENGINE_COLS)
  out$error <- NA_character_; out$firth_error <- FALSE
  r <- tryCatch(snp_tests(fit$y, fit$design, g, fit$null, maxit = control$maxit, keep_fits = keep_fits), error = function(e) e)
  if (inherits(r, "error")) { out$error <- conditionMessage(r); return(out) }
  out[SNP_COLS] <- r[SNP_COLS]
  if (keep_fits) out$fits <- list(beta = r$fits$beta, cat = r$fits$cat, firth = NULL)
  if (isTRUE(control$firth) && fit$family != "BE") {
    ff <- tryCatch(firth_snp_tests(fit$y, fit$design, g, fit$null, maxround = control$max_rounds, maxit = control$maxit,
                                   bfgs_maxit = control$bfgs_maxit, keep_fits = keep_fits,
                                   null_start = isTRUE(control$firth_null_start)), error = function(e) e)
    if (inherits(ff, "error")) { out$firth_error <- TRUE; out$nu_firth_conv <- FALSE; out$tau_firth_conv <- FALSE }
    else { out[FIRTH_COLS] <- ff[FIRTH_COLS]; if (keep_fits) out$fits$firth <- ff$fits$cat }
  }
  out$mu_F_logp <- mu_F_logp(out$mu_stat, fit$n_mid, fit$p_mean)
  out
}

## engine list -> package columns (TEST_COLS, then the engine columns when all_stats = TRUE).
## per_unit: factor from the coded predictor x' to per unit of x for the coefficients and standard errors
## (1 for a dosage, i.e. the engine's own scale; see code_predictor()).
engine_to_row <- function(e, all_stats = FALSE, firth = TRUE, per_unit = 1) {
  map <- if (firth) CLEAR_MAP else CLEAR_MAP_LRT
  row <- stats::setNames(lapply(map, function(k) e[[k]]), names(map))
  if (per_unit != 1) { for (k in c("est_mu", "est_p0", "est_p1")) row[[k]] <- row[[k]] * per_unit
    for (k in c("se_mu", "se_p0", "se_p1")) row[[k]] <- row[[k]] * abs(per_unit) }
  row$p_mu <- exp(row$log_p_mu); row$p_p0 <- exp(row$log_p_p0); row$p_p1 <- exp(row$log_p_p1)
  row$error <- e$error; row$firth_error <- e$firth_error
  row <- row[TEST_COLS]
  if (all_stats) row <- c(row, e[ENGINE_COLS])
  row
}

## ---- list of rows -> data.frame ----------------------------------------------------------------
## rows: list of named lists of length-one values. NULL entries become NA.
rows_to_df <- function(rows, cols = names(rows[[1]])) {
  if (!length(rows)) return(NULL)
  cl <- lapply(cols, function(k) unlist(lapply(rows, function(r) { v <- r[[k]]; if (is.null(v) || !length(v)) NA else v }), use.names = FALSE))
  names(cl) <- cols
  as.data.frame(cl, stringsAsFactors = FALSE, optional = TRUE)
}
## row-bind data.frames with the same columns (columns missing from one data.frame are filled with NA)
bind_df <- function(dfs) {
  dfs <- dfs[!vapply(dfs, is.null, TRUE)]
  if (!length(dfs)) return(NULL)
  cols <- unique(unlist(lapply(dfs, names)))
  cl <- lapply(cols, function(k) unlist(lapply(dfs, function(d) if (k %in% names(d)) d[[k]] else rep(NA, nrow(d))), use.names = FALSE))
  names(cl) <- cols
  as.data.frame(cl, stringsAsFactors = FALSE, optional = TRUE)
}

## ---- argument checks ---------------------------------------------------------------------------
check_activity <- function(y, what = "y") {
  if (!is.numeric(y)) stop("`", what, "` must be numeric.", call. = FALSE)
  if (any(!is.na(y) & (y < 0 | y > 1 | !is.finite(y)))) stop("`", what, "` must lie in [0, 1] (NA allowed).", call. = FALSE)
  invisible(TRUE)
}
## covariates -> numeric matrix with n rows (NULL -> n x 0)
as_covariates <- function(covariates, n, rows = NULL) {
  if (is.null(covariates)) return(matrix(numeric(0), n, 0L, dimnames = list(rows, NULL)))
  if (is.data.frame(covariates)) {
    if (!all(vapply(covariates, is.numeric, TRUE))) stop("`covariates` must be numeric; code factors as indicator columns first.", call. = FALSE)
    covariates <- as.matrix(covariates)
  }
  if (is.null(dim(covariates))) covariates <- matrix(covariates, ncol = 1L)
  if (!is.numeric(covariates)) stop("`covariates` must be a numeric matrix or data frame.", call. = FALSE)
  if (anyNA(covariates) || any(!is.finite(covariates))) stop("`covariates` must not contain missing or non-finite values.", call. = FALSE)
  if (is.null(colnames(covariates)) && ncol(covariates) > 0) colnames(covariates) <- paste0("cov", seq_len(ncol(covariates)))
  covariates
}
## the standardised design with an intercept (cov_matrix), with a clear error for constant columns
make_design <- function(Z) {
  C <- tryCatch(cov_matrix(Z), error = function(e) NULL)
  if (is.null(C)) stop("`covariates` has a constant column (after removing samples); remove it.", call. = FALSE)
  C
}
check_count <- function(x, what, min = 0) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < min || x != round(x)) stop("`", what, "` must be a single whole number >= ", min, ".", call. = FALSE)
  as.integer(x)
}
check_flag <- function(x, what) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) stop("`", what, "` must be TRUE or FALSE.", call. = FALSE)
  x
}
