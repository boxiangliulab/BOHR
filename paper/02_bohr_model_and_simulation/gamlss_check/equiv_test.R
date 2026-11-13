## Stage 1 (PILOT_PROTOCOL.md section 2): is the split-likelihood engine the same estimator as gamlss
## with sigma ~ 1? Per SNP, five models are fitted both ways -- null, full, and full minus the SNP in
## mu / nu / tau -- and the maximised log-likelihoods and the LRT statistics are compared.
## Criterion (written before running): where both converge, |dl| <= 1e-3 and |dLRT| <= 0.01; the
## engine below gamlss by > 1e-3 is a failure; gamlss below the engine is recorded (and refitted with
## a tight gamlss control to see whether gamlss was simply stopped early).
## Usage: Rscript equiv_test.R <ncores>
suppressPackageStartupMessages({ library(gamlss); library(parallel) })
source("/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/code/bohr2_data.R")
source("/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/code/bohr2_engine.R")
NC <- as.integer(commandArgs(TRUE)[1]); if (is.na(NC)) NC <- 1L
OUT <- file.path(V2, "results", if (length(commandArgs(TRUE)) >= 2) commandArgs(TRUE)[2] else "stage1_equiv"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
TISS <- c("Kidney_Cortex", "Heart_Left_Ventricle", "Skin_Sun_Exposed_Lower_leg")
PER_FAM <- 3L; N_SNP <- 10L
set.seed(20260930)

## ---- gamlss side --------------------------------------------------------------------------------
gfit <- function(df, fam, xn, mu_g, nu_g, tau_g, tight = FALSE) {
  rhs <- function(inc) paste(c(xn, if (inc) "g"), collapse = " + ")
  ctrl <- if (tight) gamlss.control(n.cyc = 2000, c.crit = 1e-8, trace = FALSE) else gamlss.control(n.cyc = 200, trace = FALSE)
  fo <- switch(fam,
    BEINF = BEINF(mu.link = "logit", sigma.link = "logit", nu.link = "log", tau.link = "log"),
    BEZI  = BEZI(mu.link = "logit", sigma.link = "log", nu.link = "logit"),
    BEOI  = BEOI(mu.link = "logit", sigma.link = "log", nu.link = "logit"),
    BE    = BE(mu.link = "logit", sigma.link = "logit"))
  a <- list(formula = as.formula(paste("y ~", rhs(mu_g))), sigma.formula = ~1, family = fo, data = df, control = ctrl)
  if (fam == "BEINF") { a$nu.formula <- as.formula(paste("~", rhs(nu_g))); a$tau.formula <- as.formula(paste("~", rhs(tau_g))) }
  if (fam == "BEZI") a$nu.formula <- as.formula(paste("~", rhs(nu_g)))
  if (fam == "BEOI") a$nu.formula <- as.formula(paste("~", rhs(tau_g)))      # gamlss BEOI nu = P(y = 1)
  tm <- system.time(m <- tryCatch(suppressWarnings(do.call(gamlss, a)), error = function(e) NULL))[["elapsed"]]
  if (is.null(m)) return(list(ll = NA_real_, conv = NA, sec = tm, m = NULL))
  list(ll = as.numeric(logLik(m)), conv = isTRUE(m$converged), sec = tm, m = m)
}
g_wald <- function(m) {                     # SNP estimate / SE per parameter, gamlss default vcov
  if (is.null(m)) return(NULL)
  est <- tryCatch(vcov(m, type = "coef"), error = function(e) NULL)
  se  <- tryCatch(vcov(m, type = "se"), error = function(e) NULL)
  if (is.null(est) || is.null(se)) return(NULL)
  k <- which(names(est) == "g"); list(est = unname(est[k]), se = unname(se[k]))
}

## ---- one SNP: both implementations ------------------------------------------------------------
one_snp <- function(task) {
  y <- task$y; C <- task$C; g <- task$g; fam <- task$family; nul <- task$nul
  t_e <- system.time(r <- snp_tests(y, C, g, nul))[["elapsed"]]
  ## engine log-likelihoods of the five joint models (parts are additive)
  lb0 <- if (!is.null(nul$beta)) nul$beta$ll else NA_real_
  lc0 <- if (!is.null(nul$cat)) nul$cat$ll else 0
  lb1 <- r$ll_beta_full; lc1 <- if (fam == "BE") 0 else r$ll_cat_full
  e <- c(null = lb0 + lc0, full = lb1 + lc1, dmu = lb0 + lc1,
         dnu = if (fam %in% c("BEINF", "BEZI")) lb1 + lc1 - r$nu_stat / 2 else NA,
         dtau = if (fam %in% c("BEINF", "BEOI")) lb1 + lc1 - r$tau_stat / 2 else NA)
  df <- data.frame(y = y, C[, -1, drop = FALSE], g = g); xn <- colnames(C)[-1]
  df <- setNames(df, c("y", paste0("x", seq_along(xn)), "g")); xn <- paste0("x", seq_along(xn))
  spec <- list(null = c(F, F, F), full = c(T, T, T), dmu = c(F, T, T), dnu = c(T, F, T), dtau = c(T, T, F))
  keep <- names(e)[is.finite(e) | names(e) %in% c("null", "full", "dmu")]
  if (fam == "BE") keep <- c("null", "full")
  G <- lapply(spec[keep], function(s) gfit(df, fam, xn, s[1], s[2], s[3]))
  wald_g <- g_wald(G$full$m)
  rows <- lapply(keep, function(k) {
    gl <- G[[k]]$ll; d <- e[[k]] - gl; tight <- NA_real_
    if (is.finite(d) && abs(d) > 1e-3) tight <- gfit(df, fam, xn, spec[[k]][1], spec[[k]][2], spec[[k]][3], tight = TRUE)$ll
    data.frame(model = k, ll_engine = e[[k]], ll_gamlss = gl, gamlss_conv = G[[k]]$conv,
               ll_gamlss_tight = tight, gamlss_sec = G[[k]]$sec)
  })
  M <- do.call(rbind, rows)
  M <- cbind(task$meta, family = fam, M, engine_sec = t_e,
             mu_conv = r$mu_conv, nu_conv = r$nu_conv, tau_conv = r$tau_conv,
             mu_sep = r$mu_sep, nu_sep = r$nu_sep, tau_sep = r$tau_sep,
             mu_est = r$mu_est, mu_se = r$mu_se, nu_est = r$nu_est, nu_se = r$nu_se,
             tau_est = r$tau_est, tau_se = r$tau_se,
             g_wald_est = if (!is.null(wald_g)) paste(signif(wald_g$est, 6), collapse = ";") else NA,
             g_wald_se = if (!is.null(wald_g)) paste(signif(wald_g$se, 6), collapse = ";") else NA,
             null_beta_conv = if (!is.null(nul$beta)) nul$beta$converged else NA,
             null_cat_conv = if (!is.null(nul$cat)) nul$cat$converged else NA)
  M
}

## ---- tasks: real promoters (all four families where present) ----------------------------------
tasks <- list(); scan_log <- list()
for (tis in TISS) {
  td <- readRDS(file.path(PATHS$work, "tissue", paste0(tis, ".rds")))
  C <- cov_matrix(td, td$peer, td$k_now)
  pr <- colnames(td$act)[td$info[colnames(td$act), "seqnames"] %in% CHRS]
  famof <- vapply(pr, function(p) { z <- prepare_y(td$act[, p], ncol(C) - 1L); if (is.na(z$skip)) z$family else NA_character_ }, "")
  pick <- unlist(lapply(c("BEINF", "BEZI", "BEOI", "BE"), function(f) { w <- names(famof)[which(famof == f)]
    if (length(w) > PER_FAM) sample(w, PER_FAM) else w }))
  for (p in pick) {
    z <- prepare_y(td$act[, p], ncol(C) - 1L); s <- names(z$y)
    chr <- td$info[p, "seqnames"]; pos <- td$info[p, "start"]
    geno <- read_geno(chr, pos - FLANK, pos + FLANK - 1L, td$U)
    if (is.null(geno)) next
    idx <- cis_index(geno, pos); Gs <- geno$G[s, idx, drop = FALSE]
    cnt <- apply(Gs, 2, function(x) min(tabulate(x + 1L, nbins = 3L)))
    Gs <- Gs[, cnt >= MIN_GENO, drop = FALSE]
    if (!ncol(Gs)) next
    Cs <- orth_design(C[s, , drop = FALSE]); y <- unname(z$y)
    nul <- promoter_null(y, Cs, z$family)
    ## N_SNP random SNPs + the SNP with the smallest engine joint p (a real-signal case)
    t_scan <- system.time(jp <- unlist(mclapply(seq_len(ncol(Gs)), function(j) snp_tests(y, Cs, Gs[, j], nul)$joint_logp,
                                                mc.cores = NC)))[["elapsed"]]
    scan_log[[length(scan_log) + 1]] <- data.frame(tissue = tis, promoter = p, family = z$family, n = length(y),
                                                   n_cov = ncol(Cs) - 1L, n_snp = ncol(Gs), sec = t_scan, cores = NC)
    js <- unique(c(sample(ncol(Gs), min(N_SNP, ncol(Gs))), which.min(jp)))
    for (j in js) tasks[[length(tasks) + 1]] <- list(y = y, C = Cs, g = Gs[, j], family = z$family, nul = nul,
      meta = data.frame(source = "real", tissue = tis, promoter = p, snp = colnames(Gs)[j], n = length(y),
                        n0 = sum(y == 0), n1 = sum(y == 1), n_cov = ncol(Cs) - 1L, n_cov_dropped = length(attr(Cs, "dropped")),
                        top = j == which.min(jp)))
  }
  cat(tis, ": promoters", length(pick), "tasks so far", length(tasks), "\n")
}

## ---- tasks: simulated data, every family, known effects ----------------------------------------
sim_one <- function(fam, n, rep) {
  X <- cbind(1, matrix(rnorm(n * 5), n)); g <- rbinom(n, 2, 0.3)
  mu <- plogis(drop(X %*% c(-0.3, 0.2, -0.1, 0.1, 0, 0.05)) + 0.25 * g); phi <- 12
  y <- rbeta(n, mu * phi, (1 - mu) * phi)
  e0 <- drop(X %*% c(-1.5, 0.3, 0, 0.1, 0, 0)) + 0.3 * g; e1 <- drop(X %*% c(-1.8, 0, 0.2, 0, -0.1, 0)) - 0.2 * g
  if (fam == "BEINF") { den <- 1 + exp(e0) + exp(e1); u <- runif(n); y[u < exp(e0) / den] <- 0; y[u > 1 - exp(e1) / den] <- 1 }
  if (fam == "BEZI") y[runif(n) < plogis(e0)] <- 0
  if (fam == "BEOI") y[runif(n) < plogis(e1)] <- 1
  list(y = y, X = X, g = g)
}
for (fam in c("BEINF", "BEZI", "BEOI", "BE")) for (rp in 1:2) {
  d <- sim_one(fam, 300, rp); colnames(d$X) <- c("(Intercept)", paste0("c", 1:5))
  nul <- promoter_null(d$y, d$X, fam)
  for (k in 1:3) { g <- if (k == 1) d$g else rbinom(300, 2, 0.3)
    tasks[[length(tasks) + 1]] <- list(y = d$y, C = d$X, g = g, family = fam, nul = nul,
      meta = data.frame(source = "sim", tissue = paste0("sim", rp), promoter = fam, snp = if (k == 1) "causal" else paste0("null", k),
                        n = 300, n0 = sum(d$y == 0), n1 = sum(d$y == 1), n_cov = 5L, n_cov_dropped = 0L, top = k == 1)) }
}
cat("total tasks:", length(tasks), "\n")
fwrite(do.call(rbind, scan_log), file.path(OUT, "engine_scan_timing.tsv"), sep = "\t")

## ---- run --------------------------------------------------------------------------------------
res <- mclapply(tasks, function(t) tryCatch(one_snp(t), error = function(e) { cat("ERR", conditionMessage(e), "\n"); NULL }),
                mc.cores = NC, mc.preschedule = FALSE)
nerr <- sum(vapply(res, is.null, TRUE))
R <- do.call(rbind, res[!vapply(res, is.null, TRUE)])
R$d_ll <- R$ll_engine - R$ll_gamlss; R$d_ll_tight <- R$ll_engine - R$ll_gamlss_tight
fwrite(R, file.path(OUT, "equiv_models.tsv"), sep = "\t")

## LRT comparison per SNP: engine statistic vs 2 * (l_full - l_reduced) from gamlss
key <- c("source", "tissue", "promoter", "snp")
W <- dcast(as.data.table(R), as.formula(paste(paste(key, collapse = "+"), "~ model")), value.var = c("ll_engine", "ll_gamlss", "gamlss_conv"))
L <- rbindlist(lapply(c("null", "dmu", "dnu", "dtau"), function(k) {
  ce <- paste0("ll_engine_", k); cg <- paste0("ll_gamlss_", k); cc <- paste0("gamlss_conv_", k)
  if (!ce %in% names(W)) return(NULL)
  W[is.finite(get(ce)), .(source, tissue, promoter, snp, test = k,
      lrt_engine = 2 * (ll_engine_full - get(ce)), lrt_gamlss = 2 * (ll_gamlss_full - get(cg)),
      both_conv = gamlss_conv_full & get(cc))]
}))
L[, d_lrt := lrt_engine - lrt_gamlss]
fwrite(L, file.path(OUT, "equiv_lrt.tsv"), sep = "\t")

sink(file.path(OUT, "equiv_summary.txt"))
cat("tasks", length(tasks), " errors", nerr, " model rows", nrow(R), "\n\n")
cat("== models where gamlss converged: |d_ll| distribution (engine - gamlss) ==\n")
Rc <- as.data.table(R)[gamlss_conv == TRUE & is.finite(d_ll)]
print(Rc[, .(n = .N, within_1e3 = sum(abs(d_ll) <= 1e-3), engine_lower = sum(d_ll < -1e-3),
             engine_higher = sum(d_ll > 1e-3), max_abs = max(abs(d_ll))), by = .(source, family)])
R2 <- as.data.table(R); R2[, sep_any := (mu_sep %in% TRUE) | (nu_sep %in% TRUE) | (tau_sep %in% TRUE) | (null_cat_conv %in% FALSE) | (null_beta_conv %in% FALSE)]
cat("\n== REGULAR (engine: no separation flag, null fits converged) & gamlss converged ==\n")
print(R2[sep_any == FALSE & gamlss_conv == TRUE & is.finite(d_ll), .(n = .N, within_1e3 = sum(abs(d_ll) <= 1e-3), engine_lower = sum(d_ll < -1e-3),
         engine_higher = sum(d_ll > 1e-3), max_abs = max(abs(d_ll))), by = .(source, family)])
cat("\n== SEPARATED (MLE does not exist): engine must not be below gamlss by > 1e-3 ==\n")
print(R2[sep_any == TRUE & is.finite(d_ll), .(n = .N, engine_lower = sum(d_ll < -1e-3), engine_higher = sum(d_ll > 1e-3),
         min_d = min(d_ll), max_d = max(d_ll)), by = .(source, tissue, family)])
cat("\n== engine lower than gamlss by > 1e-3 (FAILURE cases) ==\n")
print(Rc[d_ll < -1e-3, .(source, tissue, promoter, snp, family, model, d_ll, d_ll_tight, mu_sep, nu_sep, tau_sep)])
cat("\n== engine higher by > 1e-3: does a tight gamlss refit close the gap? ==\n")
print(Rc[d_ll > 1e-3, .(n = .N, closed_by_tight = sum(abs(d_ll_tight) <= 1e-3, na.rm = TRUE),
                        engine_still_higher = sum(d_ll_tight > 1e-3, na.rm = TRUE)), by = .(family, model)])
cat("\n== gamlss not converged / failed ==\n")
print(as.data.table(R)[gamlss_conv != TRUE | is.na(gamlss_conv), .N, by = .(source, family, model)])
cat("\n== LRT, REGULAR SNPs only: |engine - gamlss| where both converged ==\n")
regsnp <- unique(R2[sep_any == FALSE, .(source, tissue, promoter, snp)])
print(L[regsnp, on = .(source, tissue, promoter, snp), nomatch = 0][both_conv == TRUE & is.finite(d_lrt),
        .(n = .N, within_0.01 = sum(abs(d_lrt) <= 0.01), max_abs = max(abs(d_lrt)), median_abs = median(abs(d_lrt))), by = .(source, test)])
cat("\n== LRT: |engine - gamlss| where both converged (all SNPs) ==\n")
print(L[both_conv == TRUE & is.finite(d_lrt), .(n = .N, within_0.01 = sum(abs(d_lrt) <= 0.01),
                                                max_abs = max(abs(d_lrt)), median_abs = median(abs(d_lrt))), by = .(source, test)])
cat("\n== separation flags (engine) among real tasks ==\n")
print(as.data.table(R)[model == "full" & source == "real", .(n = .N, mu_sep = sum(mu_sep %in% TRUE),
      nu_sep = sum(nu_sep %in% TRUE), tau_sep = sum(tau_sep %in% TRUE)), by = .(tissue, family)])
cat("\n== time per SNP (s): engine vs gamlss (sum over fitted models) ==\n")
print(as.data.table(R)[, .(engine = engine_sec[1], gamlss = sum(gamlss_sec)), by = .(source, tissue, family, promoter, snp)][
  , .(n = .N, engine_med = median(engine), gamlss_med = median(gamlss)), by = .(source, tissue, family)])
sink()
cat(readLines(file.path(OUT, "equiv_summary.txt")), sep = "\n")
