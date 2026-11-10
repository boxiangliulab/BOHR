## BOHR Firth-penalised tests for p0/p1 and the null-model simulator.
##
## Ported verbatim from the analysis engine of the BOHR paper,
##   Figure4_puQTL_landscape_eQTL_coloc/rerun_v2/code/bohr2_addon3.R
##   (md5 057882da6606c1f6c53b312f9b72bbbe), lines 65-139 (Firth fits) and
##   212-295 (simulate_null, firth_snp_tests).
## Function bodies are identical to the original except for one change made
## for the package (author decision 2026-11-10: every setting is an argument):
## the iteration limits are arguments with the original values as defaults --
## firth_max() gains `bfgs_maxit = 1000L` (was the literal 1000 in
## optim(control = list(reltol = 1e-12, maxit = 1000))); firth_logit() and
## firth_mnl3() gain `bfgs_maxit = 1000L` and pass it to firth_max();
## firth_snp_tests() gains `maxit = 200L, bfgs_maxit = 1000L` and passes both to
## every firth_logit / firth_mnl3 call (maxround = 6L was already an argument).
## With the defaults the arithmetic is unchanged. Stage 2b (2026-11-10):
## firth_snp_tests() gains `keep_fits = FALSE`; when TRUE it also returns the
## final full penalised fit it has already computed (element `fits$cat`), for
## the coefficients of bohr(); no computation changes. Not ported: the
## variable-precision beta fit, sandwich variances and addon_snp_tests()
## (lines 8-64 and 141-210), which were sensitivity analyses and do not enter
## any published P value. Requires engine.R (nr_max, vcov_of, log1pexp,
## lrt_logp, na0).
##
## All functions in this file are internal (not exported).

## ---------------------------------------------------------------------------------------------
## Firth (Jeffreys-prior) penalised likelihood l* = l + 1/2 log|I|. `fixed` = indices held at their
## start values (the penalised LRT of logistf: the restricted model maximises the FULL model's l*
## with the tested coefficient fixed at 0).
## Firth wrapper: damped Newton (Fisher information as the Hessian proxy), then -- if that stalls,
## which happens on the flat ridges of quasi-separated data because I is near-singular -- a BFGS
## polish with the analytic gradient; converged = Newton converged or BFGS converged.
firth_max <- function(start, fixed, llf, gh, maxit, bfgs_maxit = 1000L) {
  f <- if (length(fixed)) optimise_subset(start, fixed, llf, gh, maxit) else nr_max(start, llf, gh, maxit)
  if (!f$converged) {
    free <- setdiff(seq_along(start), fixed); emb <- function(tf) { th <- f$th; th[free] <- tf; th }
    o <- tryCatch(optim(f$th[free], function(tf) { v <- llf(emb(tf)); if (is.finite(v)) -v else 1e300 },
                        function(tf) -gh(emb(tf))$g[free], method = "BFGS", control = list(reltol = 1e-12, maxit = bfgs_maxit)),
                  error = function(e) NULL)
    ## accept the polish only if the objective re-evaluated at its solution is finite and not lower
    ## (never trust the 1e300 sentinel used for non-finite values inside optim)
    if (!is.null(o)) { th2 <- emb(o$par); l2 <- llf(th2)
      if (is.finite(l2) && l2 >= f$ll - 1e-9) { f$th <- th2; f$ll <- l2; f$converged <- o$convergence == 0 } }
  }
  f
}
optimise_subset <- function(start, fixed, llf, gh, maxit) {
  free <- setdiff(seq_along(start), fixed)
  emb <- function(tf) { th <- start; th[free] <- tf; th }
  f <- nr_max(start[free], function(tf) llf(emb(tf)),
              function(tf) { G <- gh(emb(tf)); list(g = G$g[free], H = G$H[free, free, drop = FALSE]) }, maxit)
  f$th <- emb(f$th); f
}
logdet_pd <- function(I) { d <- determinant(I, logarithm = TRUE); if (d$sign <= 0) -Inf else as.numeric(d$modulus) }

firth_logit <- function(z, X, start = NULL, fixed = integer(0), maxit = 200L, bfgs_maxit = 1000L) {
  llf <- function(b) { e <- drop(X %*% b); q <- plogis(e)
    sum(z * e - log1pexp(e)) + 0.5 * logdet_pd(crossprod(X, (q * (1 - q)) * X)) }
  gh <- function(b) {
    q <- plogis(drop(X %*% b)); w <- q * (1 - q); I <- crossprod(X, w * X); Ii <- vcov_of(I)
    h <- if (is.null(Ii)) rep(0, length(z)) else w * rowSums((X %*% Ii) * X)
    list(g = drop(crossprod(X, z - q + h * (0.5 - q))), H = I)
  }
  if (is.null(start)) start <- c(qlogis(min(max(mean(z), 1e-3), 1 - 1e-3)), rep(0, ncol(X) - 1))
  f <- firth_max(start, fixed, llf, gh, maxit, bfgs_maxit)
  G <- gh(f$th); gmax <- max(abs(G$g[setdiff(seq_along(f$th), fixed)]))
  list(llpen = f$ll, coef = f$th, vcov = vcov_of(G$H), converged = f$converged, iter = f$iter, gmax = gmax)
}

firth_mnl3 <- function(y0, y1, X0, X1, start = NULL, fixed = integer(0), maxit = 200L, bfgs_maxit = 1000L) {
  p0 <- ncol(X0); p1 <- ncol(X1); n0 <- sum(y0); n1 <- sum(y1); nm <- length(y0) - n0 - n1
  probs <- function(th) {
    e0 <- drop(X0 %*% th[1:p0]); e1 <- drop(X1 %*% th[p0 + seq_len(p1)])
    m <- pmax(0, e0, e1); ld <- m + log(exp(-m) + exp(e0 - m) + exp(e1 - m))
    list(e0 = e0, e1 = e1, ld = ld, q0 = exp(e0 - ld), q1 = exp(e1 - ld))
  }
  info <- function(q0, q1) {
    H00 <- crossprod(X0, (q0 * (1 - q0)) * X0); H11 <- crossprod(X1, (q1 * (1 - q1)) * X1); H01 <- -crossprod(X0, (q0 * q1) * X1)
    rbind(cbind(H00, H01), cbind(t(H01), H11))
  }
  llf <- function(th) { u <- probs(th); sum(y0 * u$e0 + y1 * u$e1 - u$ld) + 0.5 * logdet_pd(info(u$q0, u$q1)) }
  gh <- function(th) {
    u <- probs(th); q0 <- u$q0; q1 <- u$q1; I <- info(q0, q1); Ii <- vcov_of(I)
    g <- c(crossprod(X0, y0 - q0), crossprod(X1, y1 - q1))
    if (!is.null(Ii)) {
      i0 <- 1:p0; i1 <- p0 + seq_len(p1)
      A00 <- rowSums((X0 %*% Ii[i0, i0, drop = FALSE]) * X0); A11 <- rowSums((X1 %*% Ii[i1, i1, drop = FALSE]) * X1)
      A01 <- rowSums((X0 %*% Ii[i0, i1, drop = FALSE]) * X1)
      ## dV_i/d eta_0 and dV_i/d eta_1 (V_i = diag(q) - q q'), traced against A_i = Z_i I^-1 Z_i'
      tr0 <- q0 * (1 - q0) * (1 - 2 * q0) * A00 - 2 * q0 * q1 * (1 - 2 * q0) * A01 - q0 * q1 * (1 - 2 * q1) * A11
      tr1 <- -q0 * q1 * (1 - 2 * q0) * A00 - 2 * q0 * q1 * (1 - 2 * q1) * A01 + q1 * (1 - q1) * (1 - 2 * q1) * A11
      g <- g + 0.5 * c(crossprod(X0, tr0), crossprod(X1, tr1))
    }
    list(g = g, H = I)
  }
  if (is.null(start)) start <- c(log(max(n0, 0.5) / nm), rep(0, p0 - 1), log(max(n1, 0.5) / nm), rep(0, p1 - 1))
  f <- firth_max(start, fixed, llf, gh, maxit, bfgs_maxit)
  G <- gh(f$th); gmax <- max(abs(G$g[setdiff(seq_along(f$th), fixed)]))
  list(llpen = f$ll, coef = f$th, vcov = vcov_of(G$H), converged = f$converged, iter = f$iter, gmax = gmax)
}

## ---------------------------------------------------------------------------------------------
## Simulate a phenotype from the fitted null model of a promoter (same samples, covariates and family).
## phi_mult: optional per-sample multiplier of the precision (T5: exp(delta * (g - mean(g)))).
## Categories are redrawn (up to 50 times) until every category of the family is present.
simulate_null <- function(C, nul, phi_mult = 1, y_obs = NULL) {
  ## interior values from the fitted beta part; when the promoter has too few interior samples for a
  ## beta fit (its mu test is NA in the real data too), resample the observed interior values instead --
  ## only the categorical tests are defined there.
  fam <- nul$family; n <- nrow(C); pC <- ncol(C)
  mid_obs <- if (!is.null(y_obs)) y_obs[y_obs > 0 & y_obs < 1] else numeric(0)
  if (is.null(nul$beta) && !length(mid_obs)) return(NULL)
  if (!is.null(nul$beta)) {
    mu <- plogis(drop(C %*% nul$beta$coef[1:pC])); mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
    phi <- exp(nul$beta$coef[pC + 1]) * phi_mult
  }
  need_mid <- if (!is.null(nul$beta)) pC + 3L else 1L
  for (k in 1:50) {
    y <- if (!is.null(nul$beta)) pmin(pmax(rbeta(n, mu * phi, (1 - mu) * phi), 1e-10), 1 - 1e-10) else mid_obs[sample.int(length(mid_obs), n, replace = TRUE)]
    if (fam == "BEINF") {
      e0 <- drop(C %*% nul$cat$coef[1:pC]); e1 <- drop(C %*% nul$cat$coef[pC + 1:pC])
      m <- pmax(0, e0, e1); den <- exp(-m) + exp(e0 - m) + exp(e1 - m)
      p0 <- exp(e0 - m) / den; p1 <- exp(e1 - m) / den; u <- runif(n)
      y[u < p0] <- 0; y[u >= p0 & u < p0 + p1] <- 1
    } else if (fam == "BEZI") y[runif(n) < plogis(drop(C %*% nul$cat$coef))] <- 0
    else if (fam == "BEOI") y[runif(n) < plogis(drop(C %*% nul$cat$coef))] <- 1
    ok <- switch(fam, BEINF = any(y == 0) && any(y == 1) && sum(y > 0 & y < 1) >= need_mid,
                 BEZI = any(y == 0) && sum(y > 0) >= need_mid, BEOI = any(y == 1) && sum(y < 1) >= need_mid, BE = TRUE)
    if (ok) return(y)
  }
  NULL
}

## ---------------------------------------------------------------------------------------------
## Production Firth tests (version 3). l* is not concave, so full and restricted fits are alternated
## until neither improves: full from the null / intercept-only starts; each restricted fit from the
## current full solution (SNP coefficient set to 0) and from the intercept-only start; the full fit
## again from every restricted solution. At exit the full fit is at least as high as a full fit started
## at the restricted optimum, and the restricted fit at least as high as one started from the final full
## solution (the review counter-example of 2026-10-01).
firth_snp_tests <- function(y, C, g, nul, maxround = 6L, maxit = 200L, bfgs_maxit = 1000L, keep_fits = FALSE) {
  fam <- nul$family; pC <- ncol(C); Cg <- cbind(C, g = g); gi <- pC + 1L; eps <- 1e-9
  r <- list(nu_firth_stat = NA_real_, nu_firth_logp = NA_real_, nu_firth_est = NA_real_, nu_firth_se = NA_real_, nu_firth_conv = NA,
            tau_firth_stat = NA_real_, tau_firth_logp = NA_real_, tau_firth_est = NA_real_, tau_firth_se = NA_real_, tau_firth_conv = NA,
            firth_rounds = NA_integer_)
  better <- function(a, b) if (is.null(a) || (is.finite(b$llpen) && b$llpen > a$llpen + eps)) b else a
  if (fam %in% c("BEZI", "BEOI")) {
    z <- as.numeric(if (fam == "BEZI") y == 0 else y == 1)
    s0 <- c(qlogis(min(max(mean(z), 1e-3), 1 - 1e-3)), rep(0, gi - 1))
    ff <- better(if (isTRUE(nul$cat$converged)) firth_logit(z, Cg, start = c(na0(nul$cat$coef), 0), maxit = maxit, bfgs_maxit = bfgs_maxit) else NULL, firth_logit(z, Cg, start = s0, maxit = maxit, bfgs_maxit = bfgs_maxit))
    fr <- better(firth_logit(z, Cg, start = replace(ff$coef, gi, 0), fixed = gi, maxit = maxit, bfgs_maxit = bfgs_maxit), firth_logit(z, Cg, start = s0, fixed = gi, maxit = maxit, bfgs_maxit = bfgs_maxit))
    for (k in seq_len(maxround)) {
      ch <- FALSE
      f2 <- firth_logit(z, Cg, start = fr$coef, maxit = maxit, bfgs_maxit = bfgs_maxit); if (f2$llpen > ff$llpen + eps) { ff <- f2; ch <- TRUE }
      r2 <- firth_logit(z, Cg, start = replace(ff$coef, gi, 0), fixed = gi, maxit = maxit, bfgs_maxit = bfgs_maxit); if (r2$llpen > fr$llpen + eps) { fr <- r2; ch <- TRUE }
      if (!ch) break
    }
    st <- 2 * (ff$llpen - fr$llpen); se <- if (!is.null(ff$vcov) && ff$vcov[gi, gi] > 0) sqrt(ff$vcov[gi, gi]) else NA_real_
    pre <- if (fam == "BEZI") "nu" else "tau"
    r[[paste0(pre, "_firth_stat")]] <- st; r[[paste0(pre, "_firth_logp")]] <- lrt_logp(st, 1)
    r[[paste0(pre, "_firth_est")]] <- ff$coef[gi]; r[[paste0(pre, "_firth_se")]] <- se
    r[[paste0(pre, "_firth_conv")]] <- isTRUE(ff$converged) && isTRUE(fr$converged); r$firth_rounds <- k
    if (keep_fits) r$fits <- list(cat = ff)
  } else if (fam == "BEINF") {
    y0 <- as.numeric(y == 0); y1 <- as.numeric(y == 1); inu <- gi; itau <- 2L * gi
    ff <- better(if (isTRUE(nul$cat$converged)) firth_mnl3(y0, y1, Cg, Cg, start = c(na0(nul$cat$coef[1:pC]), 0, na0(nul$cat$coef[pC + 1:pC]), 0), maxit = maxit, bfgs_maxit = bfgs_maxit) else NULL,
                 firth_mnl3(y0, y1, Cg, Cg, start = NULL, maxit = maxit, bfgs_maxit = bfgs_maxit))
    fn <- better(firth_mnl3(y0, y1, Cg, Cg, start = replace(ff$coef, inu, 0), fixed = inu, maxit = maxit, bfgs_maxit = bfgs_maxit), firth_mnl3(y0, y1, Cg, Cg, start = NULL, fixed = inu, maxit = maxit, bfgs_maxit = bfgs_maxit))
    ft <- better(firth_mnl3(y0, y1, Cg, Cg, start = replace(ff$coef, itau, 0), fixed = itau, maxit = maxit, bfgs_maxit = bfgs_maxit), firth_mnl3(y0, y1, Cg, Cg, start = NULL, fixed = itau, maxit = maxit, bfgs_maxit = bfgs_maxit))
    for (k in seq_len(maxround)) {
      ch <- FALSE
      for (rr in list(fn, ft)) { f2 <- firth_mnl3(y0, y1, Cg, Cg, start = rr$coef, maxit = maxit, bfgs_maxit = bfgs_maxit); if (f2$llpen > ff$llpen + eps) { ff <- f2; ch <- TRUE } }
      n2 <- firth_mnl3(y0, y1, Cg, Cg, start = replace(ff$coef, inu, 0), fixed = inu, maxit = maxit, bfgs_maxit = bfgs_maxit); if (n2$llpen > fn$llpen + eps) { fn <- n2; ch <- TRUE }
      t2 <- firth_mnl3(y0, y1, Cg, Cg, start = replace(ff$coef, itau, 0), fixed = itau, maxit = maxit, bfgs_maxit = bfgs_maxit); if (t2$llpen > ft$llpen + eps) { ft <- t2; ch <- TRUE }
      if (!ch) break
    }
    V <- ff$vcov
    r$nu_firth_stat <- 2 * (ff$llpen - fn$llpen); r$tau_firth_stat <- 2 * (ff$llpen - ft$llpen)
    r$nu_firth_logp <- lrt_logp(r$nu_firth_stat, 1); r$tau_firth_logp <- lrt_logp(r$tau_firth_stat, 1)
    r$nu_firth_est <- ff$coef[inu]; r$tau_firth_est <- ff$coef[itau]
    if (!is.null(V)) { if (V[inu, inu] > 0) r$nu_firth_se <- sqrt(V[inu, inu]); if (V[itau, itau] > 0) r$tau_firth_se <- sqrt(V[itau, itau]) }
    r$nu_firth_conv <- isTRUE(ff$converged) && isTRUE(fn$converged); r$tau_firth_conv <- isTRUE(ff$converged) && isTRUE(ft$converged)
    r$firth_rounds <- k
    if (keep_fits) r$fits <- list(cat = ff)
  }
  r
}
