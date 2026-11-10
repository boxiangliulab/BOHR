## BOHR model engine: maximiser, part-wise fits and per-predictor tests.
##
## Ported verbatim from the analysis engine of the BOHR paper,
##   Figure4_puQTL_landscape_eQTL_coloc/rerun_v2/code/bohr2_engine.R
##   (md5 23d725f8d7e8bfac7e6315c0c271e6d0), lines 16-224.
## Function bodies are identical to the original except for one change made
## for the package (author decision 2026-11-10: every setting is an argument):
## promoter_null() and snp_tests() gain an argument `maxit = 200L` that is passed
## to every beta_fit / logit_fit / mnl3_fit call in their bodies (the original
## calls used those functions' default, maxit = 200L). With the default the
## arithmetic is unchanged. Stage 2b (2026-11-10): snp_tests() gains
## `keep_fits = FALSE`; when TRUE it also returns the full fits it has already
## computed (element `fits`: `beta`, the beta fit with the predictor; `cat`, the
## categorical fit with the predictor, after the nesting guard), for the
## coefficients of bohr(). No computation changes; with the default the
## returned list is the original one. Original header (lines 1-14):
##
##   BEINF : l = l_mnl3( 0 / mid / 1 ; nu, tau )  +  l_beta( 0<y<1 ; mu, sigma )
##   BEZI  : l = l_logit( y==0 vs >0 ; nu )        +  l_beta( 0<y<1 ; mu, sigma )
##   BEOI  : l = l_logit( y==1 vs <1 ; tau )       +  l_beta( 0<y<1 ; mu, sigma )
##   BE    : l =                                        l_beta( all  ; mu, sigma )
##
##   gamlss parameterisation this reproduces: p0 = nu/(1+nu+tau), p1 = tau/(1+nu+tau) with log links
##   (so log nu = log p0/p_mid, log tau = log p1/p_mid: a baseline-category multinomial logit with the
##   interior as reference); BEZI/BEOI nu is a logit-linked probability; the interior is a beta with a
##   logit-linked mean and, here, a constant precision. The two parts share no parameters, so the joint
##   MLE is the pair of part-wise MLEs and every LRT is the sum of the part-wise LRTs.
##
## All functions in this file are internal (not exported).

## ---------------------------------------------------------------------------------------------
## Shared maximiser: monotone Newton / Fisher scoring with step halving and Levenberg-Marquardt
## damping. H is the (expected or observed) information, positive definite at a regular point.
## The log-likelihood never decreases. When the Newton direction fails (near-singular H, e.g. under
## separation) the damping moves the step towards the gradient. Converged = step < 1e-7 with a
## negligible gain; five consecutive negligible gains with a large step = stalled (separation /
## unbounded likelihood) -> stop, converged = FALSE.
nr_max <- function(th, llf, gh, maxit = 200L, ll_sup = Inf) {
  ll <- llf(th); conv <- FALSE; lam <- 0; stall <- 0L; it <- 0L
  for (it in seq_len(maxit)) {
    if (ll > ll_sup - 1e-8) break                   # at the supremum: complete separation (never at a regular optimum)
    G <- gh(th); dsc <- max(abs(diag(G$H)), 1e-12); ok <- FALSE
    for (tr in 1:10) {
      step <- tryCatch(solve(G$H + diag(lam * dsc, length(G$g)), G$g), error = function(e) NULL)
      if (!is.null(step) && all(is.finite(step))) {
        s <- 1
        while (s >= 1e-6) {
          thn <- th + s * step; lln <- llf(thn)
          if (is.finite(lln) && lln >= ll - 1e-10 * (1 + abs(ll))) { ok <- TRUE; break }
          s <- s / 2
        }
        if (ok) break
      }
      lam <- if (lam == 0) 1e-8 else lam * 100
    }
    if (!ok) break
    dll <- lln - ll; th <- thn; ll <- lln
    lam <- if (lam > 1e-8) lam / 100 else 0
    if (abs(dll) < 1e-11 * (1 + abs(ll))) {
      if (max(abs(s * step)) < 1e-7) { conv <- TRUE; break }
      stall <- stall + 1L; if (stall >= 5L) break
    } else stall <- 0L
  }
  list(th = th, ll = ll, converged = conv, iter = it)
}
vcov_of <- function(H) tryCatch(solve(H), error = function(e) NULL)

## ---------------------------------------------------------------------------------------------
## beta regression, logit mean, constant precision phi (optimised on log phi); Fisher scoring with the
## expected information of Ferrari & Cribari-Neto (2004).
beta_fit <- function(y, X, start = NULL, maxit = 200L) {
  p <- ncol(X); ly <- log(y); l1y <- log1p(-y); ystar <- ly - l1y
  llf <- function(th) {
    mu <- plogis(drop(X %*% th[1:p])); mu <- pmin(pmax(mu, 1e-12), 1 - 1e-12); phi <- exp(th[p + 1])
    sum(lgamma(phi) - lgamma(mu * phi) - lgamma((1 - mu) * phi) + (mu * phi - 1) * ly + ((1 - mu) * phi - 1) * l1y)
  }
  gh <- function(th) {
    mu <- plogis(drop(X %*% th[1:p])); mu <- pmin(pmax(mu, 1e-12), 1 - 1e-12); phi <- exp(th[p + 1])
    a <- mu * phi; b <- (1 - mu) * phi; t <- mu * (1 - mu)
    mustar <- digamma(a) - digamma(b); ta <- trigamma(a); tb <- trigamma(b)
    U <- c(phi * crossprod(X, t * (ystar - mustar)),
           phi * sum(mu * (ystar - mustar) + l1y - digamma(b) + digamma(phi)))
    w <- phi * (ta + tb) * t^2
    cc <- phi * (ta * mu - tb * (1 - mu))
    Kbb <- phi * crossprod(X, w * X)
    Kbp <- phi * crossprod(X, t * cc)                                  # d/d(log phi) = phi * d/d(phi)
    Kpp <- phi^2 * sum(ta * mu^2 + tb * (1 - mu)^2 - trigamma(phi))
    list(g = U, H = rbind(cbind(Kbb, Kbp), c(Kbp, Kpp)))
  }
  if (is.null(start)) {
    b0 <- tryCatch(qr.coef(qr(X), ystar), error = function(e) rep(0, p)); b0[is.na(b0)] <- 0
    mu0 <- plogis(drop(X %*% b0))
    phi0 <- max(mean(mu0 * (1 - mu0)) / max(var(y), 1e-8) - 1, 0.5)
    start <- c(b0, log(phi0))
  }
  f <- nr_max(start, llf, gh, maxit)
  G <- gh(f$th)
  list(ll = f$ll, coef = f$th, vcov = vcov_of(G$H), converged = f$converged, iter = f$iter, maxgrad = max(abs(G$g)))
}

## ---------------------------------------------------------------------------------------------
## binary logistic regression (Newton = Fisher scoring for the canonical link). Replaces glm.fit, whose
## IRLS halves the step only on non-finite deviance and can diverge under quasi-separation
## (stage 1: Kidney PromoterID_5712 ended 12 x log(eps) = 432 log-lik units below gamlss).
log1pexp <- function(e) ifelse(e > 0, e + log1p(exp(-e)), log1p(exp(e)))
logit_fit <- function(z, X, start = NULL, maxit = 200L) {
  llf <- function(b) { e <- drop(X %*% b); sum(z * e - log1pexp(e)) }
  gh <- function(b) { q <- plogis(drop(X %*% b)); list(g = drop(crossprod(X, z - q)), H = crossprod(X, (q * (1 - q)) * X)) }
  if (is.null(start)) start <- c(qlogis(min(max(mean(z), 1e-6), 1 - 1e-6)), rep(0, ncol(X) - 1))
  f <- nr_max(start, llf, gh, maxit, ll_sup = 0)
  G <- gh(f$th)
  list(ll = f$ll, coef = f$th, vcov = vcov_of(G$H), converged = f$converged, iter = f$iter, maxgrad = max(abs(G$g)))
}

## ---------------------------------------------------------------------------------------------
## 3-category baseline logit: categories zero (y0 = 1), one (y1 = 1), reference = interior.
## eq0: log(p0/pmid) = X0 %*% t0 ; eq1: log(p1/pmid) = X1 %*% t1 (designs may differ -> per-parameter tests).
## First column of X0 and X1 must be the intercept.
mnl3_fit <- function(y0, y1, X0, X1, start = NULL, maxit = 200L) {
  p0 <- ncol(X0); p1 <- ncol(X1); n0 <- sum(y0); n1 <- sum(y1); nm <- length(y0) - n0 - n1
  lden <- function(e0, e1) { m <- pmax(0, e0, e1); m + log(exp(-m) + exp(e0 - m) + exp(e1 - m)) }
  llf <- function(th) {
    e0 <- drop(X0 %*% th[1:p0]); e1 <- drop(X1 %*% th[p0 + seq_len(p1)])
    sum(y0 * e0 + y1 * e1 - lden(e0, e1))
  }
  gh <- function(th) {
    e0 <- drop(X0 %*% th[1:p0]); e1 <- drop(X1 %*% th[p0 + seq_len(p1)]); ld <- lden(e0, e1)
    q0 <- exp(e0 - ld); q1 <- exp(e1 - ld)
    g <- c(crossprod(X0, y0 - q0), crossprod(X1, y1 - q1))
    H00 <- crossprod(X0, (q0 * (1 - q0)) * X0); H11 <- crossprod(X1, (q1 * (1 - q1)) * X1)
    H01 <- -crossprod(X0, (q0 * q1) * X1)
    list(g = g, H = rbind(cbind(H00, H01), cbind(t(H01), H11)))
  }
  if (is.null(start)) start <- c(log(n0 / nm), rep(0, p0 - 1), log(n1 / nm), rep(0, p1 - 1))
  f <- nr_max(start, llf, gh, maxit, ll_sup = 0)
  G <- gh(f$th)
  list(ll = f$ll, coef = f$th, vcov = vcov_of(G$H), converged = f$converged, iter = f$iter, maxgrad = max(abs(G$g)))
}

## ---------------------------------------------------------------------------------------------
## helpers
lrt_logp <- function(stat, df) ifelse(is.finite(stat), pchisq(pmax(stat, 0), df, lower.tail = FALSE, log.p = TRUE), NA_real_)
wald_logp <- function(est, se) ifelse(is.finite(est) & is.finite(se) & se > 0,
                                      log(2) + pnorm(-abs(est / se), log.p = TRUE), NA_real_)
na0 <- function(x) { x[!is.finite(x)] <- 0; x }
sep_flag <- function(est, se, conv) !isTRUE(conv) || !is.finite(est) || abs(est) > 15 || !is.finite(se) || se > 50

## ---------------------------------------------------------------------------------------------
## Null fits for one promoter. y already has endpoint deletion applied; C has an intercept column.
promoter_null <- function(y, C, family, maxit = 200L) {
  mid <- y > 0 & y < 1
  out <- list(family = family, n = length(y), n_mid = sum(mid), n0 = sum(y == 0), n1 = sum(y == 1))
  out$beta <- if (sum(mid) > ncol(C) + 2) beta_fit(y[mid], C[mid, , drop = FALSE], maxit = maxit) else NULL
  out$cat <- switch(family,
    BEINF = mnl3_fit(as.numeric(y == 0), as.numeric(y == 1), C, C, maxit = maxit),
    BEZI  = logit_fit(as.numeric(y == 0), C, maxit = maxit),
    BEOI  = logit_fit(as.numeric(y == 1), C, maxit = maxit),
    BE    = NULL)
  out
}

## ---------------------------------------------------------------------------------------------
## All tests for one SNP. g = dosage (0/1/2) aligned with y and C. nul = promoter_null(...).
snp_tests <- function(y, C, g, nul, maxit = 200L, keep_fits = FALSE) {
  fam <- nul$family; mid <- y > 0 & y < 1; pC <- ncol(C); fits <- list()
  r <- list(mu_stat = NA_real_, mu_est = NA_real_, mu_se = NA_real_, mu_conv = NA, mu_sep = NA,
            nu_stat = NA_real_, nu_est = NA_real_, nu_se = NA_real_, nu_conv = NA, nu_sep = NA,
            tau_stat = NA_real_, tau_est = NA_real_, tau_se = NA_real_, tau_conv = NA, tau_sep = NA,
            cat_joint_stat = NA_real_, cat_df = 0L, ll_beta_full = NA_real_, ll_cat_full = NA_real_,
            mu_pe = NA_real_, p0_pe = NA_real_, p1_pe = NA_real_, mu_iter = NA_integer_, cat_iter = NA_integer_)
  ## --- beta part (mu) ---
  if (!is.null(nul$beta) && length(unique(g[mid])) > 1) {
    Xm <- cbind(C[mid, , drop = FALSE], g = g[mid])
    st <- c(nul$beta$coef[1:pC], 0, nul$beta$coef[pC + 1])
    f <- beta_fit(y[mid], Xm, start = st, maxit = maxit)
    if (keep_fits) fits$beta <- f
    r$ll_beta_full <- f$ll; r$mu_stat <- 2 * (f$ll - nul$beta$ll); r$mu_conv <- f$converged
    r$mu_est <- f$coef[pC + 1]; r$mu_se <- if (!is.null(f$vcov)) sqrt(f$vcov[pC + 1, pC + 1]) else NA_real_
    r$mu_sep <- sep_flag(r$mu_est, r$mu_se, f$converged)
    r$mu_iter <- f$iter
    b <- f$coef[1:pC]; gm <- f$coef[pC + 1]; eta <- drop(C %*% b)
    r$mu_pe <- mean(plogis(eta + 2 * gm) - plogis(eta))            # covariate-adjusted (all samples), dosage 0 -> 2
  }
  ## --- categorical part (nu / tau) ---
  ## If the covariates alone already separate the categories completely (null categorical log-lik at
  ## its supremum 0), every categorical model has sup l = 0: nu/tau tests are not assessable -> NA.
  Cg <- cbind(C, g = g)
  r$cat_complete_sep <- !is.null(nul$cat) && nul$cat$ll > -1e-6
  if (r$cat_complete_sep) {
    if (fam %in% c("BEINF", "BEZI")) { r$nu_conv <- FALSE; r$nu_sep <- TRUE }
    if (fam %in% c("BEINF", "BEOI")) { r$tau_conv <- FALSE; r$tau_sep <- TRUE }
  } else if (fam == "BEZI" || fam == "BEOI") {
    z <- as.numeric(if (fam == "BEZI") y == 0 else y == 1)
    f <- logit_fit(z, Cg, start = c(na0(nul$cat$coef), 0), maxit = maxit); r$cat_iter <- f$iter
    if (keep_fits) fits$cat <- f
    st <- 2 * (f$ll - nul$cat$ll); est <- f$coef[pC + 1]
    se <- if (!is.null(f$vcov)) sqrt(f$vcov[pC + 1, pC + 1]) else NA_real_
    eta <- drop(C %*% f$coef[1:pC]); pe <- mean(plogis(eta + 2 * est) - plogis(eta))
    if (fam == "BEZI") { r$nu_stat <- st; r$nu_est <- est; r$nu_se <- se; r$nu_conv <- f$converged
      r$nu_sep <- sep_flag(est, se, f$converged); r$p0_pe <- pe; r$p1_pe <- 0 }
    else { r$tau_stat <- st; r$tau_est <- est; r$tau_se <- se; r$tau_conv <- f$converged
      r$tau_sep <- sep_flag(est, se, f$converged); r$p1_pe <- pe; r$p0_pe <- 0 }
    r$cat_joint_stat <- st; r$cat_df <- 1L; r$ll_cat_full <- f$ll
  } else if (fam == "BEINF") {
    y0 <- as.numeric(y == 0); y1 <- as.numeric(y == 1)
    t0 <- na0(nul$cat$coef[1:pC]); t1 <- na0(nul$cat$coef[pC + 1:pC])
    full <- mnl3_fit(y0, y1, Cg, Cg, start = c(t0, 0, t1, 0), maxit = maxit)
    dnu  <- mnl3_fit(y0, y1, C,  Cg, start = c(t0, t1, 0), maxit = maxit)          # SNP dropped from nu (eq0)
    dtau <- mnl3_fit(y0, y1, Cg, C,  start = c(t0, 0, t1), maxit = maxit)          # SNP dropped from tau (eq1)
    ## nesting guard: a reduced model can never beat the full one at the optimum -> refit full from it
    best <- max(dnu$ll, dtau$ll)
    if (is.finite(best) && best > full$ll + 1e-8) {
      stt <- if (dnu$ll >= dtau$ll) c(dnu$coef[1:pC], 0, dnu$coef[pC + 1:(pC + 1)]) else c(dtau$coef[1:(pC + 1)], dtau$coef[pC + 1 + 1:pC], 0)
      f2 <- mnl3_fit(y0, y1, Cg, Cg, start = stt, maxit = maxit); if (f2$ll > full$ll) full <- f2
    }
    r$ll_cat_full <- full$ll; r$cat_iter <- full$iter
    if (keep_fits) fits$cat <- full
    r$cat_joint_stat <- 2 * (full$ll - nul$cat$ll); r$cat_df <- 2L
    V <- full$vcov
    r$nu_stat <- 2 * (full$ll - dnu$ll);  r$nu_est <- full$coef[pC + 1]
    r$tau_stat <- 2 * (full$ll - dtau$ll); r$tau_est <- full$coef[2 * (pC + 1)]
    r$nu_se  <- if (!is.null(V)) sqrt(V[pC + 1, pC + 1]) else NA_real_
    r$tau_se <- if (!is.null(V)) sqrt(V[2 * (pC + 1), 2 * (pC + 1)]) else NA_real_
    r$nu_conv <- full$converged && dnu$converged; r$tau_conv <- full$converged && dtau$converged
    r$nu_sep <- sep_flag(r$nu_est, r$nu_se, r$nu_conv); r$tau_sep <- sep_flag(r$tau_est, r$tau_se, r$tau_conv)
    th <- full$coef; e0 <- drop(C %*% th[1:pC]); e1 <- drop(C %*% th[pC + 1 + 1:pC])
    a0 <- th[pC + 1]; a1 <- th[2 * (pC + 1)]
    pr <- function(d) { l0 <- e0 + a0 * d; l1 <- e1 + a1 * d; m <- pmax(0, l0, l1)          # max-shifted softmax (review 2026-10-01)
      den <- exp(-m) + exp(l0 - m) + exp(l1 - m); cbind(exp(l0 - m) / den, exp(l1 - m) / den) }
    P0 <- pr(0); P2 <- pr(2)
    r$p0_pe <- mean(P2[, 1] - P0[, 1]); r$p1_pe <- mean(P2[, 2] - P0[, 2])
  }
  ## --- log p values; sigma-free joint = beta part + categorical part ---
  mu_df <- if (is.finite(r$mu_stat)) 1L else 0L
  joint <- sum(c(if (is.finite(r$mu_stat)) r$mu_stat, if (is.finite(r$cat_joint_stat)) r$cat_joint_stat))
  jdf <- if (isTRUE(r$cat_complete_sep)) 0L else mu_df + r$cat_df
  res <- c(r, list(joint_stat = if (jdf > 0) joint else NA_real_, joint_df = jdf,
            mu_lrt_logp = lrt_logp(r$mu_stat, 1), nu_lrt_logp = lrt_logp(r$nu_stat, 1), tau_lrt_logp = lrt_logp(r$tau_stat, 1),
            joint_logp = if (jdf > 0) lrt_logp(joint, jdf) else NA_real_,
            mu_wald_logp = wald_logp(r$mu_est, r$mu_se), nu_wald_logp = wald_logp(r$nu_est, r$nu_se),
            tau_wald_logp = wald_logp(r$tau_est, r$tau_se)))
  if (keep_fits) res$fits <- fits
  res
}
