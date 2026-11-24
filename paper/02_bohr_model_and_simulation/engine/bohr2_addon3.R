## BOHR v2 add-on tests, version 3 (Firth alternating search, finite BFGS polish).
## Version 2 note: (multi-start Firth, see PILOT_PROTOCOL.md amendment 5):
## identical to bohr2_addon.R except the Firth part of addon_snp_tests(). Original header:
## variable-precision beta, sandwich SEs,
## Firth-penalised categorical fits with penalised LRTs, and a simulator from the fitted null model.
## Requires bohr2_engine.R (nr_max, beta_fit, logit_fit, mnl3_fit, lrt_logp, wald_logp, log1pexp, na0).

## ---------------------------------------------------------------------------------------------
## Beta regression with variable precision: logit(mu) = X b, log(phi) = Z g. Fisher scoring with the
## expected information (the constant-precision engine is the special case Z = 1).
beta_fit_vp <- function(y, X, Z, start = NULL, maxit = 200L) {
  p <- ncol(X); q <- ncol(Z); ly <- log(y); l1y <- log1p(-y); ystar <- ly - l1y
  unpack <- function(th) {
    mu <- plogis(drop(X %*% th[1:p])); mu <- pmin(pmax(mu, 1e-12), 1 - 1e-12)
    phi <- exp(pmin(drop(Z %*% th[p + 1:q]), 700)); list(mu = mu, phi = phi)
  }
  llf <- function(th) { u <- unpack(th); mu <- u$mu; phi <- u$phi
    sum(lgamma(phi) - lgamma(mu * phi) - lgamma((1 - mu) * phi) + (mu * phi - 1) * ly + ((1 - mu) * phi - 1) * l1y) }
  gh <- function(th) {
    u <- unpack(th); mu <- u$mu; phi <- u$phi
    a <- mu * phi; b <- (1 - mu) * phi; t <- mu * (1 - mu)
    mustar <- digamma(a) - digamma(b); ta <- trigamma(a); tb <- trigamma(b)
    sb <- phi * t * (ystar - mustar)                                          # d l_i / d eta_mu
    sg <- phi * (mu * (ystar - mustar) + l1y - digamma(b) + digamma(phi))    # d l_i / d eta_phi
    wbb <- phi^2 * t^2 * (ta + tb)
    wbg <- t * phi * phi * (ta * mu - tb * (1 - mu))
    wgg <- phi^2 * (ta * mu^2 + tb * (1 - mu)^2 - trigamma(phi))
    Kbb <- crossprod(X, wbb * X); Kbg <- crossprod(X, wbg * Z); Kgg <- crossprod(Z, wgg * Z)
    list(g = c(crossprod(X, sb), crossprod(Z, sg)), H = rbind(cbind(Kbb, Kbg), cbind(t(Kbg), Kgg)),
         S = cbind(X * sb, Z * sg))
  }
  if (is.null(start)) {
    b0 <- tryCatch(qr.coef(qr(X), ystar), error = function(e) rep(0, p)); b0[is.na(b0)] <- 0
    mu0 <- plogis(drop(X %*% b0))
    phi0 <- max(mean(mu0 * (1 - mu0)) / max(var(y), 1e-8) - 1, 0.5)
    start <- c(b0, log(phi0), rep(0, q - 1))
  }
  f <- nr_max(start, llf, function(th) { G <- gh(th); list(g = G$g, H = G$H) }, maxit)
  G <- gh(f$th)
  list(ll = f$ll, coef = f$th, vcov = vcov_of(G$H), converged = f$converged, iter = f$iter, S = G$S)
}

## ---------------------------------------------------------------------------------------------
## Sandwich (HC0) covariance: bread = inverse expected information, meat = sum of outer products of
## the per-observation scores at the fitted values.
sandwich_V <- function(H, S) { B <- vcov_of(H); if (is.null(B)) return(NULL); B %*% crossprod(S) %*% B }
beta_scores <- function(y, X, th) {                       # constant precision, th = c(b, log phi)
  p <- ncol(X); mu <- plogis(drop(X %*% th[1:p])); mu <- pmin(pmax(mu, 1e-12), 1 - 1e-12); phi <- exp(th[p + 1])
  ystar <- log(y) - log1p(-y); a <- mu * phi; b <- (1 - mu) * phi; t <- mu * (1 - mu)
  mustar <- digamma(a) - digamma(b); ta <- trigamma(a); tb <- trigamma(b)
  sb <- phi * t * (ystar - mustar); sp <- phi * (mu * (ystar - mustar) + log1p(-y) - digamma(b) + digamma(phi))
  H <- rbind(cbind(crossprod(X, (phi^2 * t^2 * (ta + tb)) * X), crossprod(X, t * phi * phi * (ta * mu - tb * (1 - mu)))),
             c(crossprod(X, t * phi * phi * (ta * mu - tb * (1 - mu))), phi^2 * sum(ta * mu^2 + tb * (1 - mu)^2 - trigamma(phi))))
  list(H = H, S = cbind(X * sb, sp))
}
logit_scores <- function(z, X, b) { q <- plogis(drop(X %*% b)); list(H = crossprod(X, (q * (1 - q)) * X), S = X * (z - q)) }
mnl3_scores <- function(y0, y1, X0, X1, th) {
  p0 <- ncol(X0); p1 <- ncol(X1)
  e0 <- drop(X0 %*% th[1:p0]); e1 <- drop(X1 %*% th[p0 + seq_len(p1)])
  m <- pmax(0, e0, e1); ld <- m + log(exp(-m) + exp(e0 - m) + exp(e1 - m)); q0 <- exp(e0 - ld); q1 <- exp(e1 - ld)
  H00 <- crossprod(X0, (q0 * (1 - q0)) * X0); H11 <- crossprod(X1, (q1 * (1 - q1)) * X1); H01 <- -crossprod(X0, (q0 * q1) * X1)
  list(H = rbind(cbind(H00, H01), cbind(t(H01), H11)), S = cbind(X0 * (y0 - q0), X1 * (y1 - q1)))
}

## ---------------------------------------------------------------------------------------------
## Firth (Jeffreys-prior) penalised likelihood l* = l + 1/2 log|I|. `fixed` = indices held at their
## start values (the penalised LRT of logistf: the restricted model maximises the FULL model's l*
## with the tested coefficient fixed at 0).
## Firth wrapper: damped Newton (Fisher information as the Hessian proxy), then -- if that stalls,
## which happens on the flat ridges of quasi-separated data because I is near-singular -- a BFGS
## polish with the analytic gradient; converged = Newton converged or BFGS converged.
firth_max <- function(start, fixed, llf, gh, maxit) {
  f <- if (length(fixed)) optimise_subset(start, fixed, llf, gh, maxit) else nr_max(start, llf, gh, maxit)
  if (!f$converged) {
    free <- setdiff(seq_along(start), fixed); emb <- function(tf) { th <- f$th; th[free] <- tf; th }
    o <- tryCatch(optim(f$th[free], function(tf) { v <- llf(emb(tf)); if (is.finite(v)) -v else 1e300 },
                        function(tf) -gh(emb(tf))$g[free], method = "BFGS", control = list(reltol = 1e-12, maxit = 1000)),
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

firth_logit <- function(z, X, start = NULL, fixed = integer(0), maxit = 200L) {
  llf <- function(b) { e <- drop(X %*% b); q <- plogis(e)
    sum(z * e - log1pexp(e)) + 0.5 * logdet_pd(crossprod(X, (q * (1 - q)) * X)) }
  gh <- function(b) {
    q <- plogis(drop(X %*% b)); w <- q * (1 - q); I <- crossprod(X, w * X); Ii <- vcov_of(I)
    h <- if (is.null(Ii)) rep(0, length(z)) else w * rowSums((X %*% Ii) * X)
    list(g = drop(crossprod(X, z - q + h * (0.5 - q))), H = I)
  }
  if (is.null(start)) start <- c(qlogis(min(max(mean(z), 1e-3), 1 - 1e-3)), rep(0, ncol(X) - 1))
  f <- firth_max(start, fixed, llf, gh, maxit)
  G <- gh(f$th); gmax <- max(abs(G$g[setdiff(seq_along(f$th), fixed)]))
  list(llpen = f$ll, coef = f$th, vcov = vcov_of(G$H), converged = f$converged, iter = f$iter, gmax = gmax)
}

firth_mnl3 <- function(y0, y1, X0, X1, start = NULL, fixed = integer(0), maxit = 200L) {
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
  f <- firth_max(start, fixed, llf, gh, maxit)
  G <- gh(f$th); gmax <- max(abs(G$g[setdiff(seq_along(f$th), fixed)]))
  list(llpen = f$ll, coef = f$th, vcov = vcov_of(G$H), converged = f$converged, iter = f$iter, gmax = gmax)
}

## ---------------------------------------------------------------------------------------------
## Per-SNP add-on statistics (T2, T3, optional T8): standard fits as in snp_tests, then sandwich
## Wald for mu/nu/tau, Firth penalised LRT + Wald for nu/tau, and (vp = TRUE) the mu LRT with the
## precision modelled on the same covariates on both sides.
addon_snp_tests <- function(y, C, g, nul, vp_null = NULL) {
  fam <- nul$family; mid <- y > 0 & y < 1; pC <- ncol(C); Cg <- cbind(C, g = g); gi <- pC + 1L
  r <- list(mu_sw_logp = NA_real_, nu_sw_logp = NA_real_, tau_sw_logp = NA_real_,
            nu_firth_stat = NA_real_, nu_firth_logp = NA_real_, nu_firth_est = NA_real_, nu_firth_se = NA_real_, nu_firth_conv = NA,
            tau_firth_stat = NA_real_, tau_firth_logp = NA_real_, tau_firth_est = NA_real_, tau_firth_se = NA_real_, tau_firth_conv = NA,
            mu_vp_stat = NA_real_, mu_vp_logp = NA_real_, mu_vp_conv = NA)
  ## beta part: sandwich for mu (+ variable-precision LRT when requested)
  if (!is.null(nul$beta) && length(unique(g[mid])) > 1) {
    Xm <- cbind(C[mid, , drop = FALSE], g = g[mid])
    f <- beta_fit(y[mid], Xm, start = c(nul$beta$coef[1:pC], 0, nul$beta$coef[pC + 1]))
    sc <- beta_scores(y[mid], Xm, f$coef); V <- sandwich_V(sc$H, sc$S)
    if (!is.null(V) && V[gi, gi] > 0) r$mu_sw_logp <- wald_logp(f$coef[gi], sqrt(V[gi, gi]))
    if (!is.null(vp_null)) {
      Zm <- C[mid, , drop = FALSE]
      fv <- tryCatch(beta_fit_vp(y[mid], Xm, Zm, start = c(vp_null$coef[1:pC], 0, vp_null$coef[pC + seq_len(ncol(Zm))])), error = function(e) NULL)
      if (!is.null(fv)) { r$mu_vp_stat <- 2 * (fv$ll - vp_null$ll); r$mu_vp_logp <- lrt_logp(r$mu_vp_stat, 1)
        r$mu_vp_conv <- fv$converged && vp_null$converged }
    }
  }
  ## categorical part (a completely separated categorical null has no MLE: sandwich NA, Firth still defined)
  cs <- !is.null(nul$cat) && nul$cat$ll > -1e-6
  if (fam %in% c("BEZI", "BEOI")) {
    z <- as.numeric(if (fam == "BEZI") y == 0 else y == 1)
    sw <- NA_real_
    if (!cs) { f <- logit_fit(z, Cg, start = c(na0(nul$cat$coef), 0))
      sc <- logit_scores(z, Cg, f$coef); V <- sandwich_V(sc$H, sc$S)
      if (!is.null(V) && V[gi, gi] > 0) sw <- wald_logp(f$coef[gi], sqrt(V[gi, gi])) }
    ## l* is not concave in general (a stationary point below the maximum was hit in the pilot), so:
    ## full from the null start; restricted from the full fit and from the intercept-only start (best
    ## kept); full again from the best restricted fit (best kept) -> PLR >= 0 by construction.
    b <- function(a, c) if (c$llpen > a$llpen) c else a
    ff <- firth_logit(z, Cg, start = if (isTRUE(nul$cat$converged)) c(na0(nul$cat$coef), 0) else NULL)
    fr <- b(firth_logit(z, Cg, start = replace(ff$coef, gi, 0), fixed = gi),
            firth_logit(z, Cg, start = c(qlogis(min(max(mean(z), 1e-3), 1 - 1e-3)), rep(0, gi - 1)), fixed = gi))
    ff <- b(ff, firth_logit(z, Cg, start = fr$coef))
    st <- 2 * (ff$llpen - fr$llpen); se <- if (!is.null(ff$vcov)) sqrt(ff$vcov[gi, gi]) else NA_real_
    out <- list(stat = st, logp = lrt_logp(st, 1), est = ff$coef[gi], se = se, conv = ff$converged && fr$converged)
    if (fam == "BEZI") { r$nu_sw_logp <- sw; r$nu_firth_stat <- out$stat; r$nu_firth_logp <- out$logp
      r$nu_firth_est <- out$est; r$nu_firth_se <- out$se; r$nu_firth_conv <- out$conv }
    else { r$tau_sw_logp <- sw; r$tau_firth_stat <- out$stat; r$tau_firth_logp <- out$logp
      r$tau_firth_est <- out$est; r$tau_firth_se <- out$se; r$tau_firth_conv <- out$conv }
  } else if (fam == "BEINF") {
    y0 <- as.numeric(y == 0); y1 <- as.numeric(y == 1)
    t0 <- na0(nul$cat$coef[1:pC]); t1 <- na0(nul$cat$coef[pC + 1:pC])
    st0 <- c(t0, 0, t1, 0); inu <- gi; itau <- 2L * gi
    V <- NULL
    if (!cs) { f <- mnl3_fit(y0, y1, Cg, Cg, start = st0); sc <- mnl3_scores(y0, y1, Cg, Cg, f$coef); V <- sandwich_V(sc$H, sc$S) }
    if (!is.null(V)) {
      if (V[inu, inu] > 0) r$nu_sw_logp <- wald_logp(f$coef[inu], sqrt(V[inu, inu]))
      if (V[itau, itau] > 0) r$tau_sw_logp <- wald_logp(f$coef[itau], sqrt(V[itau, itau]))
    }
    b <- function(a, c) if (c$llpen > a$llpen) c else a
    ff <- firth_mnl3(y0, y1, Cg, Cg, start = if (isTRUE(nul$cat$converged)) st0 else NULL)
    fn <- b(firth_mnl3(y0, y1, Cg, Cg, start = replace(ff$coef, inu, 0), fixed = inu),
            firth_mnl3(y0, y1, Cg, Cg, start = NULL, fixed = inu))
    ft <- b(firth_mnl3(y0, y1, Cg, Cg, start = replace(ff$coef, itau, 0), fixed = itau),
            firth_mnl3(y0, y1, Cg, Cg, start = NULL, fixed = itau))
    ff <- b(b(ff, firth_mnl3(y0, y1, Cg, Cg, start = fn$coef)), firth_mnl3(y0, y1, Cg, Cg, start = ft$coef))
    r$nu_firth_stat <- 2 * (ff$llpen - fn$llpen); r$tau_firth_stat <- 2 * (ff$llpen - ft$llpen)
    r$nu_firth_logp <- lrt_logp(r$nu_firth_stat, 1); r$tau_firth_logp <- lrt_logp(r$tau_firth_stat, 1)
    r$nu_firth_est <- ff$coef[inu]; r$tau_firth_est <- ff$coef[itau]
    if (!is.null(ff$vcov)) { r$nu_firth_se <- sqrt(ff$vcov[inu, inu]); r$tau_firth_se <- sqrt(ff$vcov[itau, itau]) }
    r$nu_firth_conv <- ff$converged && fn$converged; r$tau_firth_conv <- ff$converged && ft$converged
  }
  r
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
## solution, which version 2 did not guarantee.
firth_snp_tests <- function(y, C, g, nul, maxround = 6L) {
  fam <- nul$family; pC <- ncol(C); Cg <- cbind(C, g = g); gi <- pC + 1L; eps <- 1e-9
  r <- list(nu_firth_stat = NA_real_, nu_firth_logp = NA_real_, nu_firth_est = NA_real_, nu_firth_se = NA_real_, nu_firth_conv = NA,
            tau_firth_stat = NA_real_, tau_firth_logp = NA_real_, tau_firth_est = NA_real_, tau_firth_se = NA_real_, tau_firth_conv = NA,
            firth_rounds = NA_integer_)
  better <- function(a, b) if (is.null(a) || (is.finite(b$llpen) && b$llpen > a$llpen + eps)) b else a
  if (fam %in% c("BEZI", "BEOI")) {
    z <- as.numeric(if (fam == "BEZI") y == 0 else y == 1)
    s0 <- c(qlogis(min(max(mean(z), 1e-3), 1 - 1e-3)), rep(0, gi - 1))
    ff <- better(if (isTRUE(nul$cat$converged)) firth_logit(z, Cg, start = c(na0(nul$cat$coef), 0)) else NULL, firth_logit(z, Cg, start = s0))
    fr <- better(firth_logit(z, Cg, start = replace(ff$coef, gi, 0), fixed = gi), firth_logit(z, Cg, start = s0, fixed = gi))
    for (k in seq_len(maxround)) {
      ch <- FALSE
      f2 <- firth_logit(z, Cg, start = fr$coef); if (f2$llpen > ff$llpen + eps) { ff <- f2; ch <- TRUE }
      r2 <- firth_logit(z, Cg, start = replace(ff$coef, gi, 0), fixed = gi); if (r2$llpen > fr$llpen + eps) { fr <- r2; ch <- TRUE }
      if (!ch) break
    }
    st <- 2 * (ff$llpen - fr$llpen); se <- if (!is.null(ff$vcov) && ff$vcov[gi, gi] > 0) sqrt(ff$vcov[gi, gi]) else NA_real_
    pre <- if (fam == "BEZI") "nu" else "tau"
    r[[paste0(pre, "_firth_stat")]] <- st; r[[paste0(pre, "_firth_logp")]] <- lrt_logp(st, 1)
    r[[paste0(pre, "_firth_est")]] <- ff$coef[gi]; r[[paste0(pre, "_firth_se")]] <- se
    r[[paste0(pre, "_firth_conv")]] <- isTRUE(ff$converged) && isTRUE(fr$converged); r$firth_rounds <- k
  } else if (fam == "BEINF") {
    y0 <- as.numeric(y == 0); y1 <- as.numeric(y == 1); inu <- gi; itau <- 2L * gi
    ff <- better(if (isTRUE(nul$cat$converged)) firth_mnl3(y0, y1, Cg, Cg, start = c(na0(nul$cat$coef[1:pC]), 0, na0(nul$cat$coef[pC + 1:pC]), 0)) else NULL,
                 firth_mnl3(y0, y1, Cg, Cg, start = NULL))
    fn <- better(firth_mnl3(y0, y1, Cg, Cg, start = replace(ff$coef, inu, 0), fixed = inu), firth_mnl3(y0, y1, Cg, Cg, start = NULL, fixed = inu))
    ft <- better(firth_mnl3(y0, y1, Cg, Cg, start = replace(ff$coef, itau, 0), fixed = itau), firth_mnl3(y0, y1, Cg, Cg, start = NULL, fixed = itau))
    for (k in seq_len(maxround)) {
      ch <- FALSE
      for (rr in list(fn, ft)) { f2 <- firth_mnl3(y0, y1, Cg, Cg, start = rr$coef); if (f2$llpen > ff$llpen + eps) { ff <- f2; ch <- TRUE } }
      n2 <- firth_mnl3(y0, y1, Cg, Cg, start = replace(ff$coef, inu, 0), fixed = inu); if (n2$llpen > fn$llpen + eps) { fn <- n2; ch <- TRUE }
      t2 <- firth_mnl3(y0, y1, Cg, Cg, start = replace(ff$coef, itau, 0), fixed = itau); if (t2$llpen > ft$llpen + eps) { ft <- t2; ch <- TRUE }
      if (!ch) break
    }
    V <- ff$vcov
    r$nu_firth_stat <- 2 * (ff$llpen - fn$llpen); r$tau_firth_stat <- 2 * (ff$llpen - ft$llpen)
    r$nu_firth_logp <- lrt_logp(r$nu_firth_stat, 1); r$tau_firth_logp <- lrt_logp(r$tau_firth_stat, 1)
    r$nu_firth_est <- ff$coef[inu]; r$tau_firth_est <- ff$coef[itau]
    if (!is.null(V)) { if (V[inu, inu] > 0) r$nu_firth_se <- sqrt(V[inu, inu]); if (V[itau, itau] > 0) r$tau_firth_se <- sqrt(V[itau, itau]) }
    r$nu_firth_conv <- isTRUE(ff$converged) && isTRUE(fn$converged); r$tau_firth_conv <- isTRUE(ff$converged) && isTRUE(ft$converged)
    r$firth_rounds <- k
  }
  r
}
