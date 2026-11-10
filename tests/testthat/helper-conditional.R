## Synthetic tissue with multi-signal units, and the reference conditional analysis built from the original engine's
## functions with the steps of run_conditional_v3c.R (forward rounds and leave-one-out scans) and merge_conditional_v3c.R
## (backward step), transcribed line by line (each transcribed line is grep-checked in the original scripts).

## Two chromosomes, variants every 50 kb; P201 mu with two causal variants and a perfect-LD copy of the second (an earlier
## position: ties go to the first variant in row order); P202 BEINF with p0 effects of two variants; P203/P204 an exact
## mirror pair (P203 goes to 0 with the dosage of one variant: a p0 unit at the silent promoter); P205 BEZI null; P206
## BEOI with a p1 effect. Also returns an absolute activity matrix in which only P204 (the partner) responds.
sim_cond <- function(seed = 31, n = 170, spacing = 50000, span = 1.5e6) {
  with_seed(seed, {
    s <- paste0("D", seq_len(n))
    pos <- seq(60000, span, by = spacing)
    V <- rbind(data.frame(chr = "chr1", pos = pos), data.frame(chr = "chr2", pos = pos))
    V$end <- V$pos + ifelse(stats::runif(nrow(V)) < 0.08, 2, 0)
    V$id <- sprintf("%s_%d_v", V$chr, V$pos)
    maf <- stats::runif(nrow(V), 0.2, 0.5)
    G <- t(vapply(maf, function(q) as.numeric(stats::rbinom(n, 2, q)), numeric(n))); dimnames(G) <- list(V$id, s)
    gi <- function(chr, p) { o <- which(V$chr == chr); o[which.min(abs(V$pos[o] - p))] }
    G[gi("chr1", 410000), ] <- G[gi("chr1", 510000), ]                                 # perfect LD, earlier position
    G[sample.int(length(G), round(0.005 * length(G)))] <- NA
    Z <- matrix(stats::rnorm(n * 2), n, 2, dimnames = list(s, c("pc1", "pc2")))
    g0 <- function(chr, p) { v <- G[gi(chr, p), ]; v[is.na(v)] <- 0; v }
    beta_y <- function(eta, phi = 10) stats::rbeta(n, stats::plogis(eta) * phi, (1 - stats::plogis(eta)) * phi)
    Y <- list()
    Y$PromoterID_201 <- beta_y(-0.3 + 0.55 * g0("chr1", 310000) - 0.45 * g0("chr1", 510000) + 0.2 * Z[, 1])
    y <- beta_y(0.2); u <- stats::runif(n)
    q0 <- pmin(0.02 + 0.3 * g0("chr1", 860000) + 0.3 * g0("chr1", 1110000), 0.85)
    y[u < q0] <- 0; y[u > 0.88] <- 1; Y$PromoterID_202 <- y
    y <- beta_y(0.4); u <- stats::runif(n); y[u < 0.06 + 0.2 * g0("chr1", 1310000)] <- 0; y[u > 0.9] <- 1
    Y$PromoterID_203 <- y; Y$PromoterID_204 <- 1 - y
    y <- beta_y(0); y[stats::runif(n) < 0.3] <- 0; Y$PromoterID_205 <- y
    y <- beta_y(-0.2); y[stats::runif(n) < 0.12 + 0.18 * g0("chr2", 760000)] <- 1; Y$PromoterID_206 <- y
    ph <- do.call(rbind, Y); colnames(ph) <- s
    pr <- data.frame(id = rownames(ph), chr = c(rep("chr1", 4), "chr2", "chr2"), pos = c(400000, 1000000, 1300000, 1320000, 500000, 750000),
                     gene = c("G1", "G2", "G3", "G3", "G4", "G5"), stringsAsFactors = FALSE)
    ## absolute activity: log2 counts; only P204 (the partner) moves with the driver of the P203 zeros
    lib <- stats::rnorm(n, 8, 0.3)
    A <- rbind(PromoterID_203 = lib + stats::rnorm(n, 0, 0.4), PromoterID_204 = lib + 0.6 * g0("chr1", 1310000) + stats::rnorm(n, 0, 0.4))
    colnames(A) <- s
    list(phenotype = ph, genotype = G, variants = V[, c("id", "chr", "pos", "end")], promoters = pr, covariates = Z, absolute = A)
  })
}

## Reference: run_conditional_v3c.R:92-180 for the given units of one tissue (real arm, mode "new"), with the original functions;
## window = FLANK of bohr2_data.R. units: data frame with promoter, type, leader, leader_logp and p (unit ACAT P).
## Then the backward step of merge_conditional_v3c.R:30-34. Returns units, signals (with loo_* and backward_keep), rounds, loo.
ref_conditional <- function(eng, x, units, thr, MAXSIG = 10L) {
  samples <- colnames(x$phenotype); act <- t(x$phenotype)
  Z <- x$covariates[samples, , drop = FALSE]
  C_all <- eng$cov_matrix(list(pcs = Z, tissue = "t"), Z[, 0, drop = FALSE], 0L)
  D <- x$genotype[, samples, drop = FALSE]; storage.mode(D) <- "integer"; D <- eng$impute_modal(D)   # read_geno(): per-SNP mode
  V <- x$variants; P <- x$promoters
  acat_log <- eng$acat_log
  prep_promoter <- function(p) {                                                       # :92-100
    z <- eng$prepare_y(act[, p], ncol(C_all) - 1L); s <- names(z$y); y <- unname(z$y)
    i <- match(p, P$id); on <- which(V$chr == P$chr[i])
    idx <- on[eng$cis_index(list(start = V$pos[on], end = V$end[on]), P$pos[i])]; G <- t(D[idx, s, drop = FALSE])
    cnt <- apply(G, 2, function(x) tabulate(x + 1L, nbins = 3L)); keep <- which(apply(cnt, 2, min) >= eng$MIN_GENO)
    list(promoter = p, family = z$family, y = y, C = C_all[s, , drop = FALSE], G = G[, keep, drop = FALSE])
  }
  PP <- lapply(unique(units$promoter), prep_promoter); names(PP) <- unique(units$promoter)
  mid_rank <- function(y, Cc) { m <- y > 0 & y < 1; c(n_mid = sum(m), rank = qr(Cc[m, , drop = FALSE])$rank) }
  test_one <- function(type, y, Cc, g, nul, mr = NULL) {
    if (type == "mu") { r <- eng$snp_tests(y, Cc, g, nul)
      n <- mr[["n_mid"]]; pp <- mr[["rank"]] + 1L; W <- pmax(r$mu_stat, 0)
      lF <- if (is.finite(W) && n > pp + 1L) pf((n - pp) * expm1(W / n), 1, n - pp, lower.tail = FALSE, log.p = TRUE) else NA_real_
      c(est = unname(r$mu_est), se = unname(r$mu_se), logp = lF, logp_chi = unname(r$mu_lrt_logp)) }
    else { r <- eng$firth_snp_tests(y, Cc, g, nul)
      if (type == "p0") c(est = unname(r$nu_firth_est), se = unname(r$nu_firth_se), logp = unname(r$nu_firth_logp), logp_chi = NA_real_)
      else c(est = unname(r$tau_firth_est), se = unname(r$tau_firth_se), logp = unname(r$tau_firth_logp), logp_chi = NA_real_) } }
  run_pass <- function(jobs) {                                                          # :116-134, serial; rows in test order
    nl <- lapply(jobs, function(jb) { Pp <- PP[[jb$promoter]]
      Cc <- eng$orth_design(cbind(Pp$C, Pp$G[, jb$cs, drop = FALSE]))
      list(Cc = Cc, nul = eng$promoter_null(Pp$y, Cc, Pp$family), dropped = attr(Cc, "dropped"), mr = mid_rank(Pp$y, Cc)) })
    res <- lapply(seq_along(jobs), function(b) { jb <- jobs[[b]]; Pp <- PP[[jb$promoter]]; N <- nl[[b]]
      do.call(rbind, lapply(jb$test, function(sn) {
        v <- tryCatch(test_one(jb$type, Pp$y, N$Cc, Pp$G[, sn], N$nul, N$mr), error = function(e) c(est = NA, se = NA, logp = NA, logp_chi = NA))
        data.frame(b = b, snp = sn, est = unname(v["est"]), se = unname(v["se"]), logp = unname(v["logp"]), logp_chi = unname(v["logp_chi"]),
                   stringsAsFactors = FALSE) })) })
    list(res = res, dropped = vapply(nl, function(x) length(x$dropped), 0L))
  }
  ## forward rounds (:137-164)
  st <- lapply(seq_len(nrow(units)), function(i) list(promoter = units$promoter[i], type = units$type[i], signals = units$leader[i], active = TRUE, stop = NA_character_))
  sig_rows <- lapply(seq_len(nrow(units)), function(i) data.frame(promoter = units$promoter[i], type = units$type[i], k = 1L, snp = units$leader[i],
    selected_round = 0L, logp_at_selection = units$leader_logp[i], logp_chi_at_selection = NA_real_, acat_at_selection = units$p[i], stringsAsFactors = FALSE))
  rounds <- list()
  for (r in seq_len(MAXSIG)) {
    act_i <- which(vapply(st, `[[`, TRUE, "active")); if (!length(act_i)) break
    jobs <- lapply(act_i, function(i) { s <- st[[i]]; list(i = i, promoter = s$promoter, type = s$type, cs = s$signals,
      test = setdiff(colnames(PP[[s$promoter]]$G), s$signals)) })
    ps <- run_pass(jobs)
    for (jb_i in seq_along(jobs)) { i <- jobs[[jb_i]]$i; z <- ps$res[[jb_i]]
      acat <- acat_log(z$logp); best <- if (nrow(z) && any(is.finite(z$logp))) z[which.min(z$logp), ] else NULL
      rounds[[length(rounds) + 1]] <- data.frame(promoter = st[[i]]$promoter, type = st[[i]]$type, round = r, n_cond = length(st[[i]]$signals),
        n_tested = sum(is.finite(z$logp)), n_failed = sum(!is.finite(z$logp)), acat_p = acat, acat_p_chi = acat_log(z$logp_chi),
        best_snp = if (is.null(best)) NA_character_ else best$snp, best_logp = if (is.null(best)) NA_real_ else best$logp,
        cols_dropped = ps$dropped[jb_i], stringsAsFactors = FALSE)
      go <- is.finite(acat) && acat <= thr && !is.null(best) && length(st[[i]]$signals) < MAXSIG
      if (go) { st[[i]]$signals <- c(st[[i]]$signals, best$snp)
        sig_rows[[length(sig_rows) + 1]] <- data.frame(promoter = st[[i]]$promoter, type = st[[i]]$type, k = length(st[[i]]$signals), snp = best$snp,
          selected_round = r, logp_at_selection = best$logp, logp_chi_at_selection = best$logp_chi, acat_at_selection = acat, stringsAsFactors = FALSE) }
      else { st[[i]]$active <- FALSE
        st[[i]]$stop <- if (!is.finite(acat)) "no finite conditional p" else if (acat > thr) "window ACAT above threshold" else "cap" } }
  }
  ## leave-one-out (:166-174)
  loo <- NULL
  multi <- which(vapply(st, function(s) length(s$signals) >= 2, TRUE))
  if (length(multi)) {
    jobs <- unlist(lapply(multi, function(i) { s <- st[[i]]; lapply(seq_along(s$signals), function(k) list(i = i, k = k, promoter = s$promoter, type = s$type,
      cs = s$signals[-k], test = setdiff(colnames(PP[[s$promoter]]$G), s$signals[-k]))) }), recursive = FALSE)
    ps <- run_pass(jobs)
    loo <- do.call(rbind, lapply(seq_along(jobs), function(b) { z <- ps$res[[b]]; z$promoter <- jobs[[b]]$promoter; z$type <- jobs[[b]]$type; z$k <- jobs[[b]]$k; z$b <- NULL; z }))
  }
  signals <- do.call(rbind, sig_rows)
  signals <- signals[order(match(paste(signals$promoter, signals$type), paste(units$promoter, units$type)), signals$k), , drop = FALSE]
  signals$loo_est <- signals$loo_se <- signals$loo_logp <- signals$loo_logp_chi <- signals$loo_acat <- NA_real_; signals$loo_top_snp <- NA_character_
  if (!is.null(loo)) for (j in seq_len(nrow(signals))) {
    z <- loo[loo$promoter == signals$promoter[j] & loo$type == signals$type[j] & loo$k == signals$k[j], , drop = FALSE]; if (!nrow(z)) next
    m <- match(signals$snp[j], z$snp)
    signals$loo_est[j] <- z$est[m]; signals$loo_se[j] <- z$se[m]; signals$loo_logp[j] <- z$logp[m]; signals$loo_logp_chi[j] <- z$logp_chi[m]
    signals$loo_acat[j] <- acat_log(z$logp); signals$loo_top_snp[j] <- z$snp[which.min(replace(z$logp, !is.finite(z$logp), Inf))]   # merge :30
  }
  n_sig <- stats::ave(seq_len(nrow(signals)), signals$promoter, signals$type, FUN = length)
  signals$backward_keep <- ifelse(n_sig == 1L, TRUE, ifelse(is.na(signals$loo_acat), NA, signals$loo_acat <= thr))      # merge :33
  U <- data.frame(promoter = units$promoter, type = units$type, n_signals = vapply(st, function(s) length(s$signals), 0L),
                  stop = vapply(st, `[[`, "", "stop"), stringsAsFactors = FALSE)
  list(units = U, signals = signals, rounds = do.call(rbind, rounds), loo = loo)
}
## lines of the original scripts that ref_conditional() transcribes
ref_conditional_lines <- list(
  run_conditional_v3c.R = c(
    "mid_rank <- function(y, Cc) { m <- y > 0 & y < 1; c(n_mid = sum(m), rank = qr(Cc[m, , drop = FALSE])$rank) }",
    "n <- mr[[\"n_mid\"]]; pp <- mr[[\"rank\"]] + 1L; W <- pmax(r$mu_stat, 0)",
    "lF <- if (is.finite(W) && n > pp + 1L) pf((n - pp) * expm1(W / n), 1, n - pp, lower.tail = FALSE, log.p = TRUE) else NA_real_",
    "if (type == \"p0\") c(est = unname(r$nu_firth_est), se = unname(r$nu_firth_se), logp = unname(r$nu_firth_logp), logp_chi = NA_real_)",
    "Cc <- orth_design(cbind(Pp$C, Pp$G[, jb$cs, drop = FALSE]))",
    "list(Cc = Cc, nul = promoter_null(Pp$y, Cc, Pp$family), dropped = attr(Cc, \"dropped\"), mr = mid_rank(Pp$y, Cc)) }, mc.cores = NC, mc.preschedule = FALSE)",
    "acat <- acat_log(z$logp); best <- if (nrow(z) && any(is.finite(z$logp))) z[which.min(logp)] else NULL",
    "go <- is.finite(acat) && acat <= thr && !is.null(best) && length(st[[i]]$signals) < MAXSIG && !PERM",
    "st[[i]]$stop <- if (PERM) \"perm: one round\" else if (!is.finite(acat)) \"no finite conditional p\" else if (acat > thr) \"window ACAT above threshold\" else \"cap\" } }",
    "cs = s$signals[-k], test = setdiff(colnames(PP[[s$promoter]]$G), s$signals[-k]))) }), recursive = FALSE)",
    "test = setdiff(colnames(PP[[s$promoter]]$G), s$signals)) })",
    "MAXSIG <- 10L; PERM <- arm != \"real\""),
  merge_conditional_v3c.R = c(
    "if (nrow(nl)) { bw <- nl[, .(loo_acat = acat_log(logp), loo_top_snp = snp[which.min(replace(logp, !is.finite(logp), Inf))]), by = .(promoter, type, k)]",
    "ns[, backward_keep := fifelse(n_sig == 1L, TRUE, fifelse(is.na(loo_acat), NA, loo_acat <= thr))][, n_sig := NULL]"))

## Reference labels of exact-pair boundary events (mirror_abs_labels.R:48-50, 60-66, 71-75) for one tissue.
## U: significant p0 units on boundary_exact pairs (S, O, leader); kept: kept signals (promoter, k, snp); A: absolute (promoters x samples).
ref_abs_labels <- function(eng, x, U, kept, A, ZC = 1.96) {
  samples <- colnames(x$phenotype); Z <- x$covariates[samples, , drop = FALSE]
  C <- eng$cov_matrix(list(pcs = Z, tissue = "t"), Z[, 0, drop = FALSE], 0L)
  D <- x$genotype[, samples, drop = FALSE]; storage.mode(D) <- "integer"; G <- t(eng$impute_modal(D))
  fit <- function(X, y) { f <- lm.fit(X, y); p <- ncol(X); if (f$rank < p) return(c(b = NA, se = NA, z = NA, sdy = sd(y)))
    s2 <- sum(f$residuals^2) / (length(y) - p); V <- chol2inv(qr.R(f$qr)); j <- which(f$qr$pivot == p)
    b <- unname(f$coefficients[p]); se <- sqrt(V[j, j] * s2); c(b = b, se = se, z = b / se, sdy = sd(y)) }
  res <- lapply(seq_len(nrow(U)), function(i) { r <- U[i, ]
    cond <- intersect(setdiff(kept$snp[kept$promoter == r$S], r$leader), colnames(G))
    dose <- as.numeric(G[, r$leader]); Xc <- cbind(C, G[, cond, drop = FALSE], dose); dropped <- FALSE
    if (qr(Xc)$rank < ncol(Xc)) { Xc <- cbind(C, dose); dropped <- length(cond) > 0 }
    fS <- fit(Xc, A[r$S, samples]); fO <- fit(Xc, A[r$O, samples])
    data.frame(S = r$S, O = r$O, z_S = fS[["z"]], z_O = fO[["z"]], stringsAsFactors = FALSE) })
  X <- do.call(rbind, res)
  regS <- abs(X$z_S) >= ZC; regO <- abs(X$z_O) >= ZC
  X$outcome <- ifelse(is.na(X$z_S) | is.na(X$z_O), "not estimable", ifelse(regS & !regO, "S only", ifelse(regO & !regS, "O only", ifelse(regS & regO, "both", "neither"))))
  X$label <- ifelse(X$outcome == "O only", paste0(X$O, ":p1"), paste0(X$S, ":p0"))
  X
}
ref_abs_lines <- list(mirror_abs_labels.R = c(
  "s2 <- sum(f$residuals^2) / (length(y) - p); V <- chol2inv(qr.R(f$qr)); j <- which(f$qr$pivot == p)",
  "dose <- as.numeric(G[, r$leader]); Xc <- cbind(C, G[, cond, drop = FALSE], dose); dropped <- FALSE",
  "if (qr(Xc)$rank < ncol(Xc)) { Xc <- cbind(C, dose); dropped <- length(cond) > 0 }",
  "cond <- if (isTRUE(r$sig)) intersect(setdiff(cs[S == r$S, snp], r$leader), colnames(G)) else character(0)",
  "X[, outcome := fcase(is.na(z_S) | is.na(z_O), \"not estimable\", regS & !regO, \"S only\", regO & !regS, \"O only\", regS & regO, \"both\", default = \"neither\")]"))
