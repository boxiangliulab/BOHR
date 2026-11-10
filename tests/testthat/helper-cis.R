## Synthetic tissue and the reference cis pipeline built from the original engine's functions.

## A synthetic tissue: two chromosomes, every promoter family, an exact mirror pair, a boundary-only pair,
## a three-promoter gene, a skipped promoter, missing activities, missing genotypes and indels.
sim_tissue <- function(seed = 11, n = 110, spacing = 50000, span = 1.6e6) {
  with_seed(seed, {
    s <- paste0("D", seq_len(n))
    pos <- seq(60000, span, by = spacing)
    V <- rbind(data.frame(chr = "chr1", pos = pos), data.frame(chr = "chr2", pos = pos))
    V$end <- V$pos + ifelse(stats::runif(nrow(V)) < 0.1, 3, 0)                       # some indels
    V$id <- sprintf("%s_%d_v", V$chr, V$pos)
    maf <- stats::runif(nrow(V), 0.15, 0.5)
    G <- t(vapply(maf, function(q) as.numeric(stats::rbinom(n, 2, q)), numeric(n))); dimnames(G) <- list(V$id, s)
    G[sample.int(length(G), round(0.01 * length(G)))] <- NA                            # missing genotypes
    Z <- matrix(stats::rnorm(n * 3), n, 3, dimnames = list(s, c("pc1", "pc2", "peer1")))
    zib <- function(gv, p0, p1, d0 = 0, d1 = 0, dmu = 0) {
      gv[is.na(gv)] <- 0
      mu <- stats::plogis(stats::qlogis(0.45) + dmu * gv + 0.2 * Z[, 1]); y <- stats::rbeta(n, mu * 9, (1 - mu) * 9)
      u <- stats::runif(n); q0 <- pmax(p0 + d0 * gv, 0); q1 <- pmax(p1 + d1 * gv, 0)
      y[u < q0] <- 0; y[u >= q0 & u < q0 + q1] <- 1; y }
    gi <- function(chr, p) { o <- which(V$chr == chr); o[which.min(abs(V$pos[o] - p))] }   # nearest variant
    Y <- list()
    Y$PromoterID_101 <- zib(G[gi("chr1", 310000), ], 0.2, 0.15, d0 = 0.15)                 # BEINF, p0 effect
    Y$PromoterID_102 <- zib(G[gi("chr1", 510000), ], 0.25, 0.2, dmu = 0.5)                  # exact mirror pair with 103
    Y$PromoterID_103 <- 1 - Y$PromoterID_102
    Y$PromoterID_104 <- zib(G[gi("chr1", 760000), ], 0.3, 0.1, d1 = 0.12)                   # boundary-only pair with 105
    Y$PromoterID_105 <- ifelse(Y$PromoterID_104 %in% c(0, 1), 1 - Y$PromoterID_104, 0.9 * (1 - Y$PromoterID_104))
    Y$PromoterID_106 <- zib(G[gi("chr1", 1010000), ], 0, 0, dmu = 0.6)                      # BE, n_mid >= 80
    Y$PromoterID_107 <- zib(G[gi("chr1", 1110000), ], 0.3, 0)                               # BEZI (gene with three promoters)
    Y$PromoterID_108 <- zib(G[gi("chr1", 1160000), ], 0, 0.25)                              # BEOI
    Y$PromoterID_109 <- zib(G[gi("chr1", 1210000), ], 0.15, 0.15)                           # BEINF
    Y$PromoterID_110 <- zib(G[gi("chr2", 410000), ], 0.97, 0)                               # skipped (prop_or_sd)
    y111 <- zib(G[gi("chr2", 660000), ], 0.1, 0.03, dmu = -0.4); y111[c(5, 9, 33)] <- NA     # endpoint rule + missing values
    Y$PromoterID_111 <- y111
    Y$PromoterID_112 <- zib(G[gi("chr2", 910000), ], 0.12, 0.12, d0 = -0.05)                # BEINF on chr2
    ph <- do.call(rbind, Y); colnames(ph) <- s
    pr <- data.frame(id = rownames(ph), chr = c(rep("chr1", 9), rep("chr2", 3)),
                     pos = c(300000, 500000, 520000, 750000, 770000, 1000000, 1100000, 1150000, 1200000, 400000, 650000, 900000),
                     gene = c("G1", "G2", "G2", "G3", "G3", "G4", "G5", "G5", "G5", "G6", "G7", "G8"), stringsAsFactors = FALSE)
    list(phenotype = ph, genotype = G, variants = V[, c("id", "chr", "pos", "end")], promoters = pr, covariates = Z)
  })
}

## Reference: run_final2.R prep_promoter + run_chunk (per SNP), run_muF.R (F correction), discovery_nogate.R
## (units, ACAT, leader), units_v3_estimable.R, mirror_pairs.R and units_v3m_mirror.R (v3m), with the
## original functions; window = FLANK (5e5) of bohr2_data.R.
ref_tissue <- function(eng, x) {
  samples <- colnames(x$phenotype); act <- t(x$phenotype)                      # td$act: samples x promoters
  Z <- x$covariates[samples, , drop = FALSE]
  C_all <- eng$cov_matrix(list(pcs = Z, tissue = "t"), Z[, 0, drop = FALSE], 0L)
  D <- x$genotype[, samples, drop = FALSE]; storage.mode(D) <- "integer"; D <- eng$impute_modal(D)   # read_geno: per-SNP mode
  V <- x$variants; P <- x$promoters
  snp <- list(); prom <- list()
  for (i in seq_len(nrow(P))) {
    p <- P$id[i]; pos <- P$pos[i]
    base <- list(promoter = p, gene = P$gene[i], family = NA_character_, skip = NA_character_, n = NA_integer_, n0 = NA_integer_, n1 = NA_integer_)
    z <- eng$prepare_y(act[, p], ncol(C_all) - 1L)
    if (!is.na(z$skip)) { base$skip <- z$skip; prom[[i]] <- base; next }
    if (!any(z$y > 0 & z$y < 1)) { base$skip <- "no_interior"; prom[[i]] <- base; next }
    s <- names(z$y); y <- unname(z$y); Cs <- eng$orth_design(C_all[s, , drop = FALSE])
    base$family <- z$family; base$n <- length(y); base$n0 <- sum(y == 0); base$n1 <- sum(y == 1)
    nul <- eng$promoter_null(y, Cs, z$family)
    on <- which(V$chr == P$chr[i])
    idx <- on[eng$cis_index(list(start = V$pos[on], end = V$end[on]), pos)]
    if (!length(idx)) { base$skip <- "no_snp"; prom[[i]] <- base; next }
    G <- t(D[idx, s, drop = FALSE])                                               # geno$G[s, idx]
    cnt <- apply(G, 2, function(v) tabulate(v + 1L, nbins = 3L)); keep <- which(apply(cnt, 2, min) >= eng$MIN_GENO)
    if (!length(keep)) { base$skip <- "no_snp"; prom[[i]] <- base; next }
    prom[[i]] <- base
    mid <- y > 0 & y < 1; rk <- if (sum(mid) > 0) qr(Cs[mid, , drop = FALSE])$rank else NA_integer_; pm <- rk + 1L; nmid <- sum(mid)
    for (j in keep) {
      r <- tryCatch(eng$snp_tests(y, Cs, G[, j], nul), error = function(e) NULL)
      if (is.null(r)) next
      f <- if (z$family != "BE") tryCatch(eng$firth_snp_tests(y, Cs, G[, j], nul), error = function(e) NULL) else NULL
      fc <- stats::setNames(rep(list(NA), length(BOHR:::FIRTH_COLS)), BOHR:::FIRTH_COLS); if (!is.null(f)) fc <- f[BOHR:::FIRTH_COLS]
      W <- r$mu_stat
      muF <- if (is.finite(W) && nmid > pm + 1L) stats::pf((nmid - pm) * expm1(pmax(W, 0) / nmid), 1, nmid - pm, lower.tail = FALSE, log.p = TRUE) else NA_real_
      af <- (cnt[2, j] + 2 * cnt[3, j]) / (2 * length(y))
      snp[[length(snp) + 1L]] <- c(list(promoter = p, snp = V$id[idx[j]], snp_pos = V$pos[idx[j]], alt_af = af, maf = min(af, 1 - af), min_geno = min(cnt[, j])),
                                   r[BOHR:::SNP_COLS], fc, list(mu_F_logp = muF))
    }
  }
  S <- BOHR:::rows_to_df(snp); PR <- BOHR:::rows_to_df(prom)
  S <- S[order(S$promoter, S$snp, method = "radix"), , drop = FALSE]                  # = merge(..., by = c("promoter", "snp"))
  ## units (discovery_nogate.R:34-49)
  PCOL <- c(mu = "mu_F_logp", p0 = "nu_firth_logp", p1 = "tau_firth_logp")
  EST <- c(mu = "mu_est", p0 = "nu_firth_est", p1 = "tau_firth_est"); PE <- c(mu = "mu_pe", p0 = "p0_pe", p1 = "p1_pe")
  Pf <- PR[is.na(PR$skip) & !is.na(PR$family), , drop = FALSE]
  types <- list(mu = Pf$promoter, p0 = Pf$promoter[Pf$family %in% c("BEZI", "BEINF")], p1 = Pf$promoter[Pf$family %in% c("BEOI", "BEINF")])
  U <- list()
  for (ty in names(PCOL)) for (p in unique(S$promoter[S$promoter %in% types[[ty]]])) {
    s <- S[S$promoter == p, , drop = FALSE]; v <- s[[PCOL[[ty]]]]
    j <- which.min(replace(v, !is.finite(v), Inf))
    U[[length(U) + 1L]] <- list(promoter = p, type = ty, p = eng$acat_log(v), n_snp = sum(is.finite(v)), leader = s$snp[j], leader_logp = v[j],
                                leader_est = s[[EST[[ty]]]][j], leader_pe = s[[PE[[ty]]]][j])
  }
  U <- BOHR:::rows_to_df(U); U <- U[is.finite(U$p), , drop = FALSE]
  i <- match(U$promoter, PR$promoter); U$n <- PR$n[i]; U$n0 <- PR$n0[i]; U$n1 <- PR$n1[i]
  U$n_mid <- U$n - U$n0 - U$n1; U$f0 <- U$n0 / U$n; U$f1 <- U$n1 / U$n
  ## estimability (units_v3_estimable.R:10, 15-17)
  NMID_MIN <- 80L; SHARE_MIN <- 0.10
  U$estimable_mu_rule <- !(U$type == "mu" & U$n_mid < NMID_MIN)
  U$boundary_share <- ifelse(U$type == "p0", U$f0, ifelse(U$type == "p1", U$f1, NA_real_))
  U$estimable <- U$estimable_mu_rule & (U$type == "mu" | (U$boundary_share >= SHARE_MIN & U$boundary_share <= 1 - SHARE_MIN))
  ## mirror table (mirror_pairs.R:22-50) for this one tissue; canonical by the mean over all samples
  g <- P$gene; np <- table(g); two <- names(np)[np == 2]
  num <- function(p) as.numeric(sub("PromoterID_", "", p))
  FL <- NULL
  for (xg in two) { pr <- P$id[g == xg]; pr <- pr[order(num(pr))]; a <- act[, pr[1]]; b <- act[, pr[2]]
    ok <- !is.na(a) & !is.na(b); aa <- a[ok]; bb <- b[ok]
    za <- aa == 0; oa <- aa == 1; zb <- bb == 0; ob <- bb == 1; ia <- !za & !oa; ib <- !zb & !ob
    bex <- any(ok) && all(za == ob) && all(oa == zb)
    iex <- bex && all(ia == ib) && (!any(ia) || max(abs(aa[ia] + bb[ia] - 1)) <= 1e-6)
    pm1 <- mean(a, na.rm = TRUE); pm2 <- mean(b, na.rm = TRUE)
    canon <- if (pm1 >= pm2) pr[1] else pr[2]; part <- if (canon == pr[1]) pr[2] else pr[1]
    FL <- rbind(FL, data.frame(promoter = c(canon, part), mirror_pair_id = xg, boundary_exact = bex, interior_exact = iex,
                               mirror_role = c("canonical", "partner"), mirror_partner = c(part, canon), stringsAsFactors = FALSE)) }
  FL$drop_mu <- FL$interior_exact & FL$mirror_role %in% "partner"                     # mirror_pairs.R:49
  fi <- match(U$promoter, FL$promoter)
  U$boundary_exact <- FL$boundary_exact[fi]; U$interior_exact <- FL$interior_exact[fi]; U$drop_mu <- FL$drop_mu[fi]
  U$mirror_role <- FL$mirror_role[fi]; U$mirror_partner <- FL$mirror_partner[fi]; U$mirror_pair_id <- FL$mirror_pair_id[fi]
  for (k in c("boundary_exact", "interior_exact", "drop_mu")) U[[k]][is.na(U[[k]])] <- FALSE   # discovery_nogate.R:53
  ## v3m (units_v3m_mirror.R:29, 35-37, 44)
  lab <- c(mu = "mu", p0 = "p1", p1 = "p0")
  be <- U$boundary_exact %in% TRUE
  U$dedup_drop <- ifelse(U$type == "mu", U$drop_mu, be & U$type == "p1")
  U$duplicate_of <- ifelse(U$dedup_drop, paste0(U$mirror_partner, ":", lab[U$type]), NA_character_)
  U$q <- NA_real_; k <- U$dedup_drop == FALSE & U$estimable == TRUE; U$q[k] <- stats::p.adjust(U$p[k], "BH")
  list(snp = S, prom = PR, units = U, mirror = FL)
}
## lines of the original scripts that ref_tissue() transcribes
ref_tissue_lines <- list(
  discovery_nogate.R = c('types <- rbind(P[, .(promoter, type = "mu")], P[family %in% c("BEZI", "BEINF"), .(promoter, type = "p0")], P[family %in% c("BEOI", "BEINF"), .(promoter, type = "p1")])',
                         "s[, { j <- which.min(replace(get(col), !is.finite(get(col)), Inf))",
                         ".(p = acat_log(get(col)), n_snp = sum(is.finite(get(col))), leader = snp[j], leader_logp = get(col)[j], leader_est = get(EST[[ty]])[j],",
                         "U <- U[is.finite(p)]",
                         's <- merge(s, as.data.table(F$snp)[, .(promoter, snp, mu_F_logp)], by = c("promoter", "snp"), all.x = TRUE)'),
  units_v3_estimable.R = c("NMID_MIN <- 80L; SHARE_MIN <- 0.10",
                           'U[, estimable_mu_rule := !(type == "mu" & n_mid < NMID_MIN)]',
                           'U[, boundary_share := fifelse(type == "p0", f0, fifelse(type == "p1", f1, NA_real_))]',
                           'U[, estimable := estimable_mu_rule & (type == "mu" | (boundary_share >= SHARE_MIN & boundary_share <= 1 - SHARE_MIN))]'),
  mirror_pairs.R = c("bex <- any(ok) && all(za == ob) && all(oa == zb)",
                     "iex <- bex && all(ia == ib) && (!any(ia) || max(abs(aa[ia] + bb[ia] - 1)) <= 1e-6)",
                     "M[, canonical := fifelse(pooled_mean_1 >= pooled_mean_2, promoter_1, promoter_2)]",
                     'P[, drop_mu := interior_exact & mirror_role %in% "partner"]'),
  units_v3m_mirror.R = c('U[, dedup_drop := fifelse(type == "mu", drop_mu, be & type == "p1")]',
                         'U[, duplicate_of := fifelse(dedup_drop, paste0(mirror_partner, ":", lab[type]), NA_character_)]',
                         "U[, q := NA_real_]; U[dedup_drop == FALSE & estimable == TRUE, q := p.adjust(p, \"BH\")]"))
