#!/usr/bin/env Rscript
### Linear comparators v3: all three QTLtools modes, on EXACTLY the same samples and the
### same pre-drawn permutations as the BOHR runs (permutation seed 20261001 + promoter index, as in bohr2_sim.R).
###
### A. No covariates (the simulation design). With no covariates "residualise" is only centring, so:
###   raw        : QTLtools default -- Pearson t-test of raw y on dosage, df n - 2              (= stored ols_p)
###   normal     : QTLtools v1.2 --normal -- z = qnorm((rank(y) - 0.5) / max(rank(y))), Pearson t, df n - 2
###   rintfirst  : GTEx order (RINT first) -- z = qnorm((rank(y) - 0.5) / n), Pearson t, df n - 2
###   (normal and rintfirst differ only in the denominator, which matters when many values tie at the top, e.g. y = 1)
###
### B. Covariate demonstration (calibration mode only): QTLtools --normal fails on real data because it residualises
###    the phenotype on covariates BEFORE ranking; tied 0/1 values then get ranks set by the covariates, and genotype PCs
###    correlate with genotype. Here each promoter gets one ancestry-like covariate C = rho * std(g) + sqrt(1 - rho^2) * e
###    (rho = 0.5, e ~ N(0, 1), seed 30000 + promoter index). C belongs to the genotype "slot" (like genotype PCs in the
###    real permutation design), so it is permuted TOGETHER with g. y is independent of both (null).
###   cov_raw       : QTLtools raw with --cov: residualise y on C, Pearson t with UNresidualised g, df n - 2
###   cov_normal    : QTLtools --normal with --cov: residualise y on C, THEN rank-normalise (max(rank)), Pearson t, df n - 2
###   cov_rintfirst : RINT first (n), then lm(z ~ g + C) t-test for g, df n - 3 (Figure 4's simulation recipe)
###
### Usage: Rscript baselines_v3.R <calib|power> <M1|M2> <dataset> <N> <out_rds> [K]
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); mode <- a[1]; method <- a[2]; dataset <- a[3]; N <- as.integer(a[4]); out <- a[5]
K <- if (length(a) >= 6) as.integer(a[6]) else 1000L; RHO <- 0.5
PSEED <- as.integer(Sys.getenv("PERM_SEED_BASE", "20261001"))   # permutation seed base (default = original runs)
NM <- "/projects/group/bliu/e1327989/GTEx/new_simulation/newmethod_bohr2"; CAL <- "/projects/group/bliu/e1327989/GTEx/new_simulation/calibration_perm1000"
ROOT <- "/projects/group/bliu/e1327989/GTEx/new_simulation"
src <- readLines(file.path(NM, "scripts/bohr2_sim.R"))
eval(parse(text = src[grep("^prepare_y <- function", src):(grep("^S <- load_sim", src) - 1)]))   # prepare_y, orth_design, load_sim
S <- load_sim(method, dataset, N)
ids <- colnames(S$Y); if (mode == "calib") ids <- intersect(readLines(Sys.getenv("PERM_IDS")), ids)
tp <- function(r, df) { r <- pmax(pmin(r, 1), -1); 2 * pt(-abs(r * sqrt(df / pmax(1 - r^2, 1e-300))), df) }
cor_cols <- function(v, G) { vc <- v - mean(v); Gc <- sweep(G, 2, colMeans(G)); drop(crossprod(vc, Gc)) / sqrt(sum(vc^2) * colSums(Gc^2)) }
cor_mat  <- function(V, G) { Vc <- sweep(V, 2, colMeans(V)); Gc <- sweep(G, 2, colMeans(G)); colSums(Vc * Gc) / sqrt(colSums(Vc^2) * colSums(Gc^2)) }
resid_on <- function(v, Cm) {            # residualise a fixed vector v (or matrix, column-wise) on each column of Cm (+ intercept)
  Cc <- sweep(Cm, 2, colMeans(Cm)); V <- if (is.matrix(v)) sweep(v, 2, colMeans(v)) else matrix(v - mean(v), nrow(Cm), ncol(Cm))
  sweep(Cc, 2, colSums(V * Cc) / colSums(Cc^2), `*`) * -1 + V }
rint_n   <- function(v) qnorm((rank(v) - .5) / length(v))
rint_max <- function(v) { r <- rank(v); qnorm((r - .5) / max(r)) }
one <- function(id) {
  y <- S$Y[, id]; names(y) <- rownames(S$Y); z <- prepare_y(y, 0L)
  if (!is.na(z$skip) || !any(z$y > 0 & z$y < 1)) return(NULL)
  s <- names(z$y); yy <- unname(z$y); g <- as.numeric(S$G[s, id]); n <- length(yy)
  if (anyNA(g) || min(tabulate(g + 1L, 3L)) < 5L) return(NULL)
  zn <- rint_n(yy); zm <- rint_max(yy)
  if (mode == "power") { G <- matrix(g, ncol = 1); perm <- NA_integer_; P <- matrix(seq_len(n), ncol = 1) } else {
    set.seed(PSEED + match(id, colnames(S$Y))); P <- replicate(K, sample.int(n)); G <- matrix(g[P], n, K); perm <- seq_len(K) }
  D <- data.table(id = id, perm = perm, raw = tp(cor_cols(yy, G), n - 2), normal = tp(cor_cols(zm, G), n - 2), rintfirst = tp(cor_cols(zn, G), n - 2))
  if (mode == "power") {
    b <- function(v) unname(2 * coef(lm(v ~ g))[2])
    D[, `:=`(eff_raw = b(yy), eff_normal = b(zm), eff_rintfirst = b(zn), boundary_share = mean(yy %in% c(0, 1)), n_ties_max = max(table(yy)))]
  } else {                                  # covariate demonstration under the null
    set.seed(30000L + match(id, colnames(S$Y)))
    C <- RHO * as.numeric(scale(g)) + sqrt(1 - RHO^2) * rnorm(n); Cm <- matrix(C[P], n, K)   # C moves with g
    Ry <- resid_on(yy, Cm)                                                                    # y residualised on C (per permutation)
    Zr <- apply(Ry, 2, rint_max)                                                              # --normal: rank AFTER residualising
    Rz <- resid_on(zn, Cm); Rg <- resid_on(G, Cm)                                             # RINT first, then partial out C from z and g
    D[, `:=`(cov_raw = tp(cor_mat(Ry, G), n - 2), cov_normal = tp(cor_mat(Zr, G), n - 2), cov_rintfirst = tp(cor_mat(Rz, Rg), n - 3),
             boundary_share = mean(yy %in% c(0, 1)))]
  }
  D
}
R <- rbindlist(lapply(ids, one), fill = TRUE)
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(meta = list(mode = mode, method = method, dataset = dataset, N = N, K = K, rho = RHO), res = R), out)
cat(sprintf("%s %s %s N=%d: promoters %d rows %d\n", mode, method, dataset, N, uniqueN(R$id), nrow(R)))
