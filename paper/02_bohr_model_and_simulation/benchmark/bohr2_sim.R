#!/usr/bin/env Rscript
### BOHR v2 (Figure 4 new method, 2026-10-01) applied to the simulations.
### Engine files are Figure 4's, copied VERBATIM to ../code_copy (md5 in SOURCE_MD5.txt) and sourced
### unchanged: bohr2_engine.R (promoter_null, snp_tests) + bohr2_addon3.R (firth_snp_tests), i.e. the
### same engine as Figure 4's production run (run_final2.R). From bohr2_data.R only prepare_y() and
### orth_design() are copied below (that file also loads VariantAnnotation, not needed here).
### No covariates in the simulation -> C = intercept only; the C3 covariate cap does not apply.
### Per tested genotype vector we record:
###   BOHR v2 : mu / nu / tau single-parameter LRT (log p), their Wald (log p), sigma-free joint LRT,
###             Firth penalised LRT for nu / tau (log p), convergence and separation flags
###   Linear  : (a) logit-normal mu-LRT exactly as the old Fig. 2 baseline (y_logit ~ gt, sigma ~ gt,
###                 mu-LRT drops gt from mu, gamlss NO) + its Wald
###             (b) QTLtools-style OLS: y ~ g, t-test on the slope (raw activity)
### Promoter-level rules as production run_final2.R: prepare_y (sd <= 0.05 / >= 95% boundary skip,
### <= 5% boundary trimming, family by presence of 0/1), >= 1 interior value, >= 5 samples per genotype
### class, eligibility mid >= .70 (mu), n0/n >= .35 (p0), n1/n >= .35 (p1).
### Modes:
###   calib : fixed-50 promoters (PERM_IDS file), own genotype permuted K times (all permutations drawn
###           before any fitting, seed per promoter), one row per permutation
###   power : every promoter, its actual genotype, one row per promoter (with effect_type/effect_size)
### Usage: Rscript bohr2_sim.R <calib|power> <M1|M2> <orig|realfloor|big> <N> <out_rds> [cores] [K] [chunk nchunk]
suppressPackageStartupMessages({ library(data.table); library(parallel); library(gamlss) })
a <- commandArgs(trailingOnly = TRUE)
mode <- a[1]; method <- a[2]; dataset <- a[3]; N <- as.integer(a[4]); out_rds <- a[5]
cores <- if (length(a) >= 6) as.integer(a[6]) else 8L
K <- if (length(a) >= 7) as.integer(a[7]) else 1000L
CH <- if (length(a) >= 9) as.integer(a[8]) else 1L; NCH <- if (length(a) >= 9) as.integer(a[9]) else 1L
ROOT <- "/projects/group/bliu/e1327989/GTEx/new_simulation"; NM <- file.path(ROOT, "newmethod_bohr2")
CAL <- file.path(ROOT, "calibration_perm1000")
source(file.path(NM, "code_copy", "bohr2_engine.R")); source(file.path(NM, "code_copy", "bohr2_addon3.R"))
MIN_GENO <- 5L
PSEED <- as.integer(Sys.getenv("PERM_SEED_BASE", "20261001"))   # permutation seed base (default = all original runs)
## ---- copied verbatim from code_copy/bohr2_data.R ----
prepare_y <- function(y, p_cov) {
  y <- y[!is.na(y)]
  if (length(y) < 20 || length(y) < p_cov) return(list(skip = "n"))
  p0 <- mean(y == 0); p1 <- mean(y == 1)
  if (p0 >= 0.95 || p1 >= 0.95 || sd(y) <= 0.05) return(list(skip = "prop_or_sd"))
  keep <- rep(TRUE, length(y))
  if (p0 <= 0.05) keep <- keep & y > 0
  if (p1 <= 0.05) keep <- keep & y < 1
  y <- y[keep]
  z <- any(y == 0); o <- any(y == 1)
  list(skip = NA_character_, y = y,
       family = if (z && o) "BEINF" else if (z) "BEZI" else if (o) "BEOI" else "BE")
}
orth_design <- function(C, tol = 1e-7) {
  q <- qr(C, tol = tol); r <- q$rank
  Q <- qr.Q(q)[, seq_len(r), drop = FALSE]
  s <- sign(diag(qr.R(q))[seq_len(r)]); s[s == 0] <- 1
  Q <- sweep(Q, 2, s, "*") * sqrt(nrow(C))
  Q[, 1] <- 1
  dimnames(Q) <- list(rownames(C), c("(Intercept)", if (r > 1) paste0("z", seq_len(r - 1))))
  attr(Q, "dropped") <- colnames(C)[q$pivot[-seq_len(r)]]
  Q
}
## ---- data ----
load_sim <- function(method, dataset, N) {
  e <- new.env()
  if (method == "M1") {
    d <- switch(dataset, orig = file.path(ROOT, "method1_direct", sprintf("sim_N%d_G500", N)),
                realfloor = file.path(CAL, "realfloor_pheno", sprintf("M1_N%d", N)),
                nothr = file.path(NM, "m1_nothr", sprintf("sim_N%d_G500_realfloor", N)),           # no-threshold generation, floor L = 0.005
                nothr2000 = file.path(NM, "m1_nothr_g2000", sprintf("sim_N%d_G2000_realfloor", N)), # + 2000 genes, independent genotype per gene
                file.path(NM, "m1_trunc", sub("^floor", "L", dataset), sprintf("sim_N%d_G500_realfloor", N)))   # floor0.001 -> m1_trunc/L0.001
    load(file.path(d, "sim_data.RData"), envir = e)
    Y <- e$activity_mat; rownames(Y) <- sprintf("s%04d", seq_len(nrow(Y)))
    G <- if (!is.null(e$geno_mat)) { Gm <- e$geno_mat[, colnames(Y), drop = FALSE]; dimnames(Gm) <- list(rownames(Y), colnames(Y)); Gm } else   # per-gene genotypes
         matrix(rep(e$genotypes, ncol(Y)), nrow(Y), dimnames = list(rownames(Y), colnames(Y)))                                               # one shared SNP
    des <- fread(file.path(d, "ground_truth.txt"))[, .(id = gene_id, effect_type, effect_size, true_mu_effect, true_nu_effect, true_tau_effect)]
  } else {
    d <- switch(dataset, orig = file.path(ROOT, "method2_polyester", sprintf("sim_%d", N)),
                big = file.path(CAL, "local_sims", "M2_big", sprintf("sim_%d", N)),
                nothr_b1 = file.path(NM, "m2_nothr", "b1", sprintf("sim_%d", N)),     # no-threshold generation, gene batch 1 (gene_seed 42)
                nothr_b2 = file.path(NM, "m2_nothr", "b2", sprintf("sim_%d", N)),     # no-threshold generation, gene batch 2 (gene_seed 43)
                rj_b1 = file.path(NM, "m2_nothr_realjunc", "b1", sprintf("sim_%d", N)),   # same as nothr_b1, junctions re-filtered with the
                rj_b2 = file.path(NM, "m2_nothr_realjunc", "b2", sprintf("sim_%d", N)))   # real rule (unique reads >= 5), proActiv re-run (2026-10-06)
    load(file.path(d, "recal_att_f.RData"), envir = e); X <- as.data.frame(e$recal_att_f)
    dfile <- file.path(d, "merge.txt"); if (!file.exists(dfile)) dfile <- file.path(d, "test_output", "sim_tx_info_extended.txt")
    mg <- fread(dfile)[!is.na(promoter_id)]
    mg <- unique(mg[, .(id = paste0("PromoterID_", promoter_id), snp = paste0(gene_id, "_sQTL_b38"), effect_type, effect_size, true_mu_effect, true_nu_effect, true_tau_effect)], by = "id")[id %in% colnames(X)]
    psamp <- sub("^(sample_[0-9]+).*", "\\1", rownames(X))
    vcf <- fread(cmd = paste("zcat", file.path(d, "test_output/genotypes.vcf.gz"), "| grep -v '^##'"), header = TRUE)
    setnames(vcf, 1, "CHROM"); vs <- names(vcf)[10:ncol(vcf)]; ord <- match(psamp, vs); stopifnot(!anyNA(ord))
    gd <- function(snp) { r <- vcf[ID == snp]; if (!nrow(r)) return(rep(NA_real_, length(psamp)))
      as.numeric(vapply(strsplit(as.character(unlist(r[1, ..vs])), "[|/]"), function(x) sum(as.integer(x)), integer(1))[ord]) }
    Y <- as.matrix(X[, mg$id, drop = FALSE]); rownames(Y) <- psamp
    G <- sapply(mg$snp, gd); colnames(G) <- mg$id; rownames(G) <- psamp
    des <- mg[, .(id, effect_type, effect_size, true_mu_effect, true_nu_effect, true_tau_effect)]
  }
  list(Y = Y, G = G, des = des)
}
S <- load_sim(method, dataset, N)
ids <- colnames(S$Y)
if (mode == "calib") ids <- intersect(readLines(Sys.getenv("PERM_IDS")), ids)
if (NCH > 1) { ids <- sort(ids); ids <- ids[(seq_along(ids) - 1) %% NCH == (CH - 1)] }
cat(sprintf("mode=%s method=%s dataset=%s N=%d promoters=%d K=%d chunk=%d/%d\n", mode, method, dataset, N, length(ids), K, CH, NCH))

## ---- linear baselines ----
ctrl <- gamlss.control(n.cyc = 200, trace = FALSE)
lin_logitnormal <- function(yl, g) {
  dl <- data.frame(yl = yl, gt = g)
  m1 <- tryCatch(suppressWarnings(gamlss(yl ~ gt, sigma.formula = ~gt, family = NO, data = dl, control = ctrl)), error = function(e) NULL)
  m0 <- tryCatch(suppressWarnings(gamlss(yl ~ 1,  sigma.formula = ~gt, family = NO, data = dl, control = ctrl)), error = function(e) NULL)
  if (is.null(m1) || is.null(m0)) return(c(NA, NA, 0, NA))
  p <- pchisq(as.numeric(2 * (logLik(m1) - logLik(m0))), 1, lower.tail = FALSE)
  w <- tryCatch({ invisible(capture.output(s <- as.data.frame(summary(m1)))); s[grep("^gt", rownames(s))[1], "Pr(>|t|)"] }, error = function(e) NA_real_)
  b <- coef(m1, "mu")                                   # same as refit_real_pipeline.R:156 linear_prob_effect
  c(p, w, as.numeric(isTRUE(m1$converged) && isTRUE(m0$converged)), unname(plogis(b[1] + 2 * b[2]) - plogis(b[1])))
}
ols_beta <- function(y, g) { gc <- g - mean(g); sxx <- sum(gc^2); if (sxx <= 0) NA_real_ else 2 * sum(gc * y) / sxx }   # dosage 0 -> 2 change in mean activity
ols_p <- function(y, g) { n <- length(y); gc <- g - mean(g); sxx <- sum(gc^2); if (sxx <= 0) return(NA_real_)
  b <- sum(gc * y) / sxx; res <- y - mean(y) - b * gc; s2 <- sum(res^2) / (n - 2); if (s2 <= 0) return(NA_real_)
  2 * pt(-abs(b / sqrt(s2 / sxx)), n - 2) }

COLS <- c("mu_lrt_logp", "nu_lrt_logp", "tau_lrt_logp", "mu_wald_logp", "nu_wald_logp", "tau_wald_logp", "joint_logp",
          "mu_conv", "nu_conv", "tau_conv", "mu_sep", "nu_sep", "tau_sep", "cat_complete_sep",
          "mu_est", "nu_est", "tau_est", "mu_pe", "p0_pe", "p1_pe")       # effects (added 2026-10-01 for effect recovery)
FCOLS <- c("nu_firth_logp", "tau_firth_logp", "nu_firth_conv", "tau_firth_conv", "nu_firth_est", "tau_firth_est")
test_one <- function(y, C, g, nul, yl) {
  r <- tryCatch(snp_tests(y, C, g, nul), error = function(e) NULL)
  v <- if (is.null(r)) setNames(rep(NA_real_, length(COLS)), COLS) else vapply(COLS, function(k) { x <- r[[k]]; if (is.null(x) || length(x) != 1) NA_real_ else as.numeric(x) }, 0)
  fw <- setNames(rep(NA_real_, length(FCOLS)), FCOLS)
  if (nul$family != "BE") { f <- tryCatch(firth_snp_tests(y, C, g, nul), error = function(e) NULL)
    if (!is.null(f)) fw <- vapply(FCOLS, function(k) { x <- f[[k]]; if (is.null(x) || length(x) != 1) NA_real_ else as.numeric(x) }, 0) }
  ln <- lin_logitnormal(yl, g)
  c(v, fw, lin_mu_lrt_p = ln[1], lin_mu_wald_p = ln[2], lin_conv = ln[3], ols_p = ols_p(y, g), engine_error = as.numeric(is.null(r)),
    lin_prob_effect = ln[4], ols_effect = ols_beta(y, g))
}

one_promoter <- function(j) {
  id <- ids[j]; y <- S$Y[, id]; names(y) <- rownames(S$Y)
  hdr <- data.table(id = id, skip = NA_character_, family = NA_character_, n = NA_integer_, n0 = NA_integer_, n1 = NA_integer_, min_geno = NA_integer_,
                    elig_mu = NA, elig_p0 = NA, elig_p1 = NA)
  z <- prepare_y(y, 0L)
  if (!is.na(z$skip)) { hdr$skip <- z$skip; return(list(hdr = hdr)) }
  if (!any(z$y > 0 & z$y < 1)) { hdr$skip <- "no_interior"; return(list(hdr = hdr)) }
  s <- names(z$y); yy <- unname(z$y); g <- as.numeric(S$G[s, id])
  if (anyNA(g)) { hdr$skip <- "geno_na"; return(list(hdr = hdr)) }
  n <- length(yy); n0 <- sum(yy == 0); n1 <- sum(yy == 1); mg <- min(tabulate(g + 1L, 3L))
  hdr <- data.table(id = id, skip = NA_character_, family = z$family, n = n, n0 = n0, n1 = n1, min_geno = mg,
                    elig_mu = (n - n0 - n1) / n >= .70, elig_p0 = n0 / n >= .35 && z$family %in% c("BEINF", "BEZI"),
                    elig_p1 = n1 / n >= .35 && z$family %in% c("BEINF", "BEOI"))
  if (mg < MIN_GENO) { hdr$skip <- "min_geno_lt5"; return(list(hdr = hdr)) }
  C <- orth_design(matrix(1, n, 1, dimnames = list(s, "(Intercept)")))
  nul <- tryCatch(promoter_null(yy, C, z$family), error = function(e) NULL)
  if (is.null(nul)) { hdr$skip <- "null_error"; return(list(hdr = hdr)) }
  yl <- yy; yl[yl == 0] <- .001; yl[yl == 1] <- .999; yl <- log(yl / (1 - yl))
  if (mode == "power") {
    return(list(hdr = hdr, res = as.data.table(as.list(test_one(yy, C, g, nul, yl)))[, id := id]))
  }
  set.seed(PSEED + match(id, colnames(S$Y)))                # PSEED = 20261001 unless PERM_SEED_BASE is set (seed check)
  P <- replicate(K, sample.int(n))                         # all permutations fixed before any fitting
  R <- t(vapply(seq_len(K), function(k) test_one(yy, C, g[P[, k]], nul, yl), numeric(length(COLS) + length(FCOLS) + 7)))
  list(hdr = hdr, res = data.table(id = id, perm = seq_len(K), R))
}
t0 <- proc.time()[3]
out <- mclapply(seq_along(ids), function(j) tryCatch(one_promoter(j), error = function(e) list(hdr = data.table(id = ids[j], skip = paste("ERROR:", conditionMessage(e))))),
                mc.cores = cores, mc.preschedule = FALSE)
hdr <- rbindlist(lapply(out, `[[`, "hdr"), fill = TRUE); res <- rbindlist(lapply(out, `[[`, "res"), fill = TRUE)
hdr <- merge(hdr, S$des, by = "id", all.x = TRUE)
dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(meta = list(mode = mode, method = method, dataset = dataset, N = N, K = K, chunk = CH, n_chunks = NCH, date = Sys.time(),
                         engine_md5 = tools::md5sum(file.path(NM, "code_copy", c("bohr2_engine.R", "bohr2_addon3.R")))),
             header = hdr, res = res), out_rds)
cat(sprintf("done in %.1f min | promoters ok %d / %d | rows %d | engine errors %d -> %s\n", (proc.time()[3] - t0) / 60,
            hdr[is.na(skip), .N], nrow(hdr), nrow(res), if (nrow(res)) sum(res$engine_error, na.rm = TRUE) else 0L, out_rds))
print(hdr[, .N, by = skip])
