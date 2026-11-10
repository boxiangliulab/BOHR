## Synthetic donors for the trait regression, and the reference trait analysis built from the original code: data.R
## (cov_design, trait_contrast, constants; function definitions only), run_trait.R one_prom (transcribed), mirror.R (local pair
## rule, transcribed), summarize.R to_units / final rules (transcribed) and 01_abs_labels.R / 03_apply_B.R (transcribed).

trait_code_dir <- function() {
  d <- Sys.getenv("BOHR_TRAIT_DIR", "/projects/group/bliu/e1327989/GTEx/tissue_regression")
  f <- file.path(d, c("trait_regression_20261007/code/run_trait.R", "trait_regression_20261007/code/data.R", "trait_regression_20261104/code/summarize.R",
                      "trait_regression_20261104/code/mirror_v3.R", "trait_regression_20261005/code/mirror.R",
                      "trait_regression_20261106_mirror_abs/code/01_abs_labels.R", "trait_regression_20261106_mirror_abs/code/03_apply_B.R"))
  if (nzchar(d) && all(file.exists(f))) d else NULL
}
skip_if_no_trait_code <- function() { skip_if_no_engine(); if (is.null(trait_code_dir())) testthat::skip("original trait regression code not available") }
load_trait_code <- function() {
  env <- load_engine()
  source_defs(file.path(trait_code_dir(), "trait_regression_20261007/code/data.R"), env, funs = c("cov_design", "trait_contrast", "prepare_y", "orth_design"),
              consts = c("TRAITS", "CONT_TRAITS", "MH_BINARY", "COVS", "MIN_CLASS", "MIN_HARDY"))
  env
}
expect_trait_line <- function(file, text) testthat::expect_true(any(grepl(text, readLines(file.path(trait_code_dir(), file), warn = FALSE), fixed = TRUE)),
                                                               info = paste0(file, ": ", text))

## donors with the paper's covariates (SEX, AGE, PC1-5, DTHHRDY with rare levels, MHABNWBC with missing values) and traits; features:
## an exact mirror pair (P301/P302, p0 of P301 driven by SEX, only P302's absolute activity responds), a boundary-only pair
## (P303/P304), BE / BEZI / BEOI / BEINF features with SEX, AGE and BMI effects, a feature rare in the binary trait, missing values.
sim_trait <- function(seed = 41, n = 180) {
  with_seed(seed, {
    d <- paste0("G", seq_len(n))
    meta <- data.frame(SUBJID = d, SEX = stats::rbinom(n, 1, 0.45), AGE = round(stats::runif(n, 20, 70)), BMI = stats::rnorm(n, 27, 4),
                       PC1 = stats::rnorm(n), PC2 = stats::rnorm(n), PC3 = stats::rnorm(n), PC4 = stats::rnorm(n), PC5 = stats::rnorm(n),
                       DTHHRDY = sample(c(0, 1, 2, 3, 4), n, replace = TRUE, prob = c(0.45, 0.04, 0.25, 0.2, 0.06)),
                       MHABNWBC = stats::rbinom(n, 1, 0.1), MHT2D = stats::rbinom(n, 1, 0.12), stringsAsFactors = FALSE)
    meta$MHABNWBC[sample.int(n, 6)] <- NA; meta$BMI[sample.int(n, 4)] <- NA; meta$DTHHRDY[sample.int(n, 3)] <- NA
    rownames(meta) <- d
    by <- function(eta, phi = 9) stats::rbeta(n, stats::plogis(eta) * phi, (1 - stats::plogis(eta)) * phi)
    age <- (meta$AGE - 45) / 15; bmi <- (ifelse(is.na(meta$BMI), 27, meta$BMI) - 27) / 4
    Y <- list()
    y <- by(0.3); u <- stats::runif(n); y[u < stats::plogis(-2 + 1.6 * meta$SEX)] <- 0; y[u > 0.9] <- 1
    Y$PromoterID_301 <- y; Y$PromoterID_302 <- 1 - y                                                    # exact mirror pair
    y <- by(-0.2 + 0.3 * age); u <- stats::runif(n); y[u < 0.15] <- 0; y[u > 0.8] <- 1
    Y$PromoterID_303 <- y; Y$PromoterID_304 <- ifelse(y %in% c(0, 1), 1 - y, 0.85 * (1 - y))              # boundary-only pair
    Y$PromoterID_305 <- by(0.1 + 0.35 * age + 0.2 * meta$SEX)                                            # BE, AGE mu
    y <- by(0); y[stats::runif(n) < stats::plogis(-1 + 0.6 * bmi)] <- 0; Y$PromoterID_306 <- y            # BEZI, BMI p0
    y <- by(0.4); y[stats::runif(n) < stats::plogis(-1.2 - 0.9 * meta$MHT2D)] <- 1; Y$PromoterID_307 <- y # BEOI
    y <- by(0.2 + 0.4 * meta$SEX); u <- stats::runif(n); y[u < 0.2] <- 0; y[u > 0.85] <- 1; y[c(2, 7, 30)] <- NA
    Y$PromoterID_308 <- y                                                                                # BEINF, SEX mu, missing values
    y <- by(0.1); y[meta$MHT2D == 1 & stats::runif(n) < 0.7] <- 0; y[stats::runif(n) < 0.05] <- 0; Y$PromoterID_309 <- y
    Y$PromoterID_310 <- by(0.5 * meta$SEX - 0.2)                                                         # single-promoter gene
    ph <- do.call(rbind, Y); colnames(ph) <- d
    pr <- data.frame(id = rownames(ph), gene = c("H1", "H1", "H2", "H2", "H3", "H4", "H5", "H6", "H7", "H3"), stringsAsFactors = FALSE)
    lib <- stats::rnorm(n, 8, 0.3)
    A <- rbind(PromoterID_301 = lib + stats::rnorm(n, 0, 0.3), PromoterID_302 = lib + 0.8 * meta$SEX + stats::rnorm(n, 0, 0.3),
               PromoterID_303 = lib + stats::rnorm(n, 0, 0.3), PromoterID_304 = lib + stats::rnorm(n, 0, 0.3))
    colnames(A) <- d
    list(phenotype = ph, meta = meta, promoters = pr, absolute = A)
  })
}
## the package's inputs from sim_trait(): covariates in the paper's column order, the Hardy scale as a categorical column
trait_inputs <- function(x) {
  m <- x$meta
  list(covariates = data.frame(SEX = m$SEX, AGE = m$AGE, PC1 = m$PC1, PC2 = m$PC2, PC3 = m$PC3, PC4 = m$PC4, PC5 = m$PC5,
                               DTHHRDY = as.character(m$DTHHRDY), MHABNWBC = m$MHABNWBC, row.names = rownames(m), stringsAsFactors = FALSE),
       traits = m[, c("SEX", "AGE", "BMI", "MHT2D")])
}

## Reference for one tissue: run_trait.R (one_trait / one_prom, real arm) per trait, then summarize.R units with mirror rule v3 and the
## final rules, then (with A) 01_abs_labels.R + 03_apply_B.R. Mirror flags: mirror.R's local rule on the tissue matrix, canonical by
## the pooled mean (rowMeans of the phenotype).
ref_trait <- function(env, x, traits, A = NULL) {
  meta <- x$meta; td <- list(U = colnames(x$phenotype), act = t(x$phenotype)); prs <- colnames(td$act)
  info <- data.frame(geneId = x$promoters$gene, row.names = x$promoters$id)
  MIN_CLASS <- env$MIN_CLASS; CONT_TRAITS <- env$CONT_TRAITS; COVS <- env$COVS
  ## mirror.R:35-58 (local rule), pooled mean = rowMeans
  pooled <- rowMeans(x$phenotype, na.rm = TRUE)
  inf <- data.frame(promoter = prs, gene = info[prs, "geneId"], stringsAsFactors = FALSE); np <- table(inf$gene)
  MR <- NULL
  for (g in names(np)[np == 2]) {
    pr <- inf$promoter[inf$gene == g]; Aa <- td$act[, pr[1]]; B <- td$act[, pr[2]]; ok <- !is.na(Aa) & !is.na(B)
    a <- Aa[ok]; b <- B[ok]
    bex <- length(a) > 0 && all((a == 0) == (b == 1)) && all((a == 1) == (b == 0))
    im <- a > 0 & a < 1
    iex <- bex && all((b > 0 & b < 1) == im) && (!any(im) || max(abs(a[im] + b[im] - 1)) <= 1e-6)
    can <- pr[which.max(pooled[pr])]; par <- setdiff(pr, can)
    MR <- rbind(MR, data.frame(promoter = c(can, par), boundary_exact = bex, interior_exact = iex, mirror_role = if (bex) c("canonical", "partner") else NA,
                               mirror_partner = if (bex) c(par, can) else NA, stringsAsFactors = FALSE))
  }
  MR$drop_mu <- MR$interior_exact & MR$mirror_role %in% "partner"; MR$drop_p1 <- MR$boundary_exact %in% TRUE    # mirror.R:56, mirror_v3.R:18
  out <- list(); labs <- list()
  for (trait in traits) {
    ## run_trait.R:35-47, 58-59
    cv <- setdiff(COVS, trait)
    d <- meta[td$U, c(trait, cv), drop = FALSE]
    donors <- td$U[complete.cases(d)]
    g <- meta[donors, trait]
    if (length(donors) < 20 || sd(g) == 0) next
    else if (!trait %in% CONT_TRAITS && min(table(g)) < MIN_CLASS) next
    gp <- setNames(env$trait_contrast(g, trait), donors)
    C_all <- env$cov_design(meta, donors, trait)
    rows <- lapply(prs, function(p) {                                                            # run_trait.R:62-116 (one_prom, real arm)
      base <- list(promoter = p, family = NA_character_, skip = NA_character_, n = NA, n0 = NA, n1 = NA, n_mid = NA, mu_testable = NA)
      yraw <- td$act[donors, p]
      z <- env$prepare_y(yraw, ncol(C_all) - 1L)
      if (!is.na(z$skip)) { base$skip <- z$skip; return(base) }
      if (!any(z$y > 0 & z$y < 1)) { base$skip <- "no_interior"; return(base) }
      s <- names(z$y); y <- unname(z$y); gi <- unname(gp[s])
      if (!trait %in% CONT_TRAITS) { tb <- table(gi); if (length(tb) < 2 || min(tb) < MIN_CLASS) { base$skip <- "trait_rare_in_promoter"; return(base) } }
      Cs <- env$orth_design(C_all[s, , drop = FALSE]); mid <- y > 0 & y < 1
      if (!trait %in% CONT_TRAITS) { tm <- table(factor(gi[mid], levels = sort(unique(gi))))
        base$mu_testable <- min(tm) >= MIN_CLASS } else base$mu_testable <- TRUE
      base$family <- z$family; base$n <- length(y); base$n0 <- sum(y == 0); base$n1 <- sum(y == 1); base$n_mid <- sum(mid)
      nul <- env$promoter_null(y, Cs, z$family)
      r <- env$snp_tests(y, Cs, gi, nul)
      rk <- if (sum(mid) > 0) qr(Cs[mid, , drop = FALSE])$rank else NA_integer_
      pm <- rk + 1L; W <- r$mu_stat; nm <- sum(mid)
      mu_F_logp <- if (is.finite(W) && nm > pm + 1L) pf((nm - pm) * expm1(max(W, 0) / nm), 1, nm - pm, lower.tail = FALSE, log.p = TRUE) else NA_real_
      fr <- if (z$family != "BE") tryCatch(env$firth_snp_tests(y, Cs, gi, nul), error = function(e) list(firth_error = conditionMessage(e))) else list()
      c(base, list(mu_F_logp = mu_F_logp, mu_est = 2 * r$mu_est, mu_se = 2 * r$mu_se, mu_pe = r$mu_pe, p0_pe = r$p0_pe, p1_pe = r$p1_pe,
                   cat_complete_sep = r$cat_complete_sep,
                   nu_firth_logp = if (is.null(fr$nu_firth_logp)) NA else fr$nu_firth_logp, tau_firth_logp = if (is.null(fr$tau_firth_logp)) NA else fr$tau_firth_logp,
                   nu_firth_est = if (is.null(fr$nu_firth_est)) NA else 2 * fr$nu_firth_est, nu_firth_se = if (is.null(fr$nu_firth_se)) NA else 2 * fr$nu_firth_se,
                   tau_firth_est = if (is.null(fr$tau_firth_est)) NA else 2 * fr$tau_firth_est, tau_firth_se = if (is.null(fr$tau_firth_se)) NA else 2 * fr$tau_firth_se))
    })
    R <- BOHR:::rows_to_df(lapply(rows, function(r) { for (k in c("mu_F_logp", "mu_est", "mu_se", "mu_pe", "p0_pe", "p1_pe", "cat_complete_sep", "nu_firth_logp",
                                                                   "tau_firth_logp", "nu_firth_est", "nu_firth_se", "tau_firth_est", "tau_firth_se"))
      if (is.null(r[[k]])) r[[k]] <- NA; r }))
    R$n_donors <- length(donors); R <- R[is.na(R$skip), , drop = FALSE]
    ## summarize.R:51-66 (to_units), 67-98 (mirror rule v3, final rules, BH)
    mk <- function(type, keep, testable, logp, est, se, pe) data.frame(promoter = R$promoter[keep], type = type, testable = testable[keep], logp = logp[keep],
      est = est[keep], se = se[keep], pe = pe[keep], family = R$family[keep], n = R$n[keep], n0 = R$n0[keep], n1 = R$n1[keep], n_mid = R$n_mid[keep],
      n_donors = R$n_donors[keep], cat_complete_sep = R$cat_complete_sep[keep], stringsAsFactors = FALSE)
    U <- rbind(mk("mu", rep(TRUE, nrow(R)), R$mu_testable %in% TRUE, R$mu_F_logp, R$mu_est, R$mu_se, R$mu_pe),
               mk("p0", R$family %in% c("BEINF", "BEZI"), rep(TRUE, nrow(R)), R$nu_firth_logp, R$nu_firth_est, R$nu_firth_se, R$p0_pe),
               mk("p1", R$family %in% c("BEINF", "BEOI"), rep(TRUE, nrow(R)), R$tau_firth_logp, R$tau_firth_est, R$tau_firth_se, R$p1_pe))
    U$testable <- U$testable & is.finite(U$logp)
    U$p <- ifelse(U$testable, exp(U$logp), NA_real_)
    i <- match(U$promoter, MR$promoter)
    for (v in c("boundary_exact", "interior_exact", "drop_mu", "drop_p1")) U[[v]] <- MR[[v]][i] %in% TRUE
    U$mirror_partner <- MR$mirror_partner[i]
    U$dedup_drop <- ifelse(U$type == "mu", U$drop_mu, ifelse(U$type == "p1", U$drop_p1, FALSE))
    U$primary <- U$testable & !U$dedup_drop
    U$share <- ifelse(U$type == "p0", U$n0 / U$n, ifelse(U$type == "p1", U$n1 / U$n, NA_real_))
    U$share_10_90 <- U$type == "mu" | (U$share >= 0.10 & U$share <= 0.90)
    U$estimable_mu <- !(U$type == "mu" & U$n_mid < 80)
    U$assessable <- !(U$type != "mu" & U$cat_complete_sep %in% TRUE)
    U$report_tissue <- U$n_donors >= 50
    U$final <- U$primary & U$estimable_mu & U$share_10_90 & U$assessable & U$report_tissue
    U$q_final <- NA_real_
    for (ty in c("mu", "p0", "p1")) { k <- U$final %in% TRUE & U$type == ty; U$q_final[k] <- p.adjust(U$p[k], "BH") }
    U$sig_final <- U$final %in% TRUE & U$q_final <= 0.05
    U$bh_family <- ifelse(U$final %in% TRUE, U$type, NA_character_); U$abs_relabel <- FALSE
    U <- U[U$testable, , drop = FALSE]
    if (!is.null(A)) {                                                                              # 01_abs_labels.R + 03_apply_B.R
      E <- U[U$type == "p0" & U$boundary_exact & U$sig_final, , drop = FALSE]
      if (nrow(E)) {
        X <- cbind(C_all, g = env$trait_contrast(meta[donors, trait], trait))
        ols_z <- function(X, Y) { q <- qr(X); p <- ncol(X); stopifnot(q$rank == p)
          B <- qr.coef(q, Y); R <- qr.resid(q, Y)
          s2 <- colSums(R^2) / (nrow(X) - p)
          v <- chol2inv(qr.R(q))[p, p]
          b <- B[p, ]; se <- sqrt(v * s2)
          list(b = b, se = se, z = b / se, sdy = apply(Y, 2, sd)) }
        fS <- ols_z(X, t(A[E$promoter, donors, drop = FALSE])); fO <- ols_z(X, t(A[E$mirror_partner, donors, drop = FALSE]))
        rs <- abs(fS$z) >= 1.96; ro <- abs(fO$z) >= 1.96
        outc <- ifelse(rs & ro, "both", ifelse(rs, "S_only", ifelse(ro, "O_only", "neither")))
        labs[[trait]] <- data.frame(trait = trait, S = E$promoter, O = E$mirror_partner, z_S = unname(fS$z), z_O = unname(fO$z), outcome = outc, stringsAsFactors = FALSE)
        for (j in which(outc == "O_only")) {
          iS <- which(U$promoter == E$promoter[j] & U$type == "p0"); iO <- which(U$promoter == E$mirror_partner[j] & U$type == "p1")
          for (k in c("final", "primary", "q_final", "sig_final")) U[[k]][iO] <- U[[k]][iS]
          U$dedup_drop[iO] <- FALSE; U$bh_family[iO] <- "p0"; U$abs_relabel[iO] <- TRUE
          U$final[iS] <- FALSE; U$primary[iS] <- FALSE; U$q_final[iS] <- NA; U$sig_final[iS] <- FALSE
          U$dedup_drop[iS] <- TRUE; U$bh_family[iS] <- NA; U$abs_relabel[iS] <- TRUE
        }
      }
    }
    U$trait <- trait
    out[[trait]] <- U
  }
  U <- do.call(rbind, out); U <- U[order(match(U$trait, traits), U$promoter, U$type, method = "radix"), , drop = FALSE]; rownames(U) <- NULL
  list(units = U, labels = do.call(rbind, labs))
}
## lines of the original scripts that ref_trait() transcribes
ref_trait_lines <- list(
  "trait_regression_20261007/code/run_trait.R" = c(
    "donors <- td$U[complete.cases(d)]", "if (length(donors) < 20 || sd(g) == 0) skip_trait <- \"trait_constant_or_n\"",
    "else if (!trait %in% CONT_TRAITS && min(table(g)) < MIN_CLASS) skip_trait <- \"trait_rare\"",
    "if (!trait %in% CONT_TRAITS) { tb <- table(gi); if (length(tb) < 2 || min(tb) < MIN_CLASS) { base$skip <- \"trait_rare_in_promoter\"; return(base) }",
    "if (!trait %in% CONT_TRAITS) { tm <- table(factor(gi[mid], levels = sort(unique(gi))))",
    "mu_F_logp <- if (is.finite(W) && nm > pm + 1L) pf((nm - pm) * expm1(max(W, 0) / nm), 1, nm - pm, lower.tail = FALSE, log.p = TRUE) else NA_real_",
    "mu_est = 2 * r$mu_est, mu_se = 2 * r$mu_se, nu_est = 2 * r$nu_est, nu_se = 2 * r$nu_se,",
    "o[[k]] <- if (!is.null(fr[[k]])) 2 * fr[[k]] else NA_real_"),
  "trait_regression_20261104/code/summarize.R" = c(
    "U[, testable := testable & is.finite(logp)]",
    "U[, dedup_drop := fifelse(type == \"mu\", drop_mu, fifelse(type == \"p1\", drop_p1, drop_p0))]",
    "U[, share_10_90 := type == \"mu\" | (share >= 0.10 & share <= 0.90)]",
    "U[, estimable_mu := !(type == \"mu\" & n_mid < 80)]",
    "U[, assessable := !(type != \"mu\" & cat_complete_sep %in% TRUE)]",
    "U[, report_tissue := n_donors >= 50]",
    "U[, final := primary & estimable_mu & share_10_90 & assessable & report_tissue]",
    "U[final == TRUE, q_final := p.adjust(p, \"BH\"), by = .(tissue, trait, type)]"),
  "trait_regression_20261104/code/mirror_v3.R" = "M[, `:=`(drop_p0 = FALSE, drop_p1 = boundary_exact %in% TRUE)]",
  "trait_regression_20261005/code/mirror.R" = c(
    "bex <- length(a) > 0 && all((a == 0) == (b == 1)) && all((a == 1) == (b == 0))",
    "iex <- bex && all((b > 0 & b < 1) == im) && (!any(im) || max(abs(a[im] + b[im] - 1)) <= 1e-6)",
    "can <- pr[which.max(pooled[pr])]; par <- setdiff(pr, can)"),
  "trait_regression_20261106_mirror_abs/code/01_abs_labels.R" = c(
    "X <- cbind(cov_design(meta, donors, tr), g = trait_contrast(meta[donors, tr], tr))",
    "B <- qr.coef(q, Y); R <- qr.resid(q, Y)", "v <- chol2inv(qr.R(q))[p, p]"),
  "trait_regression_20261106_mirror_abs/code/03_apply_B.R" = c(
    "set(U, iO, c(\"final\", \"primary\", \"q_final\", \"sig_final\", \"q\", \"sig\"), keepS[, .(final, primary, q_final, sig_final, q, sig)])",
    "set(U, iS, c(\"final\", \"primary\", \"q_final\", \"sig_final\", \"q\", \"sig\"), list(FALSE, FALSE, NA_real_, FALSE, NA_real_, FALSE))"))
