## =============================================================================
## run_trait.R -- regression of the 20 traits for one tissue and one arm (new engine)
## Usage: Rscript run_trait.R <tissue> <arm> <ncores> [trait ...]
##   arm = real  : real data
##         sim1  : parametric null -- for each (promoter, trait), the null model fitted to the real data (covariates only)
##                 is used to re-simulate the activity; trait and covariates keep their real values, then the identical pipeline is run
##         perm1 : permutation null -- the trait is randomly permuted within strata (stratum = SEX x age decade, omitting the variable being tested),
##                 activity and covariates are left unchanged
## For each (promoter, trait):
##   prepare_y (skip rules + <=5% endpoint removal) -> covariate orthogonalisation -> promoter_null
##   -> snp_tests (g = trait contrast coding) -> F correction of mu -> firth_snp_tests (when there are 0/1 parts)
## Output (atomic write): results/<arm>/<tissue>/<trait>.rds = list(res = data.table, meta = list)
## =============================================================================
suppressPackageStartupMessages(library(parallel))
CODE <- file.path(Sys.getenv("TRG_RUN", "/projects/group/bliu/e1327989/GTEx/tissue_regression/trait_regression_20261007"), "code")   # On NSCC this is set by TRG_RUN
source(file.path(CODE, "data.R")); source(file.path(CODE, "bohr2_engine.R")); source(file.path(CODE, "bohr2_addon3.R"))

a <- commandArgs(TRUE); tis <- a[1]; arm <- a[2]; NC <- as.integer(a[3])
traits <- if (length(a) > 3) a[-(1:3)] else TRAITS
stopifnot(arm %in% c("real", "sim1", "perm1"), all(traits %in% TRAITS))
td   <- readRDS(file.path(PATHS$work, paste0(tis, ".rds")))
meta <- readRDS(file.path(PATHS$work, "meta.rds"))
tis_idx <- match(tis, readLines(PATHS$tissues))
MAXPROM <- as.integer(Sys.getenv("TRG_MAXPROM", "0"))                 # Smoke test: run only the first N promoters
OUTROOT <- Sys.getenv("TRG_OUT", PATHS$res)
if (MAXPROM > 0 && OUTROOT == PATHS$res) stop("TRG_MAXPROM requires TRG_OUT (never into the production tree)")

age_bin <- function(a) cut(a, c(19, 29, 39, 49, 59, 69, 79))
na_res <- function() list(mu_stat = NA_real_, mu_est = NA_real_, mu_se = NA_real_)

one_trait <- function(trait) {
  out <- file.path(OUTROOT, arm, tis, paste0(trait, ".rds"))
  if (file.exists(out)) return(data.table(trait = trait, status = "exists"))
  t0 <- Sys.time()
  cv <- setdiff(COVS, trait)
  d  <- meta[td$U, c(trait, cv), drop = FALSE]
  donors <- td$U[complete.cases(d)]
  g <- meta[donors, trait]
  skip_trait <- NA_character_
  if (length(donors) < 20 || sd(g) == 0) skip_trait <- "trait_constant_or_n"
  else if (!trait %in% CONT_TRAITS && min(table(g)) < MIN_CLASS) skip_trait <- "trait_rare"
  if (!is.na(skip_trait)) {
    dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
    saveRDS(list(res = data.table(), meta = list(tissue = tis, arm = arm, trait = trait, n_donors = length(donors),
                                                 skip_trait = skip_trait)), out)
    return(data.table(trait = trait, status = skip_trait))
  }
  ## perm1: stratified permutation of the trait (fixed seed, reproducible; the mapping is saved too)
  perm_seed <- 20261001L + 1000L * tis_idx + match(trait, TRAITS)
  if (arm == "perm1") {
    st <- interaction(if (trait != "SEX") meta[donors, "SEX"] else 1,
                      if (trait != "AGE") age_bin(meta[donors, "AGE"]) else 1, drop = TRUE)
    set.seed(perm_seed)
    pi_ <- seq_along(donors)                                            # Position i receives the trait of donor pi_[i]
    for (s in split(seq_along(donors), st)) if (length(s) > 1) pi_[s] <- s[sample.int(length(s))]
    g <- g[pi_]
  }
  gp <- setNames(trait_contrast(g, trait), donors)
  C_all <- cov_design(meta, donors, trait)
  prs <- colnames(td$act); if (MAXPROM > 0) prs <- head(prs, MAXPROM)

  one_prom <- function(i) {
    p <- prs[i]; tp <- proc.time()[["elapsed"]]
    base <- list(promoter = p, gene = td$info[p, "geneId"], family = NA_character_, skip = NA_character_,
                 n = NA_integer_, n0 = NA_integer_, n1 = NA_integer_, n_mid = NA_integer_,
                 p0_raw = NA_real_, p1_raw = NA_real_, n_minor = NA_integer_, n_minor_mid = NA_integer_,
                 mu_testable = NA, n_cov = NA_integer_, rank_mid = NA_integer_)
    yraw <- td$act[donors, p]; yr <- yraw[!is.na(yraw)]
    base$p0_raw <- mean(yr == 0); base$p1_raw <- mean(yr == 1)
    z <- prepare_y(yraw, ncol(C_all) - 1L)
    if (!is.na(z$skip)) { base$skip <- z$skip; return(base) }
    if (arm == "sim1") {                                       # Fit the null model on the real y first, then simulate from it
      Cs0 <- orth_design(C_all[names(z$y), , drop = FALSE])
      nul0 <- promoter_null(unname(z$y), Cs0, z$family)
      ## The seed is set by promoter ID (not by matrix position): when the promoter set changes, the simulated data of the same promoter do not change
      set.seed(as.integer(((perm_seed %% 100000L) * 20011 + as.integer(sub("PromoterID_", "", p))) %% .Machine$integer.max))
      ys <- simulate_null(Cs0, nul0, y_obs = unname(z$y))
      if (is.null(ys)) { base$skip <- "sim_failed"; return(base) }
      z <- prepare_y(setNames(ys, names(z$y)), ncol(C_all) - 1L)
      if (!is.na(z$skip)) { base$skip <- paste0("sim_", z$skip); return(base) }
    }
    if (!any(z$y > 0 & z$y < 1)) { base$skip <- "no_interior"; return(base) }
    s <- names(z$y); y <- unname(z$y); gi <- unname(gp[s])
    if (!trait %in% CONT_TRAITS) { tb <- table(gi); if (length(tb) < 2 || min(tb) < MIN_CLASS) { base$skip <- "trait_rare_in_promoter"; return(base) }
      base$n_minor <- min(tb) }
    Cs <- orth_design(C_all[s, , drop = FALSE]); mid <- y > 0 & y < 1
    ## mu uses only the intermediate-value samples: a binary trait needs >= MIN_CLASS per class among them (otherwise leverage is extreme and the F correction does not apply)
    if (!trait %in% CONT_TRAITS) { tm <- table(factor(gi[mid], levels = sort(unique(gi))))
      base$n_minor_mid <- min(tm); base$mu_testable <- min(tm) >= MIN_CLASS } else base$mu_testable <- TRUE
    base$family <- z$family; base$n <- length(y); base$n0 <- sum(y == 0); base$n1 <- sum(y == 1); base$n_mid <- sum(mid)
    base$n_cov <- ncol(Cs) - 1L
    nul <- promoter_null(y, Cs, z$family)
    r <- tryCatch(snp_tests(y, Cs, gi, nul), error = function(e) { base$skip <<- paste("engine_error:", conditionMessage(e)); NULL })
    if (is.null(r)) return(base)
    ## Small-sample F correction of mu (same as rerun_v2/code/run_muF.R: p = rank of the covariate design on the intermediate-value samples + 1)
    rk <- if (sum(mid) > 0) qr(Cs[mid, , drop = FALSE])$rank else NA_integer_
    base$rank_mid <- rk; pm <- rk + 1L; W <- r$mu_stat; nm <- sum(mid)
    mu_F_logp <- if (is.finite(W) && nm > pm + 1L) pf((nm - pm) * expm1(max(W, 0) / nm), 1, nm - pm, lower.tail = FALSE, log.p = TRUE) else NA_real_
    ## Firth is computed as usual under complete separation of the categorical part (cat_complete_sep): the penalised likelihood remains finite
    fr <- if (z$family != "BE")
      tryCatch(firth_snp_tests(y, Cs, gi, nul), error = function(e) list(firth_error = conditionMessage(e))) else list()
    keep <- c("mu_stat", "mu_lrt_logp", "mu_wald_logp", "mu_conv", "mu_sep", "mu_pe",
              "nu_stat", "nu_lrt_logp", "nu_wald_logp", "nu_conv", "nu_sep",
              "tau_stat", "tau_lrt_logp", "tau_wald_logp", "tau_conv", "tau_sep", "p0_pe", "p1_pe", "cat_complete_sep")
    o <- c(base, r[keep], list(mu_F_logp = mu_F_logp,
           ## Estimate and SE are converted back to the link scale of "per contrast" (+1 SD / unexposed->exposed)
           mu_est = 2 * r$mu_est, mu_se = 2 * r$mu_se, nu_est = 2 * r$nu_est, nu_se = 2 * r$nu_se,
           tau_est = 2 * r$tau_est, tau_se = 2 * r$tau_se))
    for (k in c("nu_firth_stat", "nu_firth_logp", "nu_firth_conv", "tau_firth_stat", "tau_firth_logp", "tau_firth_conv", "firth_rounds"))
      o[[k]] <- if (!is.null(fr[[k]])) fr[[k]] else NA
    for (k in c("nu_firth_est", "nu_firth_se", "tau_firth_est", "tau_firth_se"))
      o[[k]] <- if (!is.null(fr[[k]])) 2 * fr[[k]] else NA_real_
    o$firth_error <- if (!is.null(fr$firth_error)) fr$firth_error else NA_character_
    o$sec <- proc.time()[["elapsed"]] - tp
    o
  }
  L <- mclapply(seq_along(prs), function(i) tryCatch(one_prom(i), error = function(e)
         list(promoter = prs[i], skip = paste("error:", conditionMessage(e)))), mc.cores = NC, mc.preschedule = FALSE)
  ## When an mclapply worker dies it returns a try-error: retry once
  bad <- which(vapply(L, function(x) inherits(x, "try-error") || is.null(x), TRUE))
  if (length(bad)) L[bad] <- lapply(bad, function(i) tryCatch(one_prom(i), error = function(e) list(promoter = prs[i], skip = paste("error:", conditionMessage(e)))))
  R <- rbindlist(L, fill = TRUE)
  R[, `:=`(tissue = tis, trait = trait)]
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(out, ".tmp")
  saveRDS(list(res = R, meta = list(tissue = tis, arm = arm, trait = trait, n_donors = length(donors),
                                    covariates = colnames(C_all), skip_trait = NA_character_,
                                    perm_seed = if (arm != "real") perm_seed else NA, date = format(Sys.time()),
                                    minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))), tmp)
  file.rename(tmp, out)
  data.table(trait = trait, status = "done", n_prom = nrow(R), tested = sum(is.na(R$skip)),
             minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
}

S <- rbindlist(lapply(traits, function(tr) { x <- one_trait(tr); print(x); x }), fill = TRUE)
print(S)
