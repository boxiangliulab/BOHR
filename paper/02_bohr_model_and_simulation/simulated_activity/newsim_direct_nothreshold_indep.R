### INDEPENDENT-GENOTYPE copy of m1_nothr/newsim_direct_nothreshold.R: no-threshold baselines,
### interior floor TRUNC_L (L = 0.005), and each gene with its own genotype vector (geno_mat); run with n_genes = 2000.
### NO-THRESHOLD copy of m1_trunc/newsim_direct_realfloor_trunc.R: baseline p0/p1 no longer depend on the
### effect type (see the block marked NO-THRESHOLD); the interior measurement floor (TRUNC_L, adopted L = 0.005) is kept.
### FLOOR copy of method1_direct/newsim_direct_realfloor.R (newmethod_bohr2/m1_trunc).
### ONLY change: interior values (0<y<1) are put on a measurement floor [L, 1-L], L = env TRUNC_L (default 0 = identical
### to the original). Rationale: real relative promoter activity = A_p / sum A over the gene's promoters with
### A = log2(DESeq2-normalised junction count + 1) (proActiv), so an active promoter cannot have y below
### log2(1+1/sf) / (log2(1+1/sf) + (k-1) log2(1+Cmax/sf)) ~ 1e-3..1e-2 (newmethod_bohr2/README.md).
### TRUNC_MODE = "clamp" (default): every BETA-PART draw y <- min(max(y, L), 1-L) -- including beta draws that qbeta
###   rounded to exactly 0/1, identified from the uniform (rBEINF == qBEINF(runif)) -- so each value moves by at most L and
###   each genotype-group interior mean by at most L (the dosage-0 vs 2 contrast by at most 2L): truth approximately preserved. "redraw": exact truncated beta by quantile mapping of the same draw -- REJECTED:
###   when a genotype group's beta has most of its mass beyond 1-L (mu near 0.99) the redraw moves that group's mean by
###   up to ~0.4 and so changes the simulated effect (gene_0029 at N = 50). No extra random numbers in either mode.
#!/usr/bin/env Rscript
### New Simulation: Direct rBEINF
### Based on v11_direct_sim.R with TWO changes only:
###   1. Data generation uses rBEINF instead of manual mixing (rep(0)+rep(1)+rnorm)
###   2. Effect boundary uses proportional scaling instead of hard clamping (pmin/pmax)
###
### Everything else (effect types, proportions, parameters, model fitting, etc.)
### is IDENTICAL to v11.
###
### Usage: Rscript newsim_direct.R <n_samples> [n_genes] [n_cores] [gene_seed] [data_seed] [effect_scale]
###   effect_scale: "prob" (default) or "link"
###     prob = effects added directly on probability scale
###     link = effects added on link scale (logit for mu, log for nu/tau)
### =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) stop("Usage: Rscript newsim_direct.R <n_samples> [n_genes] [n_cores] [gene_seed] [data_seed] [effect_scale]")

n_samples    <- as.integer(args[1])
n_genes      <- ifelse(length(args) >= 2, as.integer(args[2]), 500)
n_cores      <- ifelse(length(args) >= 3, as.integer(args[3]), 10)
gene_seed    <- ifelse(length(args) >= 4, as.integer(args[4]), 42)
data_seed    <- ifelse(length(args) >= 5, as.integer(args[5]), 123)
effect_scale <- ifelse(length(args) >= 6, args[6], "prob")
stopifnot(effect_scale %in% c("prob", "link"))

library(data.table)
library(gamlss)
library(gamlss.dist)
library(truncnorm)
library(parallel)

logit    <- function(p) log(p / (1 - p))
inv_logit <- function(x) exp(x) / (1 + exp(x))

if(effect_scale == "link") {
  out_dir <- sprintf("sim_N%d_G%d_link", n_samples, n_genes)
} else {
  out_dir <- sprintf("sim_N%d_G%d_realfloor", n_samples, n_genes)   # REALFLOOR (written in the TRUNC_L-specific cwd)
}
dir.create(out_dir, showWarnings = FALSE)

cat("=== New Simulation (rBEINF) ===\n")
cat(sprintf("Samples: %d | Genes: %d | Cores: %d\n", n_samples, n_genes, n_cores))
cat(sprintf("Gene seed: %d (fixed across all N) | Data seed: %d\n", gene_seed, data_seed))
cat(sprintf("Effect scale: %s\n", effect_scale))
cat(sprintf("Output: %s\n", out_dir))
cat("Effect type proportions: middle=0.30, zero=0.15, one=0.15, mixed=0.40\n")

### ============================================================================
### Step 1a: Generate gene parameters (FIXED seed, same for all sample sizes)
###          IDENTICAL to v11
### ============================================================================
cat(sprintf("\n--- Step 1a: Generating gene parameters (gene_seed=%d) ---\n", gene_seed))

set.seed(gene_seed)

# Effect type proportions: 30% middle, 15% zero, 15% one, 40% mixed
effect_types <- sample(c("middle_value", "zero_proportion", "one_proportion", "mixed_effect"),
                       n_genes, replace = TRUE, prob = c(0.3, 0.15, 0.15, 0.4))

ground_truth <- data.table(
  gene_id        = sprintf("gene_%04d", 1:n_genes),
  effect_type    = effect_types,
  effect_size    = NA_real_,
  p_zero         = NA_real_,
  p_one          = NA_real_,
  middle_mean    = NA_real_,
  middle_sd      = NA_real_,
  true_mu_effect = NA_real_,
  true_nu_effect = NA_real_,
  true_tau_effect = NA_real_,
  dir_mu = NA_integer_,
  dir_nu = NA_integer_,
  dir_tau = NA_integer_,
  scale_factor = NA_real_
)

for(g in 1:n_genes) {
  et <- effect_types[g]

  # Baseline parameters — IDENTICAL to v11
  ## NO-THRESHOLD generation: the baseline boundary proportions are drawn from ONE common
  ## distribution for every gene, independent of its effect type -- total boundary share p0 + p1 ~ U(0, 0.9), split
  ## uniformly. Exactly two runif() draws per gene, as in every branch of the realfloor generator, so the effect-type
  ## assignment, effect sizes, middle means/sds, genotypes and all later draws follow the same random stream.
  p01 <- runif(1, 0, 0.9)
  p0 <- p01 * runif(1)
  p1 <- p01 - p0
  mid_mean <- runif(1, 0.1, 0.9)
  mid_sd   <- rgamma(1, shape = 2, rate = 20)

  # Effect size: 10% null — IDENTICAL to v11
  if(runif(1) < 0.1) {
    eff_size <- 0
  } else {
    eff_size <- rtruncnorm(1, a = 0.001, b = 0.5, mean = 0.1, sd = 0.2)
  }

  # Determine directions — IDENTICAL to v11
  if(et == "middle_value") {
    dir_mu <- if(mid_mean < 0.5) 1L else -1L
    dir_nu <- 0L; dir_tau <- 0L
  } else if(et == "zero_proportion") {
    dir_mu <- 0L; dir_tau <- 0L
    dir_nu <- if(p0 < 0.5) 1L else -1L
  } else if(et == "one_proportion") {
    dir_mu <- 0L; dir_nu <- 0L
    dir_tau <- if(p1 < 0.5) 1L else -1L
  } else {  # mixed_effect
    dir_mu <- sample(c(-1L, 1L), 1)
    dir_nu <- sample(c(-1L, 1L), 1)
    dir_tau <- 0L
  }

  ### CHANGE #2: Proportional scaling instead of hard clamping =================
  ### Check if effect at gt=2 exceeds boundaries; if so, scale down
  ### Supports both prob and link scale effects
  p_mid_base <- max(1 - p0 - p1, 0.01)
  nu_base <- p0 / p_mid_base
  tau_base <- p1 / p_mid_base

  sf <- 1.0
  for(iter in 1:20) {
    eff_scaled <- eff_size * sf
    valid <- TRUE
    for(gt_val in c(0, 1, 2)) {
      ev <- eff_scaled * gt_val

      if(effect_scale == "prob") {
        # --- Prob scale: additive on probability ---
        if(et == "middle_value") {
          cm <- mid_mean + dir_mu * ev
          cp0 <- p0; cp1 <- p1
        } else if(et == "zero_proportion") {
          cm <- mid_mean
          cp0 <- p0 + dir_nu * ev
          cp1_max <- 1 - cp0 - 0.1
          cp1 <- min(p1, cp1_max)
        } else if(et == "one_proportion") {
          cm <- mid_mean
          cp1 <- p1 + dir_tau * ev
          cp0_max <- 1 - cp1 - 0.1
          cp0 <- min(p0, cp0_max)
        } else {  # mixed_effect
          cm <- mid_mean + dir_mu * ev
          cp0 <- p0 + dir_nu * ev
          cp1_max <- 1 - cp0 - 0.1
          cp1 <- min(p1, cp1_max)
        }
      } else {
        # --- Link scale: additive on logit(mu)/log(nu)/log(tau) ---
        if(et == "middle_value") {
          cm <- inv_logit(logit(mid_mean) + dir_mu * ev)
          cp0 <- p0; cp1 <- p1
        } else if(et == "zero_proportion") {
          cm <- mid_mean
          nu_new <- nu_base * exp(dir_nu * ev)
          cp0 <- nu_new / (1 + nu_new + tau_base)
          cp1 <- tau_base / (1 + nu_new + tau_base)
        } else if(et == "one_proportion") {
          cm <- mid_mean
          tau_new <- tau_base * exp(dir_tau * ev)
          cp0 <- nu_base / (1 + nu_base + tau_new)
          cp1 <- tau_new / (1 + nu_base + tau_new)
        } else {  # mixed_effect
          cm <- inv_logit(logit(mid_mean) + dir_mu * ev)
          nu_new <- nu_base * exp(dir_nu * ev)
          cp0 <- nu_new / (1 + nu_new + tau_base)
          cp1 <- tau_base / (1 + nu_new + tau_base)
        }
      }

      # Check bounds (same for both scales)
      if(cm < 0.01 || cm > 0.99) { valid <- FALSE; break }
      if(cp0 < 0.01 || cp0 > 0.8) { valid <- FALSE; break }
      if(cp1 < 0 || cp1 > 0.8) { valid <- FALSE; break }
      if(cp0 + cp1 > 0.9) { valid <- FALSE; break }
    }
    if(valid) break
    sf <- sf * 0.8
  }
  eff_size_adj <- eff_size * sf

  # True effects on probability scale (gt=2 vs gt=0, for correlation analysis)
  if(effect_scale == "prob") {
    true_mu  <- dir_mu  * eff_size_adj
    true_nu  <- dir_nu  * eff_size_adj
    true_tau <- dir_tau * eff_size_adj
  } else {
    # Link scale: compute prob-scale effect from link transformation
    true_mu <- inv_logit(logit(mid_mean) + dir_mu * 2 * eff_size_adj) - mid_mean
    nu_at0 <- nu_base; tau_at0 <- tau_base
    nu_at2 <- nu_base * exp(dir_nu * 2 * eff_size_adj)
    tau_at2 <- tau_base * exp(dir_tau * 2 * eff_size_adj)
    true_nu <- nu_at2/(1+nu_at2+tau_at2) - nu_at0/(1+nu_at0+tau_at0)
    true_tau <- tau_at2/(1+nu_at2+tau_at2) - tau_at0/(1+nu_at0+tau_at0)
  }

  set(ground_truth, g, "effect_size",     eff_size_adj)
  set(ground_truth, g, "p_zero",          p0)
  set(ground_truth, g, "p_one",           p1)
  set(ground_truth, g, "middle_mean",     mid_mean)
  set(ground_truth, g, "middle_sd",       mid_sd)
  set(ground_truth, g, "true_mu_effect",  true_mu)
  set(ground_truth, g, "true_nu_effect",  true_nu)
  set(ground_truth, g, "true_tau_effect", true_tau)
  set(ground_truth, g, "dir_mu",  dir_mu)
  set(ground_truth, g, "dir_nu",  dir_nu)
  set(ground_truth, g, "dir_tau", dir_tau)
  set(ground_truth, g, "scale_factor", sf)
}

cat(sprintf("Effect type distribution: %s\n",
    paste(names(table(effect_types)), table(effect_types), sep = "=", collapse = ", ")))
cat(sprintf("Null genes (effect=0): %d\n", sum(ground_truth$effect_size == 0)))
cat(sprintf("Genes needing scaling: %d\n", sum(ground_truth$scale_factor < 1)))

### ============================================================================
### Step 1b: Generate sample data (DATA seed, varies with n_samples)
### ============================================================================
cat(sprintf("\n--- Step 1b: Generating sample data (data_seed=%d, N=%d) ---\n", data_seed, n_samples))

set.seed(data_seed + n_samples)

# Fixed genotype proportions: 50:30:20 — IDENTICAL to v11
gt_probs <- c(0.5, 0.3, 0.2)
genotypes <- sample(c(0, 1, 2), n_samples, replace = TRUE, prob = gt_probs)   # kept for compatibility (not used below)
## INDEPENDENT GENOTYPES: every gene gets its OWN genotype vector, drawn from the same 50:30:20 frequencies,
## so promoters are independent replicates (the original design shared one genotype vector across all genes).
geno_mat <- matrix(NA_real_, nrow = n_samples, ncol = n_genes, dimnames = list(sample_names <- sprintf("sample_%03d", 1:n_samples), sprintf("gene_%04d", 1:n_genes)))
sample_names <- sprintf("sample_%03d", 1:n_samples)

# Storage
activity_mat <- matrix(NA, nrow = n_samples, ncol = n_genes)
colnames(activity_mat) <- sprintf("gene_%04d", 1:n_genes)
rownames(activity_mat) <- sample_names

for(g in 1:n_genes) {
  row <- ground_truth[g]
  et <- row$effect_type
  eff_size <- row$effect_size
  p0 <- row$p_zero
  p1 <- row$p_one
  mid_mean <- row$middle_mean
  mid_sd <- row$middle_sd
  dir_mu <- row$dir_mu
  dir_nu <- row$dir_nu
  dir_tau <- row$dir_tau
  gvec <- sample(c(0, 1, 2), n_samples, replace = TRUE, prob = gt_probs); geno_mat[, g] <- gvec   # INDEPENDENT GENOTYPES

  # Generate data for each genotype
  for(gt in 0:2) {
    idx <- which(gvec == gt)
    if(length(idx) == 0) next
    n_curr <- length(idx)
    effect_val <- eff_size * gt

    ### CHANGE #1: rBEINF instead of manual mixing =============================
    # Compute genotype-specific BEINF parameters
    p_mid_base <- max(1 - p0 - p1, 0.01)
    nu_base_g <- p0 / p_mid_base
    tau_base_g <- p1 / p_mid_base

    if(effect_scale == "prob") {
      # --- Prob scale: additive on probability, then convert to BEINF ---
      if(et == "middle_value") {
        curr_mid <- mid_mean + dir_mu * effect_val
        curr_p0  <- p0
        curr_p1  <- p1
      } else if(et == "zero_proportion") {
        curr_mid <- mid_mean
        curr_p0  <- p0 + dir_nu * effect_val
        max_p1   <- 1 - curr_p0 - 0.1
        curr_p1  <- max(min(p1, max_p1), 0)
      } else if(et == "one_proportion") {
        curr_mid <- mid_mean
        curr_p1  <- p1 + dir_tau * effect_val
        max_p0   <- 1 - curr_p1 - 0.1
        curr_p0  <- max(min(p0, max_p0), 0)
      } else {  # mixed_effect
        curr_mid <- mid_mean + dir_mu * effect_val
        curr_p0  <- p0 + dir_nu * effect_val
        max_p1   <- 1 - curr_p0 - 0.1
        curr_p1  <- max(min(p1, max_p1), 0)
      }
      # Safety clamp
      curr_mid <- min(max(curr_mid, 0.01), 0.99)
      curr_p0  <- min(max(curr_p0, 0.001), 0.8)
      curr_p1  <- min(max(curr_p1, 0), 0.8)
      # Convert to BEINF parameters
      p_mid <- max(1 - curr_p0 - curr_p1, 0.01)
      mu_beinf  <- curr_mid
      nu_beinf  <- curr_p0 / p_mid
      tau_beinf <- curr_p1 / p_mid
    } else {
      # --- Link scale: additive on logit(mu)/log(nu)/log(tau) ---
      if(et == "middle_value") {
        mu_beinf  <- inv_logit(logit(mid_mean) + dir_mu * effect_val)
        nu_beinf  <- nu_base_g
        tau_beinf <- tau_base_g
      } else if(et == "zero_proportion") {
        mu_beinf  <- mid_mean
        nu_beinf  <- nu_base_g * exp(dir_nu * effect_val)
        tau_beinf <- tau_base_g
      } else if(et == "one_proportion") {
        mu_beinf  <- mid_mean
        nu_beinf  <- nu_base_g
        tau_beinf <- tau_base_g * exp(dir_tau * effect_val)
      } else {  # mixed_effect
        mu_beinf  <- inv_logit(logit(mid_mean) + dir_mu * effect_val)
        nu_beinf  <- nu_base_g * exp(dir_nu * effect_val)
        tau_beinf <- tau_base_g
      }
    }

    # Clamp BEINF parameters
    sigma_beinf <- min(mid_sd / sqrt(mu_beinf * (1 - mu_beinf) + 1e-10), 0.99)
    sigma_beinf <- max(sigma_beinf, 0.01)
    nu_beinf  <- max(nu_beinf, 1e-4)
    tau_beinf <- max(tau_beinf, 1e-4)

    ## rBEINF(n, ...) is exactly qBEINF(runif(n), ...) (gamlss.dist); drawing the uniform here keeps the RNG stream identical
    ## and tells which draws come from the BETA part: qBEINF gives 0 for u <= nu/(1+nu+tau), 1 for u >= (1+nu)/(1+nu+tau).
    u_draw <- runif(n_curr)
    vals <- qBEINF(u_draw, mu = mu_beinf, sigma = sigma_beinf, nu = nu_beinf, tau = tau_beinf)
    TL <- as.numeric(Sys.getenv("TRUNC_L", "0"))                       # FLOOR on beta-part (interior) draws
    if (TL > 0) {
      from_beta <- u_draw > nu_beinf / (1 + nu_beinf + tau_beinf) & u_draw < (1 + nu_beinf) / (1 + nu_beinf + tau_beinf)
      ## includes beta draws that qbeta rounded to exactly 0 or 1 (U-shaped beta), which the original generator silently
      ## turned into boundary values (e.g. gene_0029, middle_value: 64/129 dosage-2 samples == 1 at N = 600)
      if (any(from_beta)) {
        if (Sys.getenv("TRUNC_MODE", "clamp") == "redraw") {
          sa <- mu_beinf * (1 - sigma_beinf^2) / sigma_beinf^2; sb <- sa * (1 - mu_beinf) / mu_beinf
          lo <- pbeta(TL, sa, sb); hi <- pbeta(1 - TL, sa, sb)
          vals[from_beta] <- pmin(pmax(qbeta(lo + pbeta(vals[from_beta], sa, sb) * (hi - lo), sa, sb), TL), 1 - TL)
        } else vals[from_beta] <- pmin(pmax(vals[from_beta], TL), 1 - TL)
      }
    }
    activity_mat[idx, g] <- vals
  }
}

# Save ground truth and data
fwrite(ground_truth, file.path(out_dir, "ground_truth.txt"), sep = "\t")
save(activity_mat, genotypes, geno_mat, sample_names, file = file.path(out_dir, "sim_data.RData"))   # geno_mat = per-gene genotypes (used by the loaders)
if (Sys.getenv("GEN_ONLY") == "1") { cat("GEN_ONLY: data saved, skipping self-fit\n"); quit(save = "no") }
cat(sprintf("Generated %d genes x %d samples\n", n_genes, n_samples))
cat(sprintf("Zeros: %.1f%%, Ones: %.1f%%, Middle: %.1f%%\n",
    100 * mean(activity_mat == 0, na.rm = TRUE),
    100 * mean(activity_mat == 1, na.rm = TRUE),
    100 * mean(activity_mat > 0 & activity_mat < 1, na.rm = TRUE)))

### ============================================================================
### Step 2: Fit models (ZOIB + Linear) — IDENTICAL to v11
### ============================================================================
cat("\n--- Step 2: Fitting ZOIB and Linear models ---\n")

fit_one_gene <- function(g) {
  gene_id <- colnames(activity_mat)[g]
  y <- activity_mat[, g]
  gt <- genotypes

  # Skip if too few variation
  if(sd(y) < 1e-6) return(NULL)

  ctrl <- gamlss.control(n.cyc = 200, trace = FALSE)
  dat <- data.frame(y = y, gt = gt)

  ### --- ZOIB model ---
  zoib_result <- tryCatch({
    has_zeros <- any(y == 0)
    has_ones  <- any(y == 1)

    if(has_zeros && has_ones) {
      model_type <- "BEINF"
      fam <- BEINF(mu.link = "logit", sigma.link = "logit",
                    nu.link = "log", tau.link = "log")
    } else if(has_zeros && !has_ones) {
      model_type <- "BEZI"
      fam <- BEZI(mu.link = "logit", sigma.link = "log", nu.link = "logit")
    } else if(!has_zeros && has_ones) {
      model_type <- "BEOI"
      fam <- BEOI(mu.link = "logit", sigma.link = "log", nu.link = "logit")
    } else {
      model_type <- "BE"
      fam <- BE(mu.link = "logit", sigma.link = "logit")
    }

    if(model_type %in% c("BEINF")) {
      m0 <- gamlss(y ~ 1, sigma.formula = ~1, nu.formula = ~1,
                    tau.formula = ~1, family = fam, data = dat, control = ctrl)
    } else if(model_type %in% c("BEZI", "BEOI")) {
      m0 <- gamlss(y ~ 1, sigma.formula = ~1, nu.formula = ~1,
                    family = fam, data = dat, control = ctrl)
    } else {
      m0 <- gamlss(y ~ 1, sigma.formula = ~1,
                    family = fam, data = dat, control = ctrl)
    }

    if(model_type == "BEINF") {
      m1 <- gamlss(y ~ gt, sigma.formula = ~gt, nu.formula = ~gt,
                    tau.formula = ~gt, family = fam, data = dat,
                    start.from = m0, control = ctrl)
    } else if(model_type %in% c("BEZI", "BEOI")) {
      m1 <- gamlss(y ~ gt, sigma.formula = ~gt, nu.formula = ~gt,
                    family = fam, data = dat,
                    start.from = m0, control = ctrl)
    } else {
      m1 <- gamlss(y ~ gt, sigma.formula = ~gt,
                    family = fam, data = dat,
                    start.from = m0, control = ctrl)
    }

    lr_out <- tryCatch({
      capture.output(LR.test(m0, m1))
    }, error = function(e) list(error = TRUE))
    if(is.list(lr_out) && !is.null(lr_out$error)) {
      lrt_p <- 0.9999
    } else {
      lrt_p <- as.numeric(gsub(".*p-value= (.*)", "\\1", lr_out[length(lr_out)]))
    }

    ms <- as.data.frame(summary(m1))
    gt_rows <- grep("^gt(\\.[0-9]+)?$", row.names(ms))
    pnames <- m1$parameters

    est <- ms[gt_rows, "Estimate"]
    pval <- ms[gt_rows, "Pr(>|t|)"]
    names(est) <- paste0(pnames, "_estimate")
    names(pval) <- paste0(pnames, "_pvalue")

    mu_est <- est["mu_estimate"]
    mu_pval <- pval["mu_pvalue"]
    sigma_est <- est["sigma_estimate"]
    sigma_pval <- pval["sigma_pvalue"]

    if(model_type == "BEINF") {
      nu_est  <- est["nu_estimate"]
      tau_est <- est["tau_estimate"]
      nu_pval <- pval["nu_pvalue"]
      tau_pval <- pval["tau_pvalue"]
    } else if(model_type == "BEZI") {
      nu_est  <- est["nu_estimate"]
      tau_est <- NA
      nu_pval <- pval["nu_pvalue"]
      tau_pval <- NA
    } else if(model_type == "BEOI") {
      nu_est  <- NA
      tau_est <- est["nu_estimate"]
      nu_pval <- NA
      tau_pval <- pval["nu_pvalue"]
    } else {
      nu_est <- NA; tau_est <- NA
      nu_pval <- NA; tau_pval <- NA
    }

    ### --- Probability-scale effects from model intercepts ---
    mu_coefs <- coef(m1, "mu")
    mu_int <- mu_coefs[1]; mu_beta <- mu_coefs[2]
    mu_prob_eff <- inv_logit(mu_int + 2 * mu_beta) - inv_logit(mu_int)

    p0_prob_eff <- 0; p1_prob_eff <- 0
    if(model_type == "BEINF") {
      nu_coefs <- coef(m1, "nu"); tau_coefs <- coef(m1, "tau")
      nu_int <- nu_coefs[1]; nu_b <- nu_coefs[2]
      tau_int <- tau_coefs[1]; tau_b <- tau_coefs[2]
      nu_0 <- exp(nu_int); tau_0 <- exp(tau_int)
      nu_2 <- exp(nu_int + 2*nu_b); tau_2 <- exp(tau_int + 2*tau_b)
      p0_prob_eff <- nu_2/(1+nu_2+tau_2) - nu_0/(1+nu_0+tau_0)
      p1_prob_eff <- tau_2/(1+nu_2+tau_2) - tau_0/(1+nu_0+tau_0)
    } else if(model_type == "BEZI") {
      nu_coefs <- coef(m1, "nu")
      nu_int <- nu_coefs[1]; nu_b <- nu_coefs[2]
      p0_prob_eff <- inv_logit(nu_int + 2*nu_b) - inv_logit(nu_int)
    } else if(model_type == "BEOI") {
      nu_coefs <- coef(m1, "nu")
      nu_int <- nu_coefs[1]; nu_b <- nu_coefs[2]
      p1_prob_eff <- inv_logit(nu_int + 2*nu_b) - inv_logit(nu_int)
    }

    list(model_type = model_type, lrt_pvalue = lrt_p, converged = m1$converged,
         mu_estimate = mu_est, sigma_estimate = sigma_est,
         nu_estimate = nu_est, tau_estimate = tau_est,
         mu_pvalue = mu_pval, sigma_pvalue = sigma_pval,
         nu_pvalue = nu_pval, tau_pvalue = tau_pval,
         mu_prob_effect = mu_prob_eff,
         p0_prob_effect = p0_prob_eff,
         p1_prob_effect = p1_prob_eff)
  }, error = function(e) {
    list(model_type = NA, lrt_pvalue = NA, converged = FALSE,
         mu_estimate = NA, sigma_estimate = NA,
         nu_estimate = NA, tau_estimate = NA,
         mu_pvalue = NA, sigma_pvalue = NA,
         nu_pvalue = NA, tau_pvalue = NA,
         mu_prob_effect = NA, p0_prob_effect = NA, p1_prob_effect = NA)
  })

  ### --- Linear model --- IDENTICAL to v11
  linear_result <- tryCatch({
    y_lin <- y
    y_lin[y_lin == 0] <- 0.001
    y_lin[y_lin == 1] <- 0.999
    y_logit <- log(y_lin / (1 - y_lin))
    dat_lin <- data.frame(y_logit = y_logit, gt = gt)

    m0_lin <- gamlss(y_logit ~ 1, sigma.formula = ~1,
                      family = NO, data = dat_lin, control = ctrl)
    m1_lin <- gamlss(y_logit ~ gt, sigma.formula = ~gt,
                      family = NO, data = dat_lin,
                      start.from = m0_lin, control = ctrl)

    lr_out_lin <- tryCatch({
      capture.output(LR.test(m0_lin, m1_lin))
    }, error = function(e) list(error = TRUE))
    if(is.list(lr_out_lin) && !is.null(lr_out_lin$error)) {
      lrt_p_lin <- 0.9999
    } else {
      lrt_p_lin <- as.numeric(gsub(".*p-value= (.*)", "\\1",
                                    lr_out_lin[length(lr_out_lin)]))
    }

    ms_lin <- as.data.frame(summary(m1_lin))
    gt_rows_lin <- grep("^gt(\\.[0-9]+)?$", row.names(ms_lin))
    lin_est  <- ms_lin[gt_rows_lin, "Estimate"]
    lin_pval <- ms_lin[gt_rows_lin, "Pr(>|t|)"]
    names(lin_est) <- paste0(m1_lin$parameters, "_estimate")
    names(lin_pval) <- paste0(m1_lin$parameters, "_pvalue")

    lin_coefs <- coef(m1_lin, "mu")
    lin_int <- lin_coefs[1]; lin_beta <- lin_coefs[2]
    lin_prob_eff <- inv_logit(lin_int + 2 * lin_beta) - inv_logit(lin_int)

    list(linear_lrt_pvalue = lrt_p_lin,
         linear_mu_estimate = lin_est["mu_estimate"],
         linear_sigma_estimate = lin_est["sigma_estimate"],
         linear_mu_pvalue = lin_pval["mu_pvalue"],
         linear_sigma_pvalue = lin_pval["sigma_pvalue"],
         linear_prob_effect = lin_prob_eff)
  }, error = function(e) {
    list(linear_lrt_pvalue = NA, linear_mu_estimate = NA,
         linear_sigma_estimate = NA, linear_mu_pvalue = NA,
         linear_sigma_pvalue = NA, linear_prob_effect = NA)
  })

  out <- c(list(gene_id = gene_id), zoib_result, linear_result)
  out <- lapply(out, function(x) if(is.null(x)) NA else x[1])
  as.data.table(out)
}

cat(sprintf("Fitting %d genes with %d cores...\n", n_genes, n_cores))
t0 <- proc.time()

results_list <- mclapply(1:n_genes, fit_one_gene, mc.cores = n_cores)

elapsed <- (proc.time() - t0)[3]
cat(sprintf("Done in %.1f seconds (%.2f sec/gene)\n", elapsed, elapsed / n_genes))

results_dt <- rbindlist(Filter(function(x) !is.null(x) && is.data.table(x), results_list),
                        fill = TRUE)

cat(sprintf("Successfully fitted: %d / %d genes\n", nrow(results_dt), n_genes))

### ============================================================================
### Step 3: Merge with ground truth — IDENTICAL to v11
### ============================================================================
cat("\n--- Step 3: Merging results with ground truth ---\n")

merged <- merge(ground_truth, results_dt, by = "gene_id", all = TRUE)
fwrite(merged, file.path(out_dir, "results_merged.txt"), sep = "\t")

### ============================================================================
### Step 4: Correlation analysis (signed) — IDENTICAL to v11
### ============================================================================
cat("\n=== Effect Size Recovery: signed correlation (probability scale) ===\n")

eff <- merged[effect_size > 0 & !is.na(mu_prob_effect)]
eff_lin <- merged[effect_size > 0 & !is.na(linear_prob_effect)]

summary_rows <- list()

for(et in c("middle_value", "zero_proportion", "one_proportion", "mixed_effect")) {
  sub <- eff[effect_type == et]
  sub_lin <- eff_lin[effect_type == et]
  if(nrow(sub) < 5) next

  cat(sprintf("\n%s (ZOIB n=%d, Linear n=%d):\n", et, nrow(sub), nrow(sub_lin)))

  r_zoib_mu <- tryCatch(cor(sub$mu_prob_effect, sub$true_mu_effect, use="complete"),
                         warning = function(w) NA)
  r_lin_mu  <- tryCatch(cor(sub_lin$linear_prob_effect, sub_lin$true_mu_effect, use="complete"),
                         warning = function(w) NA)
  cat(sprintf("  mu:  ZOIB r=%.4f  |  Linear r=%.4f\n", r_zoib_mu, r_lin_mu))
  summary_rows[[length(summary_rows)+1]] <- data.table(
    effect_type=et, parameter="mu", zoib_r=r_zoib_mu, linear_r=r_lin_mu,
    zoib_n=nrow(sub), linear_n=nrow(sub_lin))

  if(et %in% c("zero_proportion", "mixed_effect")) {
    r_zoib_p0 <- tryCatch(cor(sub$p0_prob_effect, sub$true_nu_effect, use="complete"),
                           warning = function(w) NA)
    r_lin_p0  <- tryCatch(cor(-sub_lin$linear_prob_effect, sub_lin$true_nu_effect, use="complete"),
                           warning = function(w) NA)
    cat(sprintf("  p0:  ZOIB r=%.4f  |  Linear r=%.4f\n", r_zoib_p0, r_lin_p0))
    summary_rows[[length(summary_rows)+1]] <- data.table(
      effect_type=et, parameter="p0", zoib_r=r_zoib_p0, linear_r=r_lin_p0,
      zoib_n=nrow(sub), linear_n=nrow(sub_lin))
  }

  if(et %in% c("one_proportion", "mixed_effect")) {
    r_zoib_p1 <- tryCatch(cor(sub$p1_prob_effect, sub$true_tau_effect, use="complete"),
                           warning = function(w) NA)
    r_lin_p1  <- tryCatch(cor(sub_lin$linear_prob_effect, sub_lin$true_tau_effect, use="complete"),
                           warning = function(w) NA)
    cat(sprintf("  p1:  ZOIB r=%.4f  |  Linear r=%.4f\n", r_zoib_p1, r_lin_p1))
    summary_rows[[length(summary_rows)+1]] <- data.table(
      effect_type=et, parameter="p1", zoib_r=r_zoib_p1, linear_r=r_lin_p1,
      zoib_n=nrow(sub), linear_n=nrow(sub_lin))
  }
}

summary_dt <- rbindlist(summary_rows)
fwrite(summary_dt, file.path(out_dir, "correlation_summary.txt"), sep = "\t")

### ============================================================================
### Step 5: Power analysis — IDENTICAL to v11
### ============================================================================
cat("\n=== Power Analysis (significance at alpha=0.05) ===\n")

for(et in c("middle_value", "zero_proportion", "one_proportion", "mixed_effect")) {
  sub <- merged[effect_type == et & effect_size > 0]
  if(nrow(sub) < 5) next

  zoib_sig <- sum(sub$lrt_pvalue < 0.05, na.rm = TRUE)
  lin_sig  <- sum(sub$linear_lrt_pvalue < 0.05, na.rm = TRUE)
  n_total  <- sum(!is.na(sub$lrt_pvalue))

  cat(sprintf("%s (n=%d): ZOIB %d/%d (%.1f%%)  |  Linear %d/%d (%.1f%%)\n",
      et, nrow(sub),
      zoib_sig, n_total, 100*zoib_sig/n_total,
      lin_sig, n_total, 100*lin_sig/n_total))
}

null_genes <- merged[effect_size == 0]
if(nrow(null_genes) > 0) {
  zoib_fp <- sum(null_genes$lrt_pvalue < 0.05, na.rm = TRUE)
  lin_fp  <- sum(null_genes$linear_lrt_pvalue < 0.05, na.rm = TRUE)
  n_null  <- sum(!is.na(null_genes$lrt_pvalue))
  cat(sprintf("\nNull genes (n=%d): ZOIB FP=%d (%.1f%%)  |  Linear FP=%d (%.1f%%)\n",
      nrow(null_genes),
      zoib_fp, 100*zoib_fp/n_null,
      lin_fp, 100*lin_fp/n_null))
}

cat(sprintf("\nAll results saved to %s/\n", out_dir))
cat("=== Done ===\n")
