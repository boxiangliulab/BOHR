### packages===============================================================
library(tidyr)
library(dplyr)
library(parallel)
library(foreach)
library(doParallel)
library(coda)
library(gamlss)
library(stringr)
library(MASS)
library(qvalue)
library(VariantAnnotation)
library(data.table)
library(GenomicRanges)
library(ACAT)

# Helper functions for probability-scale conversion ============================
.logit <- function(p) log(p / (1 - p))
.inv_logit <- function(x) exp(x) / (1 + exp(x))

# Compute proportions===================================================================
calculate_proportions <- function(data, activity_col) {
  total <- nrow(data)
  zeros <- sum(data[[activity_col]] == 0)
  ones <- sum(data[[activity_col]] == 1)

  prop_zeros <- zeros / total
  prop_ones <- ones / total

  return(list(prop_zeros = prop_zeros, prop_ones = prop_ones))
}

# Select the model===================================================================
select_model <- function(data, activity_col) {
  props <- calculate_proportions(data, activity_col)

  if (props$prop_zeros > 0 && props$prop_ones > 0) {
    return("BEINF")
  } else if (props$prop_zeros > 0 && props$prop_ones == 0) {
    return("BEZI")
  } else if (props$prop_zeros == 0 && props$prop_ones > 0) {
    return("BEOI")
  } else {
    return("BE")  # If there is no obvious inflation, use the standard Beta distribution
  }
}

# Per-parameter LRT helper ====================================================
# For one SNP and a fitted full model, fit 4 reduced models that drop the SNP
# from exactly one of mu/sigma/nu/tau and return df=1 chi-square LRT p-values
# keyed by the gamlss internal parameter name. Reduced models are built with
# explicit per-parameter formulas to mirror fit_gamlss_model() behaviour.
compute_per_param_lrt <- function(snp_data, activity_col, factors, snp,
                                  model_type, model_with) {
  ll_full <- tryCatch(as.numeric(logLik(model_with)), error = function(e) NA)
  if(!is.finite(ll_full)) return(list(mu = NA, sigma = NA, nu = NA, tau = NA))
  ctrl <- gamlss.control(n.cyc = 200, trace = FALSE)
  # mu/nu/tau formulas in the original fit include factors; sigma only the SNP
  rhs <- function(inc_snp) {
    terms <- c(factors, if(inc_snp) snp)
    if(length(terms) == 0) "1" else paste(terms, collapse = " + ")
  }
  fit_red <- function(mu_inc, sg_inc, nu_inc, tau_inc) {
    mu_f  <- as.formula(paste0(activity_col, " ~ ", rhs(mu_inc)))
    # sigma.formula matches fit_gamlss_model: ~ snp  (or ~ 1 when dropped),
    # no covariates — this preserves nesting against the full model
    sg_f  <- if(sg_inc) as.formula(paste0("~ ", snp)) else as.formula("~ 1")
    fam_obj <- switch(model_type,
      "BEINF" = BEINF(mu.link="logit", sigma.link="logit",
                      nu.link="log",   tau.link="log"),
      "BEZI"  = BEZI(mu.link="logit",  sigma.link="log", nu.link="logit"),
      "BEOI"  = BEOI(mu.link="logit",  sigma.link="log", nu.link="logit"),
      "BE"    = BE(mu.link="logit",    sigma.link="logit"))
    args <- list(formula = mu_f, sigma.formula = sg_f,
                 family = fam_obj, data = snp_data, control = ctrl)
    if(model_type %in% c("BEINF","BEZI","BEOI"))
      args$nu.formula  <- as.formula(paste0("~ ", rhs(nu_inc)))
    if(model_type == "BEINF")
      args$tau.formula <- as.formula(paste0("~ ", rhs(tau_inc)))
    tryCatch(suppressWarnings(do.call(gamlss, args)),
             error = function(e) NULL)
  }
  lrt_p <- function(m_red) {
    if(is.null(m_red) || !isTRUE(m_red$converged)) return(NA)
    ll_r <- tryCatch(as.numeric(logLik(m_red)), error = function(e) NA)
    if(!is.finite(ll_r)) return(NA)
    pchisq(2 * (ll_full - ll_r), df = 1, lower.tail = FALSE)
  }
  list(
    mu    = lrt_p(fit_red(FALSE, TRUE,  TRUE,  TRUE)),
    sigma = lrt_p(fit_red(TRUE,  FALSE, TRUE,  TRUE)),
    nu    = if(model_type %in% c("BEINF","BEZI","BEOI"))
              lrt_p(fit_red(TRUE, TRUE, FALSE, TRUE)) else NA,
    tau   = if(model_type == "BEINF")
              lrt_p(fit_red(TRUE, TRUE, TRUE, FALSE)) else NA
  )
}

# Regression model functions===================================================================
# [OPTIMIZED] Add trace=FALSE to switch off the iteration output and reduce I/O overhead
# [OPTIMIZED] Add the start.from parameter so that an existing model can initialize the fit and speed up convergence
fit_gamlss_model <- function(tissue_data, formula, model_type, start.from = NULL) {
  tryCatch({
    all_vars <- all.vars(formula)[-1]  # Drop the response variable
    # Check whether an SNP variable is present
    snp_vars <- all_vars[grepl("b38", all_vars)]
    # Construct sigma_formula
    if(length(snp_vars) > 0) {
      sigma_formula <- as.formula(paste("~", paste(snp_vars, collapse = " + ")))
    } else {
      sigma_formula <- ~ 1
    }

    # [OPTIMIZED] Use trace=FALSE to switch off the iteration output
    ctrl <- gamlss.control(n.cyc = 200, trace = FALSE)

    model <- switch(model_type,
                    "BEINF" = gamlss(formula,
                                     sigma.formula = sigma_formula,
                                     nu.formula = formula,
                                     tau.formula = formula,
                                     family = BEINF(mu.link = "logit", sigma.link = "logit", nu.link = "log", tau.link = "log"),
                                     data = tissue_data,
                                     start.from = start.from,
                                     control = ctrl),
                    "BEOI" = gamlss(formula,
                                    sigma.formula = sigma_formula,
                                    nu.formula = formula,
                                    family = BEOI(mu.link = "logit", sigma.link = "log", nu.link = "logit"),
                                    data = tissue_data,
                                    start.from = start.from,
                                    control = ctrl),
                    "BEZI" = gamlss(formula,
                                    sigma.formula = sigma_formula,
                                    nu.formula = formula,
                                    family = BEZI(mu.link = "logit", sigma.link = "log", nu.link = "logit"),
                                    data = tissue_data,
                                    start.from = start.from,
                                    control = ctrl),
                    "BE" = gamlss(formula,
                                  sigma.formula = sigma_formula,
                                  family = BE(mu.link = "logit", sigma.link = "logit"),
                                  data = tissue_data,
                                  start.from = start.from,
                                  control = ctrl)
    )
    return(model)
  }, error = function(e) {
    warning(paste("Error fitting", model_type, "model:", e$message))
    return(NULL)
  })
}


full_gamlss_model <- function(pheno, full_data, model_type) {
  tryCatch({

  # Get all variable names (except the response variable)
  all_vars <- names(full_data)[names(full_data) != as.character(substitute(pheno))]

  # Check whether an SNP variable is present (format such as chr22_19486387_A_G_b38)
  snp_vars <- all_vars[grepl("b38", all_vars)]

  # Construct sigma_formula
  if(length(snp_vars) > 0) {
    # If an SNP is present, sigma_formula contains only the SNP
    sigma_formula <- as.formula(paste("~", paste(snp_vars, collapse = " + ")))
  } else {
    # If no SNP is present, sigma_formula has no variables
    sigma_formula <- ~ 1
  }


    model <- switch(model_type,
                    "BEINF" = gamlss(pheno ~ .,
                                      sigma.formula= sigma_formula,
                                      nu.formula= ~ .,
                                      tau.formula = ~ .,
                                      family = BEINF(mu.link = "logit", sigma.link = "logit", nu.link = "log", tau.link = "log"),
                                      data = full_data,
                                      control = gamlss.control(n.cyc = 200)),
                    "BEOI" = gamlss(pheno ~ .,
                                    sigma.formula= sigma_formula,
                                    nu.formula= ~ .,
                                    data = full_data,
                                    family = BEOI(mu.link = "logit", sigma.link = "log", nu.link = "logit"),
                                    control = gamlss.control(n.cyc = 200)),
                    "BEZI" = gamlss(pheno ~ .,
                                    sigma.formula= sigma_formula,
                                    nu.formula= ~ .,
                                    data = full_data,
                                    family = BEZI(mu.link = "logit", sigma.link = "log", nu.link = "logit"),
                                    control = gamlss.control(n.cyc = 200)),
                    "BE" = gamlss(pheno ~ .,
                                  sigma.formula= sigma_formula,
                                  family = BE(mu.link = "logit", sigma.link = "logit"),
                                  data = full_data,
                                  control = gamlss.control(n.cyc = 200))
    )
    return(model)
  }, error = function(e) {
    warning(paste("Error fitting", model_type, "model:", e$message))
    return(NULL)
  })
}

# get numeric genotype===================================================================
get_geno <- function(chr, pos, flank = 500000, vcffile){
  pos <- as.numeric(pos)
  # Build a GRanges object representing the SNP position
  snp <- GRanges(chr, IRanges(pos, pos))
  # Build the flanking interval
  region <- trim(flank(snp, flank, both = TRUE), seqinfo(Hsapiens)[seqlevels(snp)])


  # Extract the interval data from the preloaded VCF data instead of re-reading the file
  #vcf_ranges <- rowRanges(vcf_data)
  #overlap_idx <- findOverlaps(vcf_ranges, region)

  #if (length(overlap_idx) == 0) {
  #  return(NULL)
  #}

  # Extract the variants overlapping the interval
  #variant_indices <- queryHits(overlap_idx)
  #vcf <- vcf_data[variant_indices]

  # Read the VCF using this interval
  vcf <- readVcf(vcffile, "hg38", param = ScanVcfParam(which = region))

  # Extract the SNP information and genotype data
  snp_info <- rowRanges(vcf)
  geno <- as.matrix(geno(vcf)$GT)
  if(nrow(geno)==0){return(NULL)}
  # Convert the genotype data to numeric coding
  geno_num <- apply(geno, 2, function(x) {
    x <- as.character(x)
    x[x == "." | x == ".|."] <- sample(names(table(x))[which(table(x) == max(table(x)))], 1)
    x[x == "0|0"] <- 0
    x[x == "0|1" | x == "1|0"] <- 1
    x[x == "1|1"] <- 2
    as.numeric(x)
  })
  if(nrow(geno) == 1){geno_num <- as.data.frame(t(geno_num))}
  row.names(geno_num) <- row.names(geno)
  geno_num <- t(geno_num)

  return(geno_num)
}



# Beta correction function===================================================================
beta_correct <- function(p, perm_pvalues, dof1, dof2) {
    perm_pvalues <- as.numeric(perm_pvalues)
    perm_pvalues <- perm_pvalues[!is.na(perm_pvalues)]

    if (length(perm_pvalues) == 0) return(NA)

    # Compute the mean and variance
    mean_perm <- mean(perm_pvalues)
    var_perm <- var(perm_pvalues)

    if (is.na(var_perm) || is.na(mean_perm)) return(NA)

    # Estimate the parameters by the method of moments
    if (var_perm != 0 && mean_perm != 1) {
      beta_mm1 <- mean_perm * (mean_perm * (1 - mean_perm) / var_perm - 1)
      beta_mm2 <- beta_mm1 * (1 / mean_perm - 1)

      if (is.na(beta_mm1) || is.na(beta_mm2) || beta_mm1 <= 0 || beta_mm2 <= 0) return(NA)

      # Optimize the parameters by maximum likelihood; if this fails, use the method of moments
      tryCatch({
        fit_beta <- fitdistr(perm_pvalues, "beta", start = list(shape1 = beta_mm1, shape2 = beta_mm2))
        beta_ml1 <- fit_beta$estimate["shape1"]
        beta_ml2 <- fit_beta$estimate["shape2"]
      }, error = function(e) {
        message("MLE failed, using MM instead!")
        beta_ml1 <- beta_mm1
        beta_ml2 <- beta_mm2
      })

      # Correct the p-value
      return(list(beta_pvalue = pbeta(p, shape1 = beta_ml1, shape2 = beta_ml2, lower.tail = FALSE),
                  beta_ml1 = beta_ml1, beta_ml2 = beta_ml2))
    } else {
      return(list(beta_pvalue = NA, beta_ml1 = NA, beta_ml2 = NA))
    }
}

# get thresholds=================================================================
get_threshold <- function(t) {
  ap <- data.frame()
  for(i in names(t)){
    if(class(t[[i]])=="character"){next}

    p <- t[[i]]$lrt_pvalue
    p[p == 0] <- .Machine$double.eps
    p[p == 1] <- 1-.Machine$double.eps

    ap <- rbind(ap, data.frame(ap=ACAT(p), promoter=i))
  }
  ap <- na.omit(ap)
  new_acat <- ap$ap
  #ap$qval <- qvalue(ap$ap)$qvalues
  ap$qval <- p.adjust(ap$ap, method="BH")
  ap <- ap[order(ap$qval, decreasing=F),]

  set1 <- subset(ap, qval <= 0.05)
  set2 <- subset(ap, qval > 0.05)

  if(nrow(set1)==0){return(NULL)}

  pid1 <- set1[which.max(set1$qval),]$promoter
  pid2 <- set2[which.min(set2$qval),]$promoter
  print(c(min(t[[pid1]]$lrt_pvalue), min(t[[pid2]]$lrt_pvalue)))
  threshold <- min(c(min(t[[pid1]]$lrt_pvalue), min(t[[pid2]]$lrt_pvalue)))
  return(threshold)
}


# fit function===============================================================================
raw_fit <- function(data, activity_col, tissue_name, tissue_col="SMTSD") {
  cat("Processing:", tissue_name, activity_col, "\n")
  tissue_data <- data[data[[tissue_col]] == tissue_name, ]
  tissue_data <- cbind(tissue_data[tissue_col],tissue_data[activity_col])
  tissue_data <- na.omit(tissue_data)

  total_count <- nrow(tissue_data)
  if(total_count < 10){
    cat("Skip:", tissue_name, activity_col, "\n")
    return(NULL)
  }
  if (length(unique(tissue_data[[activity_col]])) == 1) {
    return(list(status = "Skip_fit",
      reason = "all_same_value",
      converged = FALSE,
      model_type = NA,
      mu = mean(tissue_data[[activity_col]]),
      sigma = NA,
      p0 = mean(tissue_data[[activity_col]] == 0),
      p1 = mean(tissue_data[[activity_col]] == 1),
      overall_mean = mean(tissue_data[[activity_col]]),
      sample_size = total_count))
  }

  prop_zeros <- calculate_proportions(tissue_data, activity_col)[1]
  prop_ones <- calculate_proportions(tissue_data, activity_col)[2]
  if (prop_zeros >= 0.90 || prop_ones >= 0.90 ||
      sd(tissue_data[[activity_col]]) <= 0.05) {
      mu = mean(tissue_data[[activity_col]][tissue_data[[activity_col]] > 0 & tissue_data[[activity_col]] < 1])
      if(is.na(mu)){mu = mean(tissue_data[[activity_col]])}
      p0 = mean(tissue_data[[activity_col]] == 0)
      p1 = mean(tissue_data[[activity_col]] == 1)

      return(list(
        status = "Skip_fit",
        reason = "insufficient_diversity",
        converged = FALSE,
        model_type = NA,
        mu = mu,
        sigma = NA,
        p0 = p0,
        p1 = p1,
        overall_mean = mu*(1-p1-p0) + p1,
        sample_size = total_count
    ))
  }

  # fit model
  model_type <- select_model(tissue_data, activity_col)
  formula_base <- as.formula(paste(activity_col, " ~ 1"))

  # Try the fit
  model_base <- try(fit_gamlss_model(tissue_data, formula_base, model_type), silent = TRUE)

  # Check if model_base is NULL or has error
  if (is.null(model_base) || inherits(model_base, "try-error")) {
    mu = mean(tissue_data[[activity_col]][tissue_data[[activity_col]] > 0 & tissue_data[[activity_col]] < 1])
    if(is.na(mu)){mu = mean(tissue_data[[activity_col]])}
    p0 = mean(tissue_data[[activity_col]] == 0)
    p1 = mean(tissue_data[[activity_col]] == 1)

    return(list(
      status = "Skip_fit",
      reason = "Error or warning in model fitting",
      converged = FALSE,
      model_type = model_type,
      mu = mu,
      sigma = NA,
      p0 = p0,
      p1 = p1,
      overall_mean = mu*(1-p1-p0) + p1,
      sample_size = nrow(tissue_data)
    ))
  }

  # Extract the parameters
  fit <- model_base
  mu <- (exp(fit$mu.coefficients) / (1 + exp(fit$mu.coefficients)))[['(Intercept)']]

 if(model_type == "BEINF"){
    sigma <- (exp(fit$sigma.coefficients)/(1+exp(fit$sigma.coefficients)))[['(Intercept)']]

    nu <- exp(fit$nu.coefficients)[['(Intercept)']]
    tau <- exp(fit$tau.coefficients) [['(Intercept)']]

    p2 <- 1/(1+(nu+tau))
    p0 <- nu*p2
    p1 <- tau*p2
  }else if(model_type == "BEZI"){
    sigma <- exp(fit$sigma.coefficients)[['(Intercept)']]
    nu <- (exp(fit$nu.coefficients)/(1+exp(fit$nu.coefficients)))[['(Intercept)']]
    p0 <- nu
    p1 <- 0
  }else if(model_type == "BEOI"){
    sigma <- exp(fit$sigma.coefficients)[['(Intercept)']]
    nu <- (exp(fit$nu.coefficients)/(1+exp(fit$nu.coefficients)))[['(Intercept)']]
    p0 <- 0
    p1 <- nu
  }else if(model_type == "BE"){
    sigma <- (exp(fit$sigma.coefficients)/(1+exp(fit$sigma.coefficients)))[['(Intercept)']]
    p0 <- 0
    p1 <- 0
  }

  # Return the results
  return(
    list(
    status = "Success",
    reason = NA,
    converged = model_base$converged,
    model_type = model_type,
    mu = mu,
    sigma = sigma,
    p0 = p0,
    p1 = p1,
    overall_mean = mu*(1-p1-p0) + p1,
    para_names = model_base$parameters,
    sample_size = nrow(tissue_data)
    )
  )
}


# fit regression function===============================================================================
fit_gamlss_regression <- function(data, activity_col, tissue_name, tissue_col, factors, target) {
  cat("Processing:", tissue_name, activity_col, "\n")
  tissue_data <- data[data[[tissue_col]] == tissue_name, ]
  if(!(target %in% colnames(data))){return(NULL)}
  if(!(target %in% factors)){factors <- c(factors, target)}
  tissue_data <- cbind(tissue_data[,c(factors)], tissue_data[activity_col])

  tissue_data <- na.omit(tissue_data)
  total_count <- nrow(tissue_data)

  if (total_count < 20 || nrow(tissue_data) < length(factors)) {
    cat("Skip:", tissue_name, activity_col, "\n")
    return(NULL)
  }

  if (length(unique(tissue_data[[activity_col]])) == 1) {
    return(NULL)
  }

  # fit model
  model_type <- select_model(tissue_data, activity_col)

  # regression
  prop_zeros <- calculate_proportions(tissue_data, activity_col)[1]
  prop_ones <- calculate_proportions(tissue_data, activity_col)[2]
  if (prop_zeros >= 0.95 || prop_ones >= 0.95 ||
      sd(tissue_data[[activity_col]]) <= 0.05) {
    return(NULL)
  }


  formula_with <- as.formula(paste(activity_col, " ~ ",paste(factors, collapse = "+")))
  model_with <- fit_gamlss_model(tissue_data, formula_with, model_type)

  formula_without <- as.formula(paste(activity_col, " ~ ",paste(setdiff(factors, target), collapse = "+")))
  model_without <- fit_gamlss_model(tissue_data, formula_without, model_type)

  # LR test
    lr_result <- tryCatch({
        capture.output(LR.test(model_without, model_with))
      }, error = function(e) {
          return(list(error = e$message, pvalue = NA))
      })

      if (is.list(lr_result) && !is.null(lr_result$error)) {
        lr_result <- 0.9999
      } else {
        lr_result <- as.numeric(gsub(".*p-value= (.*)","\\1", lr_result[length(lr_result)]))
      }
  lr_result

  # effect regression
  model_summary <- as.data.frame(summary(model_with))
  para_names <- setdiff(model_with$parameters, "sigma")
  if(model_type == "BEOI"){
    para_names <- gsub("nu", "tau", para_names)
  }

  pvalues <- model_summary[grep(target, row.names(model_summary)), "Pr(>|t|)"]
  names(pvalues) <- paste0(para_names,"_",target, '_pvalue')

  estimate <- model_summary[grep(target, row.names(model_summary)), "Estimate"]
  names(estimate) <- paste0(para_names,"_",target, '_estimate')

  se <- model_summary[grep(target, row.names(model_summary)), "Std. Error"]
  names(se) <- paste0(para_names,"_",target, '_SE')

  # Return the results
  return(
    list(
      target = target,
      converged = model_with$converged,
      model_type = model_type,
      p0=mean(tissue_data[activity_col] == 0),
      p1=mean(tissue_data[activity_col] == 1),
      lr_pvalue = lr_result,
      pvalues = pvalues,
      estimate = estimate,
      se = se,
      para_names = para_names,
      sample_size = total_count
      )
  )
}


# QTL regression on the raw data=========================================================================
### [OPTIMIZED] Add the inner_cores parameter to support inner-level SNP parallelization
### [OPTIMIZED] Use data.table to optimize the data operations
### [OPTIMIZED] Initialize model_with from model_without to speed up convergence
orig_QTL_regression <- function(data, info, tissue_col=NULL, flank=50000, tissue_name=NULL, activity_col, factors, vcffile, sd_th=0.05, inner_cores=1, convert_prob=FALSE, use_lrt=FALSE, snp_start=1L, snp_end=0L) {
  cat("Processing:", tissue_name, activity_col, "\n")
  if(is.null(tissue_name)){
    tissue_data <- na.omit(data[c(activity_col, "sample")])
  }else{
    tissue_data <- data[data[[tissue_col]] == tissue_name, ]
    tissue_data <- cbind(tissue_data[,-grep("PromoterID", colnames(tissue_data))], tissue_data[activity_col])
    tissue_data <- na.omit(tissue_data)
  }

  if (nrow(tissue_data) < 20 || nrow(tissue_data) < length(factors)) {
    cat("Skip:", tissue_name, activity_col, "\n")
    return(NULL)}

    # Check the data
    prop_zeros <- calculate_proportions(tissue_data, activity_col)[1]
    prop_ones <- calculate_proportions(tissue_data, activity_col)[2]
    if (prop_zeros >= 0.95 || prop_ones >= 0.95 || sd(tissue_data[[activity_col]]) <= sd_th) {
      cat("Skip:", tissue_name, activity_col, "\n")
      return(NULL)
    }

    if (prop_zeros <= 0.05) {
      tissue_data <- subset(tissue_data, tissue_data[[activity_col]] > 0)
    }
    if (prop_ones <= 0.05) {
      tissue_data <- subset(tissue_data, tissue_data[[activity_col]] < 1)
    }

    model_type <- select_model(tissue_data, activity_col)

    .factors_rhs <- if(length(factors) == 0) "1" else paste(factors, collapse = "+")
    formula_without <- as.formula(paste(activity_col, " ~ ", .factors_rhs))
    model_without <- fit_gamlss_model(tissue_data, formula_without, model_type)

    # Pre-compute null model statistics for probability-scale conversion
    if(convert_prob) {
      obs_vals_conv <- tissue_data[[activity_col]]
      obs_p0_conv <- mean(obs_vals_conv == 0)
      obs_p1_conv <- mean(obs_vals_conv == 1)
      obs_mid_conv <- obs_vals_conv[obs_vals_conv > 0 & obs_vals_conv < 1]
      obs_mu_conv <- if(length(obs_mid_conv) > 0) mean(obs_mid_conv) else 0.5
      obs_pmid_conv <- pmax(1 - obs_p0_conv - obs_p1_conv, 0.001)
    }

    # Extract the current TSS information
    chr <- as.character(info[activity_col,]$seqnames)
    pos <- as.character(info[activity_col,]$start)

    # Extract the genotype data of the current flank
    geno_numeric <- get_geno(chr, pos, flank, vcffile)
    snplist <- colnames(geno_numeric)
    if(is.null(geno_numeric)){return(NULL)}
    geno_numeric <- geno_numeric[row.names(geno_numeric) %in% tissue_data$sample,]
    geno_numeric <- as.data.frame(geno_numeric)
    colnames(geno_numeric) <- snplist

        geno_stats <- t(apply(geno_numeric, 2, function(x) {
        # Count the three genotypes
        n0 = sum(x == 0)
        n1 = sum(x == 1)
        n2 = sum(x == 2)

        # Compute MAF
        total = n0 + n1 + n2
        maf = (n1 + 2*n2) / (2*total)

        return(c(
          n_total = total,
          n_geno0 = n0,
          n_geno1 = n1,
          n_geno2 = n2,
          maf = maf,
          min_geno_count = min(n0, n1, n2)
        ))
      }))

      # Convert to a data frame
      geno_stats <- as.data.frame(geno_stats)
      geno_stats_filter <- subset(geno_stats, min_geno_count >= 5)
    geno_numeric <- geno_numeric[colnames(geno_numeric) %in% row.names(geno_stats_filter)]


    geno_numeric$sample <- row.names(geno_numeric)

    geno_data <- merge(tissue_data, geno_numeric, by="sample", all=F)

    if(!is.null(tissue_name)){
      snplist <- colnames(geno_data)[grep("b38", colnames(geno_data))]
    }else{
      snplist <- colnames(geno_data)[grep("QTL", colnames(geno_data))]
    }
    # SNP subsetting for batch mode
    if(snp_end == 0L) snp_end <- length(snplist)
    snp_start <- max(1L, min(snp_start, length(snplist)))
    snp_end   <- max(snp_start, min(snp_end, length(snplist)))
    cat("Total SNPs:", length(snplist), "| Batch:", snp_start, "-", snp_end, "\n")
    snplist <- snplist[snp_start:snp_end]
    cat("Running", length(snplist), "SNPs for:", activity_col, "\n")

    # [OPTIMIZED] Pre-build the data.table
    snp_cols <- snplist
    non_snp_cols <- setdiff(colnames(geno_data), snp_cols)
    full_dt <- as.data.table(geno_data)

    # [OPTIMIZED] Pre-compute the formula strings
    factors_formula_str <- if(length(factors) == 0) NA_character_ else paste(factors, collapse = " + ")

    # [OPTIMIZED] Define the single-SNP test function
    test_single_snp <- function(snp) {
      tryCatch({
        # [OPTIMIZED] Select columns from the data.table
        cols_needed <- c(non_snp_cols, snp)
        snp_data <- as.data.frame(full_dt[, ..cols_needed])
        snp_data <- na.omit(snp_data)

        # If there is too little data, skip this SNP
        if (nrow(snp_data) < 10) {
          return(NULL)
        }

        formula_with <- if(is.na(factors_formula_str))
          as.formula(paste(activity_col, " ~ ", snp))
        else
          as.formula(paste(activity_col, " ~ ", factors_formula_str, " + ", snp))

        # [OPTIMIZED] Initialize from model_without to speed up convergence
        model_with <- tryCatch({
          fit_gamlss_model(snp_data, formula_with, model_type, start.from = model_without)
        }, error = function(e) {
          return(NULL)
        })

        # If the model fit fails, skip this SNP
        if (is.null(model_with)) {
          return(NULL)
        }

        #LR test
        lr_result <- tryCatch({
          capture.output(LR.test(model_without, model_with))
        }, error = function(e) {
            return(list(error = e$message, pvalue = NA))
        })

        if (is.list(lr_result) && !is.null(lr_result$error)) {
          lrt_pvalue <- 0.9999
        } else {
          lrt_pvalue <- as.numeric(gsub(".*p-value= (.*)","\\1", lr_result[length(lr_result)]))
        }

      # Extract the results
      tryCatch({
        model_summary <- as.data.frame(summary(model_with))
        para_names <- model_with$parameters
        # Get the dof information
        dof1 <- length(coef(model_with))  # Number of model parameters
        dof2 <- model_with$df.res  # Residual degrees of freedom


        snp_pvalues <- model_summary[grep("b38", row.names(model_summary)), "Pr(>|t|)"]
        names(snp_pvalues) <- paste0(para_names, '_pvalue')
        snp_estimates <- model_summary[grep("b38", row.names(model_summary)), "Estimate"]
        names(snp_estimates) <- paste0(para_names, '_estimate')
        snp_SE <- model_summary[grep("b38", row.names(model_summary)), "Std. Error"]
        names(snp_SE) <- paste0(para_names, '_SE')

        mu_pvalue <- snp_pvalues['mu_pvalue']
        sigma_pvalue <- snp_pvalues['sigma_pvalue']
        mu_estimate <- snp_estimates['mu_estimate']
        sigma_estimate <- snp_estimates['sigma_estimate']
        mu_SE <- snp_SE['mu_SE']
        sigma_SE <- snp_SE['sigma_SE']

        if(model_type == "BEINF"){
          nu_pvalue <- snp_pvalues['nu_pvalue']
          tau_pvalue <- snp_pvalues['tau_pvalue']
          nu_estimate <- snp_estimates['nu_estimate']
          tau_estimate <- snp_estimates['tau_estimate']
          nu_SE <- snp_SE['nu_SE']
          tau_SE <- snp_SE['tau_SE']
          }else if(model_type == "BEZI"){
          nu_pvalue <- snp_pvalues['nu_pvalue']
          tau_pvalue <- NA
          nu_estimate <- snp_estimates['nu_estimate']
          tau_estimate <- NA
          nu_SE <- snp_SE['nu_SE']
          tau_SE <- NA
          }else if(model_type == "BEOI"){
          nu_pvalue <- NA
          tau_pvalue <- snp_pvalues['nu_pvalue']
          nu_estimate <- NA
          tau_estimate <- snp_estimates['nu_estimate']
          nu_SE <- NA
          tau_SE <- snp_SE['nu_SE']
          }else if(model_type == "BE"){
          nu_pvalue <- NA
          tau_pvalue <- NA
          nu_estimate <- NA
          tau_estimate <- NA
          nu_SE <- NA
          tau_SE <- NA
          }

        # === Per-parameter LRT (additive: keep Wald, add *_lrt columns) ====
        mu_lrt_pvalue    <- NA_real_
        sigma_lrt_pvalue <- NA_real_
        nu_lrt_pvalue    <- NA_real_
        tau_lrt_pvalue   <- NA_real_
        if(isTRUE(use_lrt)) {
          lrt_pp <- compute_per_param_lrt(snp_data, activity_col, factors,
                                          snp, model_type, model_with)
          mu_lrt_pvalue    <- lrt_pp$mu
          sigma_lrt_pvalue <- lrt_pp$sigma
          if(model_type == "BEINF") {
            nu_lrt_pvalue  <- lrt_pp$nu
            tau_lrt_pvalue <- lrt_pp$tau
          } else if(model_type == "BEZI") {
            nu_lrt_pvalue  <- lrt_pp$nu
          } else if(model_type == "BEOI") {
            tau_lrt_pvalue <- lrt_pp$nu  # BEOI gamlss "nu" is reported as tau
          }
        }

        AIC_without <- AIC(model_without)
        AIC_with <- AIC(model_with)

        # Probability-scale conversion
        mu_prob_effect <- NA_real_
        p0_prob_effect <- NA_real_
        p1_prob_effect <- NA_real_

        if(convert_prob) {
          # Compute mean_gt for this SNP from the current data
          snp_dosage <- snp_data[[snp]]
          snp_mean_gt <- mean(snp_dosage, na.rm = TRUE)

          # mu conversion (logit link, all models) with mean_gt correction
          if(!is.na(mu_estimate)) {
            obs_mu_c <- pmin(pmax(obs_mu_conv, 0.001), 0.999)
            mu_int_obs <- .logit(obs_mu_c)
            mu_int_gt0 <- mu_int_obs - mu_estimate * snp_mean_gt
            mu_prob_effect <- .inv_logit(mu_int_gt0 + 2 * mu_estimate) - .inv_logit(mu_int_gt0)
          }

          # p0/p1 conversion (depends on model_type) with mean_gt correction
          if(model_type == "BEINF" && !is.na(nu_estimate) && !is.na(tau_estimate)) {
            nu_null <- pmax(obs_p0_conv, 0.001) / obs_pmid_conv
            tau_null <- pmax(obs_p1_conv, 0.001) / obs_pmid_conv
            nu_int_obs <- log(nu_null); tau_int_obs <- log(tau_null)
            nu_int_gt0 <- nu_int_obs - nu_estimate * snp_mean_gt
            tau_int_gt0 <- tau_int_obs - tau_estimate * snp_mean_gt
            nu_0 <- exp(nu_int_gt0); tau_0 <- exp(tau_int_gt0)
            nu_2 <- exp(nu_int_gt0 + 2 * nu_estimate); tau_2 <- exp(tau_int_gt0 + 2 * tau_estimate)
            p0_prob_effect <- nu_2 / (1 + nu_2 + tau_2) - nu_0 / (1 + nu_0 + tau_0)
            p1_prob_effect <- tau_2 / (1 + nu_2 + tau_2) - tau_0 / (1 + nu_0 + tau_0)
          } else if(model_type == "BEZI" && !is.na(nu_estimate)) {
            p0_c <- pmin(pmax(obs_p0_conv, 0.001), 0.999)
            nu_int_obs <- .logit(p0_c)
            nu_int_gt0 <- nu_int_obs - nu_estimate * snp_mean_gt
            p0_prob_effect <- .inv_logit(nu_int_gt0 + 2 * nu_estimate) - .inv_logit(nu_int_gt0)
            p1_prob_effect <- 0
          } else if(model_type == "BEOI" && !is.na(tau_estimate)) {
            p1_c <- pmin(pmax(obs_p1_conv, 0.001), 0.999)
            tau_int_obs <- .logit(p1_c)
            tau_int_gt0 <- tau_int_obs - tau_estimate * snp_mean_gt
            p1_prob_effect <- .inv_logit(tau_int_gt0 + 2 * tau_estimate) - .inv_logit(tau_int_gt0)
            p0_prob_effect <- 0
          } else if(model_type == "BE") {
            p0_prob_effect <- 0
            p1_prob_effect <- 0
          }
        }

        # return results
        return(
          list(snp=snp,
            model_type=model_type,
            converged = model_with$converged,
            lrt_pvalue=lrt_pvalue,
            mu_pvalue=mu_pvalue,
            sigma_pvalue=sigma_pvalue,
            nu_pvalue=nu_pvalue,
            tau_pvalue=tau_pvalue,
            mu_lrt_pvalue=mu_lrt_pvalue,
            sigma_lrt_pvalue=sigma_lrt_pvalue,
            nu_lrt_pvalue=nu_lrt_pvalue,
            tau_lrt_pvalue=tau_lrt_pvalue,
            mu_estimate=mu_estimate,
            mu_prob_effect=mu_prob_effect,
            sigma_estimate=sigma_estimate,
            nu_estimate=nu_estimate,
            p0_prob_effect=p0_prob_effect,
            tau_estimate=tau_estimate,
            p1_prob_effect=p1_prob_effect,
            mu_SE=mu_SE,
            sigma_SE=sigma_SE,
            nu_SE=nu_SE,
            tau_SE=tau_SE,
            AIC_without=AIC_without,
            AIC_with=AIC_with,
            AIC_diff=AIC_with-AIC_without,
            dof1=dof1,
            dof2=dof2
          ))
          }, error = function(e) {
          return(NULL)
        })
      }, error = function(e) {
        # If any error occurs, return NULL and print a warning
        warning(paste("Error processing SNP", snp, ":", e$message))
        return(NULL)
      })
    }

  # [OPTIMIZED] Choose serial or parallel execution according to the inner_cores parameter
  if (inner_cores > 1 && length(snplist) > 1) {
    # Parallel version: use mclapply
    n_cores_use <- min(inner_cores, length(snplist))
    cat("Using", n_cores_use, "cores for inner parallel\n")
    orig_results <- mclapply(snplist, test_single_snp, mc.cores = n_cores_use)
  } else {
    # Serial version: use lapply
    orig_results <- lapply(snplist, test_single_snp)
  }

  # Remove NULL results and combine the data frames
  orig_results <- Filter(Negate(is.null), orig_results)
  if (length(orig_results) > 0) {
    orig_results_df <- do.call(rbind, lapply(orig_results, as.data.frame))
  } else {
    orig_results_df <- data.frame() # Return an empty data frame
  }
  cat("Done:", activity_col, "\n")
  return(orig_results_df)
}

###conditional independent===================================================================
### [OPTIMIZED] Use data.table to optimize the data.frame operations
### [OPTIMIZED] Add the inner_cores parameter to support inner-level SNP parallelization
conditional_QTL_regression <- function(data, info, tissue_col=NULL, flank=500000, tissue_name=NULL, activity_col, factors, vcffile, orig_results_list, sd_th=0.05, inner_cores=1) {
  cat("Processing:", tissue_name, activity_col, "\n")
  if(is.null(tissue_name)){
    tissue_data <- na.omit(data)
  }else{
    tissue_data <- data[data[[tissue_col]] == tissue_name, ]
    tissue_data <- cbind(tissue_data[,-grep("PromoterID", colnames(tissue_data))], tissue_data[activity_col])
    tissue_data <- na.omit(tissue_data)
  }

  if (nrow(tissue_data) < 10 || nrow(tissue_data) < length(factors)) {
      cat("Skip:", tissue_name, activity_col, "\n")
      return(NULL)
  }


  # Check the data
    prop_zeros <- calculate_proportions(tissue_data, activity_col)[1]
    prop_ones <- calculate_proportions(tissue_data, activity_col)[2]
    if (prop_zeros >= 0.95 || prop_ones >= 0.95 || sd(tissue_data[[activity_col]]) <= sd_th) {
      cat("Skip:", tissue_name, activity_col, "\n")
      return(NULL)
    }

    if (prop_zeros <= 0.05) {
      tissue_data <- subset(tissue_data, tissue_data[[activity_col]] > 0)
    }
    if (prop_ones <= 0.05) {
      tissue_data <- subset(tissue_data, tissue_data[[activity_col]] < 1)
    }

    model_type <- select_model(tissue_data, activity_col)

    #formula_without <- as.formula(paste(activity_col, " ~ ",paste(factors, collapse = "+")))
    #model_without <- fit_gamlss_model(tissue_data, formula_without, model_type)

    # Extract the current TSS information
    chr <- as.character(info[activity_col,]$seqnames)
    pos <- as.character(info[activity_col,]$start)

    # Extract the genotype data of the current flank
    geno_numeric <- get_geno(chr, pos, flank, vcffile)
    geno_numeric <- geno_numeric[row.names(geno_numeric) %in% tissue_data$sample,]
    geno_numeric <- as.data.frame(geno_numeric)

    geno_stats <- t(apply(geno_numeric, 2, function(x) {
        # Count the three genotypes
        n0 = sum(x == 0)
        n1 = sum(x == 1)
        n2 = sum(x == 2)

        # Compute MAF
        total = n0 + n1 + n2
        maf = (n1 + 2*n2) / (2*total)

        return(c(
          n_total = total,
          n_geno0 = n0,
          n_geno1 = n1,
          n_geno2 = n2,
          maf = maf,
          min_geno_count = min(n0, n1, n2)
        ))
      }))

      # Convert to a data frame
      geno_stats <- as.data.frame(geno_stats)
      geno_stats_filter <- subset(geno_stats, min_geno_count >= 5)
    geno_numeric <- geno_numeric[colnames(geno_numeric) %in% row.names(geno_stats_filter)]

    geno_numeric$sample <- row.names(geno_numeric)

    geno_data <- merge(tissue_data, geno_numeric, by="sample", all=F)



### conditional regression
# 1. Use the existing raw results directly, without recomputing
  current_results <- orig_results_list[[paste(tissue_name, activity_col, sep="_")]]
  if(is.null(current_results) || nrow(current_results) == 0) {
    return(NULL)
  }
  current_results <- current_results[current_results$snp %in% row.names(geno_stats_filter),]
  current_results <- subset(current_results, converged == TRUE)

  # 2. Get the initial threshold
  threshold <- get_threshold(orig_results_list)
  if(is.null(threshold)){
    return(NULL)
  }

  # 3. Find the most significant SNP in the raw results
  best_snp_info <- current_results[which.min(current_results$lrt_pvalue), ]

  # Check whether it is significant
  if(best_snp_info$lrt_pvalue > threshold) {
    cat("No significant SNPs for", activity_col, "\n")
    return(list(
      tissue = tissue_name,
      promoter = activity_col,
      significant_snps = character(0)
    ))
  }

  # 4. Start the conditional analysis
  significant_snps <- c(best_snp_info$snp)

  # original threshold
  update_results_list <- orig_results_list

  # [OPTIMIZED] Extract the data of the current tissue and promoter, using data.table
  snp_cols <- grep("b38", colnames(geno_data), value = TRUE)
  genotypes <- as.matrix(geno_data[, snp_cols, drop = FALSE])
  row.names(genotypes) <- geno_data$sample

  phenotype <- geno_data[[activity_col]]
  names(phenotype) <- geno_data$sample
  covariates <- as.data.frame(geno_data[, factors, drop = FALSE])

  # [OPTIMIZED] Pre-build the complete data.table (only once!)
  full_dt <- data.table(
    pheno = phenotype,
    covariates,
    genotypes
  )
  setnames(full_dt, "pheno", activity_col)

  # [OPTIMIZED] Pre-compute the factor part of the formula (avoids repeated paste)
  factors_formula_str <- paste(factors, collapse = " + ")

  repeat {
    # [OPTIMIZED] Simplified formula construction
    formula_base <- as.formula(paste(activity_col, " ~ ",
                                  factors_formula_str, " + ",
                                  paste(significant_snps, collapse = " + ")))

    # Get the remaining SNPs (excluding those already selected)
    remaining_snps <- setdiff(snp_cols, significant_snps)
    if(length(remaining_snps) == 0) break

    # [OPTIMIZED] Select columns from the pre-built data.table
    base_cols <- c(activity_col, factors, significant_snps)
    base_data <- as.data.frame(full_dt[, ..base_cols])

    # Fit the base model (only needed once!)
    model_base <- fit_gamlss_model(base_data, formula_base, model_type)
    if(is.null(model_base)) {
      warning("Base model fitting failed")
      break
    }

    # [OPTIMIZED] Pre-compute the fixed columns of test_data (once per iteration)
    fixed_cols <- c(activity_col, factors, significant_snps)

    # [OPTIMIZED] Define the single-SNP test function
    test_single_snp <- function(snp_name) {
      tryCatch({
        # [OPTIMIZED] Select columns from the pre-built data.table (much faster than creating a data.frame each time)
        test_cols <- c(fixed_cols, snp_name)
        test_data <- as.data.frame(full_dt[, ..test_cols])

        # Full model: add the SNP currently tested
        formula_full <- as.formula(paste(activity_col, " ~ ",
                                      factors_formula_str, " + ",
                                      paste(c(significant_snps, snp_name), collapse = " + ")))

        # [OPTIMIZED] Initialize the full model from the base model to speed up convergence
        model_full <- fit_gamlss_model(test_data, formula_full, model_type, start.from = model_base)

        # LR test: compare nested models
        lr_result <- tryCatch({
          capture.output(LR.test(model_base, model_full))
        }, error = function(e) {
          return(list(error = e$message, pvalue = NA))
        })

        if (is.list(lr_result) && !is.null(lr_result$error)) {
          lr_pvalue <- 0.9999
        } else {
          lr_pvalue <- as.numeric(gsub(".*p-value= (.*)","\\1", lr_result[length(lr_result)]))
        }

        return(lr_pvalue)
      }, error = function(e) {
        warning(paste("Error testing SNP", snp_name, ":", e$message))
        return(0.9999)
      })
    }

    # [OPTIMIZED] Choose serial or parallel execution according to the inner_cores parameter
    if (inner_cores > 1 && length(remaining_snps) > 1) {
      # Parallel version: use mclapply
      n_cores_use <- min(inner_cores, length(remaining_snps))
      pvalues_list <- mclapply(remaining_snps, test_single_snp, mc.cores = n_cores_use)
      pvalues <- unlist(pvalues_list)
    } else {
      # Serial version: use sapply
      pvalues <- sapply(remaining_snps, test_single_snp)
    }

    names(pvalues) <- remaining_snps


    # Find the most significant SNP
    min_pvalue <- min(pvalues, na.rm = TRUE)
    best_snp <- names(which.min(replace(pvalues, is.na(pvalues), Inf)))

    cat("Iteration", length(significant_snps), ": Best SNP =", best_snp,
        ", p-value =", min_pvalue, "\n")

    # Decide whether it is significant
    if (min_pvalue > threshold) {
      break
    }

    # Update the threshold
    fp <- as.data.frame(pvalues)
    colnames(fp) <- "lrt_pvalue"
    fp$snp <- rownames(fp)
    update_results_list[[activity_col]] <- fp


    new_threshold <- get_threshold(update_results_list)
    if(!is.null(new_threshold)) {
      threshold <- new_threshold
    }

    # Add the most significant SNP to the list of significant SNPs
    significant_snps <- c(significant_snps, best_snp)

    # Optional: set a maximum number of SNPs
    if(length(significant_snps) >= 5) {
      cat("Reached maximum number of SNPs (5)\n")
      break
    }


  }

  ###
  return(
  list(tissue = tissue_name,
        promoter = activity_col,
        significant_snps = significant_snps,
        n_signals = length(significant_snps))
  )
}
