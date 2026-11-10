## Conditional analysis: independent signals of significant cis units.
##
## Ported from the analysis code of the BOHR paper (Figure4_puQTL_landscape_eQTL_coloc/rerun_v2/code):
##   run_conditional_v3c.R (md5 8770f799...; its code is that of run_conditional_v3m.R, md5 7db756ea..., which differs
##   only in the unit list): per-promoter inputs :92-100 (= run_final2.R), mid_rank :106, test_one :107-114, run_pass
##   :116-134, forward rounds :137-164, leave-one-out scans :166-174, signal / unit tables :175-180;
##   merge_conditional_v3c.R :30-34 (= merge_conditional_v3m.R :48-52): the backward step;
##   mirror_abs_labels.R (md5 5e8c0fe6...) :48-50 (fit), :55-66 (design), :71-75 (labels) and merge_conditional_v3c.R :48-51
##   (renaming): the optional absolute-activity labels of exact-pair boundary events.
## The 14,800 signals of the paper are the v3c tables: the v3m forward/backward results of every unit, with the two units
## whose stored conditional Firth fits had failed refitted by the same code (v3c). This file implements that code once.
## The engine functions (snp_tests, firth_snp_tests, promoter_null, orth_design, acat_log) are the package's ports.

## ---- one variant against a conditional null (test_one, run_conditional_v3c.R:107-114) -------------------
## mu: F-corrected beta-part LRT, p = rank of the conditional design on the intermediate samples + 1 (the variant);
## p0/p1: Firth-penalised LRT (unpenalised LRT with control$firth = FALSE). Returns est, se, log_p, log_p_chi.
cond_test <- function(type, y, Cc, g, nul, n_mid, p_mean, control) {
  if (type == "mu") {
    r <- snp_tests(y, Cc, g, nul, maxit = control$maxit)
    return(c(est = unname(r$mu_est), se = unname(r$mu_se), log_p = mu_F_logp(r$mu_stat, n_mid, p_mean), log_p_chi = unname(r$mu_lrt_logp)))
  }
  pre <- if (type == "p0") "nu" else "tau"
  if (isTRUE(control$firth)) {
    r <- firth_snp_tests(y, Cc, g, nul, maxround = control$max_rounds, maxit = control$maxit, bfgs_maxit = control$bfgs_maxit)
    c(est = unname(r[[paste0(pre, "_firth_est")]]), se = unname(r[[paste0(pre, "_firth_se")]]), log_p = unname(r[[paste0(pre, "_firth_logp")]]),
      log_p_chi = NA_real_)
  } else {
    r <- snp_tests(y, Cc, g, nul, maxit = control$maxit)
    c(est = unname(r[[paste0(pre, "_est")]]), se = unname(r[[paste0(pre, "_se")]]), log_p = unname(r[[paste0(pre, "_lrt_logp")]]), log_p_chi = NA_real_)
  }
}

## ---- one pass (run_pass, run_conditional_v3c.R:116-134) ------------------------------------------------------
## jobs: list of (key, type, cs = conditioning variants, test = variants to test). For each job, the null is refitted with
## the conditioning variants appended to the covariates (plan B: they enter the beta mean and both boundary equations),
## then every variant of `test` is tested. Variants run in chunks, heaviest first; the rows of a job keep the order of
## `test`. Returns res (one data frame per job, NULL when its null fit failed), dropped (collinear columns) and null_ok.
cond_pass <- function(jobs, PP, control, ncores, chunk_size) {
  nl <- par_apply(seq_along(jobs), function(b) tryCatch({
    jb <- jobs[[b]]; P <- PP[[jb$key]]
    Cc <- orth_design(cbind(P$C, P$G[, jb$cs, drop = FALSE]))
    mid <- P$y > 0 & P$y < 1
    list(Cc = Cc, nul = promoter_null(P$y, Cc, P$family, maxit = control$maxit), dropped = length(attr(Cc, "dropped")),
         n_mid = sum(mid), p_mean = rank_mid(Cc, mid) + 1L)
  }, error = function(e) NULL), ncores, preschedule = FALSE)
  tasks <- list()
  for (b in seq_along(jobs)) {
    sn <- jobs[[b]]$test
    if (!length(sn) || is.null(nl[[b]])) next
    w <- c(mu = 1, p0 = 4, p1 = 4)[[jobs[[b]]$type]] * length(PP[[jobs[[b]]$key]]$y)
    ch <- split(seq_along(sn), ceiling(seq_along(sn) / chunk_size))
    for (ci in seq_along(ch)) tasks[[length(tasks) + 1L]] <- list(b = b, ci = ci, j = ch[[ci]], w = w * length(ch[[ci]]))
  }
  if (length(tasks)) tasks <- tasks[order(-vapply(tasks, `[[`, 0, "w"))]
  one <- function(t) {
    jb <- jobs[[t$b]]; P <- PP[[jb$key]]; N <- nl[[t$b]]; sn <- jb$test[t$j]
    v <- vapply(sn, function(s) tryCatch(cond_test(jb$type, P$y, N$Cc, P$G[, s], N$nul, N$n_mid, N$p_mean, control),
                                         error = function(e) c(est = NA_real_, se = NA_real_, log_p = NA_real_, log_p_chi = NA_real_)),
                numeric(4))
    list(b = t$b, ci = t$ci, d = data.frame(variant = sn, est = v[1, ], se = v[2, ], log_p = v[3, ], log_p_chi = v[4, ],
                                            stringsAsFactors = FALSE, row.names = NULL))
  }
  out <- par_apply(tasks, one, ncores, preschedule = FALSE)
  res <- vector("list", length(jobs))
  if (length(out)) {
    bb <- vapply(out, `[[`, 0, "b"); cc <- vapply(out, `[[`, 0, "ci")
    for (b in unique(bb)) { k <- which(bb == b); k <- k[order(cc[k])]; res[[b]] <- do.call(rbind, lapply(out[k], `[[`, "d")) }
  }
  list(res = res, dropped = vapply(nl, function(x) if (is.null(x)) NA_integer_ else as.integer(x$dropped), 0L),
       null_ok = !vapply(nl, is.null, TRUE))
}

#' Independent signals of significant cis units (conditional analysis)
#'
#' Finds the conditionally independent signals of every significant unit
#' (promoter and component) of a [bohr_cis()] result, with the forward and
#' backward steps of the BOHR paper, and returns each signal's lead variant
#' with its statistics.
#'
#' The steps (Methods Section 7; SI Section 7):
#' \enumerate{
#'   \item **Forward.** The first signal is the unit's leader, the variant with
#'     the smallest P value. In each round the selected variants are appended
#'     to the covariates (they enter the mean of the intermediate values and
#'     both boundary equations), the model without a test variant is refitted,
#'     and every other variant of the cis window that [bohr_cis()] tested is
#'     tested again with the unit's own test: mu with the F-corrected
#'     likelihood-ratio test, where the parameter count is the rank of the
#'     conditional design on the intermediate samples plus one; p0 and p1 with
#'     the Firth-penalised likelihood-ratio test. The conditional P values of
#'     the window are combined with ACAT ([acat()], with the cap of the
#'     [bohr_cis()] call). The variant with the smallest conditional P value is
#'     added as a new signal when this window P value is at most `threshold`;
#'     otherwise the unit stops. A unit stops at `max_signals` signals.
#'   \item **Backward.** For a unit with two or more signals, each signal is
#'     retested conditional on the others: every variant of the window except
#'     the other signals is tested, and the signal is kept when the ACAT P value
#'     of this leave-one-out scan is at most `threshold`. A unit's only signal
#'     is always kept. The leave-one-out scans are returned in `loo`; the
#'     statistics of each signal's own variant in its scan (`loo_*`) are the
#'     per-signal conditional statistics used for colocalization.
#' }
#' The threshold of a tissue is its discovery threshold: the largest unit
#' P value with q at or below the false-discovery rate of the [bohr_cis()]
#' call. Ties of the smallest P value go to the first variant in the order of
#' the rows of `genotype` (the VCF order in the paper).
#'
#' **Exact mirror pairs.** The units are those of [bohr_cis()] after its
#' mirror rule: a boundary event of an exact mirror pair is a p0 unit at the
#' promoter that goes silent. With `absolute`, each significant p0 unit of a
#' `boundary_exact` pair is labelled by absolute activity (SI Section 2,
#' "Which promoter a boundary event regulates"): the absolute activity of each
#' promoter of the pair is regressed (least squares) on the leader's dosage,
#' with the covariates and the unit's other kept signals; a promoter changes
#' when \eqn{|z| \ge} `abs_z`. When only the other promoter changes, the
#' event is renamed p1 at that promoter (its P values are those of the p0
#' test, the same test); otherwise it stays p0 at the silent promoter. In the
#' paper no puQTL was renamed. Without `absolute`, the labels of [bohr_cis()]
#' are kept (the paper's rule "v3m"), and `settings$mirror_labels` says so.
#'
#' **Inputs.** `phenotype`, `genotype`, `variants`, `promoters`,
#' `covariates`, `tissue` and `individual` must be those of the [bohr_cis()]
#' call, as objects or files (see [bohr_cis()]). Each promoter is prepared
#' again as in [bohr_cis()] (samples, endpoint rule, covariates, genotype
#' imputation over the tissue's samples, tested variants); an error is given
#' when the result differs from the [bohr_cis()] run. For cluster jobs, run
#' [bohr_cis()] one chromosome at a time, combine the parts with
#' [bohr_combine()] (which sets the tissue-wide q values and thresholds), and
#' run `bohr_conditional()` on the combined result one chromosome at a time
#' with `chromosomes`.
#'
#' @param cis A [bohr_cis()] result (one tissue, or several tissues from
#'   `bohr_cis(tissue = )` or [bohr_combine()]). Its significant units are
#'   analysed with its settings (window, genotype rule, fitting rules and ACAT
#'   cap).
#' @param phenotype,genotype,variants,promoters,covariates,tissue,individual
#'   The inputs of the [bohr_cis()] call, as R objects or file paths.
#' @param chromosomes Optional vector of chromosome names; only units of
#'   promoters on these chromosomes are analysed.
#' @param threshold `NULL` (default: the discovery threshold of each tissue,
#'   `cis$threshold`; Methods Section 7), or a number used for every tissue, or
#'   a vector named by tissue. A signal is added, and kept by the backward
#'   step, when its window ACAT P value is at most this value.
#' @param max_signals Largest number of signals per unit. Default 10 (Methods
#'   Section 7).
#' @param absolute Optional absolute activity for the labels of exact-pair
#'   boundary events: a numeric matrix with promoters in rows (row names =
#'   promoter IDs, covering both members of every exact pair with a
#'   significant boundary unit) and samples in columns (column names = sample
#'   IDs of `phenotype`), or the path of a TSV file in the same layout. No
#'   missing values on the samples used. In the paper: proActiv's absolute
#'   promoter activity.
#' @param abs_z A promoter's absolute activity changes when \eqn{|z| \ge}
#'   `abs_z`. Default 1.96 (SI Section 2).
#' @param ncores Number of cores (forked processes; 1 on Windows).
#' @param chunk_size Number of variants per parallel task (default 50, as in
#'   the paper's code).
#' @param out_dir Optional folder (created if needed). When given, the tables
#'   `signals`, `units`, `rounds`, `loo` and (with `absolute`) `labels` are also
#'   written to it as tab-separated files (`signals.tsv` and so on).
#' @return An object of class `"bohr_conditional"`, a list with
#'   \describe{
#'     \item{`signals`}{One row per signal, in selection order `k` within each
#'       unit: `promoter`, `gene`, `component`, `k`, `variant`, `pos`,
#'       `selected_round` (0 for the leader), `log_p` (natural-log P value of
#'       the variant when it was selected: marginal for the leader,
#'       conditional on the earlier signals otherwise), `log_p_chi` (mu only:
#'       the uncorrected chi-squared log P at selection), `acat_p` (the window
#'       ACAT P value that admitted it; the unit P value for the leader),
#'       `loo_est`, `loo_se`, `loo_log_p`, `loo_log_p_chi` (the signal's
#'       variant conditional on the unit's other signals; `NA` for a unit with
#'       one signal), `loo_acat_p` (ACAT P of the leave-one-out scan),
#'       `loo_top_variant` (its smallest-P variant) and `kept` (backward step).
#'       Estimates are on the scale of the model coefficient per allele (logit
#'       of mu; log odds of 0 or 1 against an intermediate value).}
#'     \item{`units`}{One row per analysed unit: `promoter`, `gene`,
#'       `component`, `chr`, `family`, `n`, `n0`, `n1`, `p` and `q` (from
#'       [bohr_cis()]), `leader`, `n_signals` (forward), `n_kept` (after the
#'       backward step) and `stop` (why the forward step stopped).}
#'     \item{`rounds`}{One row per unit and forward round: conditioning
#'       variants, variants tested and failed, the window ACAT P value (and,
#'       for mu, that of the uncorrected P values), the best variant and its
#'       log P, collinear columns dropped, and whether the conditional null
#'       fit succeeded.}
#'     \item{`loo`}{The leave-one-out scans of units with two or more signals:
#'       one row per signal `k` and tested variant, with `est`, `se`, `log_p`
#'       and `log_p_chi`.}
#'     \item{`labels`}{With `absolute`: one row per judged exact-pair boundary
#'       event (`S` = the silent promoter, `O` = its partner), with `z_S`,
#'       `z_O`, `outcome` and `label`; `changed` marks a renamed event.}
#'     \item{`threshold`, `settings`, `seconds`}{Thresholds used, settings
#'       (including `mirror_labels`) and elapsed time.}
#'   }
#'   With several tissues, every table starts with a column `tissue`.
#' @seealso [bohr_cis()], [bohr_combine()], [acat()].
#' @examples
#' set.seed(7)
#' n <- 200; s <- paste0("S", 1:n)
#' geno <- matrix(rbinom(40 * n, 2, 0.4), 40, n, dimnames = list(paste0("v", 1:40), s))
#' vars <- data.frame(id = rownames(geno), chr = "chr1", pos = seq(1000, 40000, 1000))
#' eta <- -0.2 + 0.5 * geno[10, ] - 0.5 * geno[30, ]                    # two causal variants
#' y <- rbeta(n, 8 * plogis(eta), 8 * (1 - plogis(eta)))
#' ph <- rbind(P1 = y); colnames(ph) <- s
#' proms <- data.frame(id = "P1", chr = "chr1", pos = 20000, gene = "G1")
#' cis <- bohr_cis(ph, geno, vars, proms, window = 30000, min_n_mid = 50)
#' cond <- bohr_conditional(cis, ph, geno, vars, proms)
#' cond
#' cond$signals[, c("component", "k", "variant", "log_p", "acat_p", "loo_log_p", "kept")]
#' @export
bohr_conditional <- function(cis, phenotype, genotype, variants = NULL, promoters, covariates = NULL, tissue = NULL, individual = NULL,
                             chromosomes = NULL, threshold = NULL, max_signals = 10L, absolute = NULL, abs_z = 1.96,
                             ncores = 1L, chunk_size = 50L, out_dir = NULL) {
  t0 <- proc.time()[["elapsed"]]
  if (!inherits(cis, "bohr_cis")) stop("`cis` must be a bohr_cis() result.", call. = FALSE)
  set <- cis$settings
  if (is.null(set$window) || is.null(set$control) || is.null(set$min_genotype_count)) stop("`cis` has no settings; was it built by bohr_cis()?", call. = FALSE)
  max_signals <- check_count(max_signals, "max_signals", min = 1)
  ncores <- check_count(ncores, "ncores", min = 1); chunk_size <- check_count(chunk_size, "chunk_size", min = 1)
  if (!is.numeric(abs_z) || length(abs_z) != 1L || !is.finite(abs_z) || abs_z <= 0) stop("`abs_z` must be a single positive number.", call. = FALSE)
  if (!is.null(out_dir) && (!is.character(out_dir) || length(out_dir) != 1L || is.na(out_dir) || !nzchar(out_dir)))
    stop("`out_dir` must be the path of one folder.", call. = FALSE)
  control <- set$control; acat_cap <- if (is.null(set$acat_cap)) TRUE else set$acat_cap
  tissues <- set$tissues
  multi <- !is.null(tissues)

  ## ---- inputs, read and checked as in bohr_cis() ----
  phenotype <- check_phenotype(read_matrix(phenotype, "phenotype"))
  samples <- colnames(phenotype)
  st <- sample_info(tissue, individual, samples)
  if (multi != !is.null(st$tissue))
    stop(if (multi) "`cis` covers several tissues: give the same `tissue` as in the bohr_cis() call." else
      "`cis` is a single-tissue result: `tissue` must be NULL, as in the bohr_cis() call.", call. = FALSE)
  if (multi && !all(tissues %in% st$levels)) stop("`tissue` does not cover the tissues of `cis` (e.g. ", setdiff(tissues, st$levels)[1], ").", call. = FALSE)
  prom <- check_promoters(read_table(promoters, "promoters"), rownames(phenotype))
  if (is_vcf_path(genotype)) {
    if (!is.null(variants)) stop("`variants` must be NULL when `genotype` is a VCF file: the positions are read from the VCF.", call. = FALSE)
    vc <- read_vcf_dosage(genotype, samples = unique(st$individual)); genotype <- vc$genotype; variants <- vc$variants
  } else {
    genotype <- read_matrix(genotype, "genotype")
    if (is.null(variants)) stop("`variants` is needed (a data frame or TSV file with columns id, chr and pos) unless `genotype` is a VCF file.", call. = FALSE)
    variants <- read_table(variants, "variants")
  }
  ind_given <- !identical(st$individual, samples)
  G_all <- check_genotype(genotype, unique(st$individual), individual_given = ind_given); genotype <- NULL
  V <- check_variants(variants, rownames(G_all))
  if (multi && is_path(covariates) && !dir.exists(covariates)) covariates <- read_matrix(covariates, "covariates")
  if (!multi) covariates <- read_covariates(covariates)
  if (!is.null(absolute)) {
    absolute <- read_matrix(absolute, "absolute")
    if (is.data.frame(absolute)) absolute <- as.matrix(absolute)
    if (!is.matrix(absolute) || !is.numeric(absolute) || is.null(rownames(absolute)) || is.null(colnames(absolute)))
      stop("`absolute` must be a numeric matrix with promoter IDs as row names and sample IDs as column names.", call. = FALSE)
  }

  ## ---- thresholds and units ----
  tlev <- if (multi) tissues else NA_character_
  thr <- cis$threshold
  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || !length(threshold) || anyNA(threshold) || any(threshold <= 0 | threshold > 1))
      stop("`threshold` must be a number in (0, 1], or one per tissue.", call. = FALSE)
    if (!multi) {
      if (length(threshold) != 1L) stop("`threshold` must be a single number for a single-tissue `cis`.", call. = FALSE)
      thr <- unname(threshold)
    } else if (length(threshold) == 1L && is.null(names(threshold))) {
      thr <- stats::setNames(rep(threshold, length(tissues)), tissues)
    } else {
      if (is.null(names(threshold)) || !all(tissues %in% names(threshold))) stop("`threshold` must be one number or be named by tissue.", call. = FALSE)
      thr <- threshold[tissues]
    }
  }
  U_all <- cis$units
  if (!is.null(U_all) && nrow(U_all)) {
    U_all <- U_all[U_all$significant %in% TRUE, , drop = FALSE]
    U_all$chr <- prom$chr[match(U_all$promoter, prom$id)]
    if (anyNA(U_all$chr)) stop("`promoters` has no row for a unit of `cis` (e.g. ", U_all$promoter[is.na(U_all$chr)][1], ").", call. = FALSE)
    if (!is.null(chromosomes)) U_all <- U_all[U_all$chr %in% as.character(chromosomes), , drop = FALSE]
  }

  res <- lapply(tlev, function(tt) {
    U <- if (is.null(U_all) || !nrow(U_all)) NULL else if (multi) U_all[U_all$tissue == tt, , drop = FALSE] else U_all
    th <- if (multi) unname(thr[[tt]]) else unname(thr)
    if (is.null(U) || !nrow(U) || !is.finite(th)) return(NULL)
    w <- if (multi) st$tissue == tt else rep(TRUE, length(samples))
    s <- samples[w]
    G <- G_all[, st$individual[w], drop = FALSE]; colnames(G) <- s
    Z <- if (multi) tissue_covariates(covariates, tt) else covariates
    snp_t <- if (is.null(cis$snp)) NULL else if (multi) cis$snp[cis$snp$tissue == tt, c("promoter", "variant"), drop = FALSE] else cis$snp[, c("promoter", "variant")]
    prom_t <- if (is.null(cis$promoters)) NULL else if (multi) cis$promoters[cis$promoters$tissue == tt, , drop = FALSE] else cis$promoters
    cond_one(phenotype[, s, drop = FALSE], G, V, prom, Z, U, th, snp_t, prom_t, set, control, acat_cap, max_signals,
             absolute, abs_z, ncores, chunk_size, label = if (multi) tt else NULL)
  })
  ## ---- bind tissues ----
  tab <- function(k) {
    d <- lapply(seq_along(res), function(i) { x <- res[[i]][[k]]
      if (is.null(x) || !nrow(x)) return(NULL)
      if (multi) data.frame(tissue = rep(tlev[i], nrow(x)), x, check.names = FALSE, stringsAsFactors = FALSE) else x })
    b <- bind_df(d); if (!is.null(b)) rownames(b) <- NULL; b
  }
  out <- list(signals = tab("signals"), units = tab("units"), rounds = tab("rounds"), loo = tab("loo"),
              labels = if (is.null(absolute)) NULL else tab("labels"),
              threshold = thr,
              settings = list(tissues = tissues, chromosomes = if (is.null(chromosomes)) sort(unique(U_all$chr)) else as.character(chromosomes),
                              max_signals = max_signals, acat_cap = acat_cap, window = set$window, min_genotype_count = set$min_genotype_count,
                              fdr = set$fdr, control = control, abs_z = abs_z,
                              mirror_labels = if (is.null(absolute))
                                "v3m: exact-pair boundary events are p0 units at the promoter that goes silent (absolute activity not given)" else
                                "absolute activity: exact-pair boundary events renamed p1 at the partner when only its absolute activity changes"),
              seconds = proc.time()[["elapsed"]] - t0)
  out <- structure(out, class = "bohr_conditional")
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(out_dir)) stop("Cannot create `out_dir`: ", out_dir, call. = FALSE)
    for (k in c("signals", "units", "rounds", "loo", "labels")) if (!is.null(out[[k]]) && nrow(out[[k]])) write_tsv(out[[k]], file.path(out_dir, paste0(k, ".tsv")))
  }
  out
}

## Conditional analysis of one tissue. phenotype: promoters x the tissue's samples; G: dosages (variants x the same samples,
## unfilled); Z: the tissue's covariates (or NULL); U: the tissue's significant units (with chr); th: its threshold;
## snp_t / prom_t: the tissue's tested pairs and promoter table of the bohr_cis() run (for the consistency checks).
cond_one <- function(phenotype, G, V, prom, Z, U, th, snp_t, prom_t, set, control, acat_cap, max_signals, absolute, abs_z,
                     ncores, chunk_size, label = NULL) {
  samples <- colnames(phenotype)
  where <- if (is.null(label)) "" else paste0(" (tissue ", label, ")")
  ## covariates as in cis_one()
  Zm <- as_covariates(Z, length(samples))
  if (!is.null(Z)) {
    if (!is.null(rownames(Zm))) {
      miss <- setdiff(samples, rownames(Zm))
      if (length(miss)) stop("`covariates` has no row for ", length(miss), " sample(s)", where, ", e.g. ", miss[1], ".", call. = FALSE)
      Zm <- Zm[samples, , drop = FALSE]
    } else if (nrow(Zm) != length(samples)) stop("`covariates` must have one row per sample of `phenotype`.", call. = FALSE)
  }
  rownames(Zm) <- samples
  C_all <- make_design(Zm)

  ## variants of the windows of the analysed promoters; missing dosages filled per variant over the tissue's samples
  prs <- unique(U$promoter); pi <- match(prs, prom$id)
  vchr <- split(seq_len(nrow(V)), V$chr)
  win <- lapply(pi, function(i) { vi <- vchr[[prom$chr[i]]]; if (is.null(vi)) integer(0) else vi[cis_index(V$start[vi], V$end[vi], prom$pos[i], set$window)] })
  used <- sort(unique(unlist(win)))
  Gi <- if (length(used)) impute_modal(G[used, , drop = FALSE]) else NULL
  pos_in <- match(seq_len(nrow(V)), used)

  ## per-promoter inputs (run_conditional_v3c.R:92-100, the steps of run_final2.R)
  prep <- function(k) {
    p <- prs[k]
    fit <- fit_promoter(phenotype[p, ], C_all, control)
    if (!is.na(fit$skip)) stop("Promoter ", p, where, " is skipped (", fit$skip, ") but has a significant unit in `cis`: the inputs differ from the bohr_cis() call.", call. = FALSE)
    if (!is.null(prom_t)) { j <- match(p, prom_t$promoter)
      if (!is.na(j) && !identical(prom_t$family[j], fit$family))
        stop("Promoter ", p, where, ": family ", fit$family, " here, ", prom_t$family[j], " in `cis`; the inputs differ from the bohr_cis() call.", call. = FALSE) }
    idx <- win[[k]]
    Gp <- if (length(idx)) Gi[pos_in[idx], fit$samples, drop = FALSE] else matrix(0L, 0, length(fit$samples))
    cnt <- apply(Gp, 1, function(x) tabulate(x + 1L, nbins = 3L)); dim(cnt) <- c(3L, length(idx))
    keep <- which(apply(cnt, 2, min) >= set$min_genotype_count)
    Gt <- t(Gp[keep, , drop = FALSE]); colnames(Gt) <- V$id[idx[keep]]; rownames(Gt) <- samples[fit$samples]
    if (!is.null(snp_t)) { tested <- snp_t$variant[snp_t$promoter == p]
      if (!setequal(colnames(Gt), tested)) stop("Promoter ", p, where, ": the tested variants differ from the bohr_cis() run (", length(tested), " there, ",
                                                ncol(Gt), " here); give the inputs of that call.", call. = FALSE) }
    list(y = fit$y, family = fit$family, C = C_all[fit$samples, , drop = FALSE], G = Gt)
  }
  PP <- par_apply(seq_along(prs), prep, ncores, preschedule = TRUE)
  names(PP) <- prs

  ## ---- forward rounds (run_conditional_v3c.R:137-164), synchronous over units ----
  nu <- nrow(U)
  for (i in seq_len(nu)) if (!U$leader[i] %in% colnames(PP[[U$promoter[i]]]$G))
    stop("The leader of unit ", U$promoter[i], ":", U$component[i], where, " is not among its tested variants.", call. = FALSE)
  stt <- lapply(seq_len(nu), function(i) list(key = U$promoter[i], type = U$component[i], signals = U$leader[i], active = TRUE, stop = NA_character_))
  sig <- lapply(seq_len(nu), function(i) list(data.frame(k = 1L, variant = U$leader[i], selected_round = 0L, log_p = U$leader_log_p[i],
                                                          log_p_chi = NA_real_, acat_p = U$p[i], stringsAsFactors = FALSE)))
  rounds <- list()
  for (r in seq_len(max_signals)) {
    act <- which(vapply(stt, `[[`, TRUE, "active")); if (!length(act)) break
    jobs <- lapply(act, function(i) { s <- stt[[i]]
      list(i = i, key = s$key, type = s$type, cs = s$signals, test = setdiff(colnames(PP[[s$key]]$G), s$signals)) })
    ps <- cond_pass(jobs, PP, control, ncores, chunk_size)
    for (b in seq_along(jobs)) {
      i <- jobs[[b]]$i; z <- ps$res[[b]]
      if (is.null(z)) z <- data.frame(variant = character(0), est = numeric(0), se = numeric(0), log_p = numeric(0), log_p_chi = numeric(0))
      a <- acat_log(z$log_p, cap = acat_cap)
      bj <- if (nrow(z) && any(is.finite(z$log_p))) which.min(z$log_p) else integer(0)
      rounds[[length(rounds) + 1L]] <- data.frame(promoter = stt[[i]]$key, component = stt[[i]]$type, round = r, n_conditioned = length(stt[[i]]$signals),
        n_tested = sum(is.finite(z$log_p)), n_failed = sum(!is.finite(z$log_p)), acat_p = a, acat_p_chi = acat_log(z$log_p_chi, cap = acat_cap),
        best_variant = if (length(bj)) z$variant[bj] else NA_character_, best_log_p = if (length(bj)) z$log_p[bj] else NA_real_,
        cols_dropped = ps$dropped[b], null_ok = ps$null_ok[b], stringsAsFactors = FALSE)
      go <- is.finite(a) && a <= th && length(bj) && length(stt[[i]]$signals) < max_signals
      if (go) {
        stt[[i]]$signals <- c(stt[[i]]$signals, z$variant[bj])
        sig[[i]][[length(sig[[i]]) + 1L]] <- data.frame(k = length(stt[[i]]$signals), variant = z$variant[bj], selected_round = r, log_p = z$log_p[bj],
                                                         log_p_chi = z$log_p_chi[bj], acat_p = a, stringsAsFactors = FALSE)
      } else {
        stt[[i]]$active <- FALSE
        stt[[i]]$stop <- if (!is.finite(a)) "no finite conditional p" else if (a > th) "window ACAT above threshold" else "cap"
      }
    }
  }

  ## ---- leave-one-out scans (run_conditional_v3c.R:166-174) and the backward step (merge_conditional_v3c.R:30-34) ----
  loo <- NULL
  multi_i <- which(vapply(stt, function(s) length(s$signals) >= 2L, TRUE))
  if (length(multi_i)) {
    jobs <- unlist(lapply(multi_i, function(i) { s <- stt[[i]]
      lapply(seq_along(s$signals), function(k) list(i = i, k = k, key = s$key, type = s$type, cs = s$signals[-k],
                                                     test = setdiff(colnames(PP[[s$key]]$G), s$signals[-k]))) }), recursive = FALSE)
    ps <- cond_pass(jobs, PP, control, ncores, chunk_size)
    loo <- bind_df(lapply(seq_along(jobs), function(b) { z <- ps$res[[b]]; if (is.null(z) || !nrow(z)) return(NULL)
      data.frame(promoter = jobs[[b]]$key, component = jobs[[b]]$type, k = jobs[[b]]$k, z, stringsAsFactors = FALSE) }))
  }
  sig_tab <- bind_df(lapply(seq_len(nu), function(i) { d <- do.call(rbind, sig[[i]])
    data.frame(promoter = stt[[i]]$key, gene = U$gene[i], component = stt[[i]]$type, d, stringsAsFactors = FALSE) }))
  sig_tab$pos <- V$start[match(sig_tab$variant, V$id)]
  for (k in c("loo_est", "loo_se", "loo_log_p", "loo_log_p_chi", "loo_acat_p")) sig_tab[[k]] <- NA_real_
  sig_tab$loo_top_variant <- NA_character_
  if (!is.null(loo)) {
    key_l <- paste(loo$promoter, loo$component, loo$k, sep = "\r"); key_s <- paste(sig_tab$promoter, sig_tab$component, sig_tab$k, sep = "\r")
    for (j in which(key_s %in% key_l)) {
      rows <- which(key_l == key_s[j]); z <- loo[rows, , drop = FALSE]
      sig_tab$loo_acat_p[j] <- acat_log(z$log_p, cap = acat_cap)
      sig_tab$loo_top_variant[j] <- z$variant[which.min(replace(z$log_p, !is.finite(z$log_p), Inf))]
      m <- match(sig_tab$variant[j], z$variant)
      if (!is.na(m)) { sig_tab$loo_est[j] <- z$est[m]; sig_tab$loo_se[j] <- z$se[m]; sig_tab$loo_log_p[j] <- z$log_p[m]; sig_tab$loo_log_p_chi[j] <- z$log_p_chi[m] }
    }
  }
  n_sig <- stats::ave(seq_len(nrow(sig_tab)), sig_tab$promoter, sig_tab$component, FUN = length)
  sig_tab$kept <- ifelse(n_sig == 1L, TRUE, ifelse(is.na(sig_tab$loo_acat_p), NA, sig_tab$loo_acat_p <= th))
  sig_tab <- sig_tab[, c("promoter", "gene", "component", "k", "variant", "pos", "selected_round", "log_p", "log_p_chi", "acat_p",
                         "loo_est", "loo_se", "loo_log_p", "loo_log_p_chi", "loo_acat_p", "loo_top_variant", "kept")]
  units <- data.frame(promoter = U$promoter, gene = U$gene, component = U$component, chr = U$chr, family = U$family, n = U$n, n0 = U$n0, n1 = U$n1,
                      p = U$p, q = U$q, leader = U$leader, n_signals = vapply(stt, function(s) length(s$signals), 0L),
                      n_kept = NA_integer_, stop = vapply(stt, `[[`, "", "stop"), stringsAsFactors = FALSE)
  ukey <- paste(units$promoter, units$component, sep = "\r"); skey <- paste(sig_tab$promoter, sig_tab$component, sep = "\r")
  units$n_kept <- vapply(ukey, function(k) { v <- sig_tab$kept[skey == k]; if (anyNA(v)) NA_integer_ else sum(v) }, 0L, USE.NAMES = FALSE)
  rounds <- bind_df(rounds)

  ## ---- optional absolute-activity labels of exact-pair boundary events (mirror_abs_labels.R) ----
  labels <- NULL
  if (!is.null(absolute)) {
    lab <- abs_labels_cond(U, sig_tab, absolute, C_all, Gi, pos_in, V, abs_z, where)
    labels <- lab
    ch <- if (is.null(lab)) NULL else lab[lab$changed %in% TRUE, , drop = FALSE]
    units$relabelled_from <- NA_character_; sig_tab$relabelled_from <- NA_character_
    for (j in seq_len(if (is.null(ch)) 0L else nrow(ch))) {
      rn <- function(d, swap_n = FALSE) { w <- d$promoter == ch$S[j] & d$component == "p0"
        if (!is.null(d$relabelled_from)) d$relabelled_from[w] <- paste0(ch$S[j], ":p0")
        if (swap_n) { n0 <- d$n0[w]; d$n0[w] <- d$n1[w]; d$n1[w] <- n0 }
        d$promoter[w] <- ch$O[j]; d$component[w] <- "p1"; d }
      units <- rn(units, swap_n = TRUE); sig_tab <- rn(sig_tab)
      if (!is.null(rounds)) rounds <- rn(rounds)
      if (!is.null(loo)) loo <- rn(loo)
    }
  }
  list(signals = sig_tab, units = units, rounds = rounds, loo = loo, labels = labels)
}

## mirror_abs_labels.R:48-50: least squares of y on X; z of the last column (NA when X is rank-deficient)
abs_fit_last <- function(X, y) {
  f <- stats::lm.fit(X, y); p <- ncol(X)
  if (f$rank < p) return(c(b = NA_real_, se = NA_real_, z = NA_real_))
  s2 <- sum(f$residuals^2) / (length(y) - p); V <- chol2inv(qr.R(f$qr)); j <- which(f$qr$pivot == p)
  b <- unname(f$coefficients[p]); se <- sqrt(V[j, j] * s2)
  c(b = b, se = se, z = b / se)
}
## labels of the significant p0 units of boundary_exact pairs in one tissue (mirror_abs_labels.R:55-75): absolute activity of
## S (the silent promoter) and O (its partner) on the covariates, the unit's other kept signals and the leader's dosage, over
## all of the tissue's samples; rank-deficient -> the other signals are dropped. S only / both / neither / not estimable -> p0(S);
## O only -> p1(O).
abs_labels_cond <- function(U, sig_tab, absolute, C_all, Gi, pos_in, V, abs_z, where) {
  ev <- which(U$component == "p0" & U$boundary_exact %in% TRUE & !is.na(U$mirror_partner))
  if (!length(ev)) return(NULL)
  samples <- rownames(C_all)
  need <- unique(c(U$promoter[ev], U$mirror_partner[ev]))
  miss <- setdiff(need, rownames(absolute)); if (length(miss)) stop("`absolute` has no row for promoter ", miss[1], where, ".", call. = FALSE)
  miss <- setdiff(samples, colnames(absolute)); if (length(miss)) stop("`absolute` has no column for sample ", miss[1], where, ".", call. = FALSE)
  A <- absolute[need, samples, drop = FALSE]
  if (anyNA(A) || any(!is.finite(A))) stop("`absolute` must not contain missing values for the samples used", where, ".", call. = FALSE)
  dose_of <- function(v) { r <- pos_in[match(v, V$id)]; as.numeric(Gi[r, ]) }
  rows <- lapply(ev, function(i) {
    S <- U$promoter[i]; O <- U$mirror_partner[i]; ld <- U$leader[i]
    ks <- sig_tab[sig_tab$promoter == S & sig_tab$component == "p0" & sig_tab$kept %in% TRUE, , drop = FALSE]
    cond <- setdiff(ks$variant[order(ks$k)], ld)
    dose <- dose_of(ld)
    Xc <- cbind(C_all, if (length(cond)) vapply(cond, dose_of, numeric(length(samples))) else NULL, dose); dropped <- FALSE
    if (qr(Xc)$rank < ncol(Xc)) { Xc <- cbind(C_all, dose); dropped <- length(cond) > 0 }
    fS <- abs_fit_last(Xc, A[S, ]); fO <- abs_fit_last(Xc, A[O, ])
    data.frame(S = S, O = O, gene = U$gene[i], leader = ld, n_cond = if (dropped) 0L else length(cond), cond_dropped = dropped,
               z_S = fS[["z"]], z_O = fO[["z"]], b_S = fS[["b"]], b_O = fO[["b"]], se_S = fS[["se"]], se_O = fO[["se"]], stringsAsFactors = FALSE)
  })
  L <- do.call(rbind, rows)
  rs <- abs(L$z_S) >= abs_z; ro <- abs(L$z_O) >= abs_z
  L$outcome <- ifelse(is.na(L$z_S) | is.na(L$z_O), "not estimable", ifelse(rs & !ro, "S only", ifelse(ro & !rs, "O only", ifelse(rs & ro, "both", "neither"))))
  L$label <- ifelse(L$outcome == "O only", paste0(L$O, ":p1"), paste0(L$S, ":p0"))
  L$changed <- L$outcome == "O only"
  L
}

#' @export
print.bohr_conditional <- function(x, ...) {
  S <- x$signals; U <- x$units
  cat("BOHR conditional analysis\n")
  if (is.null(U) || !nrow(U)) { cat("  no significant unit analysed\n"); return(invisible(x)) }
  by_t <- if (!is.null(U$tissue)) U$tissue else rep("", nrow(U))
  st <- if (!is.null(S$tissue)) S$tissue else rep("", nrow(S))
  tis <- unique(by_t)
  tab <- data.frame(tissue = tis, units = vapply(tis, function(t) sum(by_t == t), 0L),
                    signals = vapply(tis, function(t) sum(st == t), 0L),
                    kept = vapply(tis, function(t) sum(S$kept[st == t] %in% TRUE), 0L),
                    units_2plus = vapply(tis, function(t) sum(U$n_kept[by_t == t] >= 2, na.rm = TRUE), 0L),
                    max_kept = vapply(tis, function(t) suppressWarnings(max(U$n_kept[by_t == t], na.rm = TRUE)), 0),
                    row.names = NULL, stringsAsFactors = FALSE)
  if (is.null(U$tissue)) tab$tissue <- NULL
  print(tab, row.names = FALSE)
  cat("Signals: forward selection at the discovery threshold; kept: after the backward (leave-one-out) step.\n")
  cat("Mirror labels:", x$settings$mirror_labels, "\n")
  if (!is.null(x$labels) && nrow(x$labels)) cat("  exact-pair boundary events judged:", nrow(x$labels), "; renamed p1 at the partner:", sum(x$labels$changed), "\n")
  invisible(x)
}
