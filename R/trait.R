## Trait regression: association of every feature with donor traits in each tissue.
##
## Ported from the trait regression of the BOHR paper (tissue_regression/):
##   trait_regression_20261007/run_trait.R (md5 c9a4cb6b...): donors and trait-level skips :35-47, coding :58-59, the
##     per-feature steps :68-116 (prepare_y, no_interior, the five-per-group rules, orthonormal design, null fit, snp_tests with
##     the trait coding, the F correction of mu, Firth p0/p1, estimates "per contrast" = 2 x the coded estimate);
##   trait_regression_20261007/data.R (md5 f090d52f...): cov_design :74-94, trait_contrast :100-102, MIN_CLASS / MIN_HARDY;
##   trait_regression_20261104/summarize.R :51-98: units, mirror rule v3 de-duplication, the final (eligibility) rules and BH
##     within tissue x trait x component; mirror_v3.R :18 and trait_regression_20261005/mirror.R :35-58 (pairs, canonical);
##   trait_regression_20261106_mirror_abs/01_abs_labels.R :39-46, 56-60, 77-81 and 03_apply_B.R :30-45 (labels of exact-pair
##     boundary events by absolute activity, option B), whose results are export_B/trait_regression_units.tsv.gz.
## The engine is the package's port (identical md5 of the engine files in the trait run and the cis run).

## ---- tables with IDs: objects or TSV files (first column = IDs); data frames keep categorical columns ------------------
read_frame <- function(x, what) {
  if (is.null(x)) return(NULL)
  if (is_path(x)) { d <- read_tsv(x, what); rn <- as.character(d[[1]]); d <- d[, -1, drop = FALSE]; rownames(d) <- rn; return(d) }
  if (is.matrix(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!is.data.frame(x)) stop("`", what, "` must be a data frame, a matrix or the path of a TSV file.", call. = FALSE)
  x
}
## rows of d for the samples: by sample ID, else by individual ID, else (no row names) in sample order
rows_for <- function(d, samples, ind, what) {
  rn <- if (is.data.frame(d) && .row_names_info(d) < 0) NULL else rownames(d)
  if (!is.null(rn) && all(samples %in% rn)) { out <- d[samples, , drop = FALSE] }
  else if (!is.null(rn) && all(ind %in% rn)) { out <- d[ind, , drop = FALSE] }
  else if (is.null(rn) && nrow(d) == length(samples)) { out <- d }
  else stop("`", what, "` needs a row for every sample (row names = sample IDs) or for every individual (row names = individual IDs).", call. = FALSE)
  rownames(out) <- samples
  out
}
## covariates of tissue tt: a list named by tissue, a folder with <tissue>.tsv files, or one table
trait_covariates <- function(covariates, tt) {
  if (is.null(covariates)) return(NULL)
  if (is.list(covariates) && !is.data.frame(covariates)) {
    if (is.null(names(covariates)) || !tt %in% names(covariates)) stop("`covariates` has no element for tissue ", tt, "; give a list named by tissue.", call. = FALSE)
    return(read_frame(covariates[[tt]], paste0("covariates$", tt)))
  }
  if (is_path(covariates) && dir.exists(covariates)) {
    f <- file.path(covariates, paste0(tt, c(".tsv", ".txt", ".tsv.gz", ".txt.gz"))); f <- f[file.exists(f)]
    if (!length(f)) stop("No covariate file for tissue ", tt, " in ", covariates, " (expected ", tt, ".tsv).", call. = FALSE)
    return(read_frame(f[1], paste0("covariates (", tt, ")")))
  }
  covariates
}

## ---- traits: numeric codes and types ----------------------------------------------------------------------------------
## Logical -> 0/1; a categorical trait (factor or character) needs exactly two levels and becomes the indicator of the second
## (factor level order, else sorted). type: "binary" (two values; contrast low -> high) or "continuous" (+1 SD).
prepare_traits <- function(tr, trait_type) {
  out <- list()
  for (v in names(tr)) {
    x <- tr[[v]]; levs <- NULL
    if (is.logical(x)) x <- as.numeric(x)
    else if (is.factor(x) || is.character(x)) {
      levs <- if (is.factor(x)) levels(droplevels(x)) else sort(unique(x[!is.na(x)]))
      if (length(levs) != 2L) stop("Trait ", v, " is categorical with ", length(levs), " levels; a categorical trait needs exactly two.", call. = FALSE)
      x <- ifelse(is.na(x), NA_real_, as.numeric(as.character(x) == levs[2]))
    } else if (!is.numeric(x)) stop("Trait ", v, " must be numeric, logical or a two-level factor.", call. = FALSE)
    x <- as.numeric(x)
    vals <- sort(unique(x[!is.na(x)]))
    ty <- if (!is.null(trait_type) && v %in% names(trait_type)) trait_type[[v]] else if (length(vals) == 2L) "binary" else "continuous"
    if (!ty %in% c("binary", "continuous")) stop("`trait_type` must be \"binary\" or \"continuous\" (trait ", v, ").", call. = FALSE)
    if (ty == "binary" && length(vals) != 2L) stop("Trait ", v, " is declared binary but has ", length(vals), " distinct values.", call. = FALSE)
    out[[v]] <- list(x = x, type = ty, values = if (ty == "binary") vals else NULL, levels = levs)
  }
  out
}

## ---- covariate design of one tissue and trait (data.R:74-94, cov_design) -----------------------------------------------
## numeric columns as they are, then each categorical column as indicators of its levels, after merging the levels with fewer
## than min_level_count donors into the most common level (the reference); columns constant on the donors dropped; then
## standardised over the donors with an intercept (= cov_matrix()). cv: data frame, rows = donors.
trait_design <- function(cv, donors, min_level_count) {
  n <- length(donors)
  if (is.null(cv) || !ncol(cv)) { Z <- matrix(numeric(0), n, 0L, dimnames = list(donors, NULL)); return(make_design(Z)) }
  isnum <- vapply(cv, function(v) is.numeric(v) || is.logical(v), TRUE)
  Z <- as.matrix(cv[, isnum, drop = FALSE]); storage.mode(Z) <- "double"
  for (v in names(cv)[!isnum]) {
    f0 <- cv[[v]]; lev <- if (is.factor(f0)) levels(f0) else sort(unique(as.character(f0)))
    x <- as.character(f0); tb <- table(factor(x, levels = lev)); ref <- names(which.max(tb))
    x[x %in% names(tb)[tb < min_level_count]] <- ref
    f <- stats::relevel(factor(x), ref = ref)
    if (nlevels(f) >= 2) { D <- stats::model.matrix(~ f)[, -1, drop = FALSE]; colnames(D) <- paste0(v, "_", levels(f)[-1]); Z <- cbind(Z, D) }
  }
  if (ncol(Z)) Z <- Z[, apply(Z, 2, function(x) stats::sd(x) > 0), drop = FALSE]
  rownames(Z) <- donors
  make_design(Z)
}

## ---- one feature against a coded trait (run_trait.R:68-116) -------------------------------------------------------------
## yv: activity over the donors (named); C_all: trait_design(); gp: coded trait over the donors (x' = 2 (x - lo)/(hi - lo) or
## 2 (x - mean)/sd); binary: apply the five-per-group rules. Returns the feature row (engine values on the coded scale).
trait_feature <- function(p, gene, yv, C_all, gp, binary, min_group, control) {
  base <- list(promoter = p, gene = gene, family = NA_character_, skip = NA_character_, n = NA_integer_, n0 = NA_integer_, n1 = NA_integer_,
               n_mid = NA_integer_, n_minor = NA_integer_, n_minor_mid = NA_integer_, mu_testable = NA, n_cov = NA_integer_,
               log_p_mu = NA_real_, log_p_p0 = NA_real_, log_p_p1 = NA_real_, est_mu = NA_real_, se_mu = NA_real_, est_p0 = NA_real_,
               se_p0 = NA_real_, est_p1 = NA_real_, se_p1 = NA_real_, effect_mu = NA_real_, effect_p0 = NA_real_, effect_p1 = NA_real_,
               complete_sep = NA, error = NA_character_, firth_error = NA)
  fit <- fit_promoter(yv, C_all, control)
  if (!is.na(fit$skip)) { base$skip <- fit$skip; return(base) }
  gi <- unname(gp[fit$samples]); mid <- fit$y > 0 & fit$y < 1
  if (binary) {
    tb <- table(gi)
    if (length(tb) < 2 || min(tb) < min_group) { base$skip <- "trait_rare_in_feature"; return(base) }
    base$n_minor <- as.integer(min(tb))
    tm <- table(factor(gi[mid], levels = sort(unique(gi))))
    base$n_minor_mid <- as.integer(min(tm)); base$mu_testable <- min(tm) >= min_group
  } else base$mu_testable <- TRUE
  base$family <- fit$family; base$n <- fit$n; base$n0 <- fit$n0; base$n1 <- fit$n1; base$n_mid <- fit$n_mid; base$n_cov <- fit$n_cov
  e <- engine_one(fit, gi, control)
  map <- if (isTRUE(control$firth)) CLEAR_MAP else CLEAR_MAP_LRT
  for (k in c("log_p_mu", "log_p_p0", "log_p_p1", "effect_mu", "effect_p0", "effect_p1", "complete_sep")) base[[k]] <- e[[map[[k]]]]
  for (k in c("est_mu", "se_mu", "est_p0", "se_p0", "est_p1", "se_p1")) base[[k]] <- 2 * e[[map[[k]]]]      # per contrast (run_trait.R:107-112)
  base$error <- e$error; base$firth_error <- e$firth_error
  base
}

#' Association of features with donor traits in each tissue
#'
#' Tests every feature (for example every promoter's relative activity) for
#' association with one or more donor traits in each tissue, with the model
#' and tests of [bohr_test()], and controls the false-discovery rate within
#' each tissue, trait and component, following the trait analysis of the
#' BOHR paper (Methods Section 6; SI Section 6).
#'
#' For each tissue and trait:
#' \enumerate{
#'   \item **Donors.** The tissue's samples with a value of the trait and of
#'     every covariate. A covariate column with the name of the tested trait is
#'     left out (in the paper, sex and age were covariates except when tested).
#'     The trait is skipped with fewer than `control$min_samples` donors or no
#'     variation (`"trait_constant_or_n"`), or, for a binary trait, fewer than
#'     `min_group` donors in a group (`"trait_rare"`).
#'   \item **Covariates.** Numeric columns are used as they are. A categorical
#'     column (factor or character) becomes indicators of its levels, after the
#'     levels with fewer than `min_level_count` donors are merged into the most
#'     common level (the reference). Columns that are constant on the donors
#'     are dropped. The columns are standardised over the donors and an
#'     intercept is added.
#'   \item **Coding.** A binary trait (two distinct values) is coded so that
#'     the contrast is the lower value to the higher (0 to 1; for a two-level
#'     factor, the first level to the second); a continuous trait so that the
#'     contrast is the donors' mean to one standard deviation above it.
#'   \item **Tests.** Each feature is fitted and tested as in [bohr_fit()] and
#'     [bohr_test()], with the fitting rules of `control` (SI Section 2): mu
#'     with the F-corrected likelihood-ratio test, p0 and p1 with the
#'     Firth-penalised likelihood-ratio test. For a binary trait, a feature
#'     needs at least `min_group` donors in each group among its samples, and
#'     mu also among its intermediate samples.
#'   \item **Units and mirror pairs.** One unit per feature and component (mu,
#'     p0 when zeros remain, p1 when ones remain) with a finite P value. With
#'     `promoters$gene`, exact mirror pairs are found in each tissue
#'     ([mirror_pairs()], among the features with data in that tissue) and
#'     de-duplicated as in [bohr_cis()]: for a `boundary_exact` pair the p1
#'     units of both members are dropped (each boundary event is tested as p0
#'     at the promoter that goes silent); for an `interior_exact` pair the mu
#'     unit of the partner is dropped.
#'   \item **Eligibility and q values.** A unit is `eligible` (enters the
#'     Benjamini-Hochberg correction) when it is not a mirror duplicate, a mu
#'     unit has at least `min_n_mid` intermediate samples, a p0 or p1 unit has
#'     a boundary share within `boundary_share`, a p0 or p1 unit's boundary
#'     part is not completely separated by the covariates
#'     (`exclude_separated`), and the tissue and trait have at least
#'     `min_donors` donors. q values are computed within each tissue, trait and
#'     component; `significant` is q at most `fdr`.
#'   \item **Labels of exact-pair boundary events.** With `absolute`, each
#'     significant p0 unit of a `boundary_exact` pair is labelled by absolute
#'     activity (SI Section 2): the absolute activity of each promoter of the
#'     pair is regressed (least squares) on the coded trait and the covariates,
#'     over the donors; a promoter changes when \eqn{|z| \ge} `abs_z`. When
#'     only the other promoter changes, the event is reported as p1 at that
#'     promoter (its own fit, the same test) with the q value of the p0 test
#'     (`bh_family = "p0"`, `abs_relabel = TRUE`), and the p0 row becomes the
#'     duplicate. Without `absolute`, every such event stays p0 at the silent
#'     promoter (the rule "v3m"), and `settings$mirror_labels` says so.
#' }
#' **Scales.** Estimates and standard errors (`est`, `se`) are for the
#' contrast, as in the paper's tables: per one standard deviation of a
#' continuous trait and per change from the lower to the higher value of a
#' binary trait. This differs from [bohr_test()], which reports coefficients
#' per unit of the predictor. Effects (`effect`) are for the same contrast:
#' the change in the mean of the intermediate values (mu) or in the
#' probability of 0 or 1 (p0, p1).
#'
#' **Inputs as objects or files.** `phenotype` (and `absolute`) are matrices
#' or TSV files as in [bohr_cis()]. `traits` and `covariates` are data frames
#' (or matrices, or TSV files whose first column holds the IDs) with one row
#' per sample or per individual: rows are matched by sample ID, else by
#' individual ID (see `individual`). Covariates can also be given per tissue
#' (a list named by tissue or a folder of `<tissue>.tsv` files).
#'
#' @param phenotype Proportions in \[0, 1\] (relative activities), features in
#'   rows (row names = feature IDs) and samples in columns (unique column
#'   names); a numeric matrix or a TSV path. A feature without any value in a
#'   tissue is not analysed there.
#' @param traits Donor traits, one column per trait: numeric (binary or
#'   continuous), logical, or a two-level factor or character column. A data
#'   frame or matrix with one row per sample or individual (row names = IDs),
#'   or a TSV path (first column = IDs). `NA` marks a missing value.
#' @param covariates `NULL`, or covariates with one row per sample or
#'   individual (numeric or categorical columns): a data frame or matrix, a
#'   TSV path, or, per tissue, a list named by tissue or a folder of
#'   `<tissue>.tsv` files.
#' @param tissue,individual As in [bohr_cis()]: `NULL` (one tissue), the tissue
#'   of each sample, or a sample table with `sample`, `tissue` and optionally
#'   `individual`; `individual` maps samples to the rows of `traits` and
#'   `covariates` given per individual.
#' @param promoters `NULL`, or a data frame (or TSV path) with columns `id` and
#'   `gene` for the features, needed for the mirror rule.
#' @param trait_type Optional character vector named by trait, `"binary"` or
#'   `"continuous"`; by default a trait with exactly two distinct values is
#'   binary and any other continuous.
#' @param min_group A binary trait needs at least this many donors in each
#'   group, overall and among a feature's samples, and for the mu test among
#'   its intermediate samples. Default 5 (SI Section 6).
#' @param min_level_count Levels of a categorical covariate with fewer donors
#'   than this are merged into the most common level. Default 10 (SI Section
#'   6, the Hardy scale).
#' @param min_donors Units of a tissue and trait with fewer donors are not
#'   eligible. Default 50 (Methods Section 6).
#' @param min_n_mid A mu unit is eligible with at least this many intermediate
#'   samples. Default 80 (Methods Section 2).
#' @param boundary_share A p0 (p1) unit is eligible when its share of samples
#'   at 0 (at 1) lies in this range, ends included. Default `c(0.10, 0.90)`
#'   (Methods Section 2).
#' @param exclude_separated Exclude p0 and p1 units whose boundary part is
#'   completely separated by the covariates from the correction. Default
#'   `TRUE` (Methods Section 6; SI Section 6).
#' @param fdr q-value threshold of `significant`. Default 0.05 (Methods
#'   Section 6).
#' @param mirror,mirror_tol,major_score The exact-mirror-pair rule, as in
#'   [bohr_cis()]; the canonical member (which keeps the mu unit of an
#'   `interior_exact` pair) has the higher `major_score`, by default the mean
#'   over all supplied samples pooled across tissues.
#' @param absolute Optional absolute activity for the labels of exact-pair
#'   boundary events: a numeric matrix (features x samples, the layout of
#'   `phenotype`) or a TSV path, with no missing values on the donors used.
#' @param abs_z A promoter's absolute activity changes when \eqn{|z| \ge}
#'   `abs_z`. Default 1.96 (SI Section 2).
#' @param control Fitting rules and optimiser settings, [bohr_control()].
#' @param ncores Number of cores (forked processes; 1 on Windows).
#' @param out_dir Optional folder; when given, the tables `results`,
#'   `features`, `traits` and `labels` are also written as TSV files.
#' @return An object of class `"bohr_trait"`, a list with
#'   \describe{
#'     \item{`results`}{One row per tissue, trait, feature and component with a
#'       finite P value: `tissue`, `trait`, `promoter`, `gene`, `component`,
#'       `p`, `log_p`, `est`, `se` (per contrast), `effect`, `eligible`, `q`,
#'       `significant`, `dedup_drop`, `duplicate_of`, the mirror columns
#'       (`mirror_role` and `mirror_partner` are given for `boundary_exact`
#'       pairs, as in the paper's tables), `family`, `n`, `n0`, `n1`, `n_mid`, `boundary_share`, `n_donors`,
#'       `complete_sep`, `bh_family` (the BH family of `q`) and `abs_relabel`.}
#'     \item{`features`}{One row per tissue, trait and feature: family or skip
#'       reason, sample and group counts, and the engine messages.}
#'     \item{`traits`}{One row per tissue and trait: type, contrast, donors,
#'       covariate columns used, skip reason.}
#'     \item{`labels`}{With `absolute`: the judged events (`S`, `O`, `z_S`,
#'       `z_O`, `outcome`, `label`, `changed`).}
#'     \item{`mirror`, `settings`, `seconds`}{The pairs per tissue, the settings
#'       (including `mirror_labels`) and the elapsed time.}
#'   }
#' @seealso [bohr_test()], [bohr_cis()], [mirror_pairs()].
#' @examples
#' set.seed(11)
#' n <- 150; s <- paste0("S", 1:n)
#' sex <- rbinom(n, 1, 0.5); age <- round(runif(n, 20, 70))
#' p0 <- plogis(-1.5 + 1.2 * sex)                                  # sex shifts p0 of P1
#' y1 <- rbeta(n, 4, 4); y1[runif(n) < p0] <- 0
#' y2 <- rbeta(n, 3 + 0.05 * age, 3)                                # age shifts mu of P2
#' ph <- rbind(P1 = y1, P2 = y2); colnames(ph) <- s
#' tr <- data.frame(SEX = sex, AGE = age, row.names = s)
#' res <- bohr_trait(ph, tr, covariates = tr, min_n_mid = 30, min_donors = 20)
#' res
#' res$results[, c("trait", "promoter", "component", "p", "q", "est", "effect")]
#' @export
bohr_trait <- function(phenotype, traits, covariates = NULL, tissue = NULL, individual = NULL, promoters = NULL, trait_type = NULL,
                       min_group = 5L, min_level_count = 10L, min_donors = 50L, min_n_mid = 80L, boundary_share = c(0.10, 0.90),
                       exclude_separated = TRUE, fdr = 0.05, mirror = TRUE, mirror_tol = 1e-6, major_score = NULL,
                       absolute = NULL, abs_z = 1.96, control = bohr_control(), ncores = 1L, out_dir = NULL) {
  t0 <- proc.time()[["elapsed"]]
  check_control(control)
  min_group <- check_count(min_group, "min_group", min = 1); min_level_count <- check_count(min_level_count, "min_level_count")
  min_donors <- check_count(min_donors, "min_donors"); min_n_mid <- check_count(min_n_mid, "min_n_mid")
  ncores <- check_count(ncores, "ncores", min = 1)
  check_flag(exclude_separated, "exclude_separated"); check_flag(mirror, "mirror")
  if (!is.numeric(boundary_share) || length(boundary_share) != 2L || anyNA(boundary_share) || boundary_share[1] > boundary_share[2])
    stop("`boundary_share` must be two increasing numbers.", call. = FALSE)
  if (!is.numeric(fdr) || length(fdr) != 1L || !(fdr > 0 && fdr < 1)) stop("`fdr` must lie strictly between 0 and 1.", call. = FALSE)
  if (!is.numeric(mirror_tol) || length(mirror_tol) != 1L || !(mirror_tol >= 0)) stop("`mirror_tol` must be a single number >= 0.", call. = FALSE)
  if (!is.numeric(abs_z) || length(abs_z) != 1L || !is.finite(abs_z) || abs_z <= 0) stop("`abs_z` must be a single positive number.", call. = FALSE)
  if (!is.null(trait_type) && (!is.character(trait_type) || is.null(names(trait_type)))) stop("`trait_type` must be a character vector named by trait.", call. = FALSE)
  if (!is.null(major_score) && (!is.numeric(major_score) || is.null(names(major_score)))) stop("`major_score` must be a numeric vector named by feature ID.", call. = FALSE)
  if (!is.null(out_dir) && (!is.character(out_dir) || length(out_dir) != 1L || is.na(out_dir) || !nzchar(out_dir)))
    stop("`out_dir` must be the path of one folder.", call. = FALSE)

  phenotype <- check_phenotype(read_matrix(phenotype, "phenotype"))
  samples <- colnames(phenotype); feats <- rownames(phenotype)
  st <- sample_info(tissue, individual, samples)
  tr_all <- read_frame(traits, "traits")
  if (!ncol(tr_all)) stop("`traits` has no columns.", call. = FALSE)
  TR <- prepare_traits(rows_for(tr_all, samples, st$individual, "traits"), trait_type)
  gene <- rep(NA_character_, length(feats))
  if (!is.null(promoters)) {
    pr <- read_table(promoters, "promoters")
    if (!is.data.frame(pr) || !all(c("id", "gene") %in% names(pr))) stop("`promoters` must be a data frame with columns id and gene.", call. = FALSE)
    i <- match(feats, as.character(pr$id))
    if (anyNA(i)) stop("`promoters` has no row for ", sum(is.na(i)), " feature(s), e.g. ", feats[is.na(i)][1], ".", call. = FALSE)
    gene <- as.character(pr$gene[i])
  }
  if (!is.null(absolute)) {
    absolute <- read_matrix(absolute, "absolute")
    if (is.data.frame(absolute)) absolute <- as.matrix(absolute)
    if (!is.matrix(absolute) || !is.numeric(absolute) || is.null(rownames(absolute)) || is.null(colnames(absolute)))
      stop("`absolute` must be a numeric matrix with feature IDs as row names and sample IDs as column names.", call. = FALSE)
  }
  multi <- !is.null(st$tissue)
  tlev <- if (multi) st$levels else NA_character_
  if (multi && is_path(covariates) && !dir.exists(covariates)) covariates <- read_frame(covariates, "covariates")
  if (!multi) covariates <- read_frame(covariates, "covariates")
  if (mirror && is.null(major_score)) major_score <- rowMeans(phenotype, na.rm = TRUE)
  opt <- list(min_group = min_group, min_level_count = min_level_count, min_donors = min_donors, min_n_mid = min_n_mid,
              boundary_share = boundary_share, exclude_separated = exclude_separated, fdr = fdr, mirror = mirror && !is.null(promoters),
              mirror_tol = mirror_tol, abs_z = abs_z, control = control, ncores = ncores)
  res <- lapply(tlev, function(tt) {
    w <- if (multi) st$tissue == tt else rep(TRUE, length(samples))
    s <- samples[w]
    cv <- trait_covariates(covariates, tt)
    if (!is.null(cv)) cv <- rows_for(read_frame(cv, "covariates"), s, st$individual[w], "covariates")
    trait_one(phenotype[, s, drop = FALSE], gene, lapply(TR, function(t) { t$x <- t$x[w]; t }), cv, major_score, absolute, opt,
              label = if (multi) tt else NULL)
  })
  tab <- function(k) {
    d <- lapply(seq_along(res), function(i) { x <- res[[i]][[k]]; if (is.null(x) || !nrow(x)) return(NULL)
      if (multi) data.frame(tissue = rep(tlev[i], nrow(x)), x, check.names = FALSE, stringsAsFactors = FALSE) else x })
    b <- bind_df(d); if (!is.null(b)) rownames(b) <- NULL; b
  }
  out <- structure(list(results = tab("results"), features = tab("features"), traits = tab("traits"),
                        labels = if (is.null(absolute)) NULL else tab("labels"), mirror = tab("mirror"),
                        settings = list(tissues = if (multi) tlev else NULL, traits = names(TR), types = vapply(TR, `[[`, "", "type"),
                                        min_group = min_group, min_level_count = min_level_count, min_donors = min_donors, min_n_mid = min_n_mid,
                                        boundary_share = boundary_share, exclude_separated = exclude_separated, fdr = fdr, mirror = opt$mirror,
                                        mirror_tol = mirror_tol, abs_z = abs_z, control = control,
                                        mirror_labels = if (is.null(absolute))
                                          "v3m: exact-pair boundary events are p0 units at the promoter that goes silent (absolute activity not given)" else
                                          "absolute activity (option B): significant exact-pair boundary events reported as p1 at the partner when only its absolute activity changes"),
                        seconds = proc.time()[["elapsed"]] - t0), class = "bohr_trait")
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(out_dir)) stop("Cannot create `out_dir`: ", out_dir, call. = FALSE)
    for (k in c("results", "features", "traits", "labels")) if (!is.null(out[[k]]) && nrow(out[[k]])) write_tsv(out[[k]], file.path(out_dir, paste0(k, ".tsv")))
  }
  out
}

## one tissue: ph = features x the tissue's samples; TR = traits over the same samples; cv = covariate frame (rows = samples)
trait_one <- function(ph, gene_all, TR, cv, major_score, absolute, opt, label = NULL) {
  control <- opt$control
  tm <- tissue_mirror(ph, gene_all, major_score, opt$mirror_tol, mirror = opt$mirror)   # pairs among the features with data here
  has <- tm$has; pairs <- tm$pairs; flags <- tm$flags
  ph <- ph[has, , drop = FALSE]; gene <- gene_all[has]; feats <- rownames(ph); s_all <- colnames(ph)
  n_prom_gene <- tm$n_prom_gene[has]
  R <- list(); Fe <- list(); Tt <- list(); L <- list()
  for (tr in names(TR)) {
    x <- TR[[tr]]$x; ty <- TR[[tr]]$type
    cvt <- if (is.null(cv)) NULL else cv[, setdiff(names(cv), tr), drop = FALSE]                   # the trait is not a covariate
    ok <- !is.na(x) & (if (is.null(cvt) || !ncol(cvt)) TRUE else stats::complete.cases(cvt))
    donors <- s_all[ok]; g <- x[ok]
    trow <- list(trait = tr, type = ty, contrast_from = NA_real_, contrast_to = NA_real_, n_donors = length(donors), n_covariates = NA_integer_,
                 covariates = NA_character_, skip = NA_character_)
    if (length(donors) < control$min_samples || stats::sd(g) == 0) trow$skip <- "trait_constant_or_n"
    else if (ty == "binary" && min(table(g)) < opt$min_group) trow$skip <- "trait_rare"
    if (!is.na(trow$skip)) { Tt[[tr]] <- trow; next }
    C_all <- trait_design(if (is.null(cvt)) NULL else cvt[donors, , drop = FALSE], donors, opt$min_level_count)
    trow$n_covariates <- ncol(C_all) - 1L; trow$covariates <- paste(colnames(C_all)[-1], collapse = ";")
    cp <- code_predictor(g, ty, if (ty == "binary") TR[[tr]]$values else NULL)
    trow$contrast_from <- cp$contrast[1]; trow$contrast_to <- cp$contrast[2]
    Tt[[tr]] <- trow
    gp <- stats::setNames(cp$x, donors)
    chunks <- split(seq_along(feats), ceiling(seq_along(feats) / max(1L, ceiling(length(feats) / (4L * opt$ncores)))))
    rows <- unlist(par_apply(chunks, function(ix) lapply(ix, function(i)
      trait_feature(feats[i], gene[i], stats::setNames(ph[i, donors], donors), C_all, gp, ty == "binary", opt$min_group, control)),
      opt$ncores, preschedule = FALSE), recursive = FALSE)
    Fd <- rows_to_df(rows); Fd <- data.frame(trait = rep(tr, nrow(Fd)), Fd, stringsAsFactors = FALSE)
    Fe[[tr]] <- Fd
    ## units (summarize.R:51-66): mu needs mu_testable; p0 for BEINF/BEZI, p1 for BEINF/BEOI; finite P only
    fit_ok <- is.na(Fd$skip)
    mk <- function(comp, keep) {
      if (!any(keep)) return(NULL)
      d <- Fd[keep, , drop = FALSE]
      data.frame(trait = tr, promoter = d$promoter, gene = d$gene, component = comp, log_p = d[[paste0("log_p_", comp)]],
                 est = d[[paste0("est_", comp)]], se = d[[paste0("se_", comp)]], effect = d[[paste0("effect_", comp)]],
                 family = d$family, n = d$n, n0 = d$n0, n1 = d$n1, n_mid = d$n_mid, n_donors = length(donors), complete_sep = d$complete_sep,
                 stringsAsFactors = FALSE)
    }
    U <- bind_df(list(mk("mu", fit_ok & Fd$mu_testable %in% TRUE), mk("p0", fit_ok & Fd$family %in% c("BEINF", "BEZI")),
                      mk("p1", fit_ok & Fd$family %in% c("BEINF", "BEOI"))))
    if (is.null(U)) next
    U <- U[is.finite(U$log_p), , drop = FALSE]
    if (!nrow(U)) next
    U$p <- exp(U$log_p)
    U$boundary_share <- ifelse(U$component == "p0", U$n0 / U$n, ifelse(U$component == "p1", U$n1 / U$n, NA_real_))
    ## mirror rule v3 (mirror_v3.R:18, summarize.R:67-74)
    fi <- if (is.null(flags)) rep(NA_integer_, nrow(U)) else match(U$promoter, flags$promoter)
    get_f <- function(k, na) { v <- if (is.null(flags)) rep(na, nrow(U)) else flags[[k]][fi]; v[is.na(fi)] <- na; v }
    U$n_prom_gene <- n_prom_gene[match(U$promoter, feats)]
    U$boundary_exact <- get_f("boundary_exact", FALSE); U$interior_exact <- get_f("interior_exact", FALSE)
    ## role and partner are named for boundary_exact pairs only, as in the trait run's mirror table (mirror.R:53-54)
    U$mirror_role <- ifelse(U$boundary_exact, get_f("mirror_role", NA_character_), NA_character_)
    U$mirror_partner <- ifelse(U$boundary_exact, get_f("mirror_partner", NA_character_), NA_character_)
    U$dedup_drop <- ifelse(U$component == "mu", get_f("drop_mu", FALSE), U$component == "p1" & U$boundary_exact)
    swap <- c(mu = "mu", p0 = "p1", p1 = "p0")
    U$duplicate_of <- ifelse(U$dedup_drop, paste0(U$mirror_partner, ":", swap[U$component]), NA_character_)
    ## eligibility and BH within component (summarize.R:84-98)
    sh <- opt$boundary_share
    U$eligible <- !U$dedup_drop & (U$component != "mu" | U$n_mid >= opt$min_n_mid) &
      (U$component == "mu" | (U$boundary_share >= sh[1] & U$boundary_share <= sh[2])) &
      (U$component == "mu" | !opt$exclude_separated | !(U$complete_sep %in% TRUE)) & length(donors) >= opt$min_donors
    U$q <- NA_real_
    for (cc in c("mu", "p0", "p1")) { k <- which(U$eligible & U$component == cc); if (length(k)) U$q[k] <- stats::p.adjust(U$p[k], "BH") }
    U$significant <- U$eligible & !is.na(U$q) & U$q <= opt$fdr
    U$bh_family <- ifelse(U$eligible, U$component, NA_character_); U$abs_relabel <- FALSE
    ## labels of exact-pair boundary events by absolute activity (01_abs_labels.R, 03_apply_B.R)
    if (!is.null(absolute)) {
      lab <- abs_labels_trait(U, absolute, C_all, gp, donors, opt$abs_z, tr, label)
      if (!is.null(lab)) {
        for (j in which(lab$outcome == "O only")) {
          iS <- which(U$promoter == lab$S[j] & U$component == "p0"); iO <- which(U$promoter == lab$O[j] & U$component == "p1")
          if (length(iS) != 1L || length(iO) != 1L) { lab$changed[j] <- FALSE; next }
          for (k in c("eligible", "q", "significant")) U[[k]][iO] <- U[[k]][iS]
          U$dedup_drop[iO] <- FALSE; U$duplicate_of[iO] <- NA_character_; U$bh_family[iO] <- "p0"; U$abs_relabel[iO] <- TRUE
          U$eligible[iS] <- FALSE; U$q[iS] <- NA_real_; U$significant[iS] <- FALSE; U$dedup_drop[iS] <- TRUE
          U$duplicate_of[iS] <- paste0(lab$O[j], ":p1"); U$bh_family[iS] <- NA_character_; U$abs_relabel[iS] <- TRUE
        }
        L[[tr]] <- lab
      }
    }
    R[[tr]] <- U[, c("trait", "promoter", "gene", "component", "p", "log_p", "est", "se", "effect", "eligible", "q", "significant",
                     "dedup_drop", "duplicate_of", "n_prom_gene", "mirror_role", "mirror_partner", "boundary_exact", "interior_exact",
                     "family", "n", "n0", "n1", "n_mid", "boundary_share", "n_donors", "complete_sep", "bh_family", "abs_relabel")]
  }
  ord <- function(d) { if (is.null(d)) return(NULL); d <- d[order(match(d$trait, names(TR)), d$promoter, if (!is.null(d$component)) d$component else 0,
                                                                 method = "radix"), , drop = FALSE]; rownames(d) <- NULL; d }
  list(results = ord(bind_df(R)), features = bind_df(Fe), traits = rows_to_df(unname(Tt)), labels = bind_df(L), mirror = pairs)
}

## 01_abs_labels.R:39-46 (ols_z) and :56-60, 77-81: absolute activity of S and O on the covariates and the coded trait, over the
## donors; |z| >= abs_z. S only / both / neither / not estimable -> p0(S); O only -> p1(O). Significant p0 units of
## boundary_exact pairs only (03_apply_B.R:22).
abs_labels_trait <- function(U, absolute, C_all, gp, donors, abs_z, tr, label) {
  ev <- which(U$component == "p0" & U$boundary_exact %in% TRUE & U$significant & !is.na(U$mirror_partner))
  if (!length(ev)) return(NULL)
  where <- if (is.null(label)) "" else paste0(" (tissue ", label, ")")
  S <- U$promoter[ev]; O <- U$mirror_partner[ev]
  miss <- setdiff(c(S, O), rownames(absolute)); if (length(miss)) stop("`absolute` has no row for feature ", miss[1], where, ".", call. = FALSE)
  miss <- setdiff(donors, colnames(absolute)); if (length(miss)) stop("`absolute` has no column for sample ", miss[1], where, ".", call. = FALSE)
  X <- cbind(C_all, g = unname(gp[donors]))
  zfit <- function(P) { Y <- t(absolute[P, donors, drop = FALSE])
    if (anyNA(Y) || any(!is.finite(Y))) stop("`absolute` must not contain missing values for the donors used", where, ".", call. = FALSE)
    q <- qr(X); p <- ncol(X); if (q$rank < p) return(rep(NA_real_, length(P)))
    B <- qr.coef(q, Y); Rr <- qr.resid(q, Y); s2 <- colSums(Rr^2) / (nrow(X) - p); v <- chol2inv(qr.R(q))[p, p]
    unname(B[p, ] / sqrt(v * s2)) }
  zS <- zfit(S); zO <- zfit(O)
  rs <- abs(zS) >= abs_z; ro <- abs(zO) >= abs_z
  out <- ifelse(is.na(zS) | is.na(zO), "not estimable", ifelse(rs & ro, "both", ifelse(rs, "S only", ifelse(ro, "O only", "neither"))))
  data.frame(trait = tr, S = S, O = O, gene = U$gene[ev], z_S = zS, z_O = zO, outcome = out,
             label = ifelse(out == "O only", paste0(O, ":p1"), paste0(S, ":p0")), changed = out == "O only", stringsAsFactors = FALSE)
}

#' @export
print.bohr_trait <- function(x, ...) {
  cat("BOHR trait regression:", length(x$settings$traits), "trait(s)", if (!is.null(x$settings$tissues)) paste0("in ", length(x$settings$tissues), " tissues"), "\n")
  R <- x$results
  if (is.null(R) || !nrow(R)) { cat("  no testable unit\n"); return(invisible(x)) }
  key <- if (is.null(R$tissue)) R$trait else paste(R$tissue, R$trait, sep = " | ")
  ks <- unique(key)
  tab <- data.frame(tissue_trait = ks, tests = vapply(ks, function(k) sum(key == k), 0L), eligible = vapply(ks, function(k) sum(R$eligible[key == k]), 0L),
                    sig_mu = vapply(ks, function(k) sum(R$significant[key == k] & R$component[key == k] == "mu"), 0L),
                    sig_p0 = vapply(ks, function(k) sum(R$significant[key == k] & R$component[key == k] == "p0"), 0L),
                    sig_p1 = vapply(ks, function(k) sum(R$significant[key == k] & R$component[key == k] == "p1"), 0L),
                    row.names = NULL, stringsAsFactors = FALSE)
  print(utils::head(tab, 20), row.names = FALSE)
  if (nrow(tab) > 20) cat("  ...", nrow(tab) - 20, "more\n")
  cat(sprintf("q within tissue x trait x component; significant at q <= %g. Estimates per contrast (+1 SD, or low -> high).\n", x$settings$fdr))
  cat("Mirror labels:", x$settings$mirror_labels, "\n")
  if (!is.null(x$labels) && nrow(x$labels)) cat("  significant exact-pair boundary events judged:", nrow(x$labels), "; reported as p1 at the partner:", sum(x$labels$changed), "\n")
  invisible(x)
}
