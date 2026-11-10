#' BOHR regression with a formula
#'
#' Fits the BOHR model of a proportion response with exact 0s and 1s, `y ~ x +
#' covariates`, and tests one term (`test`) for association with each
#' component of the response: the probability of 0 (p0), the probability of 1
#' (p1) and the mean of the values in between (mu). This is the formula
#' interface to [bohr_fit()] and [bohr_test()]: the statistics are those of
#' [bohr_test()] on the same data, and the fitted model is available through
#' [coef()][coef.bohr], [vcov()][vcov.bohr], [logLik()][logLik.bohr],
#' [predict()][predict.bohr] and [summary()][summary.bohr].
#'
#' The tested term must be a main effect: a single numeric column, or a factor
#' (or logical or character variable) with exactly two levels, which is coded
#' as an indicator of the second level. Every other term of the formula is a
#' covariate; factors among them are coded by [stats::model.matrix()]. The
#' model always has an intercept. Offsets and weights are not supported.
#'
#' The fitting rules of the BOHR paper ([bohr_fit()], set with
#' [bohr_control()]) can change the data or the model: missing values are left
#' out, a response is skipped when too few observations remain or almost all
#' lie at one boundary, the observations at a boundary holding 5% of them or
#' fewer are removed (and that component is not modelled), covariate columns
#' collinear with others are dropped, and the beta part is not fitted when too
#' few values lie between 0 and 1. `print()` and `summary()` report every rule
#' that was applied, with its counts, and notes on separation and convergence.
#'
#' **Coefficients.** The model has three linear predictors: the logit of mu
#' (beta part, with a constant precision phi on the log scale), and \eqn{\log
#' \nu} and \eqn{\log \tau}, the log odds of 0 and of 1 against a value in
#' between (family `"BEINF"`), or the logit of being at 0 (`"BEZI"`) or at 1
#' (`"BEOI"`). Coefficients are on the scale of the covariates and per unit of
#' the tested term, as in [stats::glm()]. For p0 and p1 they are the
#' Firth-penalised estimates whose tests are reported (the unpenalised ones
#' with `bohr_control(firth = FALSE)`), as in the `est_p0` and `est_p1`
#' columns of [bohr_test()]. The effects (`effect_*`) are those of
#' [bohr_test()]: the average change of each component between the two ends of
#' the contrast, from the unpenalised fit, as in the paper.
#'
#' @param formula A two-sided formula, `y ~ x + covariates`, with the response
#'   in \[0, 1\].
#' @param data A data frame (or list or environment) holding the variables of
#'   `formula`.
#' @param test Name of the term to test, as written in the formula (for
#'   example `"x"` or `"log(dose)"`). Defaults to the first term on the right
#'   of the formula.
#' @param type How the tested term is treated for the effect contrast, as in
#'   [bohr_test()]: `"auto"` (default; binary for exactly two distinct values,
#'   otherwise continuous), `"dosage"`, `"binary"` or `"continuous"`. Set
#'   `"dosage"` for allele dosages. A two-level factor is binary.
#' @param contrast Optional `c(from, to)` on the scale of the tested term;
#'   the effects compare `to` with `from`.
#' @param control Fitting rules, optimiser limits and the choice of the p0/p1
#'   test, from [bohr_control()]; the defaults are the values of the BOHR
#'   paper.
#' @param na.action How missing values are handled before fitting; default
#'   [stats::na.omit()] (observations with a missing value in any variable of
#'   the formula are left out, and reported).
#' @return An object of class `"bohr"` with, among others, `result` (the
#'   one-row data frame of [bohr_test()] with every engine statistic),
#'   `fit` (the [bohr_fit()] object), `rules` and `notes` (the rules applied
#'   and notes, as data frames with columns `rule` and `message`), `test` (the
#'   tested term, its coefficient name, type and contrast), `call`, `terms`
#'   and `model` (the model frame).
#' @seealso [bohr_test()] for many predictors of one response,
#'   [bohr_fit()], [bohr_control()], [simulate_response()].
#' @examples
#' set.seed(3)
#' d <- simulate_response(250, mu = 0.4, p0 = 0.25, p1 = 0.1,
#'                        effect_p0 = 0.08, effect_mu = 0.04)
#' d$age <- rnorm(250, 50, 10)
#' d$sex <- factor(sample(c("female", "male"), 250, replace = TRUE))
#' m <- bohr(y ~ x + age + sex, data = d, test = "x", type = "dosage")
#' m
#' summary(m)
#' coef(m, "p0")
#' head(predict(m, type = "mean"))
#' ## a two-level factor as the tested term
#' bohr(y ~ sex + x, data = d, test = "sex")
#' @export
bohr <- function(formula, data, test = NULL, type = c("auto", "dosage", "binary", "continuous"), contrast = NULL,
                 control = bohr_control(), na.action = stats::na.omit) {
  cl <- match.call()
  type <- match.arg(type); check_control(control)
  if (!inherits(formula, "formula") || length(formula) != 3L) stop("`formula` must be a two-sided formula, y ~ x + covariates.", call. = FALSE)
  if (!is.null(contrast) && (!is.numeric(contrast) || length(contrast) != 2L || anyNA(contrast) || !all(is.finite(contrast)) ||
                             contrast[1] == contrast[2])) stop("`contrast` must be two different finite numbers c(from, to).", call. = FALSE)
  if (missing(data)) data <- environment(formula)
  mf <- stats::model.frame(formula, data = data, na.action = na.action, drop.unused.levels = TRUE)
  tt <- attr(mf, "terms")
  if (attr(tt, "intercept") != 1L) stop("BOHR models always include an intercept; remove `- 1` or `0 +` from the formula.", call. = FALSE)
  if (!is.null(stats::model.offset(mf))) stop("Offsets are not supported.", call. = FALSE)
  labs <- attr(tt, "term.labels")
  if (!length(labs)) stop("The formula has no term to test.", call. = FALSE)
  if (is.null(test)) test <- labs[1]
  if (!is.character(test) || length(test) != 1L || is.na(test)) stop("`test` must name one term of the formula.", call. = FALSE)
  k <- match(test, labs)
  if (is.na(k)) stop("`test = \"", test, "\"` is not a term of the formula; its terms are: ", paste(labs, collapse = ", "), ".", call. = FALSE)
  fac <- attr(tt, "factors")
  if (attr(tt, "order")[k] != 1L) stop("The tested term must be a main effect, not an interaction.", call. = FALSE)
  vrow <- which(fac[, k] > 0)
  if (any(fac[vrow, -k, drop = FALSE] > 0)) stop("The tested term `", test, "` also appears in an interaction; BOHR tests one main effect.", call. = FALSE)
  vname <- rownames(fac)[vrow]

  ## response
  y <- stats::model.response(mf)
  if (!is.null(dim(y)) || !is.numeric(y)) stop("The response must be a numeric vector in [0, 1].", call. = FALSE)
  check_activity(y, "response")
  y <- as.vector(y)

  ## the tested term: a single numeric column or a two-level factor
  xv <- mf[[vname]]; levs <- NULL
  if (is.logical(xv)) xv <- factor(xv, levels = c(FALSE, TRUE))
  if (is.character(xv)) xv <- factor(xv)
  if (is.factor(xv)) {
    levs <- levels(droplevels(xv))
    if (length(levs) != 2L) stop("The tested term must be a single numeric column or a factor with two levels; `", test, "` has ",
                                 length(levs), " level", if (length(levs) != 1L) "s", ".", call. = FALSE)
    if (type %in% c("dosage", "continuous")) stop("A two-level factor is tested as binary; use type = \"auto\" or \"binary\".", call. = FALSE)
    x <- as.numeric(xv == levs[2]); xname <- paste0(test, levs[2])
  } else if (is.numeric(xv) && (is.null(dim(xv)) || (length(dim(xv)) == 2L && ncol(xv) == 1L))) {
    x <- as.numeric(xv); xname <- test
  } else stop("The tested term must be a single numeric column or a factor with two levels; `", test, "` is ",
              if (is.matrix(xv)) paste("a matrix with", ncol(xv), "columns") else paste0("of class ", class(xv)[1]), ".", call. = FALSE)
  if (any(is.infinite(x))) stop("The tested term must not contain infinite values.", call. = FALSE)

  ## covariates: the other terms, coded by model.matrix(), without the intercept
  cov_terms <- NULL; Z <- NULL; ctr <- NULL
  if (length(labs) > 1L) {
    cov_terms <- stats::drop.terms(tt, k, keep.response = TRUE)
    mm <- stats::model.matrix(cov_terms, mf)
    ctr <- attr(mm, "contrasts")
    Z <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
    attr(Z, "assign") <- NULL; attr(Z, "contrasts") <- NULL
  }
  n_na <- if (is.null(attr(mf, "na.action"))) 0L else length(attr(mf, "na.action"))
  yf <- y; yf[is.na(x)] <- NA                                    # as bohr_test() does for a vector x
  fit <- bohr_fit(yf, Z, control)
  tc <- test_column(fit, x, xname, type, contrast, control, all_stats = TRUE, keep_fits = TRUE)
  result <- rows_to_df(list(tc$row)); rownames(result) <- NULL

  ## coefficients on the covariates' scale, per unit of the tested term
  Zm <- if (is.null(Z)) matrix(numeric(0), length(y), 0L) else Z
  cf <- model_coefs(fit, tc, Zm, xname, control$firth)

  ## rules and notes
  R <- fit_rules(fit)
  if (n_na > 0) R <- rbind(rule_df("na_action", sprintf("%s with missing values left out by na.action (%s).",
                                                         plural(n_na, "observation"), class(attr(mf, "na.action"))[1])), R)
  TR <- test_rules(tc$row, fit$family, control$firth)
  note_codes <- c("complete_separation", "separation_p0", "separation_p1", "separation_mu", "not_converged_mu", "not_converged_p0", "not_converged_p1")
  rules <- rbind(R, TR[!TR$rule %in% note_codes, , drop = FALSE]); notes <- TR[TR$rule %in% note_codes, , drop = FALSE]
  rownames(rules) <- NULL; rownames(notes) <- NULL

  structure(list(call = cl, formula = formula, terms = tt, cov_terms = cov_terms, contrasts = ctr,
                 xlevels = stats::.getXlevels(tt, mf), model = mf, na.action = attr(mf, "na.action"),
                 test = list(term = test, variable = vname, name = xname, levels = levs, type = tc$type,
                             contrast = if (is.null(tc$coding)) c(NA_real_, NA_real_) else tc$coding$contrast,
                             offset = tc$coding$offset, per_unit = tc$coding$per_unit, user_contrast = !is.null(contrast)),
                 x = x, covariates = Zm, y = y, fit = fit, result = result,
                 coefficients = cf$coef, vcov = cf$vcov, cross = list(mu_phi = cf$cov_mu_phi, p0_p1 = cf$cov_p0_p1),
                 loglik = cf$loglik, df = cf$df,
                 rules = rules, notes = notes, control = control),
            class = "bohr")
}

## ---- coefficients on the user's scale ---------------------------------------------------------------
## The engine fits on the per-feature orthonormal design Cs (fit$design) and the coded predictor
## x' = offset + per_unit * x. Cs spans the columns kept of C_all = cbind(1, scale(Z)) on the samples used
## (orth_design), so Cs = C_kept %*% M with M = qr.solve(C_kept, Cs); undoing the standardisation
## (beta_j = theta_j / s_j; intercept - sum theta_j m_j / s_j) and the coding of x gives the linear map L from
## c(theta, gamma) to c(intercept, covariates kept, x).
coef_map <- function(fit, Z, coding) {
  Cs <- fit$design; s <- fit$samples; r <- ncol(Cs)
  C_all <- make_design(Z)
  dropped <- attr(Cs, "dropped"); kept <- setdiff(colnames(C_all), dropped)
  M <- qr.solve(C_all[s, kept, drop = FALSE], Cs)
  U <- diag(r); kc <- kept[-1]
  for (j in seq_along(kc)) {
    m <- fit$input$center[[kc[j]]]; sd <- fit$input$scale[[kc[j]]]
    U[1 + j, 1 + j] <- 1 / sd; U[1, 1 + j] <- -m / sd
  }
  L <- matrix(0, r + 1L, r + 1L)
  L[seq_len(r), seq_len(r)] <- U %*% M
  L[1, r + 1L] <- coding$offset; L[r + 1L, r + 1L] <- coding$per_unit
  list(L = L, kept = kc, dropped = dropped)
}
## map one block (coefficients of c(Cs, x') of one linear predictor, with its covariance) to named user-scale
## coefficients: (Intercept), the tested term, every covariate (NA when dropped as collinear)
map_block <- function(theta, V, mp, xname, cov_names) {
  nm <- c("(Intercept)", xname, cov_names)
  b <- stats::setNames(rep(NA_real_, length(nm)), nm)
  r1 <- nrow(mp$L); pos <- c("(Intercept)", mp$kept, xname)          # order of the rows of L
  b[pos] <- drop(mp$L %*% theta)
  W <- matrix(NA_real_, length(nm), length(nm), dimnames = list(nm, nm))
  if (!is.null(V) && all(dim(V) == r1)) W[pos, pos] <- mp$L %*% V %*% t(mp$L)
  list(coef = b, vcov = W)
}
model_coefs <- function(fit, tc, Z, xname, firth) {
  out <- list(coef = list(), vcov = list(), loglik = NA_real_, df = NA_integer_)
  if (!is.na(fit$skip) || is.null(tc$fits)) return(out)
  mp <- coef_map(fit, Z, tc$coding); r1 <- nrow(mp$L); cov_names <- colnames(Z)
  fb <- tc$fits$beta
  if (!is.null(fb)) {
    V <- fb$vcov
    blk <- map_block(fb$coef[seq_len(r1)], if (is.null(V)) NULL else V[seq_len(r1), seq_len(r1), drop = FALSE], mp, xname, cov_names)
    out$coef$mu <- blk$coef; out$vcov$mu <- blk$vcov
    out$coef$phi <- c(`(Intercept)` = fb$coef[[r1 + 1L]])
    out$vcov$phi <- matrix(if (is.null(V)) NA_real_ else V[r1 + 1L, r1 + 1L], 1, 1, dimnames = list("(Intercept)", "(Intercept)"))
    if (!is.null(V)) {                                                   # covariance of mu and log(phi)
      cvp <- drop(mp$L %*% V[seq_len(r1), r1 + 1L])
      out$cov_mu_phi <- stats::setNames(rep(NA_real_, length(blk$coef)), names(blk$coef))
      out$cov_mu_phi[c("(Intercept)", mp$kept, xname)] <- cvp
    }
  }
  fc <- if (isTRUE(firth)) tc$fits$firth else tc$fits$cat
  fam <- fit$family
  if (!is.null(fc) && fam != "BE") {
    th <- fc$coef; V <- fc$vcov
    eqs <- switch(fam, BEINF = c("p0", "p1"), BEZI = "p0", BEOI = "p1")
    for (e in seq_along(eqs)) {
      i <- (e - 1L) * r1 + seq_len(r1)
      blk <- map_block(th[i], if (is.null(V)) NULL else V[i, i, drop = FALSE], mp, xname, cov_names)
      out$coef[[eqs[e]]] <- blk$coef; out$vcov[[eqs[e]]] <- blk$vcov
    }
    if (fam == "BEINF" && !is.null(V)) {                                 # cross-covariance of the two equations
      i0 <- seq_len(r1); i1 <- r1 + seq_len(r1); pos <- c("(Intercept)", mp$kept, xname)
      nm <- names(out$coef$p0); X <- matrix(NA_real_, length(nm), length(nm), dimnames = list(nm, nm))
      X[pos, pos] <- mp$L %*% V[i0, i1, drop = FALSE] %*% t(mp$L); out$cov_p0_p1 <- X
    }
  }
  ## maximised (unpenalised) log-likelihood of the full model and its number of parameters
  rr <- tc$row
  if (!is.null(fb) && is.finite(rr$ll_beta_full)) {
    out$loglik <- rr$ll_beta_full + if (fam != "BE") rr$ll_cat_full else 0
    out$df <- length(fb$coef) + if (fam != "BE" && !is.null(tc$fits$cat)) length(tc$fits$cat$coef) else 0L
  }
  out
}

## ---- methods ---------------------------------------------------------------------------------------
COMPONENTS <- c("mu", "phi", "p0", "p1")

#' Methods for BOHR regression fits
#'
#' Extract the coefficients, their covariance, the log-likelihood, the number
#' of observations and predictions from a [bohr()] fit.
#'
#' The model has up to four linear predictors (components): `"mu"`, the logit
#' of the mean of the values between 0 and 1; `"phi"`, the log of the
#' precision of the beta part (intercept only); `"p0"` and `"p1"`, the log odds
#' of 0 and of 1 against a value in between (family `"BEINF"`), or the logit of
#' being at the boundary (`"BEZI"`: p0; `"BEOI"`: p1). Components that are not
#' in the model (removed by the boundary rule or not fitted) are left out.
#' Coefficients are per unit of the tested term and of each covariate;
#' covariate columns dropped as collinear have `NA`. The p0 and p1
#' coefficients and their covariance are those of the Firth-penalised fit
#' (the default) or of the unpenalised fit (`bohr_control(firth = FALSE)`).
#'
#' `logLik()` is the maximised (unpenalised) log-likelihood of the full model,
#' over the observations used after the fitting rules (`nobs()`), with `df`
#' the number of estimated parameters; `NA` when the beta part is not fitted.
#'
#' `predict()` returns, for each observation, the component's value: `"mu"`,
#' the mean of the values in between; `"p0"` and `"p1"`, the probabilities of
#' 0 and of 1 (0 for a component not in the model); `"mean"`, the expected
#' response \eqn{p_1 + (1 - p_0 - p_1)\mu}. Without `newdata`, predictions are
#' for the observations of the model frame (also those removed by the boundary
#' rule), padded with `NA` for those removed by `na.action = na.exclude`.
#'
#' @param object,x A `"bohr"` object.
#' @param component `NULL` (all components; names are `"component:term"`) or
#'   one of `"mu"`, `"phi"`, `"p0"` and `"p1"`.
#' @param newdata Optional data frame in which to look for the variables of
#'   the model; by default the data of the fit.
#' @param type Component to predict: `"mean"` (default), `"mu"`, `"p0"` or
#'   `"p1"`.
#' @param na.action How missing values in `newdata` are handled; default
#'   [stats::na.pass()].
#' @param digits Number of significant digits printed.
#' @param ... Not used.
#' @return `coef()`: a named numeric vector. `vcov()`: a matrix. `logLik()`: a
#'   `"logLik"` object. `nobs()`: an integer. `predict()`: a numeric vector.
#'   `summary()`: an object of class `"summary.bohr"`, a list with one
#'   coefficient table per component (`tables`), the tested term's statistics
#'   (`test_table`), the rules applied and notes.
#' @seealso [bohr()].
#' @examples
#' set.seed(4)
#' d <- simulate_response(200, p0 = 0.2, p1 = 0.15, effect_p1 = 0.06)
#' m <- bohr(y ~ x, data = d, type = "dosage")
#' coef(m)
#' sqrt(diag(vcov(m, "p1")))
#' logLik(m)
#' nobs(m)
#' predict(m, newdata = data.frame(x = 0:2), type = "p1")
#' @name bohr-methods
NULL

#' @rdname bohr-methods
#' @export
coef.bohr <- function(object, component = NULL, ...) {
  cf <- object$coefficients
  if (!is.null(component)) {
    component <- match.arg(component, COMPONENTS)
    return(if (is.null(cf[[component]])) numeric(0) else cf[[component]])
  }
  cmp <- intersect(COMPONENTS, names(cf))
  unlist(lapply(cmp, function(k) stats::setNames(cf[[k]], paste0(k, ":", names(cf[[k]])))))
}

#' @rdname bohr-methods
#' @export
vcov.bohr <- function(object, component = NULL, ...) {
  V <- object$vcov
  if (!is.null(component)) {
    component <- match.arg(component, COMPONENTS)
    return(if (is.null(V[[component]])) matrix(numeric(0), 0, 0) else V[[component]])
  }
  cmp <- intersect(COMPONENTS, names(V))
  at <- function(k) paste0(k, ":", rownames(V[[k]]))
  nm <- unlist(lapply(cmp, at))
  ## the beta part (mu, phi) and the categorical part (p0, p1) share no parameters: their cross-covariance is 0;
  ## within a part the cross blocks come from the joint fit (NA when its information matrix was singular)
  W <- matrix(0, length(nm), length(nm), dimnames = list(nm, nm))
  for (part in list(c("mu", "phi"), c("p0", "p1"))) {
    p <- intersect(part, cmp); if (length(p) == 2L) W[unlist(lapply(p, at)), unlist(lapply(p, at))] <- NA
  }
  for (k in cmp) W[at(k), at(k)] <- V[[k]]
  cr <- object$cross
  if (all(c("mu", "phi") %in% cmp) && !is.null(cr$mu_phi)) { W[at("mu"), at("phi")] <- cr$mu_phi; W[at("phi"), at("mu")] <- cr$mu_phi }
  if (all(c("p0", "p1") %in% cmp) && !is.null(cr$p0_p1)) { W[at("p0"), at("p1")] <- cr$p0_p1; W[at("p1"), at("p0")] <- t(cr$p0_p1) }
  W
}

#' @rdname bohr-methods
#' @export
logLik.bohr <- function(object, ...) {
  structure(object$loglik, df = object$df, nobs = nobs.bohr(object), class = "logLik")
}

#' @rdname bohr-methods
#' @export
nobs.bohr <- function(object, ...) if (is.na(object$fit$skip)) as.integer(object$fit$n) else 0L

#' @rdname bohr-methods
#' @export
predict.bohr <- function(object, newdata = NULL, type = c("mean", "mu", "p0", "p1"), na.action = stats::na.pass, ...) {
  type <- match.arg(type)
  if (is.null(newdata)) { x <- object$x; Z <- object$covariates; n <- length(x); rn <- rownames(object$model) }
  else {
    tt <- stats::delete.response(object$terms)
    mfn <- stats::model.frame(tt, newdata, na.action = na.action, xlev = object$xlevels)
    rn <- rownames(mfn)
    xv <- mfn[[object$test$variable]]
    if (!is.null(object$test$levels)) {
      xv <- as.character(xv)
      bad <- !is.na(xv) & !xv %in% object$test$levels
      if (any(bad)) stop("`newdata` has levels of `", object$test$term, "` not in the fit.", call. = FALSE)
      x <- as.numeric(xv == object$test$levels[2])
    } else x <- as.numeric(xv)
    n <- length(x)
    Z <- if (is.null(object$cov_terms)) matrix(numeric(0), n, 0L) else {
      mm <- stats::model.matrix(stats::delete.response(object$cov_terms), mfn, contrasts.arg = object$contrasts)
      mm[, colnames(mm) != "(Intercept)", drop = FALSE]
    }
    if (ncol(Z) != ncol(object$covariates)) stop("`newdata` gives different covariate columns from the fit.", call. = FALSE)
  }
  X <- cbind(1, x, Z)
  lp <- function(k) { b <- object$coefficients[[k]]; if (is.null(b)) return(NULL); b[is.na(b)] <- 0; drop(X %*% b) }
  fam <- object$fit$family
  na_out <- rep(NA_real_, n)
  if (!is.na(object$fit$skip) || (is.null(object$coefficients$mu) && type %in% c("mu", "mean"))) out <- na_out
  else {
    mu <- if (!is.null(object$coefficients$mu)) stats::plogis(lp("mu")) else na_out
    p0 <- p1 <- rep(0, n)
    if (fam == "BEINF") {
      e0 <- lp("p0"); e1 <- lp("p1")
      if (is.null(e0) || is.null(e1)) p0 <- p1 <- na_out
      else { m <- pmax(0, e0, e1); den <- exp(-m) + exp(e0 - m) + exp(e1 - m); p0 <- exp(e0 - m) / den; p1 <- exp(e1 - m) / den }
    } else if (fam == "BEZI") { e0 <- lp("p0"); p0 <- if (is.null(e0)) na_out else stats::plogis(e0) }
    else if (fam == "BEOI") { e1 <- lp("p1"); p1 <- if (is.null(e1)) na_out else stats::plogis(e1) }
    out <- switch(type, mu = mu, p0 = p0, p1 = p1, mean = p1 + (1 - p0 - p1) * mu)
  }
  names(out) <- rn
  if (is.null(newdata)) out <- stats::napredict(object$na.action, out)
  out
}

## ---- printing --------------------------------------------------------------------------------------
FAMILY_TEXT <- c(BEINF = "0s, 1s and values in between", BEZI = "0s and values in between; p1 not modelled",
                 BEOI = "1s and values in between; p0 not modelled", BE = "values in between only; mu only")
family_line <- function(fam) if (is.na(fam)) "none (skipped)" else paste0(fam, " (", FAMILY_TEXT[[fam]], ")")
## the tests of the modelled components, in words
tests_text <- function(fam, firth) {
  b <- switch(fam, BEINF = "p0 and p1", BEZI = "p0", BEOI = "p1", NULL)
  paste0("mu, F-corrected likelihood ratio", if (!is.null(b)) paste0("; ", b, ", ", if (firth) "Firth-penalised" else "unpenalised",
                                                                     " likelihood ratio"))
}
cat_wrapped <- function(...) cat(strwrap(paste0(...), width = max(40L, getOption("width") - 2L)), sep = "\n")
contrast_text <- function(tst) {
  ctr <- tst$contrast
  if (anyNA(ctr)) return(tst$type)
  if (!is.null(tst$levels)) return(sprintf("binary factor; effects compare %s = %s with %s = %s", tst$term, tst$levels[2], tst$term, tst$levels[1]))
  f <- function(v) format(signif(v, 4))
  what <- if (tst$type == "continuous" && !isTRUE(tst$user_contrast)) c(" (mean + 1 SD)", " (mean)") else c("", "")
  sprintf("%s; effects compare %s = %s%s with %s = %s%s", tst$type, tst$term, f(ctr[2]), what[1], tst$term, f(ctr[1]), what[2])
}
fmt_col <- function(v, digits, p = FALSE) {
  out <- rep("", length(v)); ok <- !is.na(v)
  if (any(ok)) out[ok] <- if (p) vapply(v[ok], function(q) format.pval(signif(q, digits), digits = digits), "")
                          else format(signif(v[ok], digits), digits = digits)
  out
}
## one row per modelled component for the tested term: estimate, SE, statistic, P, effect, n
test_table <- function(object) {
  r <- object$result; fam <- object$fit$family
  if (is.na(fam)) return(NULL)
  cmp <- c("mu", if (fam %in% c("BEINF", "BEZI")) "p0", if (fam %in% c("BEINF", "BEOI")) "p1")
  data.frame(estimate = vapply(cmp, function(k) r[[paste0("est_", k)]], 0), SE = vapply(cmp, function(k) r[[paste0("se_", k)]], 0),
             statistic = vapply(cmp, function(k) r[[paste0("stat_", k)]], 0), p_value = vapply(cmp, function(k) r[[paste0("p_", k)]], 0),
             effect = vapply(cmp, function(k) r[[paste0("effect_", k)]], 0),
             n = vapply(cmp, function(k) if (k == "mu") r$n_mid else r$n, 0L), row.names = cmp)
}
print_test_table <- function(tab, digits) {
  out <- data.frame(estimate = fmt_col(tab$estimate, digits), SE = fmt_col(tab$SE, digits), statistic = fmt_col(tab$statistic, digits),
                    `P value` = fmt_col(tab$p_value, digits, p = TRUE), effect = fmt_col(tab$effect, digits), n = tab$n,
                    row.names = rownames(tab), check.names = FALSE)
  print(out, right = TRUE)
}
header_lines <- function(object) {
  f <- object$fit
  cat("BOHR regression: ", paste(deparse(object$formula, width.cutoff = 500L), collapse = " "), "\n", sep = "")
  cat("Tested term: ", object$test$name, " (", contrast_text(object$test), ")\n", sep = "")
  cat("Family: ", family_line(f$family), "\n", sep = "")
  if (is.na(f$skip)) cat(sprintf("Observations used: %d (%d at 0, %d at 1, %d in between)\n", f$n, f$n0, f$n1, f$n_mid))
}
test_footer <- function(object) {
  firth <- isTRUE(object$control$firth)
  fam <- object$fit$family
  b <- switch(fam, BEINF = "; p0 and p1 as log odds against a value in between", BEZI = "; p0 on the logit scale",
              BEOI = "; p1 on the logit scale", "")
  cat_wrapped(sprintf("Estimates per unit of %s: mu on the logit scale%s.", object$test$name, b))
  cat_wrapped("Tests: ", tests_text(fam, firth), ".")
}

#' @rdname bohr-methods
#' @export
print.bohr <- function(x, digits = 4L, ...) {
  header_lines(x)
  tab <- test_table(x)
  if (!is.null(tab) && is.na(x$result$error)) { cat("\n"); print_test_table(tab, digits); test_footer(x) }
  if (nrow(x$rules) || nrow(x$notes)) cat("\n")
  print_rules(x$rules, x$notes)
  invisible(x)
}

#' @rdname bohr-methods
#' @export
summary.bohr <- function(object, ...) {
  tabs <- list(); f <- object$fit
  for (k in intersect(c("mu", "p0", "p1"), names(object$coefficients))) {
    b <- object$coefficients[[k]]; se <- sqrt(diag(object$vcov[[k]]))
    st <- pv <- ef <- rep(NA_real_, length(b)); i <- match(object$test$name, names(b))
    st[i] <- object$result[[paste0("stat_", k)]]; pv[i] <- object$result[[paste0("p_", k)]]; ef[i] <- object$result[[paste0("effect_", k)]]
    tabs[[k]] <- data.frame(Estimate = unname(b), `Std. Error` = unname(se), Statistic = st, `P value` = pv, Effect = ef,
                            row.names = names(b), check.names = FALSE)
  }
  phi <- if (!is.null(object$coefficients$phi)) c(estimate = object$coefficients$phi[[1]], SE = sqrt(object$vcov$phi[1, 1])) else NULL
  structure(list(call = object$call, formula = object$formula, test = object$test, fit = f, tables = tabs, phi = phi,
                 n = c(mu = f$n_mid, p0 = f$n, p1 = f$n), test_table = test_table(object), result = object$result,
                 loglik = logLik.bohr(object), rules = object$rules, notes = object$notes, control = object$control),
            class = "summary.bohr")
}

#' @rdname bohr-methods
#' @export
print.summary.bohr <- function(x, digits = 4L, ...) {
  header_lines(x)
  TITLE <- c(mu = "mean of the values in between (logit link)",
             p0 = if (identical(x$fit$family, "BEINF")) "probability of 0 (log odds of 0 against a value in between)" else "probability of 0 (logit link)",
             p1 = if (identical(x$fit$family, "BEINF")) "probability of 1 (log odds of 1 against a value in between)" else "probability of 1 (logit link)")
  for (k in names(x$tables)) {
    t <- x$tables[[k]]
    cat(sprintf("\n%s: %s; %s\n", k, TITLE[[k]], plural(x$n[[k]], "observation")))
    out <- data.frame(Estimate = fmt_col(t$Estimate, digits), `Std. Error` = fmt_col(t[["Std. Error"]], digits),
                      Statistic = fmt_col(t$Statistic, digits), `P value` = fmt_col(t[["P value"]], digits, p = TRUE),
                      Effect = fmt_col(t$Effect, digits), row.names = rownames(t), check.names = FALSE)
    print(out, right = TRUE)
    if (k == "mu" && !is.null(x$phi)) cat(sprintf("Precision phi: %s (log phi %s, SE %s)\n", format(signif(exp(x$phi[["estimate"]]), digits)),
                                                  format(signif(x$phi[["estimate"]], digits)), format(signif(x$phi[["SE"]], digits))))
  }
  if (length(x$tables)) {
    cat("\n")
    cat_wrapped("Statistic, P value and effect: tested term only. Tests: ", tests_text(x$fit$family, isTRUE(x$control$firth)),
                ". Effect: average change of the component for the contrast.")
    if (is.finite(x$loglik)) cat(sprintf("Log-likelihood (unpenalised, full model): %s on %d parameters\n", format(signif(as.numeric(x$loglik), 8)), attr(x$loglik, "df")))
  }
  if (nrow(x$rules) || nrow(x$notes)) cat("\n")
  print_rules(x$rules, x$notes)
  invisible(x)
}
