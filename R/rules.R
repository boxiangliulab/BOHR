## Messages for every fitting rule that changed the data or the model (stage 2b; nothing happens silently),
## and notes on conditions of the fits. Internal. Each function returns a data frame with columns
## `rule` (a short code, used by the tests) and `message`.

rule_df <- function(rule = character(0), message = character(0)) data.frame(rule = rule, message = message, stringsAsFactors = FALSE)
plural <- function(k, word, words = paste0(word, "s")) paste(k, if (k == 1) word else words)
pct <- function(x) formatC(100 * x, format = "f", digits = 1)

## Rules applied by bohr_fit(): missing responses, skips, the boundary rule, the family, collinear columns and the
## beta part. fit: a bohr_fit object (with $input).
fit_rules <- function(fit) {
  ct <- fit$control; inp <- fit$input; R <- rule_df()
  if (is.null(inp)) return(R)
  add <- function(rule, msg) R <<- rbind(R, rule_df(rule, msg))
  if (inp$n_missing > 0)
    add("missing", sprintf("%s with a missing response (or a missing value of the tested term) left out.",
                           plural(inp$n_missing, "observation")))
  nn <- inp$n_input - inp$n_missing
  if (!is.na(fit$skip)) { add("skip", skip_message(fit$skip, nn, inp, ct)); return(R) }
  for (b in c("0", "1")) {
    nb <- if (b == "0") inp$n_at0 else inp$n_at1; kept <- if (b == "0") fit$n0 else fit$n1
    if (nb > 0 && kept == 0)
      add(paste0("boundary_", b), sprintf(paste0("%s at %s (%s%% of %d) removed by the boundary rule (a boundary holding %s%% of the ",
                                                 "observations or fewer is removed); p%s is not modelled."),
                                          plural(nb, "observation"), b, pct(nb / nn), nn, format(100 * ct$endpoint), b))
  }
  if (identical(fit$family, "BE")) add("mu_only", "Fitted mu only: no observations at 0 or 1 remain (family BE).")
  if (isTRUE(fit$n_cov_dropped > 0)) {
    d <- attr(fit$design, "dropped")
    add("collinear", sprintf("%s dropped as collinear with the others (tolerance 1e-7): %s.",
                             plural(length(d), "covariate column"), paste(d, collapse = ", ")))
  }
  if (is.null(fit$null$beta))
    add("mu_not_fitted", sprintf(paste0("mu not fitted: %s between 0 and 1, but more than %d (design columns + 2) are needed; ",
                                        "mu is not tested."), plural(fit$n_mid, "observation"), ncol(fit$design) + 2L))
  R
}

## the reason of a skip, from the counts before the rules (prepare_y(): n first, then boundary share or SD)
skip_message <- function(skip, nn, inp, ct) {
  k <- length(inp$covariates)
  if (skip == "n") {
    if (nn < ct$min_samples) return(sprintf("Skipped: %s with a value, fewer than the %d required.", plural(nn, "observation"), ct$min_samples))
    return(sprintf("Skipped: %s, fewer than the %d covariate columns.", plural(nn, "observation"), k))
  }
  if (skip == "prop_or_sd") {
    why <- c(if (nn > 0 && inp$n_at0 / nn >= ct$max_boundary) sprintf("%s%% of the observations are at 0", pct(inp$n_at0 / nn)),
             if (nn > 0 && inp$n_at1 / nn >= ct$max_boundary) sprintf("%s%% of the observations are at 1", pct(inp$n_at1 / nn)),
             if (is.finite(inp$sd) && inp$sd <= ct$min_sd) sprintf("the standard deviation is %s", format(signif(inp$sd, 3))))
    return(sprintf("Skipped: %s (rule: at least %s%% at one boundary, or a standard deviation of %s or less).",
                   paste(why, collapse = "; "), format(100 * ct$max_boundary), format(ct$min_sd)))
  }
  if (skip == "no_interior") return("Skipped: no observations strictly between 0 and 1.")
  paste0("Skipped: ", skip, ".")
}

## Conditions of the tests of one predictor (a row of bohr_test()): failures of the fits, a constant predictor,
## separation and convergence. Rules (they change what is tested) and notes are told apart by `kind`.
test_rules <- function(row, family, firth = TRUE) {
  R <- rule_df()
  add <- function(rule, msg) R <<- rbind(R, rule_df(rule, msg))
  err <- row$error
  if (!is.na(err) && !startsWith(err, "skipped")) {
    if (grepl("constant", err)) add("x_constant", "The tested term is constant on the observations used; it is not tested.")
    else add("fit_error", sprintf("The fits with the tested term failed (%s); no statistics are reported.", err))
  }
  if (isTRUE(row$firth_error)) add("firth_error", "The Firth-penalised fits failed; the p0 and p1 P values are not available.")
  if (!is.na(err) || is.na(family)) return(R)
  has0 <- family %in% c("BEINF", "BEZI"); has1 <- family %in% c("BEINF", "BEOI")
  if (isTRUE(row$complete_sep))
    add("complete_separation", paste0("The covariates alone separate the boundary categories completely: the unpenalised p0/p1 tests are ",
                                      "not defined", if (firth) "; the Firth-penalised tests are reported." else " (set firth = TRUE)."))
  else for (cmp in c("p0", "p1")) if ((cmp == "p0" && has0) || (cmp == "p1" && has1)) {
    if (isTRUE(row[[paste0("sep_", cmp)]]))
      add(paste0("separation_", cmp), sprintf(paste0("%s: the unpenalised fit shows separation (no convergence, |estimate| > 15 or SE > 50 ",
                                                    "on the coded scale); %s"), cmp,
                                             if (firth) "the reported test is Firth-penalised." else "its unpenalised test (firth = FALSE) is unreliable."))
  }
  if (isTRUE(row$sep_mu) && isTRUE(row$conv_mu))
    add("separation_mu", "mu: extreme estimate in the beta fit (|estimate| > 15 or SE > 50 on the coded scale).")
  if (isFALSE(row$conv_mu)) add("not_converged_mu", "mu: the beta fit with the tested term did not converge.")
  for (cmp in c("p0", "p1")) if (((cmp == "p0" && has0) || (cmp == "p1" && has1)) && isFALSE(row[[paste0("conv_", cmp)]]) && !isTRUE(row$complete_sep))
    add(paste0("not_converged_", cmp), sprintf("%s: the %s fits did not converge.", cmp, if (firth) "penalised" else "unpenalised"))
  R
}

print_rules <- function(R, notes = NULL, indent = "  ") {
  if (!is.null(R) && nrow(R)) { cat("Fitting rules applied:\n"); cat(paste0(indent, "- ", R$message), sep = "\n") }
  if (!is.null(notes) && nrow(notes)) { cat("Notes:\n"); cat(paste0(indent, "- ", notes$message), sep = "\n") }
  invisible(NULL)
}
