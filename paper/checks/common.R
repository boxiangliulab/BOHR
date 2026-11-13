## Shared functions of the checks of run_paper.sh (authored for the repository; not part of the original runs).
## Each check script records rows (check, result, detail) with add() and writes them with finish(<status.tsv>).
## result is PASS, FAIL or SKIP. Details hold aggregate numbers only (counts, largest differences), never sample or
## donor identifiers.
suppressPackageStartupMessages(library(data.table))
.res <- list()
add <- function(check, ok, detail) {
  r <- if (identical(ok, "SKIP")) "SKIP" else if (isTRUE(ok)) "PASS" else "FAIL"
  .res[[length(.res) + 1L]] <<- data.table(check = check, result = r, detail = detail)
  cat(sprintf("%-4s  %s  %s\n", r, check, detail))
}
finish <- function(file) {
  R <- if (length(.res)) rbindlist(.res) else data.table(check = character(), result = character(), detail = character())
  fwrite(R, file, sep = "\t", append = file.exists(file))   # several check scripts may write to one sub-step's file
  invisible(R)
}
## read .gz without R.utils
fr <- function(f, ...) if (grepl("\\.gz$", f)) fread(cmd = paste("gzip -dc", shQuote(f)), ...) else fread(f, ...)
## the published Supplementary Tables (one TSV per sheet)
supp <- function(dir, sheet) fread(list.files(dir, pattern = paste0("^Supplementary_Table_", sheet, "_.*[.]tsv$"), full.names = TRUE)[1],
                                   sep = "\t", encoding = "UTF-8")
## table out of an rds: a data.frame, or the $res / $snp / ... element named by `el`
as_tab <- function(x, el = NULL) {
  if (!is.null(el)) x <- x[[el]]
  if (is.list(x) && !is.data.frame(x) && !is.null(x$res)) x <- x$res
  as.data.table(x)
}
## Compare a re-run table with the stored one.
##   keys: columns that identify a row; skip: columns not compared (run times, dates, paths)
##   numeric columns: largest |difference| must be <= tol (absolute) or <= tol relative to the stored value;
##   other columns: identical; missing values must sit in the same cells.
##   Unpenalised p0 and p1 fits under complete separation: the estimate diverges and the values where the optimiser stops
##   (estimate, SE, statistic, P, probability-scale effect, even the separation and convergence flags at the margin) are
##   arbitrary and differ between platforms. Where such columns exist (nu_* / tau_* with nu_sep / tau_sep, p0_pe / p1_pe),
##   rows flagged as separated in either run (or cat_complete_sep) are not compared for those columns, and are counted.
##   The Firth-penalised fits, which the paper uses for p0 and p1, are compared in every row.
## Returns list(ok, detail).
compare_tables <- function(new, old, keys, tol = 1e-6, skip = character()) {
  new <- as.data.table(new); old <- as.data.table(old)
  skip <- union(skip, grep("(^|_)(minutes|seconds|secs|sec|time|date|elapsed|host|runtime)$", union(names(new), names(old)), value = TRUE))
  if (nrow(new) == 0 && nrow(old) == 0) return(list(ok = TRUE, detail = "both empty"))
  miss_k <- setdiff(keys, intersect(names(new), names(old)))
  if (length(miss_k)) return(list(ok = FALSE, detail = paste("key columns missing:", paste(miss_k, collapse = ", "))))
  kn <- do.call(paste, c(new[, ..keys], sep = "\r")); ko <- do.call(paste, c(old[, ..keys], sep = "\r"))
  if (anyDuplicated(kn) || anyDuplicated(ko)) return(list(ok = FALSE, detail = "duplicated keys"))
  extra <- sum(!kn %in% ko); missing <- sum(!ko %in% kn)
  new <- new[match(intersect(ko, kn), kn)]; old <- old[match(intersect(ko, kn), ko)]
  cols <- setdiff(intersect(names(new), names(old)), c(keys, skip))
  only <- setdiff(union(names(new), names(old)), c(cols, keys, skip))
  bad <- character(); dmax <- 0; n_divg <- 0L
  sep_rows <- function(p) { s <- paste0(p, "_sep"); if (!(s %in% names(new) && s %in% names(old))) return(NULL)
    x <- new[[s]] %in% TRUE | old[[s]] %in% TRUE
    if ("cat_complete_sep" %in% names(new) && "cat_complete_sep" %in% names(old)) x <- x | new$cat_complete_sep %in% TRUE | old$cat_complete_sep %in% TRUE
    x }
  SEP <- list(nu = sep_rows("nu"), tau = sep_rows("tau"))
  n_sep <- c(nu = sum(SEP$nu), tau = sum(SEP$tau)); sep_flag_diff <- 0L
  for (p in names(SEP)) if (!is.null(SEP[[p]])) sep_flag_diff <- sep_flag_diff + sum(new[[paste0(p, "_sep")]] != old[[paste0(p, "_sep")]], na.rm = TRUE)
  unpen <- function(cl) {
    if (cl == "cat_complete_sep" && !is.null(SEP$nu) && !is.null(SEP$tau)) return(SEP$nu | SEP$tau)
    p <- if (cl == "p0_pe") "nu" else if (cl == "p1_pe") "tau" else sub("_(stat|lrt_logp|wald_logp|est|se|conv|sep)$", "", cl)
    if (p %in% names(SEP) && !is.null(SEP[[p]]) && (cl %in% c("p0_pe", "p1_pe") || grepl("^(nu|tau)_(stat|lrt_logp|wald_logp|est|se|conv|sep)$", cl))) SEP[[p]] else NULL }
  any_sep_rule <- FALSE
  for (cl in cols) {
    a <- new[[cl]]; b <- old[[cl]]
    ex <- unpen(cl)
    if (!is.null(ex)) { any_sep_rule <- TRUE; a <- a[!ex]; b <- b[!ex] }
    if (is.factor(a)) a <- as.character(a); if (is.factor(b)) b <- as.character(b)
    if (is.list(a) || is.list(b)) { if (!identical(a, b)) bad <- c(bad, cl); next }
    if (!identical(is.na(a), is.na(b))) { bad <- c(bad, paste0(cl, " (missing values)")); next }
    if (is.numeric(a) && is.numeric(b)) {
      keep <- rep(TRUE, length(a)); sp <- paste0(sub("_est$", "", cl), "_sep")
      if (is.null(ex) && grepl("_est$", cl) && sp %in% names(new) && sp %in% names(old)) {
        dv <- new[[sp]] %in% TRUE & old[[sp]] %in% TRUE; keep <- !dv; n_divg <- n_divg + sum(dv) }
      i <- keep & !is.na(a) & is.finite(a) & is.finite(b)
      if (!identical(is.finite(a[keep & !is.na(a)]), is.finite(b[keep & !is.na(b)]))) { bad <- c(bad, paste0(cl, " (infinite values)")); next }
      if (any(i)) {
        d <- abs(a[i] - b[i]); r <- d / pmax(abs(b[i]), 1e-300)
        e <- pmin(d, r); dmax <- max(dmax, max(e))
        ## Firth estimates and SEs: the penalised likelihood is flat near its optimum, so they agree to about 1e-6 while
        ## their P values agree to 1e-9 (STAGE1_2_REPORT); they are held to 1e-5
        tl <- if (grepl("_firth_(est|se)$", cl)) max(tol, 1e-5) else tol
        if (max(e) > tl) bad <- c(bad, sprintf("%s (max |diff| %.3g)", cl, max(d)))
      }
    } else if (!identical(as.character(a), as.character(b))) bad <- c(bad, paste0(cl, " (", sum(as.character(a) != as.character(b), na.rm = TRUE), " cells differ)"))
  }
  ok <- extra == 0 && missing == 0 && !length(bad)
  detail <- sprintf("%d rows (stored %d; %d extra, %d missing); %d columns compared, largest difference %.3g%s%s%s%s",
                    length(kn), length(ko), extra, missing, length(cols), dmax,
                    if (n_divg) sprintf("; %d unpenalised estimates under complete separation not compared", n_divg) else "",
                    if (any_sep_rule) sprintf("; unpenalised p0/p1 columns not compared in %d (p0) and %d (p1) rows separated in either run (flag differs in %d)",
                                              n_sep[["nu"]], n_sep[["tau"]], sep_flag_diff) else "",
                    if (length(bad)) paste0("; differ: ", paste(head(bad, 6), collapse = ", ")) else "",
                    if (length(only)) paste0("; in one table only: ", paste(head(only, 6), collapse = ", ")) else "")
  list(ok = ok, detail = detail)
}
## compare every table file (tsv, tsv.gz, csv) of a re-run folder with the file of the same relative path in the stored folder
compare_folders <- function(newdir, olddir, label, tol = 1e-6, pattern = "\\.(tsv|tsv\\.gz|csv)$") {
  fs <- list.files(newdir, pattern = pattern, recursive = TRUE)
  if (!length(fs)) { add(label, FALSE, paste("no table written under", basename(newdir))); return(invisible()) }
  for (f in fs) {
    o <- file.path(olddir, f)
    if (!file.exists(o)) { add(paste(label, f), "SKIP", "no stored file of that name to compare with"); next }
    a <- fr(file.path(newdir, f)); b <- fr(o)
    if (!identical(names(a), names(b)) || nrow(a) != nrow(b)) {
      add(paste(label, f), FALSE, sprintf("%d x %d, stored %d x %d; same columns: %s", nrow(a), ncol(a), nrow(b), ncol(b), identical(names(a), names(b))))
      next }
    a[, .row := .I]; b[, .row := .I]
    r <- compare_tables(a, b, ".row", tol)
    add(paste(label, f), r$ok, r$detail)
  }
}
