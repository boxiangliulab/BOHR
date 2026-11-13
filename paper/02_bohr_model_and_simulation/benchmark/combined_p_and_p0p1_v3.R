#!/usr/bin/env Rscript
### Author 2026-10-07: (1) null calibration of the BOHR p0 / p1 (Firth) tests: per-promoter FPR at 0.05
### (seed 1, the figure's permutations) and QQ points; (2) null tail of a promoter-level
### combined p over the promoter's existing tests (mu-F, p0, p1): ACAT (Cauchy combination) and Simes (the rule used for the
### real-data Fig 2e), against the mu test alone; five permutation seeds pooled, 95% CIs from 2,000 unit-level bootstraps.
### DATASET env (default nothr_realjunc). Outputs: data_v3_<DATASET>/diag_tail/{combined_p_*,p0p1_*}.csv
suppressMessages(library(data.table))
DS <- Sys.getenv("DATASET", "nothr_realjunc"); pd <- readLines("/projects/group/bliu/e1327989/GTEx/new_simulation/figures_bohr2/scripts/prepare_data_v3.R")
OUTD <- sprintf("/projects/group/bliu/e1327989/GTEx/new_simulation/figures_bohr2/data_v3_%s/diag_tail", DS); dir.create(OUTD, showWarnings = FALSE)
TAGS <- c("seed 1 (figure)" = "", "seed 2" = "_seed2", "seed 3" = "_seed3", "seed 4" = "_seed4", "seed 5" = "_seed5")
acat <- function(P) {                                   # rows = tests, columns = mu, p0, p1 (NA = test absent)
  P <- pmin(pmax(P, 1e-300), 1 - 1e-15); T <- tan((0.5 - P) * pi); k <- rowSums(!is.na(P))
  0.5 - atan(rowSums(T, na.rm = TRUE) / k) / pi }
simes <- function(P) {                                 # vectorised: long format, rank within row, min(k p_(i) / i)
  L <- data.table(r = rep(seq_len(nrow(P)), ncol(P)), p = as.vector(P))[!is.na(p)]; setorder(L, r, p)
  L[, `:=`(i = seq_len(.N), k = .N), by = r]; S <- L[, .(s = min(k * p / i)), by = r]
  out <- rep(NA_real_, nrow(P)); out[S$r] <- pmin(S$s, 1); out }
qq_of <- function(p) { p <- sort(p[!is.na(p)]); n <- length(p)
  d <- data.table(rank = seq_len(n), n = n, observed = -log10(pmax(p, 1e-300)), expected = -log10(seq_len(n) / (n + 1)))
  d[, `:=`(band_lo = -log10(qbeta(.975, rank, n - rank + 1)), band_hi = -log10(qbeta(.025, rank, n - rank + 1)), lambda = median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(.5, 1))]
  d[(rank <= 2000) | (rank %% pmax(1L, as.integer(n / 3000)) == 0L)] }
H <- list(); Q <- list(); P01 <- NULL
for (s in names(TAGS)) {
  Sys.setenv(SEEDTAG = TAGS[[s]], DATASET = DS); e <- new.env(); eval(parse(text = pd[1:(grep("^MODES <- ", pd) - 1)]), envir = e)
  X <- e$NUL[, .(m, N, id, perm, BOHR, BOHR_p0, BOHR_p1)]; rm(e); gc()
  M <- as.matrix(X[, .(BOHR, BOHR_p0, BOHR_p1)]); ok <- rowSums(!is.na(M)) > 0
  X[, ACAT := NA_real_][ok, ACAT := acat(M[ok, , drop = FALSE])]
  X[, SIMES := NA_real_][ok, SIMES := simes(M[ok, , drop = FALSE])]
  H[[s]] <- X[, .(K = .N, mu3 = sum(BOHR < 1e-3, na.rm = TRUE), mu4 = sum(BOHR < 1e-4, na.rm = TRUE), mu5 = sum(BOHR < 1e-5, na.rm = TRUE),
                  ac3 = sum(ACAT < 1e-3, na.rm = TRUE), ac4 = sum(ACAT < 1e-4, na.rm = TRUE), ac5 = sum(ACAT < 1e-5, na.rm = TRUE),
                  si3 = sum(SIMES < 1e-3, na.rm = TRUE), si4 = sum(SIMES < 1e-4, na.rm = TRUE), si5 = sum(SIMES < 1e-5, na.rm = TRUE),
                  K_mu = sum(!is.na(BOHR)), K_ac = sum(!is.na(ACAT))), by = .(m, N, id)][, seed := s]
  for (mm in c("M1", "M2")) for (st in c("BOHR", "ACAT", "SIMES")) Q[[length(Q) + 1]] <- qq_of(X[m == mm][[st]])[, `:=`(m = mm, seed = s, stat = st)]
  if (s == "seed 1 (figure)") {                         # p0 / p1 calibration on the figure's permutations
    P01 <- rbind(X[!is.na(BOHR_p0), .(test = "p0 (Firth)", fpr = mean(BOHR_p0 < .05), K = .N), by = .(m, N, id)],
                 X[!is.na(BOHR_p1), .(test = "p1 (Firth)", fpr = mean(BOHR_p1 < .05), K = .N), by = .(m, N, id)])
    Q01 <- rbindlist(lapply(c("M1", "M2"), function(mm) rbind(qq_of(X[m == mm, BOHR_p0])[, `:=`(m = mm, test = "p0 (Firth)")], qq_of(X[m == mm, BOHR_p1])[, `:=`(m = mm, test = "p1 (Firth)")])))
    fwrite(P01, file.path(OUTD, "p0p1_fpr_perpromoter_seed1.csv")); fwrite(Q01, file.path(OUTD, "p0p1_qq_seed1.csv"))
    cat("== p0 / p1 (Firth), seed 1: pooled FPR at 0.05 and tail ratios ==\n")
    print(rbind(X[!is.na(BOHR_p0), .(test = "p0", units = uniqueN(paste(N, id)), fpr05 = round(100 * mean(BOHR_p0 < .05), 2), r1e3 = round(mean(BOHR_p0 < 1e-3) / 1e-3, 2), r1e4 = round(mean(BOHR_p0 < 1e-4) / 1e-4, 2)), by = m],
                X[!is.na(BOHR_p1), .(test = "p1", units = uniqueN(paste(N, id)), fpr05 = round(100 * mean(BOHR_p1 < .05), 2), r1e3 = round(mean(BOHR_p1 < 1e-3) / 1e-3, 2), r1e4 = round(mean(BOHR_p1 < 1e-4) / 1e-4, 2)), by = m])) }
  rm(X, M); gc() }
H <- rbindlist(H); Q <- rbindlist(Q)
U <- H[, lapply(.SD, sum), by = .(m, N, id), .SDcols = c("K", "K_mu", "K_ac", "mu3", "mu4", "mu5", "ac3", "ac4", "ac5", "si3", "si4", "si5")]
set.seed(20261007)
ratios <- function(D) { f <- function(d) d[, c(sum(mu3) / (sum(K_mu) * 1e-3), sum(mu4) / (sum(K_mu) * 1e-4), sum(mu5) / (sum(K_mu) * 1e-5),
                                               sum(ac3) / (sum(K_ac) * 1e-3), sum(ac4) / (sum(K_ac) * 1e-4), sum(ac5) / (sum(K_ac) * 1e-5),
                                               sum(si3) / (sum(K_ac) * 1e-3), sum(si4) / (sum(K_ac) * 1e-4), sum(si5) / (sum(K_ac) * 1e-5))]
  est <- f(D); bt <- replicate(2000, f(D[sample(.N, replace = TRUE)]))
  data.table(p = rep(c("mu test alone", "ACAT (mu, p0, p1)", "Simes (mu, p0, p1)"), each = 3), threshold = rep(c("1e-3", "1e-4", "1e-5"), 3),
             ratio = round(est, 2), lo = round(apply(bt, 1, quantile, .025), 2), hi = round(apply(bt, 1, quantile, .975), 2)) }
op <- options(width = 200)
R <- rbindlist(lapply(c("M2", "M1"), function(mm) ratios(U[m == mm])[, m := mm]))
cat("\n== five seeds pooled: observed / expected [95% unit-bootstrap CI] ==\n"); print(dcast(R[, v := sprintf("%.2f [%.2f-%.2f]", ratio, lo, hi)], m + p ~ threshold, value.var = "v"))
cat("\n== lambda (pooled QQ, per seed) ==\n"); print(dcast(unique(Q[, .(m, seed, stat, lambda = round(lambda, 3))]), m + stat ~ seed, value.var = "lambda"))
options(op)
fwrite(R, file.path(OUTD, "combined_p_tail_ci.csv")); fwrite(Q, file.path(OUTD, "combined_p_qq_points.csv")); fwrite(U, file.path(OUTD, "combined_p_units.csv"))
cat("COMBINED_DONE\n")
