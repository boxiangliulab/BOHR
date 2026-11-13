## =============================================================================
## 01_abs_labels.R -- the boundary label of exact mirror pairs is decided by absolute activity (author 2026-11-06)
##   One-off step: run the real trait only and save the label table;
##   later relabelling only reads this table, routine runs no longer fit.
## Events: every testable boundary unit p0(S) on boundary_exact pairs in trait_regression_20261104 (mirror rule v3)
##   (S = the TSS with activity 0 in the event; O = its mirror partner; p0(S) and p1(O) are the same test).
## Fit: absolute activity A of each TSS (proActiv log2 normalized counts, v11/tabs_att_filter_In.rds $abs; 0 reads = 0)
##   OLS on the trait contrast coding (trait_contrast), covariates = those of the trait model (cov_design: intercept + standardized SEX/AGE/PC1-5/
##   MHABNWBC + DTHHRDY dummies; the trait being tested is not a covariate), donors = freeze donors of the tissue with complete trait and covariates
##   (same as run_trait.R; A is defined for all samples). Sensitivity: only donors with defined relative activity of S.
##   regulated = |z| >= 1.96.
## Label: only S changes -> p0(S); only O changes -> p1(O); both or neither change -> p0(S) (when undecidable, record as p0)
## Output: results/mirror_abs_labels.tsv.gz (one row per event: z, effect, outcome, label), results/01_abs_labels.log
## =============================================================================
suppressPackageStartupMessages({ library(data.table); library(Matrix) })
source("/projects/group/bliu/e1327989/GTEx/tissue_regression/trait_regression_20261104/code/data.R")   # cov_design, trait_contrast, COVS, CONT_TRAITS
OUT <- "/projects/group/bliu/e1327989/GTEx/tissue_regression/trait_regression_20261106_mirror_abs/results"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ZC <- 1.96
t0 <- Sys.time()

## ---- Events: testable p0 units on boundary_exact in the v3 unit table (all primary under v3) ----
U <- readRDS(file.path(PATHS$res, "summary/units_real.rds"))
E <- U[type == "p0" & boundary_exact %in% TRUE & testable %in% TRUE,
       .(tissue, trait, S = promoter, O = mirror_partner, gene, mirror_pair_id, mirror_role_S = mirror_role,
         n, n0, n1, p_v3 = p, est_v3 = est, primary, final, q_final_v3 = q_final, sig_final_v3 = sig_final)]
TW <- U[type == "p1" & boundary_exact %in% TRUE, .(tissue, trait, O = promoter, twin_testable = testable, p_twin = p, est_twin = est)]
rm(U); invisible(gc())
E <- merge(E, TW, by = c("tissue", "trait", "O"), all.x = TRUE)
stopifnot(!anyNA(E$O), all(E$primary))
cat("events (testable p0 on boundary_exact pairs):", nrow(E), "| final", sum(E$final), "| significant final", sum(E$sig_final_v3), "\n")

## ---- Absolute activity ----
ab <- readRDS(file.path(GTEX, "v11/tabs_att_filter_In.rds"))$abs
need <- unique(c(E$S, E$O)); stopifnot(all(need %in% colnames(ab)))
ab <- ab[, need, drop = FALSE]; invisible(gc())
meta <- readRDS(file.path(PATHS$work, "meta.rds"))

ols_z <- function(X, Y) {                                   # Multi-response OLS: b, se, z of the last column (trait)
  q <- qr(X); p <- ncol(X); stopifnot(q$rank == p)
  B <- qr.coef(q, Y); R <- qr.resid(q, Y)
  s2 <- colSums(R^2) / (nrow(X) - p)
  v <- chol2inv(qr.R(q))[p, p]
  b <- B[p, ]; se <- sqrt(v * s2)
  list(b = b, se = se, z = b / se, sdy = apply(Y, 2, sd))
}

res <- list()
for (t in sort(unique(E$tissue))) {
  td <- readRDS(file.path(PATHS$work, paste0(t, ".rds")))
  stopifnot(all(td$row_id %in% rownames(ab)))
  e_t <- E[tissue == t]
  A <- as.matrix(ab[td$row_id, unique(c(e_t$S, e_t$O)), drop = FALSE]); rownames(A) <- td$U
  for (tr in sort(unique(e_t$trait))) {
    e <- e_t[trait == tr]
    cv <- setdiff(COVS, tr)
    donors <- td$U[complete.cases(meta[td$U, c(tr, cv), drop = FALSE])]
    X <- cbind(cov_design(meta, donors, tr), g = trait_contrast(meta[donors, tr], tr))
    ## Main analysis: all donors
    fS <- ols_z(X, A[donors, e$S, drop = FALSE]); fO <- ols_z(X, A[donors, e$O, drop = FALSE])
    ## Zero-value check: among donors with defined relative activity of S, A_S == 0 exactly where relative activity of S == 0, and then A_O > 0
    rel <- td$act[donors, e$S, drop = FALSE]
    zero_ok <- vapply(seq_len(nrow(e)), function(i) { r <- rel[, i]; k <- !is.na(r)
      all((A[donors[k], e$S[i]] == 0) == (r[k] == 0)) && all(A[donors[k][r[k] == 0], e$O[i]] > 0) }, TRUE)
    ## Sensitivity: only donors with defined relative activity of S (each event's own donor set)
    sens <- t(vapply(seq_len(nrow(e)), function(i) { k <- donors[!is.na(rel[, i])]
      Xk <- cbind(cov_design(meta, k, tr), g = trait_contrast(meta[k, tr], tr))
      if (qr(Xk)$rank < ncol(Xk)) return(c(NA_real_, NA_real_, length(k)))
      c(ols_z(Xk, A[k, e$S[i], drop = FALSE])$z, ols_z(Xk, A[k, e$O[i], drop = FALSE])$z, length(k)) }, numeric(3)))
    res[[paste(t, tr)]] <- e[, `:=`(n_abs = length(donors), b_S = fS$b, se_S = fS$se, z_S = fS$z, sd_S = fS$sdy,
                                    b_O = fO$b, se_O = fO$se, z_O = fO$z, sd_O = fO$sdy, zero_pattern_ok = zero_ok,
                                    z_S_reldonors = sens[, 1], z_O_reldonors = sens[, 2], n_reldonors = sens[, 3])]
  }
  cat(t, nrow(e_t), "events\n")
}
L <- rbindlist(res)
outc <- function(zs, zo) { rs <- abs(zs) >= ZC; ro <- abs(zo) >= ZC
  fifelse(rs & ro, "both", fifelse(rs, "S_only", fifelse(ro, "O_only", "neither"))) }
L[, outcome := outc(z_S, z_O)]
L[, label_component := fifelse(outcome == "O_only", "p1", "p0")]
L[, label_promoter := fifelse(outcome == "O_only", O, S)]
L[, label := paste0(label_promoter, ":", label_component)]
L[, v3_label := paste0(S, ":p0")]
L[, outcome_reldonors := outc(z_S_reldonors, z_O_reldonors)]
L[, label_reldonors := fifelse(outcome_reldonors == "O_only", paste0(O, ":p1"), paste0(S, ":p0"))]
L[, rule := "abs_activity_2026-11-06 (|z|>=1.96; S only -> p0(S); O only -> p1(O); both/neither -> p0(S))"]
setcolorder(L, c("tissue", "trait", "gene", "S", "O", "outcome", "label", "v3_label", "z_S", "z_O"))
fwrite(L, file.path(OUT, "mirror_abs_labels.tsv.gz"), sep = "\t")

sink(file.path(OUT, "01_abs_labels.log"), split = TRUE)
cat("events", nrow(L), "| tissues", uniqueN(L$tissue), "| traits", uniqueN(L$trait), "| non-finite z", L[!is.finite(z_S) | !is.finite(z_O), .N], "\n")
cat("zero pattern (A_S == 0 exactly where relative S == 0, A_O > 0 there):", L[zero_pattern_ok == TRUE, .N], "of", nrow(L), "\n")
cat("twin p1(O) testable:", L[twin_testable %in% TRUE, .N], "; |p_twin - p_v3| > 1e-6:", L[abs(p_twin - p_v3) > 1e-6, .N], "\n")
for (s in c("all", "final", "sig_final")) {
  x <- switch(s, all = L, final = L[final == TRUE], sig_final = L[sig_final_v3 == TRUE])
  cat("\n", s, ": events", nrow(x), "\n"); print(x[, .N, by = outcome][order(outcome)])
  cat("  labels differing from v3:", x[label != v3_label, .N], "; under relative-usage donors:", x[label_reldonors != v3_label, .N], "\n")
}
cat("\nruntime", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")
sink()
