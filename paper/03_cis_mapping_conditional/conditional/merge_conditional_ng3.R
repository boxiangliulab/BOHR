## merge_conditional_ng3.R — final conditional analysis for units v3 (results/final_nogate_v3: mu n_mid >= 80 AND p0/p1
## minor class >= 10%, both before BH; author 2026-10-02). Sources in priority order:
##   results/conditional_nogate_supp3 (v3 supplement: new and continued units)  > results/conditional_nogate_supp (v2)
##   > results/conditional_nogate (v1 main). v2 units outside v3 are dropped. v3 thresholds are looser or equal (no
##   truncation). Adapted from merge_conditional_ng.R (v2).
##   unit in v2 and in the supplement (new or continued)  -> supplement result (continued units carry their v1 prefix)
##   unit in v2, not in the supplement                    -> v1 result (identical under the looser v2 threshold)
##   v1 unit not in v2 (mu with n_mid < 80)               -> dropped
## Backward (single pass, GTEx/tensorQTL): keep a signal of a multi-signal unit if its leave-one-out window ACAT <= thr_v3.
## Output: results/conditional_final_v3/ (env COND_FINAL): freeze_real_{units,signals,rounds}.tsv.gz, freeze_real_arms.tsv, loo/<tissue>.rds
suppressPackageStartupMessages(library(data.table))
V2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"; DEF <- "freeze"
OUT <- Sys.getenv("COND_FINAL", file.path(V2, "results", "conditional_final_v3")); dir.create(file.path(OUT, "loo"), recursive = TRUE, showWarnings = FALSE)
acat_log <- function(logp) { p <- exp(logp[is.finite(logp)]); if (!length(p)) return(NA_real_); if (any(p == 0)) return(0)
  if (length(p) > 1) p <- pmin(p, 1 - 1 / length(p)); sm <- p < 1e-15; t <- numeric(length(p)); t[!sm] <- tan((0.5 - p[!sm]) * pi); t[sm] <- 1 / (p[sm] * pi)
  st <- mean(t); if (st > 1e15) return(1 / (st * pi)); pcauchy(st, lower.tail = FALSE) }
rd <- function(root, t, ch) { f <- file.path(V2, "results", root, DEF, t, "real", paste0(ch, ".rds")); if (file.exists(f)) readRDS(f) else NULL }
tab <- function(x, k) if (is.null(x) || !NROW(x[[k]])) data.table() else as.data.table(x[[k]])
U_all <- list(); S_all <- list(); R_all <- list(); A_all <- list()
for (t in sort(basename(Sys.glob(file.path(V2, "results/conditional_nogate", DEF, "*"))))) {
  V <- as.data.table(readRDS(file.path(V2, "results/final_nogate_v3", DEF, paste0(t, ".real.units.rds"))))
  thr <- V[q <= 0.05, max(p)]; keep2 <- V[q <= 0.05, paste(promoter, type)]
  u <- list(); s <- list(); r <- list(); l <- list(); miss <- 0L
  for (ch in paste0("chr", 1:22)) { x1 <- rd("conditional_nogate", t, ch); x2 <- rd("conditional_nogate_supp", t, ch); x3 <- rd("conditional_nogate_supp3", t, ch)
    if (is.null(x1)) { miss <- miss + 1L; next }
    cl <- function(x) { d <- tab(x, "units"); if (nrow(d)) d[!is.na(type)] else d }
    u1 <- cl(x1); u2 <- cl(x2); u3 <- cl(x3)
    k3 <- if (nrow(u3)) paste(u3$promoter, u3$type) else character(0)
    k2 <- if (nrow(u2)) setdiff(intersect(paste(u2$promoter, u2$type), keep2), k3) else character(0)
    k1 <- if (nrow(u1)) setdiff(intersect(paste(u1$promoter, u1$type), keep2), c(k2, k3)) else character(0)
    sub <- function(d, keys) if (nrow(d)) d[paste(promoter, type) %in% keys] else d
    lab <- function(d, src) if (nrow(d)) d[, source := src] else NULL
    uu <- rbindlist(list(lab(sub(u1, k1), "main"), lab(sub(u2, k2), "supp"), lab(sub(u3, k3), "supp3")), fill = TRUE)
    if (nrow(uu)) uu[, chr := ch]; u[[ch]] <- uu
    ss <- rbindlist(list(sub(tab(x1, "signals"), k1), sub(tab(x2, "signals"), k2), sub(tab(x3, "signals"), k3)), fill = TRUE); if (nrow(ss)) ss[, chr := ch]; s[[ch]] <- ss
    rr <- rbindlist(list(lab(sub(tab(x1, "rounds"), k1), "main"), lab(sub(tab(x2, "rounds"), k2), "supp"), lab(sub(tab(x3, "rounds"), k3), "supp3")), fill = TRUE)
    r[[ch]] <- rr
    l[[ch]] <- rbindlist(list(sub(tab(x1, "loo"), k1), sub(tab(x2, "loo"), k2), sub(tab(x3, "loo"), k3)), fill = TRUE) }
  u <- rbindlist(u, fill = TRUE); s <- rbindlist(s, fill = TRUE); r <- rbindlist(r, fill = TRUE); l <- rbindlist(l, fill = TRUE)
  if (nrow(u)) u[, tissue := t]; if (nrow(s)) s[, tissue := t]; if (nrow(r)) r[, tissue := t]
  ## backward from LOO (single pass)
  if (nrow(l)) { bw <- l[, .(loo_acat = acat_log(logp), loo_top_snp = snp[which.min(replace(logp, !is.finite(logp), Inf))]), by = .(promoter, type, k)]
    s <- merge(s, bw, by = c("promoter", "type", "k"), all.x = TRUE) } else s[, `:=`(loo_acat = NA_real_, loo_top_snp = NA_character_)]
  s[, n_sig := .N, by = .(promoter, type)]
  s[, backward_keep := fifelse(n_sig == 1L, TRUE, fifelse(is.na(loo_acat), NA, loo_acat <= thr))]; s[, n_sig := NULL]
  u <- merge(u, s[, .(n_signals_backward = if (anyNA(backward_keep)) NA_integer_ else sum(backward_keep)), by = .(promoter, type)], by = c("promoter", "type"), all.x = TRUE)
  saveRDS(l, file.path(OUT, "loo", paste0(t, ".rds")))
  A_all[[t]] <- data.table(tissue = t, thr_v3 = thr, units_v3 = length(keep2), units = nrow(u), missing_chr = miss, from_supp = sum(u$source == "supp"), from_supp3 = sum(u$source == "supp3"),
    signals = nrow(s), units_ge2 = sum(u$n_signals >= 2), units_ge2_backward = sum(u$n_signals_backward >= 2, na.rm = TRUE), backward_na = sum(is.na(s$backward_keep)))
  U_all[[t]] <- u; S_all[[t]] <- s; R_all[[t]] <- r }
U <- rbindlist(U_all, fill = TRUE); S <- rbindlist(S_all, fill = TRUE); R <- rbindlist(R_all, fill = TRUE); A <- rbindlist(A_all, fill = TRUE)
## p0/p1 loo_est lost to the est-name bug of test_one(): refit by code/loo_est_patch.R (run after a first merge)
S[, loo_est_refit := FALSE]; pf <- c(file.path(V2, "results/conditional_final/loo_est_patch.tsv"), file.path(OUT, "loo_est_patch.tsv")); pf <- pf[file.exists(pf)]
if (length(pf)) { P <- unique(rbindlist(lapply(pf, fread))[match == TRUE], by = c("tissue", "promoter", "type", "k", "snp"))
  S[P, on = .(tissue, promoter, type, k, snp), `:=`(loo_est = i.loo_est, loo_est_refit = TRUE)]; cat("loo_est patch applied to", sum(S$loo_est_refit), "signals\n") }
cat("multi-signal signals with NA loo_est:", S[, n := .N, by = .(tissue, promoter, type)][n >= 2, sum(is.na(loo_est))], "\n"); S[, n := NULL]
fwrite(U, file.path(OUT, "freeze_real_units.tsv.gz"), sep = "\t"); fwrite(S, file.path(OUT, "freeze_real_signals.tsv.gz"), sep = "\t")
fwrite(R, file.path(OUT, "freeze_real_rounds.tsv.gz"), sep = "\t"); fwrite(A, file.path(OUT, "freeze_real_arms.tsv"), sep = "\t")
cat("tissues", nrow(A), "; units", nrow(U), "(v3 significant", sum(A$units_v3), "); signals", nrow(S), "; units with >= 2 signals", sum(A$units_ge2),
    "(after backward", sum(A$units_ge2_backward), "); missing chr files", sum(A$missing_chr), "; backward NA", sum(A$backward_na), "\n")
print(U[, .N, by = .(type, n_signals = pmin(n_signals, 4L))][order(type, n_signals)])
