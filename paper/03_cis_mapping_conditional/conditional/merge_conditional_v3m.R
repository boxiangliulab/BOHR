## merge_conditional_v3m.R — final conditional analysis under mirror rule v3 (author 2026-11-04: exact boundary events are
## labelled p0 on the silent promoter; results/final_nogate_v3m). Sources:
##   results/conditional_v3m (run_conditional_v3m.R; units in work/cond_v3m/units_to_run.tsv: new p0(partner) units, and
##     units whose tissue threshold moved enough to continue / truncate; all re-run from the leader)   -> used for those units
##   results/conditional_final_v3 (units v3, mirror v2; merge_conditional_ng3.R)                        -> every other unit
## Units = the v3m primary significant set (q <= 0.05). v3 units outside it (mostly v2 p1(canonical) of exact pairs, now
## held as p0(partner)) are dropped.
## Backward (single pass, as before): keep a signal of a multi-signal unit if its leave-one-out window ACAT <= thr_v3m. It is
## recomputed for EVERY unit from the stored LOO scans, because 22 tissue thresholds moved slightly.
## loo_est of carried-over units keeps the value of conditional_final_v3 (incl. its est-name patch); re-run units have the
## fixed test_one() (no patch needed).
## Output: results/conditional_final_v3m/freeze_real_{units,signals,rounds}.tsv.gz, freeze_real_arms.tsv, loo/<tissue>.rds
suppressPackageStartupMessages(library(data.table))
V2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"; DEF <- "freeze"
OUT <- file.path(V2, "results", "conditional_final_v3m"); dir.create(file.path(OUT, "loo"), recursive = TRUE, showWarnings = FALSE)
OLD <- file.path(V2, "results", "conditional_final_v3")
acat_log <- function(logp) { p <- exp(logp[is.finite(logp)]); if (!length(p)) return(NA_real_); if (any(p == 0)) return(0)
  if (length(p) > 1) p <- pmin(p, 1 - 1 / length(p)); sm <- p < 1e-15; t <- numeric(length(p)); t[!sm] <- tan((0.5 - p[!sm]) * pi); t[sm] <- 1 / (p[sm] * pi)
  st <- mean(t); if (st > 1e15) return(1 / (st * pi)); pcauchy(st, lower.tail = FALSE) }
z <- function(f) fread(cmd = paste("zcat", f))
U0 <- z(file.path(OLD, "freeze_real_units.tsv.gz")); S0 <- z(file.path(OLD, "freeze_real_signals.tsv.gz")); R0 <- z(file.path(OLD, "freeze_real_rounds.tsv.gz"))
RUN <- fread(file.path(V2, "work/cond_v3m/units_to_run.tsv"))
tab <- function(x, k) if (is.null(x) || !NROW(x[[k]])) data.table() else as.data.table(x[[k]])
U_all <- list(); S_all <- list(); R_all <- list(); A_all <- list(); miss_all <- list()
for (t in sort(unique(U0$tissue))) {
  V <- as.data.table(readRDS(file.path(V2, "results/final_nogate_v3m", DEF, paste0(t, ".real.units.rds"))))
  thr <- V[q <= 0.05, max(p)]; keep <- V[q <= 0.05, paste(promoter, type)]
  rk <- RUN[tissue == t, paste(promoter, type)]
  ## re-run units
  nu <- list(); ns <- list(); nr <- list(); nl <- list()
  for (ch in unique(RUN[tissue == t, chr])) { f <- file.path(V2, "results/conditional_v3m", DEF, t, "real", paste0(ch, ".rds"))
    if (!file.exists(f)) { miss_all[[length(miss_all) + 1]] <- data.table(tissue = t, chr = ch); next }
    x <- readRDS(f); stopifnot(isTRUE(all.equal(x$meta$thr, thr)))
    nu[[ch]] <- tab(x, "units")[!is.na(type)][, chr := ch]; ns[[ch]] <- tab(x, "signals")[, chr := ch]; nr[[ch]] <- tab(x, "rounds"); nl[[ch]] <- tab(x, "loo") }
  nu <- rbindlist(nu, fill = TRUE); ns <- rbindlist(ns, fill = TRUE); nr <- rbindlist(nr, fill = TRUE); nl <- rbindlist(nl, fill = TRUE)
  if (nrow(nu)) { nu[, source := "v3m"]; nr[, source := "v3m"] }
  ## carried-over units (not re-run, still significant)
  ck <- setdiff(keep, rk)
  ou <- U0[tissue == t & paste(promoter, type) %in% ck][, n_signals_backward := NULL]
  os <- S0[tissue == t & paste(promoter, type) %in% ck][, c("loo_acat", "loo_top_snp", "backward_keep") := NULL]
  orr <- R0[tissue == t & paste(promoter, type) %in% ck]
  ol <- readRDS(file.path(OLD, "loo", paste0(t, ".rds"))); ol <- if (NROW(ol)) as.data.table(ol)[paste(promoter, type) %in% ck] else data.table()
  if (nrow(ns)) ns[, loo_est_refit := FALSE]
  u <- rbindlist(list(ou, nu), fill = TRUE); s <- rbindlist(list(os, ns), fill = TRUE); r <- rbindlist(list(orr, nr), fill = TRUE); l <- rbindlist(list(ol, nl), fill = TRUE)
  u[, tissue := t]; s[, tissue := t]; r[, tissue := t]
  stopifnot(setequal(paste(u$promoter, u$type), keep) || nrow(miss_all) > 0)
  ## backward from LOO with the v3m threshold
  if (nrow(l)) { bw <- l[, .(loo_acat = acat_log(logp), loo_top_snp = snp[which.min(replace(logp, !is.finite(logp), Inf))]), by = .(promoter, type, k)]
    s <- merge(s, bw, by = c("promoter", "type", "k"), all.x = TRUE) } else s[, `:=`(loo_acat = NA_real_, loo_top_snp = NA_character_)]
  s[, n_sig := .N, by = .(promoter, type)]
  s[, backward_keep := fifelse(n_sig == 1L, TRUE, fifelse(is.na(loo_acat), NA, loo_acat <= thr))]; s[, n_sig := NULL]
  u <- merge(u, s[, .(n_signals_backward = if (anyNA(backward_keep)) NA_integer_ else sum(backward_keep)), by = .(promoter, type)], by = c("promoter", "type"), all.x = TRUE)
  ## carried-over units: unit_q from the v3m BH (q moves slightly because BH runs over all types of the tissue)
  u[V, on = .(promoter, type), `:=`(unit_p = i.p, unit_q = i.q)]
  saveRDS(l, file.path(OUT, "loo", paste0(t, ".rds")))
  A_all[[t]] <- data.table(tissue = t, thr_v3m = thr, units_v3m = length(keep), units = nrow(u), from_v3m = sum(u$source == "v3m", na.rm = TRUE),
    signals = nrow(s), units_ge2 = sum(u$n_signals >= 2), units_ge2_backward = sum(u$n_signals_backward >= 2, na.rm = TRUE), backward_na = sum(is.na(s$backward_keep)))
  U_all[[t]] <- u; S_all[[t]] <- s; R_all[[t]] <- r }
U <- rbindlist(U_all, fill = TRUE); S <- rbindlist(S_all, fill = TRUE); R <- rbindlist(R_all, fill = TRUE); A <- rbindlist(A_all, fill = TRUE)
M <- rbindlist(miss_all); if (nrow(M)) { cat("MISSING conditional_v3m files:", nrow(M), "\n"); print(head(M)) }
fwrite(U, file.path(OUT, "freeze_real_units.tsv.gz"), sep = "\t"); fwrite(S, file.path(OUT, "freeze_real_signals.tsv.gz"), sep = "\t")
fwrite(R, file.path(OUT, "freeze_real_rounds.tsv.gz"), sep = "\t"); fwrite(A, file.path(OUT, "freeze_real_arms.tsv"), sep = "\t")
cat("tissues", nrow(A), "; units", nrow(U), "(v3m significant", sum(A$units_v3m), "); re-run", sum(A$from_v3m), "; signals", nrow(S),
    "; units >= 2 signals", sum(A$units_ge2), "(after backward", sum(A$units_ge2_backward), "); backward NA", sum(A$backward_na), "\n")
cat("multi-signal signals with NA loo_est:", S[, n := .N, by = .(tissue, promoter, type)][n >= 2, sum(is.na(loo_est))], "\n"); S[, n := NULL]
print(U[, .N, by = .(type, n_signals = pmin(n_signals, 4L))][order(type, n_signals)])
## v3 vs v3m for the moved units: v2 p1(canonical) -> v3m p0(partner)
MP <- fread(cmd = paste("zcat", file.path(V2, "results/mirror_v3/unit_map_v2_to_v3.tsv.gz")))
o <- U0[type == "p1", .(tissue, promoter_v2 = promoter, nsig_v2 = n_signals, nbw_v2 = n_signals_backward)]
n <- U[type == "p0", .(tissue, promoter_v3 = promoter, nsig_v3 = n_signals, nbw_v3 = n_signals_backward)]
cmp <- merge(merge(MP, o, by = c("tissue", "promoter_v2")), n, by = c("tissue", "promoter_v3"))
cat("\nmoved units, signal counts v2 p1(canonical) vs v3m p0(partner): same n_signals", cmp[nsig_v2 == nsig_v3, .N], "of", nrow(cmp),
    "; same backward count", cmp[nbw_v2 == nbw_v3, .N], "\n"); print(cmp[, .N, by = .(nsig_v2 = pmin(nsig_v2, 3L), nsig_v3 = pmin(nsig_v3, 3L))][order(nsig_v2, nsig_v3)])
fwrite(cmp, file.path(OUT, "moved_units_signal_counts.tsv"), sep = "\t")
