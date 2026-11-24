## merge_conditional_v3c.R: final conditional analysis; re-run units replace their v3m rows (two in the earlier run, none in the final run).
## Base: results/conditional_final_v3m (units v3m, v3m tissue thresholds). The units in work/cond_v3c/units_to_run.tsv (if any) are re-run in
## full (forward from the leader, then backward), with the same code (earlier run: run_conditional_v3c.R) -> results/conditional_v3c; their rows
## (units, signals, rounds, LOO) replace the v3m rows. Backward as merge_conditional_v3m.R: a signal of a multi-signal unit is kept if
## its leave-one-out window ACAT <= the tissue threshold; a single signal is kept. Thresholds are unchanged (units and q are v3m).
## Then the saved mirror labels (env LABELS, default results/mirror_abs_labels_v3c/labels_real.tsv.gz) are applied as names, exactly as
## merge_conditional_v3b.R (relabelled event: promoter S -> O, type p0 -> p1, n0/n1 swapped); with 0 changes nothing is renamed.
## Output: results/conditional_final_v3c/freeze_real_{units,signals,rounds}.tsv.gz, freeze_real_arms.tsv, loo/<tissue>.rds
suppressPackageStartupMessages(library(data.table))
V2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"; DEF <- "freeze"
OLD <- file.path(V2, "results/conditional_final_v3m"); OUT <- file.path(V2, "results/conditional_final_v3c"); dir.create(file.path(OUT, "loo"), recursive = TRUE, showWarnings = FALSE)
LABELS <- Sys.getenv("LABELS", file.path(V2, "results/mirror_abs_labels_v3c/labels_real.tsv.gz"))
acat_log <- function(logp) { p <- exp(logp[is.finite(logp)]); if (!length(p)) return(NA_real_); if (any(p == 0)) return(0)
  if (length(p) > 1) p <- pmin(p, 1 - 1 / length(p)); sm <- p < 1e-15; t <- numeric(length(p)); t[!sm] <- tan((0.5 - p[!sm]) * pi); t[sm] <- 1 / (p[sm] * pi)
  st <- mean(t); if (st > 1e15) return(1 / (st * pi)); pcauchy(st, lower.tail = FALSE) }
z <- function(f) fread(cmd = paste("zcat", f))
tab <- function(x, k) if (is.null(x) || !NROW(x[[k]])) data.table() else as.data.table(x[[k]])
U <- z(file.path(OLD, "freeze_real_units.tsv.gz")); S <- z(file.path(OLD, "freeze_real_signals.tsv.gz")); R <- z(file.path(OLD, "freeze_real_rounds.tsv.gz"))
A0 <- fread(file.path(OLD, "freeze_real_arms.tsv"))
RUN <- fread(file.path(V2, "work/cond_v3c/units_to_run.tsv"))
U[, refit_v3c := FALSE]
for (t in unique(RUN$tissue)) {
  rk <- RUN[tissue == t, paste(promoter, type)]; thr <- A0[tissue == t, thr_v3m]
  nu <- list(); ns <- list(); nr <- list(); nl <- list()
  for (ch in unique(RUN[tissue == t, chr])) { x <- readRDS(file.path(V2, "results/conditional_v3c", DEF, t, "real", paste0(ch, ".rds")))
    stopifnot(isTRUE(all.equal(x$meta$thr, thr)))
    nu[[ch]] <- tab(x, "units")[!is.na(type)][, chr := ch]; ns[[ch]] <- tab(x, "signals")[, chr := ch]; nr[[ch]] <- tab(x, "rounds"); nl[[ch]] <- tab(x, "loo") }
  nu <- rbindlist(nu, fill = TRUE)[, source := "v3c refit"]; ns <- rbindlist(ns, fill = TRUE)[, loo_est_refit := FALSE]; nr <- rbindlist(nr, fill = TRUE)[, source := "v3c refit"]
  nl <- rbindlist(nl, fill = TRUE); stopifnot(setequal(paste(nu$promoter, nu$type), rk))
  if (nrow(nl)) { bw <- nl[, .(loo_acat = acat_log(logp), loo_top_snp = snp[which.min(replace(logp, !is.finite(logp), Inf))]), by = .(promoter, type, k)]
    ns <- merge(ns, bw, by = c("promoter", "type", "k"), all.x = TRUE) } else ns[, `:=`(loo_acat = NA_real_, loo_top_snp = NA_character_)]
  ns[, n_sig := .N, by = .(promoter, type)]
  ns[, backward_keep := fifelse(n_sig == 1L, TRUE, fifelse(is.na(loo_acat), NA, loo_acat <= thr))][, n_sig := NULL]
  nu <- merge(nu, ns[, .(n_signals_backward = if (anyNA(backward_keep)) NA_integer_ else sum(backward_keep)), by = .(promoter, type)], by = c("promoter", "type"), all.x = TRUE)
  old_u <- U[tissue == t & paste(promoter, type) %in% rk]
  nu[old_u, on = .(promoter, type), `:=`(unit_p = i.unit_p, unit_q = i.unit_q)]; nu[, `:=`(tissue = t, refit_v3c = TRUE)]; ns[, tissue := t]; nr[, tissue := t]
  stopifnot(all(nu$n_signals_backward >= 1))                              # a significant puQTL keeps at least one signal
  U <- rbind(U[!(tissue == t & paste(promoter, type) %in% rk)], nu, fill = TRUE)
  S <- rbind(S[!(tissue == t & paste(promoter, type) %in% rk)], ns, fill = TRUE)
  R <- rbind(R[!(tissue == t & paste(promoter, type) %in% rk)], nr, fill = TRUE)
  l <- readRDS(file.path(OLD, "loo", paste0(t, ".rds"))); l <- if (NROW(l)) as.data.table(l)[!paste(promoter, type) %in% rk] else data.table()
  saveRDS(rbind(l, nl, fill = TRUE), file.path(OUT, "loo", paste0(t, ".rds"))) }
for (f in list.files(file.path(OLD, "loo"), full.names = TRUE)) if (!file.exists(file.path(OUT, "loo", basename(f)))) file.copy(f, file.path(OUT, "loo", basename(f)))
## mirror labels as names (merge_conditional_v3b.R logic)
## LABELS=none: no renames (first pass; the label step then reads this table's signals, and the merge is re-run with the labels
## only if the label step changes a name).
CH <- if (LABELS == "none") data.table(tissue = character(0), S = character(0), O = character(0)) else fread(cmd = paste("zcat", LABELS))[changed == TRUE, .(tissue, S, O)]
ren <- function(d) { d <- copy(d); d[, relabelled_from := NA_character_]
  for (i in seq_len(nrow(CH))) { r <- CH[i]; w <- d$tissue == r$tissue & d$promoter == r$S & d$type == "p0"
    d[w, `:=`(promoter = r$O, type = "p1", relabelled_from = paste0(r$S, ":p0"))]
    if (all(c("n0", "n1") %in% names(d))) d[w, c("n0", "n1") := .(n1, n0)] }
  d }
U <- ren(U); S <- ren(S); R <- ren(R)
for (i in seq_len(nrow(CH))) { r <- CH[i]; f <- file.path(OUT, "loo", paste0(r$tissue, ".rds")); l <- as.data.table(readRDS(f))
  w <- l$promoter == r$S & l$type == "p0"; l[w, `:=`(promoter = r$O, type = "p1")]; saveRDS(l, f) }
setorder(U, tissue, promoter, type); setorder(S, tissue, promoter, type, k); setorder(R, tissue, promoter, type, round)
fwrite(U, file.path(OUT, "freeze_real_units.tsv.gz"), sep = "\t"); fwrite(S, file.path(OUT, "freeze_real_signals.tsv.gz"), sep = "\t")
fwrite(R, file.path(OUT, "freeze_real_rounds.tsv.gz"), sep = "\t"); file.copy(file.path(OLD, "freeze_real_arms.tsv"), file.path(OUT, "freeze_real_arms.tsv"), overwrite = TRUE)
U0 <- z(file.path(OLD, "freeze_real_units.tsv.gz")); S0 <- z(file.path(OLD, "freeze_real_signals.tsv.gz"))
cat("mirror relabels applied:", nrow(CH), "\n")
cat("units", nrow(U), "(v3m", nrow(U0), "); signals", nrow(S), "(v3m", nrow(S0), "); kept", S[backward_keep %in% TRUE, .N], "(v3m", S0[backward_keep %in% TRUE, .N], ")\n")
cat("units with >= 2 signals: forward", U[n_signals >= 2, .N], "(v3m", U0[n_signals >= 2, .N], "); after backward", U[n_signals_backward >= 2, .N], "(v3m", U0[n_signals_backward >= 2, .N], ")\n")
print(merge(U0[, .(v3m = .N, ge2_bw_v3m = sum(n_signals_backward >= 2, na.rm = TRUE), kept_v3m = sum(n_signals_backward, na.rm = TRUE)), by = type],
            U[, .(v3c = .N, ge2_bw_v3c = sum(n_signals_backward >= 2, na.rm = TRUE), kept_v3c = sum(n_signals_backward, na.rm = TRUE)), by = type], by = "type"))
cat("\nrefitted units: old vs new signal sets\n")
for (i in seq_len(nrow(RUN))) { r <- RUN[i]
  cat(sprintf("%s %s %s\n  v3m: %s\n  v3c: %s\n", r$tissue, r$promoter, r$type,
      S0[tissue == r$tissue & promoter == r$promoter & type == r$type, paste(sprintf("k%d %s logp %.2f kept %s", k, snp, logp_at_selection, backward_keep), collapse = " | ")],
      S[tissue == r$tissue & ((promoter == r$promoter & type == r$type) | relabelled_from %in% paste0(r$promoter, ":", r$type)), paste(sprintf("k%d %s logp %.2f kept %s", k, snp, logp_at_selection, backward_keep), collapse = " | ")])) }
