## Patch for a labelling bug in the conditional drivers' test_one() (run_conditional*.R, found 2026-10-02):
## for p0/p1, c(est = r$nu_firth_est, ...) inherits the coefficient's name "g" whenever the winning Firth fit started
## from the default start, so the element is named "est.g" and v["est"] returns NA. Only the link-scale `est` of the
## p0/p1 tests is affected; se, logp, selection, stopping, the LOO window ACAT and the backward step use logp/se and are
## not. The main run and the Firth add-on store estimates through data.table columns and are not affected.
## This script refits the leave-one-out (joint) model at the signal SNP itself for every signal of a multi-signal
## p0/p1 unit whose loo_est is NA, checks that se and logp reproduce the stored loo_se / loo_logp, and writes
## <COND_FINAL>/loo_est_patch.tsv. merge_conditional_ng.R and merge_conditional_ng3.R fill loo_est from it (rows with match == TRUE) and set
## loo_est_refit = TRUE. The per-SNP LOO window tables (loo/<tissue>.rds) keep est = NA for those scans.
## Usage: NC=4 Rscript code/loo_est_patch.R
suppressPackageStartupMessages({ library(data.table); library(parallel) })
CODE <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/code"
source(file.path(CODE, "bohr2_data.R")); source(file.path(CODE, "bohr2_engine.R")); source(file.path(CODE, "bohr2_addon3.R"))
FIN <- Sys.getenv("COND_FINAL", file.path(V2, "results/conditional_final")); NC <- as.integer(Sys.getenv("NC", "4"))
SF <- fread(cmd = paste("zcat", file.path(FIN, "freeze_real_signals.tsv.gz")))
SF[, n_sig := .N, by = .(tissue, promoter, type)]
todo <- SF[n_sig >= 2 & type %in% c("p0", "p1") & is.na(loo_est)]
cat("signals to refit:", nrow(todo), "in", uniqueN(todo$tissue), "tissues\n")
pos_of <- function(id) as.integer(vapply(strsplit(id, "_", fixed = TRUE), `[`, "", 2))
one_tissue <- function(t) {
  td <- readRDS(file.path(V2, "work/final/freeze", paste0(t, ".rds"))); C_all <- cov_matrix(td, td$peer, td$K)
  x <- todo[tissue == t]; out <- list()
  for (p in unique(x$promoter)) {
    z <- prepare_y(td$act[, p], ncol(C_all) - 1L); if (!is.na(z$skip)) stop(p, ": skip differs from the main run: ", z$skip)
    s <- names(z$y); y <- unname(z$y); ch <- as.character(td$info[p, "seqnames"])
    for (ty in unique(x[promoter == p]$type)) {
      sig <- SF[tissue == t & promoter == p & type == ty][order(k)]; ps <- pos_of(sig$snp)
      g <- read_geno(ch, min(ps), max(ps), td$U); stopifnot(all(sig$snp %in% colnames(g$G)))
      G <- g$G[s, sig$snp, drop = FALSE]; pre <- if (ty == "p0") "nu" else "tau"
      for (kk in x[promoter == p & type == ty]$k) {
        cs <- sig[k != kk]$snp; sn <- sig[k == kk]$snp            # conditioning set in selection order, as in the driver
        Cc <- orth_design(cbind(C_all[s, , drop = FALSE], G[, cs, drop = FALSE]))
        r <- firth_snp_tests(y, Cc, G[, sn], promoter_null(y, Cc, z$family))
        out[[length(out) + 1]] <- data.table(tissue = t, promoter = p, type = ty, k = kk, snp = sn,
          loo_est = unname(r[[paste0(pre, "_firth_est")]]), se_refit = unname(r[[paste0(pre, "_firth_se")]]),
          logp_refit = unname(r[[paste0(pre, "_firth_logp")]]), se_stored = sig[k == kk]$loo_se, logp_stored = sig[k == kk]$loo_logp)
      } } }
  rbindlist(out) }
TIS <- sort(unique(todo$tissue))
RES <- mclapply(TIS, function(t) tryCatch(one_tissue(t), error = function(e) { message(t, ": ", conditionMessage(e)); NULL }), mc.cores = NC, mc.preschedule = FALSE)
ok <- !vapply(RES, is.null, TRUE); if (any(!ok)) cat("FAILED tissues:", TIS[!ok], "\n")
P <- rbindlist(RES[ok])
P[, match := is.finite(loo_est) & abs(logp_refit - logp_stored) <= 1e-6 * pmax(1, abs(logp_stored)) & abs(se_refit - se_stored) <= 1e-6 * pmax(1, se_stored)]
fwrite(P, file.path(FIN, "loo_est_patch.tsv"), sep = "\t")
cat("refit", nrow(P), "of", nrow(todo), "; match (se and logp reproduce the stored values):", sum(P$match), "\n")
print(P[, .(n = .N, match = sum(match), max_dlogp = max(abs(logp_refit - logp_stored)), max_dse = max(abs(se_refit - se_stored))), by = type])
