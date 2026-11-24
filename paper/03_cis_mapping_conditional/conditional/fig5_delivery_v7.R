## Copy of rerun_v2/code/fig5_delivery_v7.R as run in the final run: the earlier run's counts (13,811 leads, 23 with k != 1) are reported, not asserted.
## fig5_delivery_v7.R: one definition of a puQTL's
## lead signal = the first signal kept after the backward step, i.e. the row with the smallest k among backward_keep == TRUE of
## the unit (tissue x promoter x type). Delivery v7 = delivery v6 with two columns appended to signals.tsv.gz:
##   lead               TRUE on the lead-signal row of each unit, FALSE elsewhere
##   lead_logp_marginal on the lead row: the lead variant's own-type discovery log P (natural log) from the main per-SNP results
##                      (mu: mu_F_logp, F-corrected LRT; p0: nu_firth_logp; p1: tau_firth_logp; signals_pertype.tsv.gz, which
##                      looks them up in the main-run SNP tables, code/fig5_pertype.R). Equals logp_at_selection when k == 1
##                      (copied as text); for the units whose k = 1 was dropped it is the MARGINAL P of the new lead, not its
##                      conditional P. Empty (NA) on non-lead rows.
## The first 31 columns of every row are written byte-identically (text-level append); every other v6 file is copied unchanged.
suppressPackageStartupMessages(library(data.table))
V2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"
D6 <- file.path(V2, "results/fig5_delivery_v6"); D7 <- file.path(V2, "results/fig5_delivery_v7"); dir.create(D7, showWarnings = FALSE)
for (f in setdiff(list.files(D6), c("signals.tsv.gz", "README.txt"))) stopifnot(file.copy(file.path(D6, f), file.path(D7, f), overwrite = TRUE, copy.date = TRUE))
raw <- readLines(gzcon(file(file.path(D6, "signals.tsv.gz"), "rb")))
S <- fread(cmd = paste("zcat", file.path(D6, "signals.tsv.gz")), colClasses = list(character = c("logp_at_selection")), na.strings = NULL)
stopifnot(nrow(S) == length(raw) - 1L, identical(strsplit(raw[1], "\t")[[1]], names(S)))
P <- fread(cmd = paste("zcat", file.path(D6, "signals_pertype.tsv.gz")), colClasses = list(character = c("mu_F_logp", "nu_firth_logp", "tau_firth_logp")), na.strings = NULL)
S[, row := .I]
S[, bk := backward_keep %in% c("TRUE", TRUE)]
S[, lead_k := min(k[bk]), by = .(tissue, promoter, type)]
stopifnot(S[, all(is.finite(lead_k))])
S[, lead := k == lead_k]
S <- merge(S, P[, .(tissue, promoter, type, k, snp, mu_F_logp, nu_firth_logp, tau_firth_logp)], by = c("tissue", "promoter", "type", "k", "snp"), all.x = TRUE, sort = FALSE)
setorder(S, row); stopifnot(S[, identical(row, seq_len(.N))])
S[, own := fcase(type == "mu", mu_F_logp, type == "p0", nu_firth_logp, type == "p1", tau_firth_logp)]
stopifnot(S[lead == TRUE, !anyNA(own) && all(nzchar(own))])
## k == 1 leads: the pertype own-type value must equal logp_at_selection; the text of logp_at_selection is used
stopifnot(S[lead == TRUE & k == 1, all(abs(as.numeric(own) - as.numeric(logp_at_selection)) <= 1e-9 * pmax(1, abs(as.numeric(own))))])
S[, lead_logp_marginal := fifelse(lead & k == 1, logp_at_selection, fifelse(lead, own, ""))]
n_units <- uniqueN(S[, .(tissue, promoter, type)])
cat("signals", nrow(S), "; units", n_units, "; lead TRUE", S[lead == TRUE, .N], "; lead with k != 1", S[lead == TRUE & k != 1, .N], "\n")
stopifnot(S[lead == TRUE, .N] == n_units, S[lead == TRUE, all(bk)])
out <- c(paste0(raw[1], "\tlead\tlead_logp_marginal"), paste0(raw[-1], "\t", ifelse(S$lead, "TRUE", "FALSE"), "\t", S$lead_logp_marginal))
con <- gzfile(file.path(D7, "signals.tsv.gz"), "w"); writeLines(out, con); close(con)
## check: the first 31 columns are byte-identical to v6
chk <- readLines(gzcon(file(file.path(D7, "signals.tsv.gz"), "rb")))
stopifnot(identical(sub("\t[^\t]*\t[^\t]*$", "", chk), raw))
X <- S[lead == TRUE & k != 1, .(tissue, promoter, type, k_lead = k, lead_snp = snp, lead_logp_at_selection = as.numeric(logp_at_selection),
                                 lead_logp_marginal = as.numeric(lead_logp_marginal))]
K1 <- S[k == 1 & paste(tissue, promoter, type) %in% X[, paste(tissue, promoter, type)], .(tissue, promoter, type, old_k1_snp = snp, old_k1_logp = as.numeric(logp_at_selection))]
X <- merge(X, K1, by = c("tissue", "promoter", "type"))
X <- merge(X, S[, .(n_kept = sum(bk), n_forward = .N), by = .(tissue, promoter, type)], by = c("tissue", "promoter", "type"))
setorder(X, tissue, promoter, type)
fwrite(X, file.path(D7, "lead_changes_v6_to_v7.tsv"), sep = "\t")
print(X[, .(tissue, promoter, type, k_lead, old_k1_snp, lead_snp, old_k1_log10p = round(old_k1_logp / log(10), 2), new_marg_log10p = round(lead_logp_marginal / log(10), 2), n_kept, n_forward)])
cat("DONE\n")
