## fig5_pe_moments_v6.R: copy of fig5_pe_moments_v5.R for Figure 5 delivery v6 (conditional refit v3c). Only the delivery path differs.
## fig5_pe_moments_v5.R: copy of fig5_pe_moments_v4.R for Figure 5 delivery v5 (abs-activity mirror labels, units v3a). Only the delivery path differs.
## fig5_pe_moments_v4.R: copy of fig5_pe_moments.R writing Figure 5 delivery v4 (mirror rule v3: exact switches = p0 on the silent TSS;
## in boundary_only_partner_pheno, columns canonical / partner mean the promoter holding the boundary unit / the other promoter).
## Items 1 and 6, units v3.
## For every signal of the v3 delivery (signals.tsv.gz: all k, including the 25 that failed backward), on the unit's
## promoter, MARGINAL (no conditioning), with the main-run design (cov_matrix -> prepare_y -> orth_design -> promoter_null):
##   item 1: per-type effects on the probability scale (dosage 0 -> 2) with delta-method SEs, from the same fits as the
##           tests (code/pe_delta.R: mu = beta regression; p0/p1 = the Firth full fit). The link-scale est / se of the
##           same fits are repeated and must equal signals_pertype (checked below).
##   item 6: interior-only moments by genotype: donors with 0 < y < 1 of the RAW usage, modal-imputed genotypes, as in
##           genotype_pheno.
## Computed once per (tissue, promoter, snp) and joined to every signal row of that SNP.
## Usage: Rscript fig5_pe_moments.R [tissue ...]   (env NC; OUT_FILE overrides the output, e.g. for tests)
## Output: results/fig5_delivery_v6/signals_pe_moments.tsv.gz
suppressPackageStartupMessages({ library(data.table); library(parallel) })
CODE <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/code"
for (f in c("bohr2_data.R", "bohr2_engine.R", "bohr2_addon3.R", "pe_delta.R")) source(file.path(CODE, f))
OUT <- file.path(V2, "results/fig5_delivery_v6"); NC <- as.integer(Sys.getenv("NC", "16"))
OUT_FILE <- Sys.getenv("OUT_FILE", file.path(OUT, "signals_pe_moments.tsv.gz"))
S <- fread(cmd = paste("zcat", file.path(OUT, "signals.tsv.gz")))[, .(tissue, promoter, type, k, snp, backward_keep)]
args <- commandArgs(TRUE); TIS <- if (length(args)) args else sort(unique(S$tissue))
pos_of <- function(id) as.integer(vapply(strsplit(id, "_", fixed = TRUE), `[`, "", 2))
nz <- function(x) if (is.null(x)) NA else unname(x)
one_tissue <- function(t) {
  td <- readRDS(file.path(V2, "work/final/freeze", paste0(t, ".rds"))); C_all <- cov_matrix(td, td$peer, td$K); s <- S[tissue == t]; out <- list()
  for (p in unique(s$promoter)) {
    z <- prepare_y(td$act[, p], ncol(C_all) - 1L); if (!is.na(z$skip)) { message(t, " ", p, " skip ", z$skip); next }
    smp <- names(z$y); y <- unname(z$y); fam <- z$family; Cs <- orth_design(C_all[smp, , drop = FALSE]); nul <- promoter_null(y, Cs, fam)
    ids <- unique(s[promoter == p]$snp); ps <- pos_of(ids); ch <- as.character(td$info[p, "seqnames"])
    g <- read_geno(ch, min(ps), max(ps), td$U); yr <- td$act[td$U, p]
    for (sn in ids) {
      row <- data.table(tissue = t, promoter = p, snp = sn, family = fam, n = length(y), n_mid = sum(y > 0 & y < 1), geno_ok = !is.null(g) && sn %in% colnames(g$G))
      if (row$geno_ok) { gg <- g$G[smp, sn]
        mu <- tryCatch(pe_mu_fit(y, Cs, gg, nul), error = function(e) NULL)
        ft <- tryCatch(firth_snp_tests_ff(y, Cs, gg, nul), error = function(e) NULL)
        pb <- if (!is.null(ft)) tryCatch(pe_firth(ft$ff, Cs, fam), error = function(e) NULL) else NULL
        row[, `:=`(mu_est = nz(mu$est), mu_se = nz(mu$se), mu_pe = nz(mu$pe), mu_pe_se = nz(mu$pe_se), mu_conv = nz(mu$conv),
                   nu_firth_est = nz(ft$nu_firth_est), nu_firth_se = nz(ft$nu_firth_se), p0_pe_firth = nz(pb$p0_pe), p0_pe_firth_se = nz(pb$p0_pe_se), nu_firth_conv = nz(ft$nu_firth_conv),
                   tau_firth_est = nz(ft$tau_firth_est), tau_firth_se = nz(ft$tau_firth_se), p1_pe_firth = nz(pb$p1_pe), p1_pe_firth_se = nz(pb$p1_pe_se), tau_firth_conv = nz(ft$tau_firth_conv),
                   fit_error = is.null(mu) || is.null(ft) || is.null(pb))]
        row <- cbind(row, as.data.table(as.list(interior_moments(yr, g$G[td$U, sn])))) }
      out[[length(out) + 1]] <- row } }
  rbindlist(out, fill = TRUE) }
R <- rbindlist(mclapply(TIS, function(t) tryCatch(one_tissue(t), error = function(e) { message(t, ": ", conditionMessage(e)); NULL }), mc.cores = NC, mc.preschedule = FALSE), fill = TRUE)
X <- merge(S[tissue %in% TIS], R, by = c("tissue", "promoter", "snp"), all.x = TRUE)
setcolorder(X, c("tissue", "promoter", "type", "k", "snp", "backward_keep")); setorder(X, tissue, promoter, type, k)
fwrite(X, OUT_FILE, sep = "\t")
## checks against signals_pertype (same fits -> identical link-scale est / se; mu_pe identical; p0/p1 pe = Firth vs ML)
SP <- fread(cmd = paste("zcat", file.path(OUT, "signals_pertype.tsv.gz")))[tissue %in% TIS]
C <- merge(X, SP[, .(tissue, promoter, type, k, snp, sp_mu_est = mu_est, sp_mu_se = mu_se, sp_mu_pe = mu_pe, sp_nu_est = nu_firth_est, sp_nu_se = nu_firth_se,
                     sp_tau_est = tau_firth_est, sp_tau_se = tau_firth_se, sp_p0_pe = p0_pe, sp_p1_pe = p1_pe)], by = c("tissue", "promoter", "type", "k", "snp"))
md <- function(a, b) { d <- abs(a - b); c(max = suppressWarnings(max(d, na.rm = TRUE)), na_mismatch = sum(is.na(a) != is.na(b))) }
cat("rows", nrow(X), "of", S[tissue %in% TIS, .N], "; matched to signals_pertype", nrow(C), "; geno ok", X[geno_ok == TRUE, .N], "; fit errors", X[fit_error == TRUE, .N], "\n")
print(rbind(mu_est = md(C$mu_est, C$sp_mu_est), mu_se = md(C$mu_se, C$sp_mu_se), mu_pe = md(C$mu_pe, C$sp_mu_pe), nu_est = md(C$nu_firth_est, C$sp_nu_est),
            nu_se = md(C$nu_firth_se, C$sp_nu_se), tau_est = md(C$tau_firth_est, C$sp_tau_est), tau_se = md(C$tau_firth_se, C$sp_tau_se)))
cat("Firth vs ML pe: cor p0", round(cor(C$p0_pe_firth, C$sp_p0_pe, use = "complete.obs"), 4), "p1", round(cor(C$p1_pe_firth, C$sp_p1_pe, use = "complete.obs"), 4),
    "; median |Firth| / |ML| p0", round(median(abs(C$p0_pe_firth / C$sp_p0_pe), na.rm = TRUE), 3), "p1", round(median(abs(C$p1_pe_firth / C$sp_p1_pe), na.rm = TRUE), 3), "\n")
cat("pe present but SE NA: mu", X[!is.na(mu_pe) & is.na(mu_pe_se), .N], "p0", X[!is.na(p0_pe_firth) & is.na(p0_pe_firth_se), .N], "p1", X[!is.na(p1_pe_firth) & is.na(p1_pe_firth_se), .N],
    "; unconverged mu/nu/tau", X[mu_conv == FALSE, .N], X[nu_firth_conv == FALSE, .N], X[tau_firth_conv == FALSE, .N], "\n")
cat("NA pattern by family (rows with NA pe):\n"); print(X[, .(rows = .N, mu_na = sum(is.na(mu_pe)), p0_na = sum(is.na(p0_pe_firth)), p1_na = sum(is.na(p1_pe_firth)), mid_na = sum(is.na(m0_mid) & is.na(m1_mid) & is.na(m2_mid))), by = family])
