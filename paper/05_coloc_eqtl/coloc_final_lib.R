# Shared library: puQTL (rerun_v2 FINAL ungated unit set v3: mirror v2, mu n_mid >= 80, p0/p1 minor class 10-90%) x GWAS coloc, with per-signal
# (LOO-conditional) statistics for multi-signal units (v3).
# Input rules (unchanged from 10-01):
#   mu : p = exp(mu_F_logp) (F-corrected LRT), type="quant", sdY=1, N = n_mid
#   p0 : p = exp(nu_firth_logp),  type="cc", s = n0/n, N = n      (Firth from main run, NAs filled from the add-on)
#   p1 : p = exp(tau_firth_logp), type="cc", s = n1/n, N = n
#   multi-signal units: signal k uses the LOO scan conditional on the unit's other signals (`logp`: mu F-corrected,
#   p0/p1 Firth). Single-signal units use the marginal scan.
# MAF = main-run `maf` on both arms. GWAS arm: beta/varbeta, type="quant", N from config (as canonical).
# chr:pos match; multi-allelic positions keep the min-p puQTL row.
# Priors p1 = p2 = 1e-4; p12 = 1e-5 MAIN (pre-specified, as published), 5e-6 sensitivity (same lABFs).
suppressPackageStartupMessages({ library(data.table); library(coloc) })

V2   <- "/home/e1327989/project/GTEx/QTL_figure/rerun_v2/results"
NDIR <- "/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3"
## 2026-11-04 mirror rule v3 (author; Figure 4 final_nogate_v3m): exact two-promoter switches are p0
## on the silent TSS. UNITSET = "v3" (default) reproduces the 2026-10-12 build and its paths exactly; "v3m" reads
## Figure 4's v3m units / conditional analysis and writes to parallel folders (results_v3m/, coloc_ldunion_20261012/run_v3m/).
## The GWAS side (significant SNPs, clumping, GWAS windows; LDU0) does not depend on the unit set and is shared.
UNITSET   <- Sys.getenv("UNITSET", "v3"); stopifnot(UNITSET %in% c("v3", "v3m"))
FN_ROOT   <- paste0("final_nogate_", UNITSET)
COND_ROOT <- paste0("conditional_final_", UNITSET)
RES_DIR   <- if (UNITSET == "v3") file.path(NDIR, "results") else file.path(NDIR, paste0("results_", UNITSET))
LDU0      <- file.path(NDIR, "coloc_ldunion_20261012")
LDU       <- if (UNITSET == "v3") LDU0 else file.path(LDU0, paste0("run_", UNITSET))
if (UNITSET != "v3") { dir.create(RES_DIR, showWarnings = FALSE); dir.create(LDU, showWarnings = FALSE) }
GWAS_DIR <- "/home/e1327989/project/GTEx/GWAS_sumstats/prepare_required_columns"
GWAS_SS  <- "/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/scripts/config_gwas_sample_sizes_v3eur.txt"   # 2026-10-09: module copy; eczema 40,835 (EAGLE European, no 23andMe), + Dastani adiponectin 29,347
# 2026-10-09: six GWAS inputs corrected to the intended European results. Four keep their file
# names but are read from gwas_eur_20261009/prepared/ (MAGIC FG / FI / 2-h glucose European; Okada RA European), adiponectin
# is a new file (Dastani 2012, European), eczema keeps its file (already EAGLE's European meta) with N corrected in the config.
# The pre-2026-10-09 library is in scripts/other/pre_eur_20261009/.
GWAS_DIR_EUR <- "/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/gwas_eur_20261009/prepared"
TRAITS   <- readLines("/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/scripts/traits_51_v3eur.txt")
gwas_path <- function(tr) if (file.exists(file.path(GWAS_DIR_EUR, tr))) file.path(GWAS_DIR_EUR, tr) else file.path(GWAS_DIR, tr)
MRAS     <- paste0("PromoterID_", 45111:45114)
WIN <- 200000L; P1 <- 1e-4; P2 <- 1e-4; P12_MAIN <- 1e-5; P12_SENS <- 5e-6
pos_of <- function(snp) as.integer(sub("^chr[^_]+_([0-9]+)_.*$", "\\1", snp))
chr_of <- function(snp) sub("^(chr[^_]+)_.*$", "\\1", snp)

read_units_final <- function() {
  U <- rbindlist(lapply(Sys.glob(file.path(V2, FN_ROOT, "freeze/*.real.units.rds")), function(f) {
    a <- readRDS(f); setDT(a)
    a[, .(tissue, promoter, type, gene, q, leader, n, n0, n1, n_mid, dedup_drop, estimable,
          n_prom_gene, mirror_role, mirror_partner, twin, eq_class, merged_component)] }))
  U[, primary := dedup_drop == FALSE & estimable == TRUE]
  U[, sig := primary & !is.na(q) & q <= 0.05]
  U[]
}

gwas_N <- local({
  g <- tryCatch(fread(GWAS_SS, skip = 5), error = function(e) data.table())
  if (nrow(g)) setnames(g, c("tn", "ss"))
  function(tr) { x <- g[tn == tr, ss]; if (length(x) == 1 && !is.na(x)) as.integer(x) else 50000L }
})

read_gwas <- function(tr) {
  gw <- fread(gwas_path(tr), select = c("SNPID_hg38", "BETA", "SE", "Pvalue"), showProgress = FALSE)
  # 2026-10-09 (supp-table review, Finding 1): four files (Alzheimer IGAP 2013, dementia Bellenguez 2022, HDL/LDL GLGC
  # 2013) carry P values below the double range (e.g. 1.15e-802), so fread reads the whole column as TEXT and min() then
  # returned the alphabetical minimum (gwas_min_p wrong in 155 of their 160 HC rows; posteriors and gwas_lead unaffected).
  # Always numeric: values below ~4.9e-324 become 0 (= underflow, P smaller than representable).
  if (!is.numeric(gw$Pvalue)) gw[, Pvalue := suppressWarnings(as.numeric(Pvalue))]
  stopifnot(is.numeric(gw$Pvalue))
  gw[, c("gchr", "gpos") := tstrsplit(SNPID_hg38, ":", fixed = TRUE)]
  gw[, gpos := suppressWarnings(as.integer(gpos))]
  gw <- gw[!is.na(gpos) & !is.na(BETA) & !is.na(SE) & SE > 0]
  gw[, gwas_id := paste0(gchr, ":", gpos)]
  gw[!duplicated(gwas_id), .(gwas_id, BETA, SE, Pvalue)]
}

# Marginal per-SNP stats for a set of promoters on one chromosome of one tissue, merged exactly as
# code/discovery_nogate.R: main run + muF + Firth add-on (fills NA Firth columns, by promoter x snp).
FC <- c("nu_firth_logp", "nu_firth_est", "tau_firth_logp", "tau_firth_est")
load_marginal <- function(tissue, cs, proms) {
  f <- function(root) file.path(V2, root, "freeze", tissue, "real", paste0(cs, ".rds"))
  M <- readRDS(f("final")); s <- as.data.table(M$snp)[promoter %in% proms, c("promoter", "snp", "snp_pos", "maf", "mu_est", FC), with = FALSE]; rm(M)
  Fm <- readRDS(f("final_muF")); s <- merge(s, as.data.table(Fm$snp)[promoter %in% proms, .(promoter, snp, mu_F_logp)], by = c("promoter", "snp"), all.x = TRUE); rm(Fm)
  if (file.exists(f("final_firth_addon"))) {
    A <- readRDS(f("final_firth_addon"))
    if (NROW(A$snp)) { ad <- as.data.table(A$snp)[promoter %in% proms, c("promoter", "snp", FC), with = FALSE]; setnames(ad, FC, paste0("a_", FC))
      s <- merge(s, ad, by = c("promoter", "snp"), all.x = TRUE)
      for (k in FC) s[is.na(get(k)) & !is.na(get(paste0("a_", k))), (k) := get(paste0("a_", k))]
      s[, paste0("a_", FC) := NULL] }
    rm(A) }
  gc(verbose = FALSE); s[]
}
PCOL <- c(mu = "mu_F_logp", p0 = "nu_firth_logp", p1 = "tau_firth_logp")

# SNP table for one analysis row r (one signal of one unit): marginal or LOO statistic, signal window, min-p dedup.
signal_snps <- function(r, marg, loo = NULL) {
  if (r$stat == "loo") {
    x <- loo[promoter == r$promoter & type == r$type & k == r$k, .(snp, logp)]
    x <- merge(x, marg[promoter == r$promoter, .(snp, snp_pos, maf)], by = "snp")
  } else {
    x <- marg[promoter == r$promoter, .(snp, snp_pos, maf, logp = get(PCOL[[r$type]]))]
  }
  x <- x[snp_pos >= r$signal_pos - WIN & snp_pos <= r$signal_pos + WIN & is.finite(logp) & !is.na(maf) & maf > 0 & maf < 1]
  x[, p := exp(logp)]; x <- x[p > 0]
  x[maf > 0.5, maf := 1 - maf]
  x[, gwas_id := paste0(sub("^chr", "", r$chr), ":", snp_pos)]
  x[order(p)][!duplicated(gwas_id), .(snp, gwas_id, snp_pos, maf, p)]
}

# One coloc for one analysis row; gw must already be the rows for this SNP set (keyed join done by caller).
coloc_signal <- function(r, x, gw, N_gwas, keep_res = FALSE) {
  sh <- intersect(x$gwas_id, gw$gwas_id)
  if (length(sh) < 10) return(NULL)
  q <- x[match(sh, gwas_id)]; g <- gw[match(sh, gwas_id)]
  d1 <- if (r$type == "mu")
    list(pvalues = pmax(q$p, .Machine$double.xmin), N = r$n_mid, MAF = q$maf, snp = sh, type = "quant", sdY = 1)
  else
    list(pvalues = pmax(q$p, .Machine$double.xmin), N = r$n, MAF = q$maf, snp = sh, type = "cc",
         s = if (r$type == "p0") r$n0 / r$n else r$n1 / r$n)
  d2 <- list(beta = g$BETA, varbeta = g$SE^2, N = N_gwas, MAF = q$maf, snp = sh, type = "quant")
  res <- tryCatch({ invisible(capture.output(o <- coloc.abf(d1, d2, p1 = P1, p2 = P2, p12 = P12_MAIN))); o }, error = function(e) NULL)
  if (is.null(res)) return(NULL)
  pp <- res$summary
  sens <- tryCatch({ invisible(capture.output(cb <- coloc:::combine.abf(res$results$lABF.df1, res$results$lABF.df2,
            P1, P2, P12_SENS, quiet = TRUE))); cb }, error = function(e) c(PP.H4.abf = NA_real_))
  rs <- res$results
  row <- data.table(tissue = r$tissue, promoter = r$promoter, type = r$type, k = r$k, n_signals = r$n_signals,
                    stat = r$stat, rule = r$rule, signal_snp = r$signal_snp, gene = r$gene, unit_q = r$q, sig_unit = r$sig,
                    n_snps = length(sh),
                    PP.H0 = pp[["PP.H0.abf"]], PP.H1 = pp[["PP.H1.abf"]], PP.H2 = pp[["PP.H2.abf"]],
                    PP.H3 = pp[["PP.H3.abf"]], PP.H4 = pp[["PP.H4.abf"]], PP.H4_p12_5e6 = unname(sens[["PP.H4.abf"]]),
                    pqtl_min_p = min(q$p), gwas_min_p = min(g$Pvalue, na.rm = TRUE),
                    top_snp_H4 = q$snp[match(rs$snp[which.max(rs$SNP.PP.H4)], sh)], top_snp_PP = max(rs$SNP.PP.H4),
                    gwas_lead = q$snp[which.min(g$Pvalue)])
  list(row = row, res = if (keep_res) res else NULL, q = if (keep_res) q else NULL)
}

cred_set <- function(res, q, level = 0.95) {
  rs <- as.data.table(res$results)[, .(gwas_id = snp, pp = SNP.PP.H4)][order(-pp)]
  rs[, cum := cumsum(pp)]; k <- which(rs$cum >= level)[1]
  cs <- rs[seq_len(k)]; cs[, snp := q$snp[match(gwas_id, q$gwas_id)]]; cs[]
}
