## Per-type MARGINAL tests at every signal SNP (part 1) and at every LD proxy (part 2), and
## the boundary-only partner phenotype at proxy level (part 3). Statistics are looked up from the main-run SNP tables,
## exactly as the ungated discovery assembles them (code/discovery_nogate.R, one_chr):
##   mu : mu_F_logp (F-corrected LRT, primary) and mu_chi_logp (uncorrected), mu_est / mu_se (logit mean), mu_pe
##   p0 : nu_firth_logp / nu_firth_est / nu_firth_se (Firth; main run for old-gate promoters, add-on otherwise), p0_pe
##   p1 : tau_firth_logp / tau_firth_est / tau_firth_se, p1_pe
## Testability (family of the promoter after the <= 5% endpoint rule): p0 needs zeros (BEZI / BEINF), p1 needs ones
## (BEOI / BEINF); statistics of an untestable type are set to NA. mu_estimable = n_mid >= 80 (the unit-set rule).
## All logp are natural-log p. *_pe = marginal probability-scale effect, dosage 0 -> 2. Conditional versions: not
## computed (the LOO scans store only the unit's own type).
## Output (results/fig5_delivery/): signals_pertype.tsv.gz, proxies_pertype.tsv.gz, boundary_only_partner_pheno_proxies.tsv.gz
## Usage: NC=8 Rscript code/fig5_pertype.R [tissue ...]
suppressPackageStartupMessages({ library(data.table); library(parallel) })
CODE <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/code"; source(file.path(CODE, "bohr2_data.R")); source(file.path(CODE, "ld_tool.R"))
OUT <- Sys.getenv("FIG5_OUT", file.path(V2, "results/fig5_delivery")); NC <- as.integer(Sys.getenv("NC", "8"))
rd <- function(f) fread(cmd = paste("zcat", file.path(OUT, f)))
SF <- rd("signals.tsv.gz"); PX <- rd("genotype_pheno_proxies.tsv.gz"); UF <- rd("units.tsv.gz"); BP <- rd("boundary_only_partner_pheno.tsv.gz")
args <- commandArgs(TRUE); TIS <- if (length(args)) args else sort(unique(SF$tissue))
FC <- c("nu_firth_logp", "nu_firth_est", "nu_firth_se", "tau_firth_logp", "tau_firth_est", "tau_firth_se")
KEEP <- c("mu_F_logp", "mu_chi_logp", "mu_est", "mu_se", "mu_pe", "nu_firth_logp", "nu_firth_est", "nu_firth_se", "p0_pe",
          "tau_firth_logp", "tau_firth_est", "tau_firth_se", "p1_pe", "alt_af", "maf")
pos_of <- function(id) as.integer(vapply(strsplit(id, "_", fixed = TRUE), `[`, "", 2))
gstat <- function(y, g) { ok <- !is.na(y) & !is.na(g); y <- y[ok]; g <- g[ok]
  f <- function(d, fun) { v <- y[g == d]; if (length(v)) fun(v) else NA_real_ }
  c(n0 = sum(g == 0), n1 = sum(g == 1), n2 = sum(g == 2), m0 = f(0, mean), m1 = f(1, mean), m2 = f(2, mean), sd0 = f(0, sd), sd1 = f(1, sd), sd2 = f(2, sd)) }
one_chr <- function(t, ch, proms) {   # = discovery_nogate.R one_chr, restricted to the promoters needed
  f <- function(root) file.path(V2, "results", root, "freeze", t, "real", paste0(ch, ".rds"))
  M <- readRDS(f("final")); F <- readRDS(f("final_muF")); A <- readRDS(f("final_firth_addon"))
  s <- as.data.table(M$snp)[promoter %in% proms, c("promoter", "snp", "mu_est", "mu_se", "mu_pe", "p0_pe", "p1_pe", "alt_af", "maf", FC), with = FALSE]
  s <- merge(s, as.data.table(F$snp)[promoter %in% proms, .(promoter, snp, mu_F_logp, mu_chi_logp)], by = c("promoter", "snp"), all.x = TRUE)
  if (NROW(A$snp)) { ad <- as.data.table(A$snp)[promoter %in% proms, c("promoter", "snp", FC), with = FALSE]; setnames(ad, FC, paste0("a_", FC))
    s <- merge(s, ad, by = c("promoter", "snp"), all.x = TRUE)
    for (k in FC) set(s, j = k, value = as.numeric(s[[k]]))   # numeric first: an all-NA column is logical, and := would store the add-on's values as TRUE
    for (k in FC) s[is.na(get(k)) & !is.na(get(paste0("a_", k))), (k) := get(paste0("a_", k))]; s[, paste0("a_", FC) := NULL] }
  pr <- as.data.table(M$prom)[promoter %in% proms, .(promoter, family, n_mid = n - n0 - n1)]
  s <- merge(s, pr, by = "promoter")
  s[!family %in% c("BEZI", "BEINF"), c("nu_firth_logp", "nu_firth_est", "nu_firth_se", "p0_pe") := NA_real_]
  s[!family %in% c("BEOI", "BEINF"), c("tau_firth_logp", "tau_firth_est", "tau_firth_se", "p1_pe") := NA_real_]
  s[, mu_estimable := n_mid >= 80]; s }
one_tissue <- function(t) {
  s <- SF[tissue == t]; px <- PX[tissue == t]; out1 <- list(); out2 <- list()
  for (ch in unique(s$chr)) { pr <- unique(s[chr == ch]$promoter); S <- one_chr(t, ch, pr)
    out1[[ch]] <- merge(s[chr == ch, .(tissue, promoter, type, k, snp)], S, by = c("promoter", "snp"), all.x = TRUE)
    out2[[ch]] <- merge(px[promoter %in% pr, .(tissue, promoter, type, k, signal_snp, proxy_snp, r2)], setnames(copy(S), "snp", "proxy_snp"),
                        by = c("promoter", "proxy_snp"), all.x = TRUE) }
  ## part 3: partner phenotype at every proxy of the canonical's boundary signals (pairs of boundary_only_partner_pheno)
  bp <- unique(BP[tissue == t, .(canonical, partner, canonical_type, partner_type, gene, eq_class)]); out3 <- list()
  if (nrow(bp)) { td <- readRDS(file.path(V2, "work/final/freeze", paste0(t, ".rds")))
    for (i in seq_len(nrow(bp))) { A <- bp$canonical[i]; B <- bp$partner[i]
      x <- px[promoter == A & type == bp$canonical_type[i]]; mus <- s[promoter == B & type == "mu"]
      if (!nrow(x)) next
      ids <- unique(c(x$proxy_snp, mus$snp)); ps <- pos_of(ids); ch <- as.character(td$info[A, "seqnames"])
      g <- read_geno(ch, min(ps), max(ps), td$U); if (is.null(g)) next; G <- g$G[, intersect(ids, colnames(g$G)), drop = FALSE]
      Gr <- ld_raw_region(ch, min(ps), max(ps), td$U)   # LD = PLINK rule (ld_tool.R)
      yB <- td$act[rownames(G), B]
      out3[[i]] <- rbindlist(lapply(seq_len(nrow(x)), function(j) { pxs <- x$proxy_snp[j]; if (!pxs %in% colnames(G)) return(NULL)
        st <- as.data.table(as.list(gstat(yB, G[, pxs])))
        base <- cbind(data.table(tissue = t, gene = bp$gene[i], eq_class = bp$eq_class[i], canonical = A, partner = B, canonical_type = bp$canonical_type[i],
                                 partner_type = bp$partner_type[i], boundary_k = x$k[j], boundary_snp = x$signal_snp[j], proxy_snp = pxs, r2_proxy = x$r2[j]), st)
        if (!nrow(mus)) return(base[, `:=`(mu_k = NA_integer_, mu_snp = NA_character_, r2_proxy_mu = NA_real_)])
        rbindlist(lapply(seq_len(nrow(mus)), function(m) copy(base)[, `:=`(mu_k = mus$k[m], mu_snp = mus$snp[m],
          r2_proxy_mu = if (!is.null(Gr) && all(c(pxs, mus$snp[m]) %in% colnames(Gr))) ld_r(Gr[, pxs], Gr[, mus$snp[m]])^2 else NA_real_)])) })) } }
  list(sig = rbindlist(out1, fill = TRUE), prox = rbindlist(out2, fill = TRUE), bp = rbindlist(out3, fill = TRUE)) }
RES <- mclapply(TIS, function(t) tryCatch(one_tissue(t), error = function(e) { message(t, ": ", conditionMessage(e)); NULL }), mc.cores = NC, mc.preschedule = FALSE)
ok <- !vapply(RES, is.null, TRUE); if (any(!ok)) cat("FAILED tissues:", TIS[!ok], "\n")
pick <- function(k) rbindlist(lapply(RES[ok], `[[`, k), fill = TRUE)
S1 <- pick("sig"); S2 <- pick("prox"); S3 <- pick("bp")
setcolorder(S1, c("tissue", "promoter", "type", "k", "snp", "family", "n_mid", "mu_estimable", KEEP))
setcolorder(S2, c("tissue", "promoter", "type", "k", "signal_snp", "proxy_snp", "r2", "family", "n_mid", "mu_estimable", KEEP))
fwrite(S1, file.path(OUT, "signals_pertype.tsv.gz"), sep = "\t"); fwrite(S2, file.path(OUT, "proxies_pertype.tsv.gz"), sep = "\t")
fwrite(S3, file.path(OUT, "boundary_only_partner_pheno_proxies.tsv.gz"), sep = "\t")
cat("tissues", sum(ok), "; signal rows", nrow(S1), "(missing stats:", S1[is.na(mu_est) & is.na(nu_firth_est) & is.na(tau_firth_est), .N], "); proxy rows", nrow(S2),
    "(missing:", S2[is.na(mu_est) & is.na(nu_firth_est) & is.na(tau_firth_est), .N], "); boundary-partner proxy rows", nrow(S3), "\n")
## own-type consistency: the signal's own-type logp at k = 1 must equal the leader's discovery logp
chk <- merge(S1[k == 1], SF[k == 1, .(tissue, promoter, type, logp_at_selection)], by = c("tissue", "promoter", "type"))
chk[, own := fifelse(type == "mu", mu_F_logp, fifelse(type == "p0", nu_firth_logp, tau_firth_logp))]
cat("k = 1 own-type logp vs leader logp at discovery: max |diff| =", chk[, max(abs(own - logp_at_selection), na.rm = TRUE)], "; NA own:", chk[is.na(own), .N], "\n")
