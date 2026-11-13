## Per-genotype phenotype at the LD proxies of every delivered signal, for
## panels e/f tested on lead + proxy SNPs (as the published _leader / _proxy variants).
## Proxy = a SNP the main run tested for the unit's promoter (+-500 kb window, MIN_GENO rule) with r2 >= 0.8 to the
## signal SNP; r2 over the freeze donors (td$U, the genotypes of crosstype_ld). The signal SNP itself is a row (r2 = 1).
## Phenotype convention = genotype_pheno: raw relative usage, freeze donors with a defined value, before the <= 5%
## endpoint deletion, not covariate-adjusted; genotypes modal-imputed (the model's); dosage = alt-allele count.
## Output: results/fig5_delivery/genotype_pheno_proxies.tsv.gz
## Usage: NC=16 Rscript code/fig5_proxies.R [tissue ...]
suppressPackageStartupMessages({ library(data.table); library(parallel) })
CODE <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/code"; source(file.path(CODE, "bohr2_data.R")); source(file.path(CODE, "ld_tool.R"))
## 2026-10-03 (PLINK as the LD standard): proxies are selected by the PLINK rule (raw calls, donors missing either
## SNP dropped per pair; ld_tool.R) and carry the signed r; the phenotype moments keep the model's genotypes.
OUT <- Sys.getenv("FIG5_OUT", file.path(V2, "results/fig5_delivery")); NC <- as.integer(Sys.getenv("NC", "8")); R2MIN <- 0.8
SF <- fread(cmd = paste("zcat", file.path(OUT, "signals.tsv.gz")))
args <- commandArgs(TRUE); TIS <- if (length(args)) args else sort(unique(SF$tissue))
gstat <- function(y, g) { ok <- !is.na(y) & !is.na(g); y <- y[ok]; g <- g[ok]
  f <- function(d, fun) { v <- y[g == d]; if (length(v)) fun(v) else NA_real_ }
  c(n0 = sum(g == 0), n1 = sum(g == 1), n2 = sum(g == 2), m0 = f(0, mean), m1 = f(1, mean), m2 = f(2, mean),
    sd0 = f(0, sd), sd1 = f(1, sd), sd2 = f(2, sd)) }
one_promoter <- function(p, s, td, tested) {
  ch <- as.character(td$info[p, "seqnames"]); pos <- td$info[p, "start"]
  geno <- read_geno(ch, pos - FLANK, pos + FLANK - 1L, td$U); if (is.null(geno)) return(NULL)
  G <- geno$G[, cis_index(geno, pos), drop = FALSE]; G <- G[, colnames(G) %in% tested$snp, drop = FALSE]
  Gr <- ld_raw_region(ch, pos - FLANK, pos + FLANK - 1L, td$U); Gr <- Gr[, colnames(G), drop = FALSE]   # LD = PLINK rule (ld_tool.R)
  y <- td$act[rownames(G), p]
  rbindlist(lapply(seq_len(nrow(s)), function(i) { sn <- s$snp[i]; if (!sn %in% colnames(G)) return(NULL)
    rr <- ld_r(Gr[, sn], Gr); r2 <- rr^2; px <- names(r2)[is.finite(r2) & r2 >= R2MIN]
    px <- unique(c(sn, px[order(-r2[px])]))
    rbindlist(lapply(px, function(x) as.data.table(as.list(gstat(y, G[, x])))[, `:=`(proxy_snp = x, r = if (x == sn) 1 else unname(rr[x]), r2 = if (x == sn) 1 else unname(r2[x]))]))[,
      `:=`(tissue = s$tissue[i], promoter = p, type = s$type[i], k = s$k[i], signal_snp = sn)] })) }
one_tissue <- function(t) {
  td <- readRDS(file.path(V2, "work/final/freeze", paste0(t, ".rds"))); s <- SF[tissue == t]
  out <- list()
  for (ch in unique(s$chr)) {
    M <- as.data.table(readRDS(file.path(V2, "results/final/freeze", t, "real", paste0(ch, ".rds")))$snp)[promoter %in% s$promoter, .(promoter, snp, maf)]
    for (p in unique(s[chr == ch]$promoter)) {
      r <- tryCatch(one_promoter(p, s[promoter == p], td, M[promoter == p]), error = function(e) { message(t, " ", p, ": ", conditionMessage(e)); NULL })
      if (NROW(r)) out[[p]] <- merge(r, M[promoter == p, .(proxy_snp = snp, maf)], by = "proxy_snp", all.x = TRUE) } }
  rbindlist(out) }
RES <- mclapply(TIS, function(t) tryCatch(one_tissue(t), error = function(e) { message(t, ": ", conditionMessage(e)); NULL }), mc.cores = NC, mc.preschedule = FALSE)
ok <- !vapply(RES, is.null, TRUE); if (any(!ok)) cat("FAILED tissues:", TIS[!ok], "\n")
X <- rbindlist(RES[ok])
setcolorder(X, c("tissue", "promoter", "type", "k", "signal_snp", "proxy_snp", "r", "r2", "n0", "n1", "n2", "m0", "m1", "m2", "sd0", "sd1", "sd2", "maf"))
setorder(X, tissue, promoter, type, k, -r2, proxy_snp)
fwrite(X, file.path(OUT, "genotype_pheno_proxies.tsv.gz"), sep = "\t")
S2 <- SF[tissue %in% TIS[ok]]; got <- unique(X[r2 == 1 & proxy_snp == signal_snp, .(tissue, promoter, type, k)])
cat("tissues", sum(ok), "; signals", nrow(S2), "; signals with a self row", nrow(got), "; rows", nrow(X),
    "; proxies per signal (incl. self): median", median(X[, .N, by = .(tissue, promoter, type, k)]$N), "max", max(X[, .N, by = .(tissue, promoter, type, k)]$N), "\n")
