## Figure 4 panel e, step 2: puQTL x eQTL coloc for one tissue, with Figure 6's library (coloc_final_lib.R) unchanged:
##   analysis rows (step 0), puQTL statistics via signal_snps() (marginal or LOO scan, signal +-200 kb, the min-p row
##   per position), and coloc_signal() (coloc.abf; p1 = p2 = 1e-4; p12 = 1e-5 main, 5e-6 sensitivity; mu: quant,
##   sdY = 1, N = n_mid; p0/p1: cc, s = n0/n or n1/n, N = n; MAF from the main run).
## Trait arm = the eQTL of the unit's own gene in the same tissue (GTEx v11 all-pairs): beta = slope, varbeta =
## slope_se^2, N = the tissue's freeze donors, type quant (as the GWAS arm). Variants are matched on the full ID
## (chr_pos_ref_alt) and then carry Figure 6's gwas_id key, so coloc sees exactly the puQTL rows of the GWAS analysis.
## Every analysis row gets a status: tested / not assessable (< 10 shared SNPs) / no eQTL rows / coloc error.
## Usage: Rscript 02_coloc_tissue.R <tissue> [n_cores]   (idempotent)
suppressPackageStartupMessages(library(parallel))
source("/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/scripts/coloc_final_lib.R")
OUT <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/results/coloc_eqtl_v3"
a <- commandArgs(TRUE); TIS <- a[1]; NC <- if (length(a) > 1) as.integer(a[2]) else 8L
OF <- file.path(OUT, "per_tissue", paste0(TIS, ".rds"))
if (file.exists(OF)) { cat("exists, skip\n"); quit(save = "no") }
A <- readRDS(file.path(OUT, "analysis_set.rds"))[tissue == TIS]
loo <- if (any(A$stat == "loo")) { l <- readRDS(file.path(V2, "conditional_final_v3/loo", paste0(TIS, ".rds"))); setDT(l)
         l[paste(promoter, type) %in% A[stat == "loo", paste(promoter, type)], .(snp, logp, promoter, type, k)] } else NULL
E <- fread(cmd = paste("zcat", file.path(OUT, "eqtl_subset", paste0(TIS, ".tsv.gz"))))
setkey(E, gene_base, variant_id)
X <- list()
for (cs in unique(A$chr)) {
  Ac <- A[chr == cs]; marg <- load_marginal(TIS, cs, unique(Ac$promoter))
  for (i in seq_len(nrow(Ac))) { r <- Ac[i]; X[[length(X) + 1]] <- list(r = r, x = signal_snps(r, marg, loo)) }
  rm(marg); gc(verbose = FALSE)
}
cat(TIS, ":", length(X), "analysis signals loaded\n")
base_row <- function(r, n_shared, status) data.table(tissue = r$tissue, promoter = r$promoter, type = r$type, k = r$k, n_signals = r$n_signals,
  stat = r$stat, rule = r$rule, signal_snp = r$signal_snp, gene = r$gene, unit_q = r$q, sig_unit = r$sig, n_snps = n_shared, status = status)
run1 <- function(z) { r <- z$r; x <- z$x
  e <- E[J(r$gene_base, x$snp), nomatch = 0L]
  if (!nrow(e)) return(base_row(r, 0L, "no eQTL rows"))
  gw <- data.table(gwas_id = x$gwas_id[match(e$variant_id, x$snp)], BETA = e$slope, SE = e$slope_se, Pvalue = e$pval_nominal)[is.finite(BETA) & is.finite(SE) & SE > 0]
  n_sh <- length(intersect(x$gwas_id, gw$gwas_id))
  o <- coloc_signal(r, x, gw, r$n_eqtl)
  if (is.null(o)) return(base_row(r, n_sh, if (n_sh < 10) "not assessable (<10 shared SNPs)" else "coloc error"))
  o$row[, status := "tested"][] }
res <- rbindlist(mclapply(X, run1, mc.cores = NC), fill = TRUE)
res <- merge(res, A[, .(promoter, type, k, gene_base, egene, eqtl_qval, n_eqtl)], by = c("promoter", "type", "k"), all.x = TRUE)
setnames(res, c("gwas_min_p", "gwas_lead"), c("eqtl_min_p", "eqtl_lead"), skip_absent = TRUE)
saveRDS(res, OF)
cat(sprintf("[%s] %d rows; tested %d; H4>=0.8 %d; H3>=0.8 %d; H1>=0.8 %d\n", TIS, nrow(res), res[status == "tested", .N],
            res[status == "tested" & PP.H4 >= 0.8, .N], res[status == "tested" & PP.H3 >= 0.8, .N], res[status == "tested" & PP.H1 >= 0.8, .N]))
