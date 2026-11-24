## Figure 4 panel e, LD-window coloc, step 6: puQTL x eQTL coloc for one tissue under the lab LD-window rule.
## Identical to coloc_eqtl_v3/02_coloc_tissue.R (Figure 6's library sourced unchanged: load_marginal(),
## coloc_signal(); per-class P, priors, N, MAF unchanged) except the region:
##   puQTL statistics: signal_snps() of the library, but over the pair's union [uni_start, uni_end] (step 4) instead of
##     signal +-200 kb -- copied below as signal_snps_region(); everything else in it is unchanged (finite P, MAF in
##     (0, 1), folded MAF, min-p row per position).
##   validity: SNPs with statistics in both datasets inside the intersection [int_start, int_end] >= 5; otherwise the
##     signal is not tested and gets status "different locus (LD windows disjoint)".
##   coloc input: SNPs with statistics in both datasets inside the union; >= 10 required (coloc_signal()).
## Diagnostic only, never used for a call: for "different locus" signals coloc is still run on the union and stored in
## the diag_* columns.
## Status per signal: tested / different locus (LD windows disjoint) / not assessable (<10 shared SNPs) / no eQTL rows /
## coloc error.
## Usage: Rscript 06_coloc_tissue.R <tissue> [n_cores]   (idempotent)
suppressPackageStartupMessages(library(parallel))
source("/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/scripts/coloc_final_lib.R")
V2L <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/results"
OUT <- file.path(V2L, "coloc_eqtl_ldwin")
MIN_INT <- 5L
a <- commandArgs(TRUE); TIS <- a[1]; NC <- if (length(a) > 1) as.integer(a[2]) else 8L
dir.create(file.path(OUT, "per_tissue"), showWarnings = FALSE)
OF <- file.path(OUT, "per_tissue", paste0(TIS, ".rds"))
if (file.exists(OF)) { cat("exists, skip\n"); quit(save = "no") }

signal_snps_region <- function(r, marg, loo, start, end) {
  if (r$stat == "loo") {
    x <- loo[promoter == r$promoter & type == r$type & k == r$k, .(snp, logp)]
    x <- merge(x, marg[promoter == r$promoter, .(snp, snp_pos, maf)], by = "snp")
  } else {
    x <- marg[promoter == r$promoter, .(snp, snp_pos, maf, logp = get(PCOL[[r$type]]))]
  }
  x <- x[snp_pos >= start & snp_pos <= end & is.finite(logp) & !is.na(maf) & maf > 0 & maf < 1]
  x[, p := exp(logp)]; x <- x[p > 0]
  x[maf > 0.5, maf := 1 - maf]
  x[, gwas_id := paste0(sub("^chr", "", r$chr), ":", snp_pos)]
  x[order(p)][!duplicated(gwas_id), .(snp, gwas_id, snp_pos, maf, p)]
}

A <- readRDS(file.path(V2L, "coloc_eqtl_v3/analysis_set.rds"))[tissue == TIS]
P <- fread(cmd = paste("zcat", file.path(OUT, "pairs.tsv.gz")))[tissue == TIS]
A <- merge(A, P[, .(promoter, type, k, eqtl_top, e_rare, p_start, p_end, e_start, e_end, int_start, int_end, uni_start, uni_end, lead_dist)],
           by = c("promoter", "type", "k"))
stopifnot(nrow(A) == nrow(P))
loo <- if (any(A$stat == "loo")) { l <- readRDS(file.path(V2, "conditional_final_v3/loo", paste0(TIS, ".rds"))); setDT(l)
         l[paste(promoter, type) %in% A[stat == "loo", paste(promoter, type)], .(snp, logp, promoter, type, k)] } else NULL
E <- fread(cmd = paste("zcat", file.path(OUT, "eqtl_subset", paste0(TIS, ".tsv.gz"))))
setkey(E, gene_base, variant_id)
X <- list()
for (cs in unique(A$chr)) {
  Ac <- A[chr == cs]; marg <- load_marginal(TIS, cs, unique(Ac$promoter))
  for (i in seq_len(nrow(Ac))) { r <- Ac[i]; X[[length(X) + 1]] <- list(r = r, x = signal_snps_region(r, marg, loo, r$uni_start, r$uni_end)) }
  rm(marg); gc(verbose = FALSE)
}
cat(TIS, ":", length(X), "analysis signals loaded\n")
base_row <- function(r, n_shared, n_int, status) data.table(tissue = r$tissue, promoter = r$promoter, type = r$type, k = r$k, n_signals = r$n_signals,
  stat = r$stat, rule = r$rule, signal_snp = r$signal_snp, gene = r$gene, unit_q = r$q, sig_unit = r$sig, n_snps = n_shared, n_int = n_int, status = status)
run1 <- function(z) { r <- z$r; x <- z$x
  e <- E[J(r$gene_base, x$snp), nomatch = 0L]
  if (!nrow(e)) return(base_row(r, 0L, 0L, "no eQTL rows"))
  gw <- data.table(gwas_id = x$gwas_id[match(e$variant_id, x$snp)], BETA = e$slope, SE = e$slope_se, Pvalue = e$pval_nominal)[is.finite(BETA) & is.finite(SE) & SE > 0]
  sh <- intersect(x$gwas_id, gw$gwas_id); n_sh <- length(sh)
  pos <- x$snp_pos[match(sh, x$gwas_id)]
  n_int <- if (r$int_start > r$int_end) 0L else sum(pos >= r$int_start & pos <= r$int_end)
  o <- coloc_signal(r, x, gw, r$n_eqtl)
  if (n_int < MIN_INT) {
    b <- base_row(r, n_sh, n_int, "different locus (LD windows disjoint)")
    if (!is.null(o)) b[, `:=`(diag_PP.H0 = o$row$PP.H0, diag_PP.H1 = o$row$PP.H1, diag_PP.H2 = o$row$PP.H2, diag_PP.H3 = o$row$PP.H3, diag_PP.H4 = o$row$PP.H4)]
    return(b)
  }
  if (is.null(o)) return(base_row(r, n_sh, n_int, if (n_sh < 10) "not assessable (<10 shared SNPs)" else "coloc error"))
  o$row[, `:=`(n_int = n_int, status = "tested")][] }
res <- rbindlist(mclapply(X, run1, mc.cores = NC), fill = TRUE)
stopifnot(nrow(res) == nrow(A))
res <- merge(res, A[, .(promoter, type, k, gene_base, egene, eqtl_qval, n_eqtl, eqtl_top, e_rare, lead_dist, p_start, p_end, e_start, e_end,
                        int_start, int_end, uni_start, uni_end)], by = c("promoter", "type", "k"), all.x = TRUE)
setnames(res, c("gwas_min_p", "gwas_lead"), c("eqtl_min_p", "eqtl_lead"), skip_absent = TRUE)
saveRDS(res, OF)
cat(sprintf("[%s] %d rows; tested %d; different locus %d; H4>=0.8 %d; H3>=0.8 %d; H1>=0.8 %d\n", TIS, nrow(res), res[status == "tested", .N],
            res[status == "different locus (LD windows disjoint)", .N], res[status == "tested" & PP.H4 >= 0.8, .N],
            res[status == "tested" & PP.H3 >= 0.8, .N], res[status == "tested" & PP.H1 >= 0.8, .N]))
