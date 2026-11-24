## Mirror rule v3 coloc, step 5: coloc for the NEW rows of one tissue (analysis_set.rds, carried == FALSE), both lead sets:
##   top   : the gene's GTEx eGenes-file top variant (as coloc_eqtl_ldwin/06_coloc_tissue.R; used for non-eGenes and as the
##           "GTEx-top" comparison) -> per_tissue_top/<T>.rds
##   clump : every clumped eQTL lead of an eGene (as coloc_eqtl_ldwin/sens_clump/s5_coloc_tissue_v2.R) -> per_tissue_clump/<T>.rds
## Figure 6's library sourced unchanged (load_marginal, coloc_signal; per-class P, priors, N, MAF). LOO scans from
## conditional_final_v3c. eQTL rows: the existing subsets (coloc_eqtl_ldwin/eqtl_subset, sens_clump/eqtl_subset) plus the
## v3m subsets of the new genes (results/coloc_eqtl_v3m/eqtl_subset_{top,clump}).
## Usage: Rscript 05_coloc_tissue_v3m.R <tissue> [n_cores]   (idempotent)
suppressPackageStartupMessages(library(parallel))
source("/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/scripts/coloc_final_lib.R")
V2L <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/results"; O <- file.path(V2L, "coloc_eqtl_v3m"); LW <- file.path(V2L, "coloc_eqtl_ldwin")
MIN_INT <- 5L
a <- commandArgs(TRUE); TIS <- a[1]; NC <- if (length(a) > 1) as.integer(a[2]) else 8L
for (d in c("per_tissue_top", "per_tissue_clump")) dir.create(file.path(O, d), showWarnings = FALSE)
OT <- file.path(O, "per_tissue_top", paste0(TIS, ".rds")); OC <- file.path(O, "per_tissue_clump", paste0(TIS, ".rds"))
if (file.exists(OT) && file.exists(OC)) { cat("exists, skip\n"); quit(save = "no") }
signal_snps_region <- function(r, marg, loo, start, end) {            # as in coloc_eqtl_ldwin/06_coloc_tissue.R
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
rdE <- function(fs) { fs <- fs[file.exists(fs)]; E <- unique(rbindlist(lapply(fs, function(f) fread(cmd = paste("zcat", f)))), by = c("gene_base", "variant_id")); setkey(E, gene_base, variant_id); E }
A <- readRDS(file.path(O, "analysis_set.rds"))[tissue == TIS & carried == FALSE]
PT <- fread(cmd = paste("zcat", file.path(O, "pairs_top_new_all.tsv.gz")))[tissue == TIS]
PC <- fread(cmd = paste("zcat", file.path(O, "pairs_clump_new_all.tsv.gz")))[tissue == TIS]
stopifnot(nrow(PT) == nrow(A))
A <- merge(A, PT[, .(promoter, type, k, eqtl_top, e_rare, p_start, p_end, e_start, e_end, int_start, int_end, uni_start, uni_end, lead_dist)], by = c("promoter", "type", "k"))
SC <- PC[, .(span_start = min(uni_start), span_end = max(uni_end)), by = .(promoter, type, k)]
A <- merge(A, SC, by = c("promoter", "type", "k"), all.x = TRUE)
A[, `:=`(load_start = pmin(uni_start, span_start, na.rm = TRUE), load_end = pmax(uni_end, span_end, na.rm = TRUE))]
loo <- if (any(A$stat == "loo")) { l <- readRDS(file.path(V2L, "conditional_final_v3c/loo", paste0(TIS, ".rds"))); setDT(l)
         l[paste(promoter, type) %in% A[stat == "loo", paste(promoter, type)], .(snp, logp, promoter, type, k)] } else NULL
ET <- rdE(c(file.path(LW, "eqtl_subset", paste0(TIS, ".tsv.gz")), file.path(O, "eqtl_subset_top", paste0(TIS, ".tsv.gz"))))
EC <- rdE(c(file.path(LW, "sens_clump/eqtl_subset", paste0(TIS, ".tsv.gz")), file.path(O, "eqtl_subset_clump", paste0(TIS, ".tsv.gz"))))
X <- list()
for (cs in unique(A$chr)) {
  Ac <- A[chr == cs]; marg <- load_marginal(TIS, cs, unique(Ac$promoter))
  for (i in seq_len(nrow(Ac))) { r <- Ac[i]; X[[length(X) + 1]] <- list(r = r, x = signal_snps_region(r, marg, loo, r$load_start, r$load_end),
                                                                      pr = PC[promoter == r$promoter & type == r$type & k == r$k]) }
  rm(marg); gc(verbose = FALSE)
}
cat(TIS, ":", length(X), "new analysis signals;", nrow(PC), "clumped signal x lead pairs\n")
## top-variant lead (06 logic)
base_row <- function(r, n_shared, n_int, status) data.table(tissue = r$tissue, promoter = r$promoter, type = r$type, k = r$k, n_signals = r$n_signals,
  stat = r$stat, rule = r$rule, signal_snp = r$signal_snp, gene = r$gene, unit_q = r$q, sig_unit = r$sig, n_snps = n_shared, n_int = n_int, status = status)
run_top <- function(z) { r <- z$r; x <- z$x[snp_pos >= z$r$uni_start & snp_pos <= z$r$uni_end]
  e <- ET[J(r$gene_base, x$snp), nomatch = 0L]
  if (!nrow(e)) return(base_row(r, 0L, 0L, "no eQTL rows"))
  gw <- data.table(gwas_id = x$gwas_id[match(e$variant_id, x$snp)], BETA = e$slope, SE = e$slope_se, Pvalue = e$pval_nominal)[is.finite(BETA) & is.finite(SE) & SE > 0]
  sh <- intersect(x$gwas_id, gw$gwas_id); n_sh <- length(sh); pos <- x$snp_pos[match(sh, x$gwas_id)]
  n_int <- if (r$int_start > r$int_end) 0L else sum(pos >= r$int_start & pos <= r$int_end)
  o <- coloc_signal(r, x, gw, r$n_eqtl)
  if (n_int < MIN_INT) { b <- base_row(r, n_sh, n_int, "different locus (LD windows disjoint)")
    if (!is.null(o)) b[, `:=`(diag_PP.H0 = o$row$PP.H0, diag_PP.H1 = o$row$PP.H1, diag_PP.H2 = o$row$PP.H2, diag_PP.H3 = o$row$PP.H3, diag_PP.H4 = o$row$PP.H4)]
    return(b) }
  if (is.null(o)) return(base_row(r, n_sh, n_int, if (n_sh < 10) "not assessable (<10 shared SNPs)" else "coloc error"))
  o$row[, `:=`(n_int = n_int, status = "tested")][] }
RT <- rbindlist(mclapply(X, run_top, mc.cores = NC), fill = TRUE); stopifnot(nrow(RT) == nrow(A))
RT <- merge(RT, A[, .(promoter, type, k, gene_base, egene, eqtl_qval, n_eqtl, eqtl_top, e_rare, lead_dist, p_start, p_end, e_start, e_end,
                      int_start, int_end, uni_start, uni_end)], by = c("promoter", "type", "k"), all.x = TRUE)
setnames(RT, c("gwas_min_p", "gwas_lead"), c("eqtl_min_p", "eqtl_lead"), skip_absent = TRUE)
saveRDS(RT, OT)
## clumped leads (s5 v2 logic), eGene rows only
run_cl <- function(z) { r <- z$r; xa <- z$x; if (!nrow(z$pr)) return(NULL)
  rbindlist(lapply(seq_len(nrow(z$pr)), function(j) { pr <- z$pr[j]
    b <- data.table(tissue = r$tissue, promoter = r$promoter, type = r$type, k = r$k, gene_base = r$gene_base, lead = pr$lead, rank = pr$rank,
                    lead_dist = pr$lead_dist, n_snps = 0L, n_int = 0L, status = NA_character_)
    x <- xa[snp_pos >= pr$uni_start & snp_pos <= pr$uni_end]
    e <- EC[J(r$gene_base, x$snp), nomatch = 0L]
    if (!nrow(e)) return(b[, status := "no eQTL rows"])
    gw <- data.table(gwas_id = x$gwas_id[match(e$variant_id, x$snp)], BETA = e$slope, SE = e$slope_se, Pvalue = e$pval_nominal)[is.finite(BETA) & is.finite(SE) & SE > 0]
    sh <- intersect(x$gwas_id, gw$gwas_id); pos <- x$snp_pos[match(sh, x$gwas_id)]
    n_int <- if (pr$int_start > pr$int_end) 0L else sum(pos >= pr$int_start & pos <= pr$int_end)
    ni <- n_int; b[, `:=`(n_snps = length(sh), n_int = ni)]   # 2026-11-21: was n_int = n_int, which data.table resolved to b's own 0L column; the filter below always used the local value
    if (n_int < MIN_INT) return(b[, status := "invalid (LD windows disjoint)"])
    o <- coloc_signal(r, x, gw, r$n_eqtl)
    if (is.null(o)) return(b[, status := if (length(sh) < 10) "not assessable (<10 shared SNPs)" else "coloc error"])
    cbind(b[, status := "tested"], o$row[, .(PP.H0, PP.H1, PP.H2, PP.H3, PP.H4, PP.H4_p12_5e6, pqtl_min_p, eqtl_min_p = gwas_min_p,
                                         top_snp_H4, top_snp_PP, eqtl_lead = gwas_lead)])
  }), fill = TRUE) }
RC <- rbindlist(mclapply(X[vapply(X, function(z) isTRUE(z$r$egene), TRUE)], run_cl, mc.cores = NC), fill = TRUE)
stopifnot(nrow(RC) == nrow(PC))
saveRDS(RC, OC)
cat(sprintf("[%s] top: %d rows, tested %d, DL %d, H4 %d, H3 %d | clump: %d pairs, tested %d, H4 %d, H3 %d\n", TIS, nrow(RT), RT[status == "tested", .N],
            RT[status == "different locus (LD windows disjoint)", .N], RT[status == "tested" & PP.H4 >= 0.8, .N], RT[status == "tested" & PP.H3 >= 0.8, .N],
            nrow(RC), RC[status == "tested", .N], RC[status == "tested" & PP.H4 >= 0.8, .N], RC[status == "tested" & PP.H3 >= 0.8, .N]))
