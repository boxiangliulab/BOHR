## v2 (2026-10-12, author: the clumped eQTL leads become the main analysis): identical to
## s5_coloc_tissue.R except that every coloc_signal() column is kept (pqtl_min_p, top_snp_H4, top_snp_PP for S5d) and the
## output goes to sens_clump/per_tissue_v2/. s5_coloc_tissue.R and its per_tissue/ outputs are unchanged.
## Figure 4 panel e, LD-window coloc, eQTL leads by recursive clumping, step s5: coloc for one tissue, one
## row per eGene signal x clumped eQTL lead (s3). Identical to 06_coloc_tissue.R (Figure 6's library sourced unchanged;
## per-class P, priors, N, MAF unchanged; puQTL rows over the pair's union; validity >= 5 SNPs with statistics in both
## datasets inside the intersection; >= 10 shared SNPs in the union), applied to each signal x lead pair.
## Pair status: tested / invalid (LD windows disjoint) / not assessable (<10 shared SNPs) / no eQTL rows / coloc error.
## Usage: Rscript s5_coloc_tissue.R <tissue> [n_cores]   (idempotent)
suppressPackageStartupMessages(library(parallel))
source("/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/scripts/coloc_final_lib.R")
V2L <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/results"
SO  <- file.path(V2L, "coloc_eqtl_ldwin", "sens_clump")
MIN_INT <- 5L
a <- commandArgs(TRUE); TIS <- a[1]; NC <- if (length(a) > 1) as.integer(a[2]) else 8L
dir.create(file.path(SO, "per_tissue_v2"), showWarnings = FALSE)
OF <- file.path(SO, "per_tissue_v2", paste0(TIS, ".rds"))
if (file.exists(OF)) { cat("exists, skip\n"); quit(save = "no") }

signal_snps_region <- function(r, marg, loo, start, end) {            # as in 06_coloc_tissue.R
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

A <- readRDS(file.path(V2L, "coloc_eqtl_v3/analysis_set.rds"))[tissue == TIS & egene == TRUE]
P <- fread(cmd = paste("zcat", file.path(SO, "pairs.tsv.gz")))[tissue == TIS]
S <- P[, .(span_start = min(uni_start), span_end = max(uni_end)), by = .(promoter, type, k)]
A <- merge(A, S, by = c("promoter", "type", "k"))
loo <- if (any(A$stat == "loo")) { l <- readRDS(file.path(V2, "conditional_final_v3/loo", paste0(TIS, ".rds"))); setDT(l)
         l[paste(promoter, type) %in% A[stat == "loo", paste(promoter, type)], .(snp, logp, promoter, type, k)] } else NULL
E <- fread(cmd = paste("zcat", file.path(SO, "eqtl_subset", paste0(TIS, ".tsv.gz"))))
setkey(E, gene_base, variant_id)
X <- list()
for (cs in unique(A$chr)) {
  Ac <- A[chr == cs]; marg <- load_marginal(TIS, cs, unique(Ac$promoter))
  for (i in seq_len(nrow(Ac))) { r <- Ac[i]
    X[[length(X) + 1]] <- list(r = r, x = signal_snps_region(r, marg, loo, r$span_start, r$span_end),
                               pr = P[promoter == r$promoter & type == r$type & k == r$k]) }
  rm(marg); gc(verbose = FALSE)
}
cat(TIS, ":", length(X), "eGene signals loaded;", nrow(P), "signal x lead pairs\n")
run1 <- function(z) { r <- z$r; xa <- z$x
  rbindlist(lapply(seq_len(nrow(z$pr)), function(j) { pr <- z$pr[j]
    b <- data.table(tissue = r$tissue, promoter = r$promoter, type = r$type, k = r$k, gene_base = r$gene_base, lead = pr$lead, rank = pr$rank,
                    lead_dist = pr$lead_dist, n_snps = 0L, n_int = 0L, status = NA_character_)
    x <- xa[snp_pos >= pr$uni_start & snp_pos <= pr$uni_end]
    e <- E[J(r$gene_base, x$snp), nomatch = 0L]
    if (!nrow(e)) return(b[, status := "no eQTL rows"])
    gw <- data.table(gwas_id = x$gwas_id[match(e$variant_id, x$snp)], BETA = e$slope, SE = e$slope_se, Pvalue = e$pval_nominal)[is.finite(BETA) & is.finite(SE) & SE > 0]
    sh <- intersect(x$gwas_id, gw$gwas_id); pos <- x$snp_pos[match(sh, x$gwas_id)]
    n_int <- if (pr$int_start > pr$int_end) 0L else sum(pos >= pr$int_start & pos <= pr$int_end)
    b[, `:=`(n_snps = length(sh), n_int = n_int)]
    if (n_int < MIN_INT) return(b[, status := "invalid (LD windows disjoint)"])
    o <- coloc_signal(r, x, gw, r$n_eqtl)
    if (is.null(o)) return(b[, status := if (length(sh) < 10) "not assessable (<10 shared SNPs)" else "coloc error"])
    cbind(b[, status := "tested"], o$row[, .(PP.H0, PP.H1, PP.H2, PP.H3, PP.H4, PP.H4_p12_5e6, pqtl_min_p, eqtl_min_p = gwas_min_p,
                                         top_snp_H4, top_snp_PP, eqtl_lead = gwas_lead)])
  }), fill = TRUE) }
res <- rbindlist(mclapply(X, run1, mc.cores = NC), fill = TRUE)
stopifnot(nrow(res) == nrow(P))
saveRDS(res, OF)
cat(sprintf("[%s] %d pairs; tested %d; invalid %d; H4>=0.8 %d; H3>=0.8 %d\n", TIS, nrow(res), res[status == "tested", .N],
            res[status == "invalid (LD windows disjoint)", .N], res[status == "tested" & PP.H4 >= 0.8, .N], res[status == "tested" & PP.H3 >= 0.8, .N]))
