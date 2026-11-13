#!/usr/bin/env Rscript
# LD-window GWAS coloc (author 2026-10-12), step 14: puQTL x GWAS coloc for one tissue, one row per candidate pair
# (puQTL analysis signal x GWAS locus; step 12). Identical to scripts/03_coloc_tissue.R (library sourced unchanged:
# load_marginal(), coloc_signal(); per-class P, N, MAF, priors) except the region, which follows the lab LD-window rule
# exactly as Figure 4's puQTL x eQTL version (QTL_figure/rerun_v2/code/coloc_eqtl_ldwin/06_coloc_tissue.R):
#   puQTL statistics: signal_snps() of the library over the pair's union [uni_start, uni_end] (signal_snps_region below;
#     finite P, MAF in (0, 1), folded MAF, min-p row per position -- unchanged);
#   validity: SNPs with statistics in both arms inside the intersection [int_start, int_end] >= 5, otherwise the pair is
#     not tested (status "different locus (LD windows disjoint)");
#   coloc input: SNPs with statistics in both arms inside the union; >= 10 required (coloc_signal()).
# Every candidate pair is run, whatever its GWAS lead P (< 1e-5); the threshold tiers are applied in step 15.
# Usage: Rscript 14_coloc_tissue.R <tissue> [n_cores]   (idempotent)
suppressPackageStartupMessages(library(parallel))
source("/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/scripts/coloc_final_lib.R")
D <- LDU; MIN_INT <- 5L
a <- commandArgs(TRUE); TIS <- a[1]; NC <- if (length(a) > 1) as.integer(a[2]) else 8L
dir.create(file.path(D, "per_tissue"), showWarnings = FALSE)
OF <- file.path(D, "per_tissue", paste0("coloc_ldwin_", TIS, ".rds"))
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

A <- readRDS(file.path(RES_DIR, "analysis_set.rds"))[tissue == TIS]
P <- readRDS(file.path(D, "pairs", paste0(TIS, ".rds")))
if (!nrow(P)) { saveRDS(data.table(), OF); quit(save = "no") }
P[, sid := paste(promoter, type, k)]; A[, sid := paste(promoter, type, k)]
stopifnot(all(P$sid %in% A$sid))
A <- A[sid %in% P$sid]
loo <- if (any(A$stat == "loo")) { l <- readRDS(file.path(V2, COND_ROOT, "loo", paste0(TIS, ".rds"))); setDT(l)
         l[paste(promoter, type) %in% A[stat == "loo", paste(promoter, type)], .(snp, logp, promoter, type, k)] } else NULL
span <- P[, .(s = min(uni_start), e = max(uni_end)), by = sid]
X <- list()
for (cs in unique(A$chr)) {
  Ac <- A[chr == cs]; marg <- load_marginal(TIS, cs, unique(Ac$promoter))
  for (i in seq_len(nrow(Ac))) { r <- Ac[i]; sp <- span[sid == r$sid]
    X[[r$sid]] <- list(r = r, x = signal_snps_region(r, marg, loo, sp$s, sp$e)) }
  rm(marg); gc(verbose = FALSE)
}
cat(TIS, ":", length(X), "analysis signals with candidate pairs;", nrow(P), "pairs\n")
base_row <- function(r, n_shared, n_int, status) data.table(tissue = r$tissue, promoter = r$promoter, type = r$type, k = r$k,
  n_signals = r$n_signals, stat = r$stat, rule = r$rule, signal_snp = r$signal_snp, gene = r$gene, unit_q = r$q, sig_unit = r$sig,
  n_snps = n_shared, n_int = n_int, status = status)
PK <- c("trait", "gwas_locus_lead", "gwas_locus_id", "lead_p", "lead_absz", "gwas_locus_rank", "gwas_lead_in_panel", "gwas_lead_partners",
        "p_start", "p_end", "g_start", "g_end", "int_start", "int_end", "uni_start", "uni_end")
res <- list()
for (tr in intersect(TRAITS, unique(P$trait))) {
  Pt <- P[trait == tr]
  gw <- readRDS(file.path(D, "gwas_subset", sub("\\.txt$", ".rds", tr)))
  ids <- unique(unlist(lapply(X[unique(Pt$sid)], function(z) z$x$gwas_id)))
  gw <- gw[gwas_id %chin% ids]; setkey(gw, gwas_id)
  Ng <- gwas_N(tr)
  # NB: the signal SNP table has a column "p" (P value), so the pair row is "pr" and its bounds are plain scalars
  rr <- mclapply(seq_len(nrow(Pt)), function(j) { pr <- Pt[j]; z <- X[[pr$sid]]; r <- z$r
    us <- pr$uni_start; ue <- pr$uni_end; i0 <- pr$int_start; i1 <- pr$int_end
    x <- z$x[snp_pos >= us & snp_pos <= ue]
    g <- gw[J(unique(x$gwas_id)), nomatch = 0L]
    sh <- intersect(x$gwas_id, g$gwas_id); pos <- x$snp_pos[match(sh, x$gwas_id)]
    n_int <- sum(pos >= i0 & pos <= i1)
    if (n_int < MIN_INT) return(cbind(base_row(r, length(sh), n_int, "different locus (LD windows disjoint)"), pr[, ..PK]))
    o <- coloc_signal(r, x, g, Ng)
    if (is.null(o)) return(cbind(base_row(r, length(sh), n_int, if (length(sh) < 10) "not assessable (<10 shared SNPs)" else "coloc error"), pr[, ..PK]))
    cbind(o$row[, `:=`(n_int = n_int, status = "tested")], pr[, ..PK]) }, mc.cores = NC)
  bad <- vapply(rr, function(z) inherits(z, "try-error") || is.null(z), TRUE)
  if (any(bad)) stop(TIS, " ", tr, ": ", sum(bad), " pairs failed, e.g. ", as.character(rr[[which(bad)[1]]]))
  res[[tr]] <- rbindlist(rr, fill = TRUE)
}
out <- rbindlist(res, fill = TRUE)
stopifnot(nrow(out) == nrow(P))
saveRDS(out, OF)
cat(sprintf("[%s] %d pairs; tested %d; different locus %d; not assessable %d; HC (tested, PP.H4 >= 0.8) %d, of which lead P < 5e-8: %d\n", TIS,
            nrow(out), out[status == "tested", .N], out[status == "different locus (LD windows disjoint)", .N], out[grepl("^not assessable", status), .N],
            out[status == "tested" & PP.H4 >= 0.8, .N], out[status == "tested" & PP.H4 >= 0.8 & lead_p < 5e-8, .N]))
