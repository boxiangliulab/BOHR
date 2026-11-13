## =============================================================================
## mirror.R -- mirror-pair de-duplication table of run 20261005; rule = Figure 4 shared table v2 (author decision)
##   Authoritative table: QTL_figure/rerun_v2/results/mirror/{promoter_flags,mirror_pairs}.tsv.gz (v2, decided per group)
##   Definitions (genes with exactly two promoters in the within-tissue matrix):
##     boundary_exact : on every donor with values for both, I(y_A=0)==I(y_B=1) and I(y_A=1)==I(y_B=0) (indicator variables, no tolerance)
##                      -> p0(A) == p1(B), p1(A) == p0(B); drop the partner's p0/p1 (drop_boundary)
##     interior_exact : boundary_exact, and on donors with intermediate values |y_A + y_B - 1| <= 1e-6
##                      -> mu(A) == mu(B) (opposite sign); drop the partner's mu (drop_mu)
##     canonical      : the one with the higher mean usage over all pooled freeze samples (the pooled mean in the shared table is used first, so it agrees with every figure)
##   The promoter set of run 20261005 (the TPM >= 1 set of v13) overlaps the shared table by about 80% and has 6 more small tissues:
##     pairs (same partner) found in the shared table use its flags directly; the other pairs are completed on the run 20261005 matrix by the same rule; the overlap is re-checked per group
##   Author 2026-10-01: mirroring is the biology of promoter switching itself; de-duplication only affects counts and BH, not the meaning of the results; both promoters are drawn in figures.
## Output: work/mirror.tsv, work/mirror_check.tsv
## =============================================================================
source("/projects/group/bliu/e1327989/GTEx/tissue_regression/trait_regression_20261005/code/data.R")
MD <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/results/mirror"
SH <- fread(cmd = paste0("zcat ", MD, "/promoter_flags.tsv.gz"))
SP <- fread(cmd = paste0("zcat ", MD, "/mirror_pairs.tsv.gz"))
stopifnot(all(c("boundary_exact", "interior_exact", "drop_mu", "drop_boundary") %in% names(SH)))   # must be v2
TT <- fread(file.path(PATHS$work, "tissues.tsv"))[donors > 0]
TD <- lapply(TT$tissue, function(t) readRDS(file.path(PATHS$work, paste0(t, ".rds")))); names(TD) <- TT$tissue

## Pooled mean usage for the canonical rule: shared table first, otherwise all freeze samples of run 20261005
sums <- list(); cnts <- list()
for (t in names(TD)) { A <- TD[[t]]$act; s <- colSums(A, na.rm = TRUE); k <- colSums(!is.na(A))
  for (p in names(s)) { sums[[p]] <- c(sums[[p]], s[[p]]); cnts[[p]] <- c(cnts[[p]], k[[p]]) } }
pooled <- vapply(names(sums), function(p) sum(sums[[p]]) / sum(cnts[[p]]), 0)
PM <- unique(rbind(SP[, .(promoter = promoter_1, pm = pooled_mean_1)], SP[, .(promoter = promoter_2, pm = pooled_mean_2)]))[, .(pm = pm[1]), by = promoter]
k <- PM$promoter %in% names(pooled); pooled[PM$promoter[k]] <- PM$pm[k]

rows <- list()
for (t in names(TD)) {
  td <- TD[[t]]; inf <- data.table(promoter = rownames(td$info), gene = td$info$geneId)
  inf[, n_prom_gene := .N, by = gene]
  for (g in inf[n_prom_gene == 2, unique(gene)]) {
    pr <- inf[gene == g, promoter]; A <- td$act[, pr[1]]; B <- td$act[, pr[2]]; ok <- !is.na(A) & !is.na(B)
    a <- A[ok]; b <- B[ok]
    bex <- length(a) > 0 && all((a == 0) == (b == 1)) && all((a == 1) == (b == 0))
    im <- a > 0 & a < 1
    iex <- bex && all((b > 0 & b < 1) == im) && (!any(im) || max(abs(a[im] + b[im] - 1)) <= 1e-6)
    can <- pr[which.max(pooled[pr])]; par <- setdiff(pr, can)
    rows[[length(rows) + 1]] <- data.table(tissue = t, promoter = c(can, par), gene = g, n_prom_gene = 2L,
                                           b_local = bex, i_local = iex, role_local = c("canonical", "partner"),
                                           partner_local = c(par, can))
  }
}
L <- rbindlist(rows)
M <- merge(L, SH[, .(tissue, promoter, sh_b = boundary_exact, sh_i = interior_exact, sh_role = mirror_role, sh_partner = mirror_partner)],
           by = c("tissue", "promoter"), all.x = TRUE)
M[, covered := !is.na(sh_b) & !is.na(sh_partner) & sh_partner == partner_local]
M[, `:=`(boundary_exact = fifelse(covered, sh_b, b_local), interior_exact = fifelse(covered, sh_i, i_local),
         source = fifelse(covered, "shared_v2", "local_v2"))]
M[, mirror_role := fifelse(boundary_exact, fifelse(covered & !is.na(sh_role), sh_role, role_local), NA_character_)]
M[, mirror_partner := fifelse(boundary_exact, partner_local, NA_character_)]
M[, `:=`(exact_mirror = boundary_exact & interior_exact,
         drop_mu = interior_exact & mirror_role %in% "partner",
         drop_boundary = boundary_exact & mirror_role %in% "partner",
         mirror_pair_id = paste(tissue, gene, sep = "|"))]
fwrite(M[, .(tissue, promoter, gene, n_prom_gene, mirror_pair_id, boundary_exact, interior_exact, exact_mirror,
             mirror_role, mirror_partner, drop_mu, drop_boundary, source)], file.path(PATHS$work, "mirror.tsv"), sep = "\t")

CK <- M[, .(pairs = uniqueN(gene), covered = uniqueN(gene[covered]),
            boundary_exact = uniqueN(gene[boundary_exact]), interior_exact = uniqueN(gene[interior_exact]),
            disagree_boundary = uniqueN(gene[covered & sh_b != b_local]), disagree_interior = uniqueN(gene[covered & sh_i != i_local]),
            disagree_canonical = uniqueN(gene[covered & sh_b & b_local & sh_role != role_local])), by = tissue]
fwrite(CK, file.path(PATHS$work, "mirror_check.tsv"), sep = "\t")
tot <- CK[, lapply(.SD, sum), .SDcols = -1]
print(tot)
cat(sprintf("boundary_exact %.1f%%, interior_exact %.1f%% of %d two-promoter pairs; %d pairs covered by the shared v2 table\n",
            100 * tot$boundary_exact / tot$pairs, 100 * tot$interior_exact / tot$pairs, tot$pairs, tot$covered))
