## =============================================================================
## mirror_v3.R -- mirror rule v3 (author, 2026-11-04, final: exact two-promoter pairs are tested as p0)
##   Pair definition is identical to v2 (boundary_exact / interior_exact / canonical / partner); only which boundary unit is kept changes:
##     v2: for boundary_exact keep p0 and p1 of the canonical and drop p0/p1 of the partner (drop_boundary)
##     v3: for boundary_exact keep p0 of both members and drop p1 of both (p1(A) == p0(B))
##         -> each exact switch event is recorded once, as "the p0 of the TSS that loses all its reads"
##         v2 p1(canonical) -> v3 p0(partner); v2 p0(canonical) is unchanged
##   mu (interior_exact) is unchanged: drop the mu of the partner, keep the mu of the canonical. Non-exact pairs and genes with >= 3 promoters are unchanged.
##   Same columns as the Figure 4 v3 table (QTL_figure/rerun_v2/results/mirror_v3/promoter_flags.tsv.gz):
##     drop_p0 / drop_p1 / dropped_boundary_member / kept_as
## Input: work/mirror_v2.tsv (= trait_regression_20261007/work/mirror.tsv, v2 rule, unchanged)
## Output: work/mirror_v3.tsv, work/mirror_v3_check.tsv
## =============================================================================
source("/projects/group/bliu/e1327989/GTEx/tissue_regression/trait_regression_20261104/code/data.R")
M <- fread(file.path(PATHS$work, "mirror_v2.tsv"))
stopifnot(!anyDuplicated(M[, .(tissue, promoter)]), all(c("drop_mu", "drop_boundary") %in% names(M)))
setnames(M, "drop_boundary", "drop_boundary_v2")
M[, `:=`(drop_p0 = FALSE, drop_p1 = boundary_exact %in% TRUE)]
M[, dropped_boundary_member := fifelse(drop_p1, paste0(promoter, ":p1"), NA_character_)]
M[, kept_as := fifelse(drop_p1, paste0(mirror_partner, ":p0"), NA_character_)]
M[, mirror_rule := "v3"]
fwrite(M, file.path(PATHS$work, "mirror_v3.tsv"), sep = "\t")

## ---- comparison with the Figure 4 v3 table (shared tissue x promoter pairs; our pairs = shared v2 table + pairs completed by this run's matrix) ----
F4 <- fread(cmd = "zcat /projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/results/mirror_v3/promoter_flags.tsv.gz")
X <- merge(M[, .(tissue, promoter, mirror_partner, boundary_exact, drop_p1, kept_as)],
           F4[, .(tissue, promoter, f4_partner = mirror_partner, f4_b = boundary_exact, f4_drop_p1 = drop_p1, f4_kept = kept_as)],
           by = c("tissue", "promoter"))
B <- X[boundary_exact & f4_b %in% TRUE]
CK <- data.table(
  item = c("promoters in mirror_v3.tsv", "boundary_exact promoters (p1 dropped)", "boundary_exact pairs",
           "shared with Figure 4 v3 (same tissue, promoter)", "  same boundary_exact flag", "  same drop_p1",
           "  boundary_exact in both", "    same partner and same kept_as"),
  value = c(nrow(M), sum(M$drop_p1), uniqueN(M[drop_p1 == TRUE, mirror_pair_id]),
            nrow(X), sum(X$boundary_exact == X$f4_b, na.rm = TRUE), sum(X$drop_p1 == X$f4_drop_p1, na.rm = TRUE),
            nrow(B), sum(B$mirror_partner == B$f4_partner & B$kept_as == B$f4_kept, na.rm = TRUE)))
fwrite(CK, file.path(PATHS$work, "mirror_v3_check.tsv"), sep = "\t")
print(CK)
