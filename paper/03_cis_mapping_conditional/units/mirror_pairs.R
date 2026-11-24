## Authoritative mirror table (de-duplicate EXACT mirror pairs only, keep the types).
## v2: exactness is judged PER COMPONENT.
##   boundary_exact: I(y_A = 0) == I(y_B = 1) and I(y_A = 1) == I(y_B = 0) in every observed donor (exact indicators)
##                   -> the categorical data are swapped exactly, so p0(A) == p1(B) and p1(A) == p0(B).
##   interior_exact: boundary_exact, same interior donor set, and |y_A + y_B - 1| <= 1e-6 on it -> mu(A) == mu(B).
##   A third promoter that takes a share only in interior samples breaks interior_exact but not boundary_exact.
## (v1 used one strict flag = interior_exact; it missed 5,062 boundary-only twins. v1 kept in results/mirror/v1_strict.)
## A tissue x gene is an exact mirror pair when the gene has exactly two promoters in that tissue's freeze activity
## matrix and y_A + y_B = 1 (|dev| < 1e-6) in every sample where both are observed (no sample with both = 0).
## Canonical member: the promoter with the higher mean relative usage pooled over all freeze samples of all 50
## tissues (fixed, genotype-independent; ties -> smaller PromoterID number). The partner's units duplicate the
## canonical's:  mu(partner) = mu(canonical), sign flipped;  p0(partner) = p1(canonical);  p1(partner) = p0(canonical).
## Outputs (results/mirror/): mirror_pairs.tsv.gz (one row per tissue x two-promoter gene),
##   promoter_flags.tsv.gz (one row per tissue x promoter of the freeze matrices)
suppressPackageStartupMessages(library(data.table))
V2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"; OUT <- file.path(V2, "results", "mirror")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
fs <- list.files(file.path(V2, "work/final/freeze"), pattern = "\\.rds$", full.names = TRUE)
TD <- lapply(fs, readRDS); TD <- TD[vapply(TD, function(x) is.list(x) && !is.null(x$act) && !is.null(x$info), TRUE)]
names(TD) <- vapply(TD, function(x) x$tissue, "")
cat("tissues:", length(TD), "\n")
pool <- rbindlist(lapply(TD, function(td) data.table(promoter = colnames(td$act), s = colSums(td$act, na.rm = TRUE), k = colSums(!is.na(td$act)))))[
  , .(pooled_mean = sum(s) / sum(k)), by = promoter]
pm <- setNames(pool$pooled_mean, pool$promoter)
num <- function(p) as.numeric(sub("PromoterID_", "", p))
L <- lapply(names(TD), function(t) { td <- TD[[t]]; g <- td$info[colnames(td$act), "geneId"]; np <- table(g)
  P <- data.table(tissue = t, promoter = colnames(td$act), gene = g, n_prom_gene = as.integer(np[g]))
  two <- names(np)[np == 2]
  M <- rbindlist(lapply(two, function(x) { pr <- colnames(td$act)[g == x]; pr <- pr[order(num(pr))]; a <- td$act[, pr[1]]; b <- td$act[, pr[2]]
    ok <- !is.na(a) & !is.na(b); d <- abs(a[ok] + b[ok] - 1); aa <- a[ok]; bb <- b[ok]
    za <- aa == 0; oa <- aa == 1; zb <- bb == 0; ob <- bb == 1; ia <- !za & !oa; ib <- !zb & !ob
    bex <- any(ok) && all(za == ob) && all(oa == zb)
    iex <- bex && all(ia == ib) && (!any(ia) || max(abs(aa[ia] + bb[ia] - 1)) <= 1e-6)
    data.table(tissue = t, gene = x, promoter_1 = pr[1], promoter_2 = pr[2], n_samples = sum(ok), boundary_exact = bex, interior_exact = iex,
               max_dev = if (any(ok)) max(d) else NA_real_,
               n_both0 = sum(a[ok] == 0 & b[ok] == 0), tissue_mean_1 = mean(a, na.rm = TRUE), tissue_mean_2 = mean(b, na.rm = TRUE)) }))
  list(P = P, M = M) })
P <- rbindlist(lapply(L, `[[`, "P")); M <- rbindlist(lapply(L, `[[`, "M"))
M[, exact_mirror := interior_exact & boundary_exact]          # = v1 'strict' flag (both components)
M[, `:=`(pooled_mean_1 = pm[promoter_1], pooled_mean_2 = pm[promoter_2])]
M[, canonical := fifelse(pooled_mean_1 >= pooled_mean_2, promoter_1, promoter_2)]
M[, partner := fifelse(canonical == promoter_1, promoter_2, promoter_1)]
M[, mirror_pair_id := paste(tissue, gene, sep = "|")]
## per tissue x promoter flags
P <- merge(P, rbind(M[, .(tissue, promoter = canonical, mirror_pair_id, boundary_exact, interior_exact, exact_mirror, mirror_role = "canonical", mirror_partner = partner)],
                    M[, .(tissue, promoter = partner, mirror_pair_id, boundary_exact, interior_exact, exact_mirror, mirror_role = "partner", mirror_partner = canonical)]),
           by = c("tissue", "promoter"), all.x = TRUE)
for (k in c("boundary_exact", "interior_exact", "exact_mirror")) P[is.na(get(k)), (k) := FALSE]
P[, drop_mu := interior_exact & mirror_role %in% "partner"]          # mu(partner) == mu(canonical), sign flipped
P[, drop_boundary := boundary_exact & mirror_role %in% "partner"]    # p0(partner) == p1(canonical); p1(partner) == p0(canonical)
P[, dedup_drop := drop_mu & drop_boundary]                           # kept for v1 compatibility: all three units dropped
fwrite(M, file.path(OUT, "mirror_pairs.tsv.gz"), sep = "\t"); fwrite(P, file.path(OUT, "promoter_flags.tsv.gz"), sep = "\t")
cat("tissue x gene with two promoters:", nrow(M), "; boundary_exact:", sum(M$boundary_exact), sprintf("(%.1f%%)", 100 * mean(M$boundary_exact)), "; interior_exact:", sum(M$interior_exact), sprintf("(%.1f%%)", 100 * mean(M$interior_exact)), "\n")
cat("tissue x promoter rows:", nrow(P), "; partner mu units dropped:", sum(P$drop_mu), "; partner boundary units dropped:", sum(P$drop_boundary), "\n")
M[, canonical_is_tissue_major := fifelse(canonical == promoter_1, tissue_mean_1 >= tissue_mean_2, tissue_mean_2 >= tissue_mean_1)]
cat("exact pairs whose canonical (pooled major) is also the major promoter in that tissue:", sprintf("%.1f%%", 100 * M[exact_mirror == TRUE, mean(canonical_is_tissue_major)]), "\n")
fwrite(M, file.path(OUT, "mirror_pairs.tsv.gz"), sep = "\t")
