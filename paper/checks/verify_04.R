## run_paper.sh 04.3 (authored for the repository): the stored final characteristic tests and ROKU calls against
## Supplementary Tables 3b, 3g and 3h. Usage: Rscript verify_04.R <tissue_regression folder> <fit/roku_20261007> <Supplementary Tables folder> <status.tsv>
a <- commandArgs(TRUE)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
TR <- fr(file.path(a[1], "trait_regression_20261106_mirror_abs/results/export_B/trait_regression_units.tsv.gz"))
S3b <- supp(a[3], "S3b")
add("significant characteristic tests (4,018)", TR[sig_final == TRUE, .N] == nrow(S3b) && nrow(S3b) == 4018,
    sprintf("%d stored, %d in Supplementary Table 3b", TR[sig_final == TRUE, .N], nrow(S3b)))
k1 <- TR[sig_final == TRUE, paste(tissue, trait, promoter, type)]
k2 <- S3b[, paste(SMTSD, characteristic, promoter, sub("μ", "mu", component, fixed = TRUE))]
add("the same tests as Supplementary Table 3b (tissue, characteristic, promoter, component)", setequal(k1, k2),
    sprintf("%d in both", length(intersect(k1, k2))))
m <- merge(TR[sig_final == TRUE, .(k = paste(tissue, trait, promoter, type), q_final)], data.table(k = k2, q = S3b$q), by = "k")
add("q values of Supplementary Table 3b", nrow(m) == nrow(S3b) && max(abs(m$q_final - m$q) / m$q) < 1e-6,
    sprintf("max relative difference %.3g", max(abs(m$q_final - m$q) / m$q)))
R <- readRDS(file.path(a[2], "v11_tissue_specific_ROKU_results_delta010.rds")); K <- R$true_outliers; K[is.na(K)] <- 0
S3g <- supp(a[3], "S3g"); S3h <- supp(a[3], "S3h")
add("ROKU calls after the effect filter (Supplementary Table 3g)", sum(K != 0) == nrow(S3g),
    sprintf("%d calls stored (%d elevated, %d reduced), %d rows published", sum(K != 0), sum(K == 1), sum(K == -1), nrow(S3g)))
pt <- data.table(SMTSD = colnames(K), elevated = colSums(K == 1), reduced = colSums(K == -1))
h <- merge(pt, S3h[, .(SMTSD, e = elevated, r = reduced)], by = "SMTSD")
add("ROKU calls per tissue (Supplementary Table 3h)", nrow(h) == nrow(S3h) && all(h$elevated == h$e & h$reduced == h$r),
    sprintf("%d of %d tissues equal", sum(h$elevated == h$e & h$reduced == h$r), nrow(S3h)))
finish(a[4])
