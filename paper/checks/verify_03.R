## run_paper.sh 03.3 (authored for the repository): the stored units (v3c) and independent signals of all 50 tissues
## against Supplementary Tables 4a, 4d and 5d.
## Usage: Rscript verify_03.R <run root (rerun_v2)> <Supplementary Tables folder> <status.tsv>
a <- commandArgs(TRUE)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
V <- file.path(a[1], "results")
S4a <- supp(a[2], "S4a"); S4d <- supp(a[2], "S4d"); S5d <- supp(a[2], "S5d")
fs <- list.files(file.path(V, "final_nogate_v3c/freeze"), "\\.real\\.units\\.rds$", full.names = TRUE)
U <- rbindlist(lapply(fs, function(f) as.data.table(readRDS(f))[, .(tissue, promoter, type, gene, q)]))
U[, sig := !is.na(q) & q <= 0.05]
add("unit tables of the 50 tissues (final_nogate_v3c)", length(fs) == 50 && setequal(unique(U$tissue), S4a$SMTSD),
    sprintf("%d files, %d tissues, the tissues of Supplementary Table 4a: %s", length(fs), uniqueN(U$tissue), setequal(unique(U$tissue), S4a$SMTSD)))
A <- U[, .(units_tested = sum(!is.na(q)), mu_t = sum(!is.na(q) & type == "mu"), b_t = sum(!is.na(q) & type != "mu"),
           genes_tested = uniqueN(gene[!is.na(q)]), units_significant = sum(sig), mu_s = sum(sig & type == "mu"),
           b_s = sum(sig & type != "mu"), genes_significant = uniqueN(gene[sig])), by = .(SMTSD = tissue)]
M <- merge(A, S4a, by = "SMTSD", suffixes = c("", ".pub"))
cmp <- list("units tested" = c("units_tested", "units_tested.pub"), "μ units tested" = c("mu_t", "units_tested_μ"),
            "boundary units tested" = c("b_t", "units_tested_boundary"), "genes tested" = c("genes_tested", "genes_tested.pub"),
            "significant units" = c("units_significant", "units_significant.pub"), "significant μ units" = c("mu_s", "units_significant_μ"),
            "significant boundary units" = c("b_s", "units_significant_boundary"), "genes with a puQTL" = c("genes_significant", "genes_significant.pub"))
for (k in names(cmp)) { x <- M[[cmp[[k]][1]]]; y <- M[[cmp[[k]][2]]]
  add(paste0("Supplementary Table 4a, ", k, " per tissue"), nrow(M) == 50 && all(x == y),
      sprintf("%d of 50 tissues equal; total %d (published %d)", sum(x == y), sum(x), sum(y))) }
add("significant units, all tissues (13,811)", sum(U$sig) == 13811, sprintf("%d", sum(U$sig)))
SG <- fr(file.path(V, "conditional_final_v3c/freeze_real_signals.tsv.gz"))[backward_keep == TRUE]
add("independent signals kept after the backward step (14,800)", nrow(SG) == 14800, sprintf("%d", nrow(SG)))
k1 <- SG[, paste(tissue, promoter, type, k, snp)]
k2 <- S5d[, paste(SMTSD, promoter, sub("μ", "mu", class, fixed = TRUE), k, signal_variant)]
add("signals equal to Supplementary Table 5d (tissue, promoter, component, signal, variant)", setequal(k1, k2),
    sprintf("%d stored, %d published, %d in both", length(k1), length(k2), length(intersect(k1, k2))))
cov <- merge(U[sig == TRUE, .(tissue, promoter, type)], unique(SG[, .(tissue, promoter, type)]))
add("every significant unit has at least one signal", nrow(cov) == sum(U$sig), sprintf("%d of %d", nrow(cov), sum(U$sig)))
pt <- SG[, .N, by = .(SMTSD = tissue)]; d4 <- merge(pt, S4d[class == "any class", .(SMTSD, independent_signals)], by = "SMTSD")
add("Supplementary Table 4d, independent signals per tissue", nrow(d4) == 50 && all(d4$N == d4$independent_signals),
    sprintf("%d of 50 tissues equal", sum(d4$N == d4$independent_signals)))
finish(a[3])
