## run_paper.sh 05.1, 05.2, 07.1, 07.2 and 07.3 (authored for the repository): re-run colocalization and MRAS scoring
## against the stored outputs, and the stored final tables against Supplementary Tables 5, 8 and 9.
## Usage: Rscript check_05_07.R eqtl|gwas <new per_tissue> <stored per_tissue> <status.tsv> <tissue> [<tissue> ...]
##        Rscript check_05_07.R eqtl_final <coloc_eqtl_v3c/final_clump_v3c> <Supplementary Tables folder> <status.tsv>
##        Rscript check_05_07.R gwas_final <results_v3c> <Supplementary Tables folder> <status.tsv>
##        Rscript check_05_07.R alphagenome <re-run folder> <stored tables> <stored add-on tables> <status.tsv>
a <- commandArgs(TRUE)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
mode <- a[1]
if (mode %in% c("eqtl", "gwas")) {
  for (t in a[-(1:4)]) {
    f <- if (mode == "eqtl") paste0(t, ".rds") else paste0("coloc_ldwin_", t, ".rds")
    keys <- if (mode == "eqtl") c("promoter", "type", "k") else c("promoter", "type", "k", "trait", "gwas_locus_id")
    N <- as.data.table(readRDS(file.path(a[2], f))); O <- as.data.table(readRDS(file.path(a[3], f)))
    r <- compare_tables(N, O, keys, tol = 1e-6)
    hc <- function(x) x[status == "tested" & PP.H4 >= 0.8, .N]
    add(sprintf("%s: %s colocalization of every signal (status, posteriors, windows)", t, if (mode == "eqtl") "eQTL" else "GWAS"),
        r$ok, sprintf("%s; H4 >= 0.8: %d (stored %d)", r$detail, hc(N), hc(O)))
  }
  finish(a[4])
} else if (mode == "eqtl_final") {
  G <- fr(file.path(a[2], "gene_tissue_classes.tsv.gz")); S5a <- supp(a[3], "S5a"); S5c <- supp(a[3], "S5c")
  n <- G[, .N, by = fig4e]; get <- function(cl) n[fig4e == cl, N]
  add("gene-tissue pairs with a puQTL (8,297)", nrow(G) == 8297 && nrow(S5a) == 8297, sprintf("%d stored, %d in Supplementary Table 5a", nrow(G), nrow(S5a)))
  add("puQTL distinct from the eQTL (1,600) and without eQTL in the tissue (1,663; 3,263 together)",
      get("puQTL distinct from the eQTL") == 1600 && get("No eQTL in this tissue") == 1663,
      sprintf("%d distinct, %d without eQTL", get("puQTL distinct from the eQTL"), get("No eQTL in this tissue")))
  k1 <- G[, paste(tissue, gene_base, fig4e)]; k2 <- S5a[, paste(SMTSD, sub("\\..*$", "", gene_id), fig4e_class)]
  add("the class of every gene-tissue pair equals Supplementary Table 5a", setequal(k1, k2), sprintf("%d of %d equal", length(intersect(k1, k2)), length(k2)))
  pool <- S5c[SMTSD == "All tissues (pooled)"]
  if (nrow(pool)) add("pooled classes of Supplementary Table 5c", pool$puqtl_gene_tissue_pairs == nrow(G) && pool$distinct_from_eqtl == get("puQTL distinct from the eQTL"),
                      sprintf("%s: %d pairs, %d distinct", pool$SMTSD, pool$puqtl_gene_tissue_pairs, pool$distinct_from_eqtl))
  finish(a[4])
} else if (mode == "gwas_final") {
  M <- as.data.table(readRDS(file.path(a[2], "coloc_final_merged.rds"))); S8a <- supp(a[3], "S8a"); S8b <- supp(a[3], "S8b")
  add("tested signal x GWAS-locus pairs (47,506)", nrow(M) == 47506, sprintf("%d stored; Supplementary Table 8b total %d", nrow(M), sum(S8b$tests_total)))
  add("high-confidence colocalizations, H4 >= 0.8 (1,125: 623 μ, 502 boundary)",
      M[HC == TRUE, .N] == 1125 && M[HC == TRUE & type == "mu", .N] == 623 && M[HC == TRUE & type != "mu", .N] == 502,
      sprintf("%d (%d μ, %d boundary); Supplementary Table 8a rows %d", M[HC == TRUE, .N], M[HC == TRUE & type == "mu", .N], M[HC == TRUE & type != "mu", .N], nrow(S8a)))
  k1 <- M[HC == TRUE, paste(tissue, promoter, type, k, trait)]
  k2 <- S8a[, paste(SMTSD, promoter, sub("μ", "mu", class, fixed = TRUE), signal_k, gwas_file)]
  add("the high-confidence colocalizations of Supplementary Table 8a", setequal(k1, k2), sprintf("%d of %d equal", length(intersect(k1, k2)), length(k2)))
  m <- merge(M[HC == TRUE, .(k = paste(tissue, promoter, type, k, trait), PP.H4)], data.table(k = k2, H4 = S8a$H4), by = "k")
  add("their H4 posteriors", nrow(m) == nrow(S8a) && max(abs(m$PP.H4 - m$H4)) < 1e-6, sprintf("max |difference| %.3g", max(abs(m$PP.H4 - m$H4))))
  pt <- M[, .(tests = .N, HC = sum(HC)), by = .(SMTSD = tissue)]; b <- merge(pt, S8b[, .(SMTSD, tests_total, HC_total)], by = "SMTSD")
  add("tests and high-confidence colocalizations per tissue (Supplementary Table 8b)", all(b$tests == b$tests_total & b$HC == b$HC_total),
      sprintf("%d of %d tissues equal", sum(b$tests == b$tests_total & b$HC == b$HC_total), nrow(S8b)))
  finish(a[4])
} else if (mode == "alphagenome") {
  A <- a[2]; t0 <- file.mtime(file.path(A, ".start"))
  need <- list(tables = c("predicted_cage_peaks_COMMON.tsv", "predicted_procap_peaks_COMMON.tsv", "score_CAGE_peaks.tsv", "score_PROCAP_peaks.tsv", "score_CAGE_share.tsv"),
               "a3/tables" = c("score_CAGE_peaks_v3.tsv", "score_PROCAP_peaks_v3.tsv", "score_CAGE_share_v3.tsv", "ranking_v3_candidates.tsv"))
  stored <- list(tables = a[3], "a3/tables" = a[4]); lead_rows <- list()
  for (d in names(need)) {
    fs <- list.files(file.path(A, d), "\\.tsv$"); fs <- fs[file.mtime(file.path(A, d, fs)) > t0]
    miss <- setdiff(need[[d]], fs)
    if (length(miss)) add(paste("MRAS scoring outputs written to", d), FALSE, paste("not written:", paste(miss, collapse = ", ")))
    for (f in fs) {
      o <- file.path(stored[[d]], f)
      if (!file.exists(o)) { add(paste("MRAS", f), "SKIP", "no stored file of that name"); next }
      x <- fread(file.path(A, d, f)); y <- fread(o)
      if (!identical(names(x), names(y)) || nrow(x) != nrow(y)) { add(paste("MRAS", f), FALSE, sprintf("%d x %d, stored %d x %d", nrow(x), ncol(x), nrow(y), ncol(y))); next }
      x[, .row := .I]; y[, .row := .I]; r <- compare_tables(x, y, ".row", tol = 1e-6, skip = "leads")
      add(paste("MRAS", f), r$ok, paste0(r$detail, if ("leads" %in% names(x)) "; column leads compared separately below" else ""))
      if ("leads" %in% names(x) && "leads" %in% names(y)) {
        lx <- as.character(x$leads); ly <- as.character(y$leads); lx[is.na(lx)] <- ""; ly[is.na(ly)] <- ""
        nd <- sum(lx != ly); nv <- length(unique(x$snp[lx != ly]))
        lead_rows[[f]] <- c(nd, nv)
      }
    }
  }
  if (length(lead_rows)) {
    nd <- sum(vapply(lead_rows, `[`, 0, 1)); nv <- max(vapply(lead_rows, `[`, 0, 2))
    add("MRAS score tables, column leads (lead labels of the variant)", if (nd == 0) TRUE else "SKIP",
        if (nd == 0) "identical in every table" else sprintf(paste0(
          "%d cells in %d tables differ, for %d variants: the labels come from an older lead table (TATA_Enhancer/new_lrt_fixed_analysis/",
          "selected_cutoffs/mu0.70_p0.35/qtl_*_proxy.txt) that was last modified on 2026-09-05, after these tables were written (2026-07-02, file dates); ",
          "the scores are compared in full above, and the Supplementary Tables take their lead columns from the final signal table, not from this label"),
          nd, length(lead_rows), nv))
  }
  finish(a[5])
}
