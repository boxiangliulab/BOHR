## run_paper.sh 01.1 and 01.3 (authored for the repository): promoter counts recomputed from the stored proActiv
## annotations and filter outputs, against the published numbers.
## Usage: Rscript verify_01.R filter  <GTEx v11 folder> <Supplementary Tables folder> <status.tsv>
##        Rscript verify_01.R release <GTEx folder>     <Supplementary Tables folder> <status.tsv>
a <- commandArgs(TRUE)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
S1a <- supp(a[3], "S1a"); S1d <- supp(a[3], "S1d")
counts <- function(info) { g <- table(info$geneId)
  c(promoters = nrow(info), genes = length(g), single = sum(g == 1), multi = sum(g > 1), external = sum(!info$internalPromoter)) }
if (a[1] == "filter") {
  v11 <- a[2]
  e <- new.env(); load(file.path(v11, "Promoter_info_all.RData"), envir = e); n <- counts(e$info_all)
  w <- S1a[gtex_release == "v11"]; wd <- S1d[gtex_release == "v11"]
  add("v11 annotation: promoters, genes, multi-promoter genes (Supplementary Table 1a)",
      n[["promoters"]] == w$total_promoters && n[["genes"]] == w$total_genes && n[["multi"]] == w$multi_promoter_genes,
      sprintf("%d promoters, %d genes, %d multi-promoter (published %d, %d, %d)", n[["promoters"]], n[["genes"]], n[["multi"]],
              w$total_promoters, w$total_genes, w$multi_promoter_genes))
  ext <- e$info_all[!e$info_all$internalPromoter, ]
  add("v11 external promoters (Supplementary Table 1d; SI §1: 92,557 in 52,037 genes)",
      nrow(ext) == wd$promoters_external && length(unique(ext$geneId)) == 52037,
      sprintf("%d external promoters in %d genes (published %d in 52,037)", nrow(ext), length(unique(ext$geneId)), wd$promoters_external))
  e2 <- new.env(); load(file.path(v11, "Promoter_info_filter_5.RData"), envir = e2); f <- e2$info
  add("promoters analysed (Methods §1: 10,826 in 4,744 genes)", nrow(f) == 10826 && length(unique(f$geneId)) == 4744,
      sprintf("%d promoters in %d genes", nrow(f), length(unique(f$geneId))))
  add("analysed promoters are external promoters of the annotation", all(rownames(f) %in% rownames(ext)),
      sprintf("%d of %d among the external promoters", sum(rownames(f) %in% rownames(ext)), nrow(f)))
  e3 <- new.env(); load(file.path(v11, "filter_5_recal_att_filter.RData"), envir = e3); X <- e3$recal_att_f
  pc <- grep("^PromoterID_", colnames(X), value = TRUE)          # the other columns are sample, SMTS and SMTSD
  add("relative activity matrix: 22,865 samples x the analysed promoters", nrow(X) == 22865 && setequal(pc, rownames(f)),
      sprintf("%d samples x %d promoter columns (plus %s); same promoters as the filter: %s", nrow(X), length(pc),
              paste(setdiff(colnames(X), pc), collapse = ", "), setequal(pc, rownames(f))))
  rg <- vapply(pc, function(p) { v <- X[[p]]; v <- v[is.finite(v)]; if (length(v)) c(length(v), min(v), max(v)) else c(0, NA, NA) }, numeric(3))
  add("relative activity values within [0, 1]", all(rg[2, ] >= 0 & rg[3, ] <= 1, na.rm = TRUE),
      sprintf("%.0f finite values, range %.3g to %.3g", sum(rg[1, ]), min(rg[2, ], na.rm = TRUE), max(rg[3, ], na.rm = TRUE)))
} else {
  G <- a[2]; vc <- file.path(G, "version_comparison")
  ld <- function(f) { e <- new.env(); load(f, envir = e); e$info_all }
  ## the figure's v6p is the GENCODE v19 catalogue stored as _v7 (GTEx v6p and v7 both use v19; the _v6 file holds the same promoters)
  rel <- list(v6p = ld(file.path(vc, "Promoter_info_all_v7.RData")), v8 = ld(file.path(vc, "Promoter_info_all_v8.RData")),
              v11 = ld(file.path(G, "v11/Promoter_info_all.RData")))
  for (r in names(rel)) {
    n <- counts(rel[[r]]); w <- S1a[gtex_release == r]; wd <- S1d[gtex_release == r]
    ok <- n[["single"]] == w$single_promoter_genes && n[["multi"]] == w$multi_promoter_genes && n[["promoters"]] == w$total_promoters &&
          n[["external"]] == wd$promoters_external
    add(sprintf("Fig. 1A %s: single- and multi-promoter genes, promoters, external promoters", r), ok,
        sprintf("%d single, %d multi, %d promoters, %d external (published %d, %d, %d, %d)", n[["single"]], n[["multi"]], n[["promoters"]],
                n[["external"]], w$single_promoter_genes, w$multi_promoter_genes, w$total_promoters, wd$promoters_external))
  }
  v6 <- ld(file.path(vc, "Promoter_info_all_v6.RData"))
  add("v6 and v7 annotation files hold the same promoters", nrow(v6) == nrow(rel$v6p) && identical(v6$geneId, rel$v6p$geneId) && identical(v6$start, rel$v6p$start),
      sprintf("%d and %d promoters", nrow(v6), nrow(rel$v6p)))
}
finish(a[4])
