## Figure 4 panel e, step 4: reclassify gene x tissue (and promoter x tissue) with the author's rule of 2026-10-07:
##   "At the gene level, a gene counts as fully shared only when all its signals are shared; the supplementary figure can show the breakdown."
## The signal-level coloc calls of step 3 are unchanged (PP >= 0.8, p12 = 1e-5). Only the aggregation changes:
##   step 3 took the highest-ranked call over a gene's signals (H4 > H3 > H1 > H2 > H0 > ambiguous), so one shared signal
##   made the whole gene "shared". Here, over the gene's TESTED signals (all promoters, all types, all kept signals):
##     shared, all signals   every tested signal is H4
##     shared, some signals  >= 1 H4 but not all   (detail: + distinct [>= 1 H3] / + ambiguous [rest ambiguous] / + other)
##     distinct (H3)         >= 1 H3 and no H4     (detail: all distinct / distinct + other)
##     H1 / H2 / H0          no H4 or H3; the highest of H1 > H2 > H0
##     ambiguous             every tested signal ambiguous
##     not assessable        no tested signal
##   puQTL-only / eQTL-only as in step 3. The same rule is applied to promoter x tissue (question 1 of 2026-10-07).
## Second rule for the MAIN figure (author 2026-10-07, later: "any difference counts as distinct"): any H3 signal makes the pair
## distinct, whatever the other signals are. class_main_anyd:
##     distinct (any H3)          >= 1 H3                      (detail: distinct only / + ambiguous or other / + shared)
##     shared, all signals         every tested signal is H4
##     shared, rest unresolved     >= 1 H4, no H3, not all H4   (rest ambiguous, H1, H2 or H0)
##     H1 / H2 / H0 / ambiguous / not assessable as above.
## Read-only on step 3 outputs. Output: results/coloc_eqtl_v3/all_shared_rule/ (+ figure_v3/data/panel_e_allshared.tsv
## and the supplementary-figure tables).
suppressPackageStartupMessages(library(data.table))
V2  <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"; IN <- file.path(V2, "results/coloc_eqtl_v3")
OUT <- file.path(IN, "all_shared_rule"); dir.create(OUT, showWarnings = FALSE)
FIGD <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/figure_v3/data"
S  <- fread(cmd = paste("zcat", file.path(IN, "coloc_eqtl_signals.tsv.gz")))
G  <- fread(cmd = paste("zcat", file.path(IN, "coloc_eqtl_gene_tissue.tsv.gz")))
AB <- c("shared (H4)" = "H4", "distinct (H3)" = "H3", "puQTL-specific (H1)" = "H1", "eQTL-specific (H2)" = "H2",
        "neither (H0)" = "H0", "ambiguous" = "amb", "not assessable" = "NA")
S[, cc := AB[as.character(call)]]
stopifnot(!anyNA(S$cc))
classify <- function(cc) {
  t <- cc[cc != "NA"]
  if (!length(t)) return(c("not assessable", "not assessable"))
  if (all(t == "H4")) return(c("shared, all signals", if (length(t) == 1) "shared, single signal" else "shared, all of several signals"))
  if (any(t == "H4")) return(c("shared, some signals",
    if (any(t == "H3")) "shared + distinct" else if (all(t %in% c("H4", "amb"))) "shared + ambiguous" else "shared + other"))
  if (any(t == "H3")) return(c("distinct (H3)", if (all(t == "H3")) "distinct, all signals" else "distinct + other"))
  for (h in c("H1", "H2", "H0")) if (any(t == h)) {
    lab <- c(H1 = "puQTL-specific (H1)", H2 = "eQTL-specific (H2)", H0 = "neither (H0)")[[h]]; return(c(lab, lab)) }
  c("ambiguous", "ambiguous")
}
cnt <- function(cc, h) sum(cc == h)
summ <- function(D, by) D[, { k <- classify(cc)
  .(n_signals = .N, n_tested = sum(cc != "NA"), n_H4 = cnt(cc, "H4"), n_H3 = cnt(cc, "H3"), n_H1 = cnt(cc, "H1"),
    n_H2 = cnt(cc, "H2"), n_H0 = cnt(cc, "H0"), n_amb = cnt(cc, "amb"), types = paste(sort(unique(type)), collapse = "+"),
    n_promoters = uniqueN(promoter), class_main = k[1], class_detail = k[2]) }, by = by]

## ---- gene x tissue (the universe of step 3: genes tested for puQTL in the tissue) ----
GS <- summ(S[egene == TRUE], c("tissue", "gene_base"))
GT <- merge(G[, .(tissue, gene_base, puqtl, egene, class_old = class)], GS, by = c("tissue", "gene_base"), all.x = TRUE)
GT[puqtl & !egene, `:=`(class_main = "puQTL-only", class_detail = "puQTL-only")]
GT[!puqtl & egene, `:=`(class_main = "eQTL-only", class_detail = "eQTL-only")]
GT <- GT[puqtl | egene]
stopifnot(!anyNA(GT$class_main))
GT[, class_main_anyd := fcase(class_main %in% c("puQTL-only", "eQTL-only"), class_main,
                              n_H3 > 0, "distinct (any H3)",
                              class_main == "shared, all signals", "shared, all signals",
                              class_main == "shared, some signals", "shared, rest unresolved",
                              default = class_main)]
GT[, class_detail_anyd := fcase(class_main_anyd != "distinct (any H3)", class_detail,
                                n_H4 > 0, "distinct + shared", n_H3 == n_tested, "distinct only", default = "distinct + ambiguous or other")]
fwrite(GT, file.path(OUT, "gene_tissue_classes.tsv.gz"), sep = "\t")
MAIN <- c("eQTL-only", "puQTL-only", "shared, all signals", "shared, some signals", "distinct (H3)", "puQTL-specific (H1)",
          "eQTL-specific (H2)", "neither (H0)", "ambiguous", "not assessable")
PM <- GT[, .N, by = .(class = class_main)][, class := factor(class, levels = MAIN)][order(class)]
fwrite(PM, file.path(OUT, "panel_e_counts_allshared.tsv"), sep = "\t"); fwrite(PM, file.path(FIGD, "panel_e_allshared.tsv"), sep = "\t")
DET <- GT[puqtl & egene, .N, by = .(class_main, class_detail)][order(-N)]
fwrite(DET, file.path(OUT, "supp_gene_detail.tsv"), sep = "\t"); fwrite(DET, file.path(FIGD, "supp_coloc_gene_detail.tsv"), sep = "\t")
MAIN2 <- c("puQTL-only", "eQTL-only", "shared, all signals", "shared, rest unresolved", "distinct (any H3)", "puQTL-specific (H1)",
           "eQTL-specific (H2)", "neither (H0)", "ambiguous", "not assessable")
PM2 <- GT[, .N, by = .(class = class_main_anyd)][, class := factor(class, levels = MAIN2)][order(class)]
PM2 <- rbind(PM2, data.table(class = "distinct (any H3): of which also a shared signal", N = GT[class_detail_anyd == "distinct + shared", .N]))
fwrite(PM2, file.path(OUT, "panel_e_counts_anydistinct.tsv"), sep = "\t"); fwrite(PM2, file.path(FIGD, "panel_e_anydistinct.tsv"), sep = "\t")
DET2 <- GT[puqtl & egene, .N, by = .(class_main_anyd, class_detail_anyd)][order(-N)]
fwrite(DET2, file.path(OUT, "supp_gene_detail_anydistinct.tsv"), sep = "\t"); fwrite(DET2, file.path(FIGD, "supp_coloc_gene_detail_anydistinct.tsv"), sep = "\t")
cat("gene x tissue shown in panel e:", nrow(GT), "\n"); print(PM)
cat("\nany-distinct rule (main figure):\n"); print(PM2); print(DET2)
cat("\nold (highest-ranked call) vs new main class, puQTL & eGene:\n"); print(dcast(GT[puqtl & egene, .N, by = .(class_old, class_main)], class_main ~ class_old, value.var = "N", fill = 0))
cat("\ndetail:\n"); print(DET)

## ---- promoter x tissue (question 1: one promoter, several types) ----
PT <- summ(S, c("tissue", "promoter", "gene_base", "egene"))
PT[, class_main_anyd := fcase(n_H3 > 0, "distinct (any H3)", class_main == "shared, all signals", "shared, all signals",
                              class_main == "shared, some signals", "shared, rest unresolved", default = class_main)]
PT[, class_detail_anyd := fcase(class_main_anyd != "distinct (any H3)", class_detail,
                                n_H4 > 0, "distinct + shared", n_H3 == n_tested, "distinct only", default = "distinct + ambiguous or other")]
fwrite(PT, file.path(OUT, "promoter_tissue_classes.tsv.gz"), sep = "\t")
PT[, combo := fcase(types == "mu", "mu", types == "p0", "p0", types == "p1", "p1", types == "mu+p1", "mu+p1", types == "mu+p0", "mu+p0",
                    types == "p0+p1", "p0+p1", default = "mu+p0+p1")]
PC <- PT[egene == TRUE, .N, by = .(combo, class_main, class_detail, class_main_anyd, class_detail_anyd)][order(combo, -N)]
fwrite(PC, file.path(OUT, "supp_promoter_by_types.tsv"), sep = "\t"); fwrite(PC, file.path(FIGD, "supp_coloc_promoter_by_types.tsv"), sep = "\t")
cat("\npromoter x tissue on eGenes:", PT[egene == TRUE, .N], "; with >= 2 types:", PT[egene == TRUE & grepl("\\+", types), .N], "\n")
print(dcast(PT[egene == TRUE & grepl("\\+", types), .N, by = .(combo, class_main)], class_main ~ combo, value.var = "N", fill = 0))
cat("any-distinct rule, promoters with >= 2 types:\n")
print(dcast(PT[egene == TRUE & grepl("\\+", types), .N, by = .(combo, class_main_anyd)], class_main_anyd ~ combo, value.var = "N", fill = 0))
print(PT[egene == TRUE & grepl("\\+", types), .N, by = class_detail_anyd][order(-N)])

## ---- unit level: each type on its own (lead signal k = 1), eGenes ----
UT <- S[k == 1 & egene == TRUE & cc != "NA", .N, by = .(type, cc)]
fwrite(UT, file.path(OUT, "supp_unit_by_type.tsv"), sep = "\t"); fwrite(UT, file.path(FIGD, "supp_coloc_unit_by_type.tsv"), sep = "\t")
cat("\nunit lead-signal calls by type (eGenes):\n"); print(dcast(UT, cc ~ type, value.var = "N", fill = 0))

## ---- cross-type agreement by LD between the mu and p0/p1 lead signals of one promoter ----
L <- S[k == 1 & cc != "NA", .(tissue, promoter, type, cc, signal_snp)]
W <- merge(L[type == "mu"], L[type != "mu"], by = c("tissue", "promoter"), suffixes = c(".mu", ".sw"))
P <- fread(cmd = paste("zcat", file.path(V2, "results/ld_v3/signal_pairs.tsv.gz")))[same_promoter == TRUE & k_a == 1 & k_b == 1]
P <- rbind(P[type_a == "mu" & type_b != "mu", .(tissue, promoter = promoter_a, type.sw = type_b, r2)],
           P[type_b == "mu" & type_a != "mu", .(tissue, promoter = promoter_a, type.sw = type_a, r2)])[, .(r2 = max(r2)), by = .(tissue, promoter, type.sw)]
W <- merge(W, P, by = c("tissue", "promoter", "type.sw"), all.x = TRUE)
W[, ld := fcase(signal_snp.mu == signal_snp.sw | r2 >= 0.8, "same locus (r2 >= 0.8)", r2 >= 0.2, "partial LD (0.2-0.8)", r2 < 0.2, "different loci (r2 < 0.2)",
                default = "LD not available")]
W[, pair := fcase(cc.mu == "H4" & cc.sw == "H4", "both H4", (cc.mu == "H4") != (cc.sw == "H4") & (cc.mu == "H3" | cc.sw == "H3"), "H4 vs H3",
                  (cc.mu == "H4") != (cc.sw == "H4"), "H4 vs other", cc.mu == cc.sw, "same, not H4", default = "different, no H4")]
CT <- W[, .N, by = .(ld, pair)][order(ld, pair)]
fwrite(CT, file.path(OUT, "supp_crosstype_ld.tsv"), sep = "\t"); fwrite(CT, file.path(FIGD, "supp_coloc_crosstype_ld.tsv"), sep = "\t")
cat("\ncross-type pairs:", nrow(W), "\n"); print(dcast(CT, ld ~ pair, value.var = "N", fill = 0))
cat("same call by LD class:\n"); print(W[, .(n = .N, pct_same = round(100 * mean(cc.mu == cc.sw), 1)), by = ld])
cat("DONE\n")
