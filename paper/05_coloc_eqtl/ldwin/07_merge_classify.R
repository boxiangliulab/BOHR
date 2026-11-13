## Figure 4 panel e, LD-window coloc, step 7: merge the per-tissue coloc (step 6), classify, compare with the +-200 kb run.
## Signal call (PP >= 0.8, p12 = 1e-5, precedence H4 > H3 > H1 > H2 > H0, else ambiguous) as in coloc_eqtl_v3/03, plus
##   "different locus (LD windows disjoint)" (DL) for signals whose LD window and the eQTL lead's LD window share < 5 SNPs
##   with statistics in both datasets (author 2026-10-12: not tested, counted under Distinct, count shown separately).
## Gene x tissue (Fig. 4e; the rule of coloc_eqtl_v3/04 class_main_anyd and figure_v3 panel_e_v3e.py), over the pair's
## signals, DL treated like H3 ("any difference counts as different"):
##   distinct      >= 1 H3 or DL, or H1 with no H3 / DL / H4   (shown separately: any H3 / DL, no H3 / H1 only)
##   no eQTL       gene not a GTEx v11 eGene in the tissue (q > 0.05)
##   shared, all   every tested signal H4
##   partly        >= 1 H4, no H3 or DL, the others unresolved (ambiguous / H1 / H2 / H0 / not assessable)
##   unresolved    otherwise (ambiguous, H2, H0, not assessable)
## Read-only on coloc_eqtl_v3. Output: results/coloc_eqtl_ldwin/ (signals, flags for Figure 6, gene x tissue, counts,
## per-tissue table, comparison with v3).
suppressPackageStartupMessages(library(data.table))
V2  <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"; OUT <- file.path(V2, "results/coloc_eqtl_ldwin")
OLD <- file.path(V2, "results/coloc_eqtl_v3")
A <- readRDS(file.path(OLD, "analysis_set.rds"))
R <- rbindlist(lapply(Sys.glob(file.path(OUT, "per_tissue/*.rds")), readRDS), fill = TRUE)
stopifnot(nrow(R) == nrow(A), !anyDuplicated(R[, .(tissue, promoter, type, k)]))
DLS <- "different locus (LD windows disjoint)"
LEV <- c("shared (H4)", "distinct (H3)", DLS, "puQTL-specific (H1)", "eQTL-specific (H2)", "neither (H0)", "ambiguous", "not assessable")
R[, call := fifelse(status == DLS, DLS, fifelse(status != "tested", "not assessable",
  fifelse(PP.H4 >= 0.8, "shared (H4)", fifelse(PP.H3 >= 0.8, "distinct (H3)", fifelse(PP.H1 >= 0.8, "puQTL-specific (H1)",
  fifelse(PP.H2 >= 0.8, "eQTL-specific (H2)", fifelse(PP.H0 >= 0.8, "neither (H0)", "ambiguous")))))))]
R[, call := factor(call, levels = LEV)]
AB <- setNames(c("H4", "H3", "DL", "H1", "H2", "H0", "amb", "NA"), LEV)
R[, cc := AB[as.character(call)]]
front <- c("tissue", "gene", "gene_base", "promoter", "type", "k", "signal_snp", "stat", "rule", "n_signals", "unit_q", "egene", "eqtl_qval",
           "n_eqtl", "status", "n_snps", "n_int", "call", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "PP.H4_p12_5e6")
setcolorder(R, c(front, setdiff(names(R), front)))
setorder(R, tissue, promoter, type, k)

## ---- gene x tissue ----
classify <- function(cc) {
  t <- cc[cc != "NA"]
  if (any(cc == "H3")) return("distinct: any H3")
  if (any(cc == "DL")) return("distinct: different locus (no H3)")
  if (length(t) && all(t == "H4")) return("shared, all signals")
  if (any(t == "H4")) return("shared, rest unresolved")
  if (any(t == "H1")) return("distinct: H1 only")
  if (any(t == "H2")) return("unresolved: H2")
  if (any(t == "H0")) return("unresolved: H0")
  if (any(t == "amb")) return("unresolved: ambiguous")
  "unresolved: not assessable"
}
G <- R[, .(gene = gene[1], egene = egene[1], eqtl_qval = eqtl_qval[1], n_signals = .N, n_tested = sum(status == "tested"),
           n_H4 = sum(cc == "H4"), n_H3 = sum(cc == "H3"), n_DL = sum(cc == "DL"), n_H1 = sum(cc == "H1"), n_H2 = sum(cc == "H2"),
           n_H0 = sum(cc == "H0"), n_amb = sum(cc == "amb"), n_NA = sum(cc == "NA"), detail = classify(cc)), by = .(tissue, gene_base)]
G[egene == FALSE, detail := "no eQTL in this tissue"]
G[, fig4e := fcase(detail == "no eQTL in this tissue", "No eQTL in this tissue", startsWith(detail, "distinct"), "puQTL distinct from the eQTL",
                   detail == "shared, all signals", "puQTL shared with the eQTL", detail == "shared, rest unresolved", "Partly shared with the eQTL",
                   default = "Unresolved")]
F4 <- c("puQTL distinct from the eQTL", "No eQTL in this tissue", "puQTL shared with the eQTL", "Partly shared with the eQTL", "Unresolved")
G[, fig4e := factor(fig4e, levels = F4)]
stopifnot(nrow(G) == 8283L)

## ---- the +-200 kb run (coloc_eqtl_v3), same rule ----
S0 <- fread(cmd = paste("zcat", file.path(OLD, "all_shared_rule/coloc_eqtl_signals_flags.tsv.gz")))
G0 <- fread(cmd = paste("zcat", file.path(OLD, "all_shared_rule/gene_tissue_classes.tsv.gz")))[puqtl == TRUE]
MAP0 <- c("distinct (any H3)" = F4[1], "puQTL-specific (H1)" = F4[1], "puQTL-only" = F4[2], "shared, all signals" = F4[3],
          "shared, rest unresolved" = F4[4], "eQTL-specific (H2)" = F4[5], "neither (H0)" = F4[5], "ambiguous" = F4[5], "not assessable" = F4[5])
G0[, fig4e_200kb := factor(MAP0[class_main_anyd], levels = F4)]
stopifnot(!anyNA(G0$fig4e_200kb), nrow(G0) == 8283L)
G <- merge(G, G0[, .(tissue, gene_base, fig4e_200kb, class_200kb = class_main_anyd)], by = c("tissue", "gene_base"), all.x = TRUE)
stopifnot(!anyNA(G$fig4e_200kb))
fwrite(G, file.path(OUT, "gene_tissue_classes.tsv.gz"), sep = "\t")

## ---- signal tables (same rows and core columns as coloc_eqtl_v3) ----
##   coloc_eqtl_signals.tsv.gz        : the signal results (+ the +-200 kb call and PPs for comparison)
##   coloc_eqtl_signals_flags.tsv.gz  : the same rows and columns + the flags Figure 6 reads for S8d:
##                                      puqtl_specific_H1, different_locus, gene_class_fig4e, gene_class_detail
R <- merge(R, S0[, .(tissue, promoter, type, k, call_200kb = call, PP.H4_200kb = PP.H4, PP.H3_200kb = PP.H3, PP.H1_200kb = PP.H1)],
           by = c("tissue", "promoter", "type", "k"), all.x = TRUE)
setcolorder(R, c(front, setdiff(names(R), front))); setorder(R, tissue, promoter, type, k)
fwrite(R[, !"cc"], file.path(OUT, "coloc_eqtl_signals.tsv.gz"), sep = "\t")
R[, `:=`(puqtl_specific_H1 = call == "puQTL-specific (H1)", different_locus = call == DLS)]
R <- merge(R, G[, .(tissue, gene_base, gene_class_fig4e = as.character(fig4e), gene_class_detail = detail)], by = c("tissue", "gene_base"))
setcolorder(R, c(front, setdiff(names(R), front))); setorder(R, tissue, promoter, type, k)
fwrite(R[, !"cc"], file.path(OUT, "coloc_eqtl_signals_flags.tsv.gz"), sep = "\t")

## ---- counts ----
cnt  <- G[, .N, by = fig4e][order(fig4e)]
cnt0 <- G[, .N, by = .(fig4e = fig4e_200kb)][order(fig4e)]
CT <- merge(cnt0[, .(class = fig4e, N_200kb = N)], cnt[, .(class = fig4e, N_ldwin = N)], by = "class", all = TRUE)[order(factor(class, levels = F4))]
CT[, `:=`(pct_200kb = round(100 * N_200kb / sum(N_200kb), 1), pct_ldwin = round(100 * N_ldwin / sum(N_ldwin), 1))]
fwrite(CT, file.path(OUT, "panel_e_counts.tsv"), sep = "\t")
DET <- G[, .N, by = .(fig4e, detail)][order(fig4e, -N)]
fwrite(DET, file.path(OUT, "panel_e_detail.tsv"), sep = "\t")
key <- c("puQTL distinct from the eQTL" = "independent", "No eQTL in this tissue" = "no_eqtl", "puQTL shared with the eQTL" = "shared_all",
         "Partly shared with the eQTL" = "shared_rest", "Unresolved" = "unresolved")
TT <- dcast(G[, .N, by = .(tissue, k = key[as.character(fig4e)])], tissue ~ k, value.var = "N", fill = 0)
TT <- merge(TT, unique(R[, .(tissue, donors = n_eqtl)]), by = "tissue"); TT[, pairs := independent + no_eqtl + shared_all + shared_rest + unresolved]
setcolorder(TT, c("tissue", "independent", "no_eqtl", "shared_all", "shared_rest", "unresolved", "donors", "pairs"))
fwrite(TT, file.path(OUT, "panel_e_puqtl_by_tissue.tsv"), sep = "\t")
X <- dcast(G[, .N, by = .(fig4e_200kb, fig4e)], fig4e_200kb ~ fig4e, value.var = "N", fill = 0)
fwrite(X, file.path(OUT, "compare_fig4e_200kb_vs_ldwin.tsv"), sep = "\t")
XS <- dcast(R[, .N, by = .(call_200kb, call)], call_200kb ~ call, value.var = "N", fill = 0)
fwrite(XS, file.path(OUT, "compare_signal_calls_200kb_vs_ldwin.tsv"), sep = "\t")

## ---- report ----
cat("signals:", nrow(R), "; status:\n"); print(R[, .N, by = status][order(-N)])
cat("signal calls (all / eGenes):\n"); print(dcast(R[, .N, by = .(call, egene)], call ~ egene, value.var = "N", fill = 0))
cat("\nFig. 4e, +-200 kb vs LD window:\n"); print(CT)
cat("\ndetail:\n"); print(DET)
cat("\ngene x tissue, rows = +-200 kb, columns = LD window:\n"); print(X)
cat("\nsignal calls, rows = +-200 kb, columns = LD window:\n"); print(XS)
cat("\nDL signals: rare eQTL lead", R[call == DLS & e_rare == TRUE, .N], "of", R[call == DLS, .N], "; lead distance (kb) quantiles:\n")
print(round(quantile(R[call == DLS]$lead_dist / 1e3, c(.1, .25, .5, .75, .9))))
cat("DL diagnostic coloc (union; not used): PP >= 0.8 calls\n")
print(R[call == DLS, .N, by = .(diag = fifelse(is.na(diag_PP.H4), "not run (<10 shared)", fifelse(diag_PP.H4 >= 0.8, "H4", fifelse(diag_PP.H3 >= 0.8, "H3",
  fifelse(diag_PP.H1 >= 0.8, "H1", fifelse(diag_PP.H2 >= 0.8, "H2", fifelse(diag_PP.H0 >= 0.8, "H0", "ambiguous")))))))][order(-N)])
cat("distinct pairs with an H4 signal as well:", G[startsWith(detail, "distinct") & n_H4 > 0, .N], "(any H3:", G[detail == "distinct: any H3" & n_H4 > 0, .N],
    "; DL:", G[detail == "distinct: different locus (no H3)" & n_H4 > 0, .N], ")\n")
cat("coloc region (SNPs in the union with both statistics), tested signals: median", R[status == "tested", median(n_snps)], "; v3 median",
    S0[status == "tested", median(n_snps)], "\n")
cat("DONE\n")
