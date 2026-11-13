## Figure 4 panel e, LD-window coloc, sensitivity "eQTL leads by recursive clumping", step s6: signal calls over all of a
## signal's eQTL leads, gene x tissue classes, and the comparison with the main analysis (one eQTL lead = GTEx top variant)
## and with the +-200 kb run.
## Signal call (eGene signals): over the signal x lead pairs that were tested, the highest call in the order
##   H4 > H3 > H1 > H2 > H0 > ambiguous (a puQTL signal that colocalizes with any eQTL signal of the gene is shared);
##   no tested pair and >= 1 pair with disjoint windows -> different locus; otherwise not assessable.
## Non-eGene signals keep the main-analysis call (their pairs are "no eQTL in this tissue" whatever the region).
## Gene x tissue: the rule of 07_merge_classify.R (different locus treated like H3).
## Output: sens_clump/signals_sens.tsv.gz, gene_tissue_sens.tsv.gz, compare_*.tsv; report on stdout.
suppressPackageStartupMessages(library(data.table))
O <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/results/coloc_eqtl_ldwin"; SO <- file.path(O, "sens_clump")
R <- rbindlist(lapply(Sys.glob(file.path(SO, "per_tissue/*.rds")), readRDS), fill = TRUE)
P <- fread(cmd = paste("zcat", file.path(SO, "pairs.tsv.gz")))
stopifnot(nrow(R) == nrow(P))
M <- fread(cmd = paste("zcat", file.path(O, "coloc_eqtl_signals.tsv.gz")))      # main analysis (LD window, GTEx top)
G0 <- fread(cmd = paste("zcat", file.path(O, "gene_tissue_classes.tsv.gz")))
DLS <- "different locus (LD windows disjoint)"
LEV <- c("shared (H4)", "distinct (H3)", DLS, "puQTL-specific (H1)", "eQTL-specific (H2)", "neither (H0)", "ambiguous", "not assessable")
pc <- function(h0, h1, h2, h3, h4) fifelse(h4 >= 0.8, 1L, fifelse(h3 >= 0.8, 2L, fifelse(h1 >= 0.8, 4L, fifelse(h2 >= 0.8, 5L, fifelse(h0 >= 0.8, 6L, 7L)))))
R[status == "tested", code := pc(PP.H0, PP.H1, PP.H2, PP.H3, PP.H4)]
cat("signal x lead pairs:", nrow(R), "; status:\n"); print(R[, .N, by = status][order(-N)])
cat("tested pairs by lead rank (1 = GTEx-top-like) and call:\n")
print(dcast(R[status == "tested", .N, by = .(rank = fifelse(rank >= 3, "3+", as.character(rank)), call = LEV[code])], call ~ rank, value.var = "N", fill = 0))
S <- R[, .(n_leads = .N, n_tested = sum(status == "tested"), n_invalid = sum(status == "invalid (LD windows disjoint)"),
           best = if (any(status == "tested")) min(code, na.rm = TRUE) else NA_integer_,
           best_lead = if (any(status == "tested")) lead[status == "tested"][which.min(code[status == "tested"])] else NA_character_,
           best_rank = if (any(status == "tested")) rank[status == "tested"][which.min(code[status == "tested"])] else NA_integer_),
       by = .(tissue, promoter, type, k)]
S[, call_sens := fifelse(!is.na(best), LEV[best], fifelse(n_invalid > 0 & n_invalid == n_leads, DLS, "not assessable"))]
S <- merge(M[, .(tissue, gene_base, promoter, type, k, egene, call_main = call, call_200kb)], S, by = c("tissue", "promoter", "type", "k"), all.x = TRUE)
S[egene == FALSE, call_sens := call_main]
S[egene == TRUE & is.na(call_sens), call_sens := "not assessable"]          # eGene without a typed significant variant
fwrite(S, file.path(SO, "signals_sens.tsv.gz"), sep = "\t")
cat("\neGene signals, rows = main (GTEx top), columns = sensitivity (all clumped leads):\n")
XS <- dcast(S[egene == TRUE, .N, by = .(call_main, call_sens)], call_main ~ call_sens, value.var = "N", fill = 0); print(XS)
fwrite(XS, file.path(SO, "compare_signal_main_vs_sens.tsv"), sep = "\t")
cat("main H4 -> sens:\n"); print(S[egene == TRUE & call_main == "shared (H4)", .N, by = call_sens])
cat("main different locus -> sens (eGene):\n"); print(S[egene == TRUE & call_main == DLS, .N, by = call_sens][order(-N)])
cat("the 11 signals H4 at +-200 kb and different locus in main -> sens:\n"); print(S[call_200kb == "shared (H4)" & call_main == DLS, .N, by = call_sens])
## gene x tissue
AB <- setNames(c("H4", "H3", "DL", "H1", "H2", "H0", "amb", "NA"), LEV)
S[, cc := AB[call_sens]]; stopifnot(!anyNA(S$cc))
classify <- function(cc) {                                                   # as in 07_merge_classify.R
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
G <- S[, .(egene = egene[1], n_H4 = sum(cc == "H4"), n_H3 = sum(cc == "H3"), n_DL = sum(cc == "DL"), detail = classify(cc)), by = .(tissue, gene_base)]
G[egene == FALSE, detail := "no eQTL in this tissue"]
F4 <- c("puQTL distinct from the eQTL", "No eQTL in this tissue", "puQTL shared with the eQTL", "Partly shared with the eQTL", "Unresolved")
G[, fig4e_sens := factor(fcase(detail == "no eQTL in this tissue", F4[2], startsWith(detail, "distinct"), F4[1], detail == "shared, all signals", F4[3],
                               detail == "shared, rest unresolved", F4[4], default = F4[5]), levels = F4)]
G <- merge(G, G0[, .(tissue, gene_base, fig4e_main = factor(fig4e, levels = F4), detail_main = detail, fig4e_200kb = factor(fig4e_200kb, levels = F4))],
           by = c("tissue", "gene_base"))
stopifnot(nrow(G) == 8283L)
fwrite(G, file.path(SO, "gene_tissue_sens.tsv.gz"), sep = "\t")
CT <- Reduce(function(a, b) merge(a, b, by = "class", all = TRUE), list(
  G[, .(N_200kb = .N), by = .(class = fig4e_200kb)], G[, .(N_main = .N), by = .(class = fig4e_main)], G[, .(N_sens = .N), by = .(class = fig4e_sens)]))
CT <- CT[order(factor(class, levels = F4))]
fwrite(CT, file.path(SO, "compare_fig4e_counts.tsv"), sep = "\t")
cat("\nFig. 4e: +-200 kb / main (LD window, GTEx top) / sensitivity (LD window, all clumped eQTL leads):\n"); print(CT)
cat("detail (sensitivity):\n"); print(G[, .N, by = .(fig4e_sens, detail)][order(fig4e_sens, -N)])
cat("\ngene x tissue, rows = main, columns = sensitivity:\n")
XG <- dcast(G[, .N, by = .(fig4e_main, fig4e_sens)], fig4e_main ~ fig4e_sens, value.var = "N", fill = 0); print(XG)
fwrite(XG, file.path(SO, "compare_fig4e_main_vs_sens.tsv"), sep = "\t")
cat("same class in main and sensitivity:", G[fig4e_main == fig4e_sens, .N], "of", nrow(G), sprintf("(%.1f%%)", 100 * G[, mean(fig4e_main == fig4e_sens)]), "\n")
cat("DONE\n")
