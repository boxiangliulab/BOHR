## Figure 4 panel e, LD-window coloc, step 8 (FINAL; author 2026-10-12: "OK, use the clump lead throughout, switch"; "no overlap, the positions are far apart,
## it should not be one causal variant, recording it as distinct is fine with me"). eQTL leads come from recursive clumping, as the GWAS leads of Figure 6.
##   eGene signals: paired with every clumped eQTL lead of the gene (sens_clump/ s1-s4, s5 v2). The signal takes the
##     highest-ranked call over its tested pairs (H4 > H3 > H1 > H2 > H0 > ambiguous). Its reported pair is the pair of
##     that call; ties go to the strongest eQTL lead (lowest rank). With no tested pair the signal is a different locus
##     when every pair has disjoint windows, otherwise not assessable; the reported pair is then the rank-1 lead.
##   non-eGene signals: the gene has no significant eQTL variant, so its one lead is the GTEx top variant (the LD-window
##     run of step 7, coloc_eqtl_signals.tsv.gz).
## Gene x tissue: the rule of 07_merge_classify.R; a different locus counts as distinct.
## Output: results/coloc_eqtl_ldwin/final_clump/ (signals, flags for Figure 6 S8d, gene x tissue, Fig. 4e tables,
## comparisons with the +-200 kb run and with the GTEx-top-only LD-window run).
suppressPackageStartupMessages(library(data.table))
O  <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/results/coloc_eqtl_ldwin"
SO <- file.path(O, "sens_clump"); FO <- file.path(O, "final_clump"); dir.create(FO, showWarnings = FALSE)
DLS <- "different locus (LD windows disjoint)"
LEV <- c("shared (H4)", "distinct (H3)", DLS, "puQTL-specific (H1)", "eQTL-specific (H2)", "neither (H0)", "ambiguous", "not assessable")
M <- fread(cmd = paste("zcat", file.path(O, "coloc_eqtl_signals.tsv.gz")))       # step 7: all signals, GTEx top variant as lead
stopifnot(nrow(M) == 14783L)
R <- rbindlist(lapply(Sys.glob(file.path(SO, "per_tissue_v2/*.rds")), readRDS), fill = TRUE)
P <- fread(cmd = paste("zcat", file.path(SO, "pairs.tsv.gz")))
stopifnot(nrow(R) == nrow(P), uniqueN(R$tissue) == 50L)
R <- merge(R, P[, .(tissue, promoter, type, k, lead, lead_p, n_absorbed, n_partners, p_start, p_end, l_start, l_end, int_start, int_end, uni_start, uni_end)],
           by = c("tissue", "promoter", "type", "k", "lead"))
stopifnot(nrow(R) == nrow(P))
PAN <- fread(cmd = paste("zcat", file.path(O, "panel/panel_snps.tsv.gz")), header = FALSE, select = 2)$V2
R[, lead_rare := !(lead %in% PAN)]
pc <- function(h0, h1, h2, h3, h4) fifelse(h4 >= 0.8, 1L, fifelse(h3 >= 0.8, 2L, fifelse(h1 >= 0.8, 4L, fifelse(h2 >= 0.8, 5L, fifelse(h0 >= 0.8, 6L, 7L)))))
R[, code := fifelse(status == "tested", pc(PP.H0, PP.H1, PP.H2, PP.H3, PP.H4), NA_integer_)]
setorder(R, tissue, promoter, type, k, code, rank, na.last = TRUE)                # best call first; ties -> strongest lead
N <- R[, .(n_eqtl_leads = .N, n_pairs_tested = sum(status == "tested"), n_pairs_disjoint = sum(status == "invalid (LD windows disjoint)")),
       by = .(tissue, promoter, type, k)]
B <- R[, .SD[if (any(status == "tested")) 1L else which.min(rank)], by = .(tissue, promoter, type, k)]
B <- merge(B, N, by = c("tissue", "promoter", "type", "k"))
B[, status_f := fifelse(n_pairs_tested > 0, "tested", fifelse(n_pairs_disjoint == n_eqtl_leads, DLS,
                 fifelse(status == "not assessable (<10 shared SNPs)", "not assessable (<10 shared SNPs)", "no eQTL rows")))]
## eGene rows: metadata from step 7, results from the reported pair
keep_meta <- c("tissue", "gene", "gene_base", "promoter", "type", "k", "signal_snp", "stat", "rule", "n_signals", "unit_q", "sig_unit", "egene",
               "eqtl_qval", "n_eqtl", "eqtl_top", "call_200kb", "PP.H4_200kb", "PP.H3_200kb", "PP.H1_200kb")
E <- merge(M[egene == TRUE, ..keep_meta], B[, .(tissue, promoter, type, k, status = status_f, n_snps, n_int, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4, PP.H4_p12_5e6,
                                               pqtl_min_p, eqtl_min_p, top_snp_H4, top_snp_PP, eqtl_lead, eqtl_coloc_lead = lead, eqtl_coloc_lead_rank = rank,
                                               eqtl_coloc_lead_p = lead_p, e_rare = lead_rare, lead_dist, n_eqtl_leads, n_pairs_tested, n_pairs_disjoint,
                                               p_start, p_end, e_start = l_start, e_end = l_end, int_start, int_end, uni_start, uni_end)],
           by = c("tissue", "promoter", "type", "k"))
stopifnot(nrow(E) == M[egene == TRUE, .N])
## non-eGene rows: the GTEx top variant is the gene's only lead (step 7 results)
X <- M[egene == FALSE, c(keep_meta, "status", "n_snps", "n_int", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "PP.H4_p12_5e6", "pqtl_min_p", "eqtl_min_p",
                         "top_snp_H4", "top_snp_PP", "eqtl_lead", "e_rare", "lead_dist", "p_start", "p_end", "e_start", "e_end", "int_start", "int_end",
                         "uni_start", "uni_end"), with = FALSE]
X[, `:=`(eqtl_coloc_lead = eqtl_top, eqtl_coloc_lead_rank = 1L, eqtl_coloc_lead_p = NA_real_, n_eqtl_leads = 1L,
         n_pairs_tested = as.integer(status == "tested"), n_pairs_disjoint = as.integer(status == DLS))]
F <- rbind(E, X, use.names = TRUE)
stopifnot(nrow(F) == 14783L, !anyDuplicated(F[, .(tissue, promoter, type, k)]))
F[, call := fifelse(status == DLS, DLS, fifelse(status != "tested", "not assessable",
  fifelse(PP.H4 >= 0.8, "shared (H4)", fifelse(PP.H3 >= 0.8, "distinct (H3)", fifelse(PP.H1 >= 0.8, "puQTL-specific (H1)",
  fifelse(PP.H2 >= 0.8, "eQTL-specific (H2)", fifelse(PP.H0 >= 0.8, "neither (H0)", "ambiguous")))))))]
F[, call := factor(call, levels = LEV)]
AB <- setNames(c("H4", "H3", "DL", "H1", "H2", "H0", "amb", "NA"), LEV)
F[, cc := AB[as.character(call)]]
## agreement with the sensitivity merge of s6 (same rule; s5 v1 and v2 differ only in the columns kept)
S6 <- fread(cmd = paste("zcat", file.path(SO, "signals_sens.tsv.gz")))
chk <- merge(F[, .(tissue, promoter, type, k, call = as.character(call))], S6[, .(tissue, promoter, type, k, call_sens)], by = c("tissue", "promoter", "type", "k"))
cat("signal calls equal to the s6 sensitivity merge:", chk[call == call_sens, .N], "of", nrow(chk), "\n")
## gene x tissue
classify <- function(cc) {                                                     # as in 07_merge_classify.R
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
G <- F[, .(gene = gene[1], egene = egene[1], eqtl_qval = eqtl_qval[1], eqtl_top = eqtl_top[1], n_eqtl_leads = max(n_eqtl_leads), n_signals = .N,
           n_tested = sum(status == "tested"), n_H4 = sum(cc == "H4"), n_H3 = sum(cc == "H3"), n_DL = sum(cc == "DL"), n_H1 = sum(cc == "H1"),
           n_H2 = sum(cc == "H2"), n_H0 = sum(cc == "H0"), n_amb = sum(cc == "amb"), n_NA = sum(cc == "NA"), detail = classify(cc)), by = .(tissue, gene_base)]
G[egene == FALSE, detail := "no eQTL in this tissue"]
F4 <- c("puQTL distinct from the eQTL", "No eQTL in this tissue", "puQTL shared with the eQTL", "Partly shared with the eQTL", "Unresolved")
G[, fig4e := factor(fcase(detail == "no eQTL in this tissue", F4[2], startsWith(detail, "distinct"), F4[1], detail == "shared, all signals", F4[3],
                          detail == "shared, rest unresolved", F4[4], default = F4[5]), levels = F4)]
G0 <- fread(cmd = paste("zcat", file.path(O, "gene_tissue_classes.tsv.gz")))
G <- merge(G, G0[, .(tissue, gene_base, fig4e_ldwin_top = factor(fig4e, levels = F4), fig4e_200kb = factor(fig4e_200kb, levels = F4))], by = c("tissue", "gene_base"))
stopifnot(nrow(G) == 8283L)
fwrite(G, file.path(FO, "gene_tissue_classes.tsv.gz"), sep = "\t")
## signal tables
front <- c("tissue", "gene", "gene_base", "promoter", "type", "k", "signal_snp", "stat", "rule", "n_signals", "unit_q", "egene", "eqtl_qval", "n_eqtl",
           "status", "n_snps", "n_int", "call", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "PP.H4_p12_5e6")
F <- merge(F, M[, .(tissue, promoter, type, k, call_ldwin_top = call)], by = c("tissue", "promoter", "type", "k"))
setcolorder(F, c(front, setdiff(names(F), front))); setorder(F, tissue, promoter, type, k)
fwrite(F[, !"cc"], file.path(FO, "coloc_eqtl_signals.tsv.gz"), sep = "\t")
FL <- merge(F, G[, .(tissue, gene_base, gene_class_fig4e = as.character(fig4e), gene_class_detail = detail)], by = c("tissue", "gene_base"))
FL[, `:=`(puqtl_specific_H1 = call == "puQTL-specific (H1)", different_locus = call == DLS)]
setcolorder(FL, c(front, setdiff(names(FL), front))); setorder(FL, tissue, promoter, type, k)
fwrite(FL[, !"cc"], file.path(FO, "coloc_eqtl_signals_flags.tsv.gz"), sep = "\t")
## Fig. 4e tables
CT <- Reduce(function(a, b) merge(a, b, by = "class", all = TRUE), list(G[, .(N_200kb = .N), by = .(class = fig4e_200kb)],
  G[, .(N_ldwin_top = .N), by = .(class = fig4e_ldwin_top)], G[, .(N_final_clump = .N), by = .(class = fig4e)]))
CT <- CT[order(factor(class, levels = F4))]; CT[, pct_final_clump := round(100 * N_final_clump / sum(N_final_clump), 1)]
fwrite(CT, file.path(FO, "panel_e_counts.tsv"), sep = "\t")
DET <- G[, .N, by = .(fig4e, detail)][order(fig4e, -N)]; fwrite(DET, file.path(FO, "panel_e_detail.tsv"), sep = "\t")
key <- c("puQTL distinct from the eQTL" = "independent", "No eQTL in this tissue" = "no_eqtl", "puQTL shared with the eQTL" = "shared_all",
         "Partly shared with the eQTL" = "shared_rest", "Unresolved" = "unresolved")
TT <- dcast(G[, .N, by = .(tissue, k = key[as.character(fig4e)])], tissue ~ k, value.var = "N", fill = 0)
TT <- merge(TT, unique(F[, .(tissue, donors = n_eqtl)]), by = "tissue"); TT[, pairs := independent + no_eqtl + shared_all + shared_rest + unresolved]
setcolorder(TT, c("tissue", "independent", "no_eqtl", "shared_all", "shared_rest", "unresolved", "donors", "pairs"))
fwrite(TT, file.path(FO, "panel_e_puqtl_by_tissue.tsv"), sep = "\t")
fwrite(dcast(G[, .N, by = .(fig4e_ldwin_top, fig4e)], fig4e_ldwin_top ~ fig4e, value.var = "N", fill = 0), file.path(FO, "compare_fig4e_ldwin_top_vs_final.tsv"), sep = "\t")
fwrite(dcast(G[, .N, by = .(fig4e_200kb, fig4e)], fig4e_200kb ~ fig4e, value.var = "N", fill = 0), file.path(FO, "compare_fig4e_200kb_vs_final.tsv"), sep = "\t")
fwrite(dcast(F[, .N, by = .(call_200kb, call)], call_200kb ~ call, value.var = "N", fill = 0), file.path(FO, "compare_signal_calls_200kb_vs_final.tsv"), sep = "\t")
fwrite(dcast(F[, .N, by = .(call_ldwin_top, call)], call_ldwin_top ~ call, value.var = "N", fill = 0), file.path(FO, "compare_signal_calls_ldwin_top_vs_final.tsv"), sep = "\t")
## report
cat("signals:", nrow(F), "; status:\n"); print(F[, .N, by = status][order(-N)])
cat("signal calls (all / eGenes):\n"); print(dcast(F[, .N, by = .(call, egene)], call ~ egene, value.var = "N", fill = 0))
cat("\nFig. 4e:\n"); print(CT); print(DET)
cat("distinct pairs that also have an H4 signal:", G[startsWith(detail, "distinct") & n_H4 > 0, .N], "(any H3:", G[detail == "distinct: any H3" & n_H4 > 0, .N],
    "; different locus:", G[detail == "distinct: different locus (no H3)" & n_H4 > 0, .N], ")\n")
cat("unresolved H2 / H0 / not assessable:", G[detail == "unresolved: H2", .N], "/", G[detail == "unresolved: H0", .N], "/", G[detail == "unresolved: not assessable", .N], "\n")
cat("eGene signals: reported eQTL lead rank 1 for", F[egene == TRUE & eqtl_coloc_lead_rank == 1, .N], "of", F[egene == TRUE, .N],
    "; reported lead absent from the panel (MAF < 0.05):", F[egene == TRUE & e_rare == TRUE, .N], "\n")
cat("puqtl_specific_H1 (all / eGenes):", FL[puqtl_specific_H1 == TRUE, .N], "/", FL[puqtl_specific_H1 == TRUE & egene == TRUE, .N], "\n")
cat("same Fig. 4e class as the GTEx-top LD-window run:", G[fig4e == fig4e_ldwin_top, .N], "of", nrow(G), "\n")
cat("DONE\n")
