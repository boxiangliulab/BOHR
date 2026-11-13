## Mirror rule v3 coloc, step 6: the final per-signal eQTL coloc table re-keyed to units v3m, with the rules of
## coloc_eqtl_ldwin/08_final_clump.R (unchanged):
##   carried rows: the final_clump row as it is (same unit, signal, scan and windows); unit_q updated to the v3m BH value
##   new rows    : eGene -> best call over its clumped eQTL-lead pairs (H4 > H3 > H1 > H2 > H0 > ambiguous; ties -> strongest
##                 lead); no tested pair -> different locus if every pair is disjoint, else not assessable / no eQTL rows.
##                 non-eGene -> the GTEx top-variant run. call_200kb is NA for new rows (no +-200 kb run).
## Gene x tissue classes and the Fig. 4e tables as in step 8. Comparisons with final_clump (mirror v2).
## Output: results/coloc_eqtl_v3m/final_clump_v3m/
suppressPackageStartupMessages(library(data.table))
R2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"; O <- file.path(R2, "results/coloc_eqtl_v3m"); LW <- file.path(R2, "results/coloc_eqtl_ldwin")
FO <- file.path(O, "final_clump_v3m"); dir.create(FO, showWarnings = FALSE)
DLS <- "different locus (LD windows disjoint)"
LEV <- c("shared (H4)", "distinct (H3)", DLS, "puQTL-specific (H1)", "eQTL-specific (H2)", "neither (H0)", "ambiguous", "not assessable")
z <- function(f) fread(cmd = paste("zcat", f))
A <- readRDS(file.path(O, "analysis_set.rds"))
F0 <- z(file.path(LW, "final_clump/coloc_eqtl_signals.tsv.gz"))
## carried rows
C <- merge(A[carried == TRUE, .(tissue, promoter, type, k, signal_snp, q_v3m = q)], F0, by = c("tissue", "promoter", "type", "k", "signal_snp"))
stopifnot(nrow(C) == A[carried == TRUE, .N])
C[, unit_q := q_v3m][, q_v3m := NULL][, source := "carried (final_clump, mirror v2)"]
## new rows: GTEx-top run (step-7 analogue)
RT <- rbindlist(lapply(Sys.glob(file.path(O, "per_tissue_top/*.rds")), readRDS), fill = TRUE)
N <- A[carried == FALSE]; stopifnot(nrow(RT) == nrow(N))
callf <- function(status, h0, h1, h2, h3, h4) fifelse(status == DLS, DLS, fifelse(status != "tested", "not assessable",
  fifelse(h4 >= 0.8, "shared (H4)", fifelse(h3 >= 0.8, "distinct (H3)", fifelse(h1 >= 0.8, "puQTL-specific (H1)",
  fifelse(h2 >= 0.8, "eQTL-specific (H2)", fifelse(h0 >= 0.8, "neither (H0)", "ambiguous")))))))
RT[, call_ldwin_top := callf(status, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4)]
meta <- N[, .(tissue, promoter, type, k, gene, gene_base, signal_snp, stat, rule, n_signals, unit_q = q, sig_unit = sig, egene, eqtl_qval, n_eqtl)]
TOPm <- merge(meta, RT[, .(tissue, promoter, type, k, eqtl_top, status, n_snps, n_int, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4, PP.H4_p12_5e6, pqtl_min_p, eqtl_min_p,
                           top_snp_H4, top_snp_PP, eqtl_lead, e_rare, lead_dist, p_start, p_end, e_start, e_end, int_start, int_end, uni_start, uni_end, call_ldwin_top)],
              by = c("tissue", "promoter", "type", "k"))
## new eGene rows: clumped leads (step-8 analogue)
RC <- rbindlist(lapply(Sys.glob(file.path(O, "per_tissue_clump/*.rds")), readRDS), fill = TRUE)
P <- z(file.path(O, "pairs_clump_new_all.tsv.gz")); stopifnot(nrow(RC) == nrow(P))
RC <- merge(RC, P[, .(tissue, promoter, type, k, lead, lead_p, n_absorbed, n_partners, p_start, p_end, l_start, l_end, int_start, int_end, uni_start, uni_end)],
            by = c("tissue", "promoter", "type", "k", "lead")); stopifnot(nrow(RC) == nrow(P))
PAN <- fread(cmd = paste("zcat", file.path(LW, "panel/panel_snps.tsv.gz")), header = FALSE, select = 2)$V2
RC[, lead_rare := !(lead %in% PAN)]
pc <- function(h0, h1, h2, h3, h4) fifelse(h4 >= 0.8, 1L, fifelse(h3 >= 0.8, 2L, fifelse(h1 >= 0.8, 4L, fifelse(h2 >= 0.8, 5L, fifelse(h0 >= 0.8, 6L, 7L)))))
RC[, code := fifelse(status == "tested", pc(PP.H0, PP.H1, PP.H2, PP.H3, PP.H4), NA_integer_)]
setorder(RC, tissue, promoter, type, k, code, rank, na.last = TRUE)
NN <- RC[, .(n_eqtl_leads = .N, n_pairs_tested = sum(status == "tested"), n_pairs_disjoint = sum(status == "invalid (LD windows disjoint)")), by = .(tissue, promoter, type, k)]
B <- RC[, .SD[if (any(status == "tested")) 1L else which.min(rank)], by = .(tissue, promoter, type, k)]
B <- merge(B, NN, by = c("tissue", "promoter", "type", "k"))
B[, status_f := fifelse(n_pairs_tested > 0, "tested", fifelse(n_pairs_disjoint == n_eqtl_leads, DLS,
                 fifelse(status == "not assessable (<10 shared SNPs)", "not assessable (<10 shared SNPs)", "no eQTL rows")))]
km <- c(names(meta), "eqtl_top", "call_ldwin_top")
E <- merge(TOPm[egene == TRUE, ..km], B[, .(tissue, promoter, type, k, status = status_f, n_snps, n_int, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4, PP.H4_p12_5e6,
                                            pqtl_min_p, eqtl_min_p, top_snp_H4, top_snp_PP, eqtl_lead, eqtl_coloc_lead = lead, eqtl_coloc_lead_rank = rank,
                                            eqtl_coloc_lead_p = lead_p, e_rare = lead_rare, lead_dist, n_eqtl_leads, n_pairs_tested, n_pairs_disjoint,
                                            p_start, p_end, e_start = l_start, e_end = l_end, int_start, int_end, uni_start, uni_end)],
           by = c("tissue", "promoter", "type", "k"))
stopifnot(nrow(E) == TOPm[egene == TRUE, .N])
X <- copy(TOPm[egene == FALSE])
X[, `:=`(eqtl_coloc_lead = eqtl_top, eqtl_coloc_lead_rank = 1L, eqtl_coloc_lead_p = NA_real_, n_eqtl_leads = 1L,
         n_pairs_tested = as.integer(status == "tested"), n_pairs_disjoint = as.integer(status == DLS))]
NEW <- rbind(E, X, use.names = TRUE, fill = TRUE)
NEW[, call := callf(status, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4)]
NEW[, `:=`(call_200kb = NA_character_, PP.H4_200kb = NA_real_, PP.H3_200kb = NA_real_, PP.H1_200kb = NA_real_, source = "new (mirror v3 re-run)")]
F <- rbind(C, NEW, use.names = TRUE, fill = TRUE)
stopifnot(nrow(F) == nrow(A), !anyDuplicated(F[, .(tissue, promoter, type, k)]))
F[, call := factor(call, levels = LEV)]
AB <- setNames(c("H4", "H3", "DL", "H1", "H2", "H0", "amb", "NA"), LEV); F[, cc := AB[as.character(call)]]
classify <- function(cc) { t <- cc[cc != "NA"]
  if (any(cc == "H3")) return("distinct: any H3"); if (any(cc == "DL")) return("distinct: different locus (no H3)")
  if (length(t) && all(t == "H4")) return("shared, all signals"); if (any(t == "H4")) return("shared, rest unresolved")
  if (any(t == "H1")) return("distinct: H1 only"); if (any(t == "H2")) return("unresolved: H2"); if (any(t == "H0")) return("unresolved: H0")
  if (any(t == "amb")) return("unresolved: ambiguous"); "unresolved: not assessable" }
G <- F[, .(gene = gene[1], egene = egene[1], eqtl_qval = eqtl_qval[1], eqtl_top = eqtl_top[1], n_eqtl_leads = max(n_eqtl_leads), n_signals = .N,
           n_tested = sum(status == "tested"), n_H4 = sum(cc == "H4"), n_H3 = sum(cc == "H3"), n_DL = sum(cc == "DL"), n_H1 = sum(cc == "H1"),
           n_H2 = sum(cc == "H2"), n_H0 = sum(cc == "H0"), n_amb = sum(cc == "amb"), n_NA = sum(cc == "NA"), detail = classify(cc)), by = .(tissue, gene_base)]
G[egene == FALSE, detail := "no eQTL in this tissue"]
F4 <- c("puQTL distinct from the eQTL", "No eQTL in this tissue", "puQTL shared with the eQTL", "Partly shared with the eQTL", "Unresolved")
G[, fig4e := factor(fcase(detail == "no eQTL in this tissue", F4[2], startsWith(detail, "distinct"), F4[1], detail == "shared, all signals", F4[3],
                          detail == "shared, rest unresolved", F4[4], default = F4[5]), levels = F4)]
G0 <- z(file.path(LW, "final_clump/gene_tissue_classes.tsv.gz"))
G <- merge(G, G0[, .(tissue, gene_base, fig4e_mirror_v2 = factor(fig4e, levels = F4), detail_mirror_v2 = detail)], by = c("tissue", "gene_base"), all.x = TRUE)
fwrite(G, file.path(FO, "gene_tissue_classes.tsv.gz"), sep = "\t")
front <- c("tissue", "gene", "gene_base", "promoter", "type", "k", "signal_snp", "stat", "rule", "n_signals", "unit_q", "egene", "eqtl_qval", "n_eqtl",
           "status", "n_snps", "n_int", "call", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "PP.H4_p12_5e6", "source")
setcolorder(F, c(front, setdiff(names(F), front))); setorder(F, tissue, promoter, type, k)
fwrite(F[, !"cc"], file.path(FO, "coloc_eqtl_signals.tsv.gz"), sep = "\t")
FL <- merge(F, G[, .(tissue, gene_base, gene_class_fig4e = as.character(fig4e), gene_class_detail = detail)], by = c("tissue", "gene_base"))
FL[, `:=`(puqtl_specific_H1 = call == "puQTL-specific (H1)", different_locus = call == DLS)]
setcolorder(FL, c(front, setdiff(names(FL), front))); setorder(FL, tissue, promoter, type, k)
fwrite(FL[, !"cc"], file.path(FO, "coloc_eqtl_signals_flags.tsv.gz"), sep = "\t")
## Fig. 4e tables
CT <- merge(G0[, .(N_mirror_v2 = .N), by = .(class = fig4e)], G[, .(N_mirror_v3 = .N), by = .(class = fig4e)], by = "class", all = TRUE)
CT <- CT[order(factor(class, levels = F4))]; CT[, pct_mirror_v3 := round(100 * N_mirror_v3 / sum(N_mirror_v3), 1)]
fwrite(CT, file.path(FO, "panel_e_counts.tsv"), sep = "\t")
DET <- G[, .N, by = .(fig4e, detail)][order(fig4e, -N)]; fwrite(DET, file.path(FO, "panel_e_detail.tsv"), sep = "\t")
key <- c("puQTL distinct from the eQTL" = "independent", "No eQTL in this tissue" = "no_eqtl", "puQTL shared with the eQTL" = "shared_all",
         "Partly shared with the eQTL" = "shared_rest", "Unresolved" = "unresolved")
TT <- dcast(G[, .N, by = .(tissue, k = key[as.character(fig4e)])], tissue ~ k, value.var = "N", fill = 0)
TT <- merge(TT, unique(F[, .(tissue, donors = n_eqtl)]), by = "tissue"); TT[, pairs := independent + no_eqtl + shared_all + shared_rest + unresolved]
setcolorder(TT, c("tissue", "independent", "no_eqtl", "shared_all", "shared_rest", "unresolved", "donors", "pairs"))
fwrite(TT, file.path(FO, "panel_e_puqtl_by_tissue.tsv"), sep = "\t")
GX <- merge(G0[, .(tissue, gene_base, v2 = fig4e)], G[, .(tissue, gene_base, v3 = as.character(fig4e))], by = c("tissue", "gene_base"), all = TRUE)
GX[is.na(v2), v2 := "(not in v2)"]; GX[is.na(v3), v3 := "(not in v3)"]
fwrite(dcast(GX[, .N, by = .(v2, v3)], v2 ~ v3, value.var = "N", fill = 0), file.path(FO, "compare_fig4e_mirror_v2_vs_v3.tsv"), sep = "\t")
## moved units: call of v2 p1(O) signal k vs v3 p0(S) signal k
MP <- z(file.path(R2, "results/mirror_v3/unit_map_v2_to_v3.tsv.gz"))
o <- F0[type == "p1", .(tissue, promoter_v2 = promoter, k, call_v2 = call)]; n <- F[type == "p0", .(tissue, promoter_v3 = promoter, k, call_v3 = as.character(call))]
MV <- merge(merge(MP[, .(tissue, promoter_v2, promoter_v3)], o, by = c("tissue", "promoter_v2")), n, by = c("tissue", "promoter_v3", "k"))
fwrite(dcast(MV[, .N, by = .(call_v2, call_v3)], call_v2 ~ call_v3, value.var = "N", fill = 0), file.path(FO, "compare_moved_signal_calls.tsv"), sep = "\t")
sink(file.path(FO, "log.txt"), split = TRUE)
cat("signals:", nrow(F), "(carried", F[source != "new (mirror v3 re-run)", .N], "; new", F[source == "new (mirror v3 re-run)", .N], ")\n")
cat("status:\n"); print(F[, .N, by = status][order(-N)])
cat("signal calls (all / eGenes):\n"); print(dcast(F[, .N, by = .(call, egene)], call ~ egene, value.var = "N", fill = 0))
cat("\nmoved signals, same call v2 p1(O) -> v3 p0(S):", MV[call_v2 == call_v3, .N], "of", nrow(MV), "\n")
print(MV[call_v2 != call_v3, .N, by = .(call_v2, call_v3)][order(-N)])
cat("\nFig. 4e gene x tissue (mirror v2 -> v3):\n"); print(CT); print(DET)
cat("\ngene x tissue pairs v2", nrow(G0), "; v3", nrow(G), "; same class", GX[v2 == v3, .N], "\n")
print(GX[v2 != v3, .N, by = .(v2, v3)][order(-N)])
sink(); cat("DONE\n")
