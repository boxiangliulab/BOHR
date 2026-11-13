## Figure 4 panel e, step 3: merge the per-tissue puQTL x eQTL coloc and classify.
## Signal level (every analysis signal of step 0): status tested / not assessable (< 10 shared SNPs) / no eQTL rows /
##   coloc error; for tested rows the call at PP >= 0.8 (p12 = 1e-5): shared (H4), distinct (H3), puQTL-specific (H1),
##   eQTL-specific (H2), neither (H0), otherwise ambiguous.
## Unit level: the call of the unit's lead signal (k = 1), plus any_shared / any_distinct over all its signals.
## Gene x tissue level (the unit of the old panel e), over genes TESTED for puQTL in that tissue (>= 1 estimable primary
##   unit in units v3):
##   puQTL and eGene  -> the strongest call over the gene's signals, precedence H4 > H3 > H1 > H2 > H0 > ambiguous >
##                       not assessable;
##   puQTL, not eGene -> puQTL-only;   eGene, no significant puQTL unit -> eQTL-only.
## eGene = GTEx v11 eGenes qval <= 0.05. Output: results/coloc_eqtl_v3/coloc_eqtl_{signals,units,gene_tissue}.tsv.gz,
## panel_e_counts.tsv (+ figure_v3/data/panel_e.tsv).
suppressPackageStartupMessages(library(data.table))
V2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"; OUT <- file.path(V2, "results/coloc_eqtl_v3")
EQ <- "/home/e1327989/project/GTEx/Rawdata/gtex_v11/GTEx_Analysis_v11_eQTL"
A <- readRDS(file.path(OUT, "analysis_set.rds"))
R <- rbindlist(lapply(Sys.glob(file.path(OUT, "per_tissue/*.rds")), readRDS), fill = TRUE)
stopifnot(nrow(R) == nrow(A), !anyDuplicated(R[, .(tissue, promoter, type, k)]))
call_of <- function(st, h0, h1, h2, h3, h4) fifelse(st != "tested", "not assessable",
  fifelse(h4 >= 0.8, "shared (H4)", fifelse(h3 >= 0.8, "distinct (H3)", fifelse(h1 >= 0.8, "puQTL-specific (H1)",
  fifelse(h2 >= 0.8, "eQTL-specific (H2)", fifelse(h0 >= 0.8, "neither (H0)", "ambiguous"))))))
R[, call := call_of(status, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4)]
LEV <- c("shared (H4)", "distinct (H3)", "puQTL-specific (H1)", "eQTL-specific (H2)", "neither (H0)", "ambiguous", "not assessable")
R[, call := factor(call, levels = LEV)]
setcolorder(R, c("tissue", "gene", "gene_base", "promoter", "type", "k", "signal_snp", "stat", "rule", "n_signals", "unit_q", "egene", "eqtl_qval",
                 "n_eqtl", "status", "n_snps", "call", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "PP.H4_p12_5e6"))
setorder(R, tissue, promoter, type, k)
fwrite(R, file.path(OUT, "coloc_eqtl_signals.tsv.gz"), sep = "\t")
## unit level
Un <- R[, .(gene = gene[1], gene_base = gene_base[1], egene = egene[1], n_signals_kept = .N, lead_call = call[k == min(k)][1],
            any_shared = any(call == "shared (H4)"), any_distinct = any(call == "distinct (H3)"), any_tested = any(status == "tested")),
        by = .(tissue, promoter, type)]
fwrite(Un, file.path(OUT, "coloc_eqtl_units.tsv.gz"), sep = "\t")
## gene x tissue level over puQTL-tested genes
U <- rbindlist(lapply(Sys.glob(file.path(V2, "results/final_nogate_v3/freeze/*.real.units.rds")), function(f) {
  a <- as.data.table(readRDS(f)); a[dedup_drop == FALSE & estimable == TRUE, .(tissue, gene_base = sub("\\..*$", "", gene), sig = !is.na(q) & q <= 0.05)] }))
GT <- U[, .(puqtl = any(sig)), by = .(tissue, gene_base)]
EG <- rbindlist(lapply(unique(GT$tissue), function(t) { e <- fread(cmd = paste("zcat", file.path(EQ, paste0(t, ".v11.eGenes.txt.gz"))), select = c("gene_id", "qval"))
  e[, .(tissue = t, gene_base = sub("\\..*$", "", gene_id), egene = qval <= 0.05)] }))
GT <- merge(GT, EG, by = c("tissue", "gene_base"), all.x = TRUE)
GT[is.na(egene), egene := FALSE]
best <- R[, .(best = LEV[min(as.integer(call))]), by = .(tissue, gene_base)]
GT <- merge(GT, best, by = c("tissue", "gene_base"), all.x = TRUE)
GT[, class := fcase(puqtl & egene, best, puqtl & !egene, "puQTL-only", !puqtl & egene, "eQTL-only", default = "neither puQTL nor eGene")]
fwrite(GT, file.path(OUT, "coloc_eqtl_gene_tissue.tsv.gz"), sep = "\t")
ORD <- c("eQTL-only", "puQTL-only", LEV)
cnt <- GT[class != "neither puQTL nor eGene", .N, by = class][, class := factor(class, levels = ORD)][order(class)]
fwrite(cnt, file.path(OUT, "panel_e_counts.tsv"), sep = "\t")
fwrite(cnt, "/projects/group/bliu/e1327989/GTEx/QTL_figure/figure_v3/data/panel_e.tsv", sep = "\t")
## report
cat("signals", nrow(R), "; status:\n"); print(R[, .N, by = status])
cat("signal-level calls (tested):\n"); print(R[status == "tested", .N, by = call][order(call)])
cat("signal-level calls on eGenes / non-eGenes:\n"); print(dcast(R[, .N, by = .(egene, call)], call ~ egene, value.var = "N", fill = 0))
cat("unit-level lead calls:\n"); print(Un[, .N, by = .(egene, lead_call)][order(egene, lead_call)])
cat("gene x tissue (puQTL-tested genes): ", nrow(GT), "; puQTL genes", GT[puqtl == TRUE, .N], "; eGenes", GT[egene == TRUE, .N], "\n"); print(cnt)
cat("by type, lead-signal calls on eGenes:\n"); print(dcast(Un[egene == TRUE, .N, by = .(type, lead_call)], lead_call ~ type, value.var = "N", fill = 0))
