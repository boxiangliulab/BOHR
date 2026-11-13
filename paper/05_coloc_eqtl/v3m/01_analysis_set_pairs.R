## Mirror rule v3 (author 2026-11-04 final: exact switches are p0 on the silent TSS): eQTL coloc re-keyed to units v3m.
## Step 1: the v3m analysis set (rules of coloc_eqtl_v3/00_build_analysis_set.R, unchanged: every significant primary unit x its
## backward-kept signals; n_kept >= 2 -> LOO scan, else marginal) and the rows that need a new coloc run.
##   carried : unit not re-run in the v3m conditional analysis, same (tissue, promoter, type, k, signal_snp, stat) as before
##             -> its final_clump result is reused unchanged (same SNPs, same statistics, same windows)
##   new     : every signal of a re-run unit (mostly p0(partner) twins of v2 p1(canonical)); these get new pairs:
##             puQTL LD window = the window of (tissue, signal SNP) in lead_windows.tsv (it depends only on the SNP);
##             eGene leads = the gene's clumped eQTL leads of sens_clump/pairs.tsv.gz (they depend only on tissue x gene);
##             non-eGene lead = the GTEx top variant of pairs.tsv.gz; intersection / union recomputed as in 04_pairs.R.
## Output: results/coloc_eqtl_v3m/{analysis_set.rds, new_rows.tsv, pairs_clump_new.tsv.gz, pairs_top_new.tsv.gz, missing.tsv}
source("/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/scripts/coloc_final_lib.R")
R2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"; OUT <- file.path(R2, "results/coloc_eqtl_v3m"); LW <- file.path(R2, "results/coloc_eqtl_ldwin")
EQ <- "/home/e1327989/project/GTEx/Rawdata/gtex_v11/GTEx_Analysis_v11_eQTL"
U <- rbindlist(lapply(Sys.glob(file.path(R2, "results/final_nogate_v3m/freeze/*.real.units.rds")), function(f) { a <- readRDS(f); setDT(a)
  a[, .(tissue, promoter, type, gene, q, leader, n, n0, n1, n_mid, dedup_drop, estimable, n_prom_gene, mirror_role, mirror_partner, twin, eq_class, merged_component)] }))
U[, primary := dedup_drop == FALSE & estimable == TRUE]; U[, sig := primary & !is.na(q) & q <= 0.05]
S <- fread(cmd = paste("zcat", file.path(R2, "results/conditional_final_v3m/freeze_real_signals.tsv.gz")))[, .(tissue, promoter, type, k, snp, backward_keep)]
A <- merge(U[sig == TRUE], S, by = c("tissue", "promoter", "type")); stopifnot(uniqueN(A[, .(tissue, promoter, type)]) == U[sig == TRUE, .N])
A[, n_signals := .N, by = .(tissue, promoter, type)]; A[, n_kept := sum(backward_keep), by = .(tissue, promoter, type)]
A <- A[backward_keep == TRUE]
A[, stat := fifelse(n_kept >= 2, "loo", "marginal")]
A[, rule := fcase(n_signals == 1, "single", n_kept == 1, "kept1_marginal", n_kept < n_signals, "loo_conditions_on_failed", default = "loo")]
A[, `:=`(signal_snp = snp, signal_pos = pos_of(snp), chr = chr_of(snp), gene_base = sub("\\..*$", "", gene))]
EG <- rbindlist(lapply(unique(A$tissue), function(t) { e <- fread(cmd = paste("zcat", file.path(EQ, paste0(t, ".v11.eGenes.txt.gz"))), select = c("gene_id", "qval"))
  e[, .(tissue = t, gene_base = sub("\\..*$", "", gene_id), eqtl_qval = qval)] }))
A <- merge(A, EG, by = c("tissue", "gene_base"), all.x = TRUE); A[, `:=`(egene = !is.na(eqtl_qval) & eqtl_qval <= 0.05, eqtl_tested = !is.na(eqtl_qval))]
NE <- data.table(tissue = unique(A$tissue))[, n_eqtl := vapply(tissue, function(t) length(readLines(file.path(V2, "ld_v3/donors", paste0(t, ".txt")))), 1L)]
A <- merge(A, NE, by = "tissue"); setorder(A, tissue, promoter, type, k)
## carried vs new
RUN <- fread(file.path(R2, "work/cond_v3m/units_to_run.tsv"))
A0 <- readRDS(file.path(R2, "results/coloc_eqtl_v3/analysis_set.rds"))
key <- function(d) paste(d$tissue, d$promoter, d$type, d$k, d$signal_snp, d$stat)
A[, rerun := paste(tissue, promoter, type) %in% RUN[, paste(tissue, promoter, type)]]
A[, carried := !rerun & key(A) %in% key(A0)]
cat("analysis rows v3m:", nrow(A), "; carried", A[carried == TRUE, .N], "; new", A[carried == FALSE, .N], "(re-run units", A[rerun == TRUE, .N], ")\n")
cat("non-carried rows of units not re-run (should be 0):", A[carried == FALSE & rerun == FALSE, .N], "\n")
saveRDS(A, file.path(OUT, "analysis_set.rds"))
N <- A[carried == FALSE]
## puQTL windows (depend on tissue x SNP only)
L <- fread(file.path(LW, "lead_windows.tsv"))
PW <- unique(L[role %in% c("puQTL", "both") | TRUE, .(tissue, signal_snp = snp, absent, n_partners, win_start, win_end)], by = c("tissue", "signal_snp"))
N <- merge(N, PW, by = c("tissue", "signal_snp"), all.x = TRUE)
cat("new rows: window found", N[!is.na(win_start), .N], "; SNP in lead_windows but absent from the panel", N[absent == TRUE, .N], "; SNP not in lead_windows", N[is.na(absent), .N], "\n")
fwrite(N[is.na(win_start), .(tissue, promoter, type, k, signal_snp, absent)], file.path(OUT, "missing.tsv"), sep = "\t")
setnames(N, c("win_start", "win_end", "n_partners"), c("p_start", "p_end", "p_partners"))
## eGene pairs: the gene's clumped leads
CL <- unique(fread(cmd = paste("zcat", file.path(LW, "sens_clump/pairs.tsv.gz")))[, .(tissue, gene_base, eqtl_top, lead, rank, lead_p, n_absorbed, n_partners, l_start, l_end)])
PC <- merge(N[egene == TRUE & !is.na(p_start), .(tissue, gene_base, promoter, type, k, chr, signal_snp, signal_pos, p_start, p_end)], CL, by = c("tissue", "gene_base"), allow.cartesian = TRUE)
cat("eGene new rows:", N[egene == TRUE, .N], "; with clumped leads", uniqueN(PC[, .(tissue, promoter, type, k)]), "; genes without clumped leads", uniqueN(N[egene == TRUE & !paste(tissue, gene_base) %in% CL[, paste(tissue, gene_base)], .(tissue, gene_base)]), "\n")
PC[, `:=`(int_start = pmax(p_start, l_start), int_end = pmin(p_end, l_end), uni_start = pmin(p_start, l_start), uni_end = pmax(p_end, l_end))]
PC[, overlap_bp := pmax(0L, int_end - int_start)]; PC[, lead_dist := abs(signal_pos - pos_of(lead))]
fwrite(PC, file.path(OUT, "pairs_clump_new.tsv.gz"), sep = "\t")
## non-eGene pairs: GTEx top variant of the gene (pairs.tsv.gz of the LD-window run)
TP <- unique(fread(cmd = paste("zcat", file.path(LW, "pairs.tsv.gz")))[, .(tissue, gene_base, eqtl_top, eqtl_top_pos, e_rare, e_partners, e_start, e_end)], by = c("tissue", "gene_base"))
PT <- merge(N[!is.na(p_start), .(tissue, gene_base, promoter, type, k, chr, signal_snp, signal_pos, egene, p_partners, p_start, p_end)], TP, by = c("tissue", "gene_base"))
PT[, `:=`(int_start = pmax(p_start, e_start), int_end = pmin(p_end, e_end), uni_start = pmin(p_start, e_start), uni_end = pmax(p_end, e_end))]
PT[, overlap_bp := pmax(0L, int_end - int_start)]; PT[, lead_dist := abs(signal_pos - eqtl_top_pos)]
cat("top-variant pairs for new rows:", nrow(PT), "of", N[!is.na(p_start), .N], "rows with a window\n")
fwrite(PT, file.path(OUT, "pairs_top_new.tsv.gz"), sep = "\t")
## check: for old rows, the lead_windows window of the signal SNP equals the p window stored in the old pairs
OP <- unique(fread(cmd = paste("zcat", file.path(LW, "pairs.tsv.gz")))[, .(tissue, signal_snp, p_start, p_end)])
ck <- merge(OP, PW, by = c("tissue", "signal_snp"))
cat("check old pairs vs lead_windows: rows", nrow(ck), "; identical windows", ck[p_start == win_start & p_end == win_end, .N], "\n")
