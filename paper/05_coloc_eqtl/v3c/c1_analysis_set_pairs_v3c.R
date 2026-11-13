## eQTL coloc on units v3c (author 2026-11-06 "correction": conditional refit of two units), step c1: the analysis set and the pairs of
## its new rows. Rules of coloc_eqtl_v3m/01_analysis_set_pairs.R and 04_pairs_new_all.R, unchanged (every significant primary unit x
## its backward-kept signals; n_kept >= 2 -> LOO scan, else marginal; LD windows = farthest r2 >= 0.1 partner +/- 50 kb; eGenes:
## the gene's clumped leads; non-eGenes: the GTEx top variant; intersection / union as before).
##   carried : same (tissue, promoter, type, k, signal_snp, stat) as the v3m analysis set and the unit not re-run -> v3m result reused
##   new     : the rows of the re-run units (work/cond_v3c/units_to_run.tsv)
## Windows and leads are looked up in the existing tables (coloc_eqtl_ldwin lead/rarelead windows, the v3m NSCC windows, eqtl
## leads, clumped leads); a missing one stops the script (it would need a PLINK run).
## Output: results/coloc_eqtl_v3c/{analysis_set.rds, pairs_top_new_all.tsv.gz, pairs_clump_new_all.tsv.gz}
source("/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/scripts/coloc_final_lib.R")
R2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"; OUT <- file.path(R2, "results/coloc_eqtl_v3c"); dir.create(OUT, showWarnings = FALSE)
LW <- file.path(R2, "results/coloc_eqtl_ldwin"); O3 <- file.path(R2, "results/coloc_eqtl_v3m"); PAD <- 50000L
EQ <- "/home/e1327989/project/GTEx/Rawdata/gtex_v11/GTEx_Analysis_v11_eQTL"
U <- rbindlist(lapply(Sys.glob(file.path(R2, "results/final_nogate_v3c/freeze/*.real.units.rds")), function(f) { a <- readRDS(f); setDT(a)
  a[, .(tissue, promoter, type, gene, q, leader, n, n0, n1, n_mid, dedup_drop, estimable, n_prom_gene, mirror_role, mirror_partner, twin, eq_class, merged_component)] }))
U[, primary := dedup_drop == FALSE & estimable == TRUE]; U[, sig := primary & !is.na(q) & q <= 0.05]
S <- fread(cmd = paste("zcat", file.path(R2, "results/conditional_final_v3c/freeze_real_signals.tsv.gz")))[, .(tissue, promoter, type, k, snp, backward_keep)]
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
RUN <- fread(file.path(R2, "work/cond_v3c/units_to_run.tsv"))
A0 <- readRDS(file.path(O3, "analysis_set.rds"))
key <- function(d) paste(d$tissue, d$promoter, d$type, d$k, d$signal_snp, d$stat)
A[, rerun := paste(tissue, promoter, type) %in% RUN[, paste(tissue, promoter, type)]]
A[, carried := !rerun & key(A) %in% key(A0)]
cat("analysis rows v3c:", nrow(A), "(v3m", nrow(A0), "); carried", A[carried == TRUE, .N], "; new", A[carried == FALSE, .N], "\n")
cat("non-carried rows of units not re-run (must be 0):", A[carried == FALSE & rerun == FALSE, .N], "\n"); stopifnot(A[carried == FALSE & rerun == FALSE, .N] == 0)
cat("v3m rows dropped (re-run units):\n"); print(A0[paste(tissue, promoter, type) %in% RUN[, paste(tissue, promoter, type)], .(tissue, promoter, type, k, signal_snp, stat, rule)])
cat("new rows:\n"); print(A[carried == FALSE, .(tissue, promoter, type, k, signal_snp, stat, rule, egene, gene_base)])
saveRDS(A, file.path(OUT, "analysis_set.rds"))
N <- A[carried == FALSE]
## windows: lead_windows (panel leads), rarelead A', v3m NSCC windows (04_pairs_new_all.R construction)
L <- fread(file.path(LW, "lead_windows.tsv")); RW <- fread(file.path(LW, "rarelead_windows.tsv"))
W0 <- rbind(L[absent == FALSE, .(tissue, snp, n_partners, w_start = win_start, w_end = win_end)], RW[, .(tissue, snp, n_partners = nAp, w_start = sAp, w_end = eAp)])
nw <- rbindlist(lapply(Sys.glob(file.path(O3, "nscc/out/*.win.tsv")), function(f) { b <- basename(f); x <- fread(f, colClasses = list(character = "snp"))
  if (!nrow(x)) return(NULL); x[, tissue := sub("\\.chr.*$", "", b)] }))
nw[, pos := pos_of(snp)]; nw[n_partners == 0, `:=`(min_partner_pos = pos, max_partner_pos = pos)]
W1 <- nw[, .(tissue, snp, n_partners, w_start = pmin(pos, min_partner_pos) - PAD, w_end = pmax(pos, max_partner_pos) + PAD)]
W <- rbind(W0, W1[!paste(tissue, snp) %in% W0[, paste(tissue, snp)]])
N <- merge(N, W[, .(tissue, signal_snp = snp, p_partners = n_partners, p_start = w_start, p_end = w_end)], by = c("tissue", "signal_snp"), all.x = TRUE)
stopifnot(!anyNA(N$p_start))                                                   # every new signal SNP already has its LD window
## GTEx top-variant pairs (all new rows) and clumped-lead pairs (eGene rows): gene-level inputs from the existing pair tables
TP <- unique(rbind(fread(cmd = paste("zcat", file.path(LW, "pairs.tsv.gz")))[, .(tissue, gene_base, eqtl_top, eqtl_top_pos, e_rare, e_partners, e_start, e_end)],
                   fread(cmd = paste("zcat", file.path(O3, "pairs_top_new_all.tsv.gz")))[, .(tissue, gene_base, eqtl_top, eqtl_top_pos, e_rare, e_partners, e_start, e_end)]),
             by = c("tissue", "gene_base"))
PT <- merge(N[, .(tissue, gene_base, promoter, type, k, chr, signal_snp, signal_pos, egene, p_partners, p_start, p_end)], TP, by = c("tissue", "gene_base"))
stopifnot(nrow(PT) == nrow(N))
PT[, `:=`(int_start = pmax(p_start, e_start), int_end = pmin(p_end, e_end), uni_start = pmin(p_start, e_start), uni_end = pmax(p_end, e_end))]
PT[, overlap_bp := pmax(0L, int_end - int_start)]; PT[, lead_dist := abs(signal_pos - eqtl_top_pos)]
CLc <- c("tissue", "gene_base", "eqtl_top", "lead", "rank", "lead_p", "n_absorbed", "n_partners", "l_start", "l_end")
CL <- unique(rbind(fread(cmd = paste("zcat", file.path(LW, "sens_clump/pairs.tsv.gz")))[, ..CLc], fread(cmd = paste("zcat", file.path(O3, "pairs_clump_new_all.tsv.gz")))[, ..CLc]),
             by = c("tissue", "gene_base", "lead"))
PC <- merge(N[egene == TRUE, .(tissue, gene_base, promoter, type, k, chr, signal_snp, signal_pos, p_start, p_end)], CL, by = c("tissue", "gene_base"), allow.cartesian = TRUE)
stopifnot(all(N[egene == TRUE, paste(tissue, promoter, type, k)] %in% PC[, paste(tissue, promoter, type, k)]))
PC[, `:=`(int_start = pmax(p_start, l_start), int_end = pmin(p_end, l_end), uni_start = pmin(p_start, l_start), uni_end = pmax(p_end, l_end))]
PC[, overlap_bp := pmax(0L, int_end - int_start)]; PC[, lead_dist := abs(signal_pos - pos_of(lead))]
fwrite(PT, file.path(OUT, "pairs_top_new_all.tsv.gz"), sep = "\t"); fwrite(PC, file.path(OUT, "pairs_clump_new_all.tsv.gz"), sep = "\t")
cat("pairs: top", nrow(PT), "; clumped", nrow(PC), "\n"); print(PT[, .(tissue, promoter, k, signal_snp, p_start, p_end, eqtl_top, e_start, e_end, uni_start, uni_end)])
print(PC[, .(tissue, promoter, k, lead, rank, l_start, l_end, uni_start, uni_end)])
## eQTL rows available for the new regions (05 reads the coloc_eqtl_ldwin and v3m subsets)
for (i in seq_len(nrow(PT))) { r <- PT[i]
  fs <- c(file.path(LW, "eqtl_subset", paste0(r$tissue, ".tsv.gz")), file.path(O3, "eqtl_subset_top", paste0(r$tissue, ".tsv.gz")),
          file.path(LW, "sens_clump/eqtl_subset", paste0(r$tissue, ".tsv.gz")), file.path(O3, "eqtl_subset_clump", paste0(r$tissue, ".tsv.gz")))
  E <- unique(rbindlist(lapply(fs[file.exists(fs)], function(f) fread(cmd = paste("zcat", f))))[gene_base == r$gene_base], by = "variant_id")
  sp <- range(pos_of(E$variant_id)); need <- c(r$uni_start, max(r$uni_end, PC[promoter == r$promoter & k == r$k, max(uni_end, -Inf)]))
  cat(sprintf("eQTL rows for %s %s: %d variants spanning %d..%d; coloc region needs %d..%d -> %s\n", r$tissue, r$gene_base, nrow(E), sp[1], sp[2],
              min(r$uni_start, PC[promoter == r$promoter & k == r$k, min(uni_start, Inf)]), need[2],
              if (sp[1] <= min(r$uni_start, PC[promoter == r$promoter & k == r$k, min(uni_start, Inf)]) + PAD && sp[2] >= need[2] - PAD) "covered" else "CHECK")) }
