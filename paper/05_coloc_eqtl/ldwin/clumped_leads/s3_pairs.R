## Figure 4 panel e, LD-window coloc, sensitivity "eQTL leads by recursive clumping", step s3: signal x eQTL-lead pairs.
## Every eGene signal of the main analysis (pairs.tsv; puQTL window unchanged) is paired with EVERY clumped eQTL lead of its
## gene in that tissue (s2, NSCC; fetched to sens_clump/nscc_out/<tissue>/). Lead window as in the main analysis: r2 >= 0.1
## partners among the panel variants, farthest partner +/- 50 kb (lead +/- 50 kb without a partner). Intersection and
## union per pair as in 04_pairs.R. Genes that are not eGenes are not part of the sensitivity (no significant eQTL variant;
## Fig. 4e class "no eQTL in this tissue" whatever the region).
## Output: sens_clump/pairs.tsv.gz (one row per signal x lead), sens_clump/eqtl_windows.tsv (per tissue x gene).
suppressPackageStartupMessages(library(data.table))
O <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/results/coloc_eqtl_ldwin"; SO <- file.path(O, "sens_clump")
PAD <- 50000L
tis <- readLines(file.path(O, "tissues.txt"))
L <- rbindlist(lapply(tis, function(t) { x <- fread(file.path(SO, "nscc_out", t, "leads.tsv")); x[, tissue := t] }))
Q <- rbindlist(lapply(tis, function(t) { f <- file.path(SO, "nscc_out", t, "lead_partners.tsv"); if (!file.exists(f) || file.size(f) == 0) return(NULL)
  x <- fread(f); x[, tissue := t] }))
cat("eGene x tissue:", uniqueN(L[, .(tissue, gene_base)]), "; without a typed significant variant:", L[is.na(lead), .N], "; leads:", L[!is.na(lead), .N], "\n")
L <- L[!is.na(lead)]
L <- merge(L, Q[, .(tissue, lead, n_partners, min_partner_pos, max_partner_pos)], by = c("tissue", "lead"), all.x = TRUE)
stopifnot(!anyNA(L$n_partners))
L[, pos := as.integer(sub("^chr[0-9]+_([0-9]+)_.*$", "\\1", lead))]
L[n_partners == 0, `:=`(min_partner_pos = pos, max_partner_pos = pos)]
L[, `:=`(l_start = pmin(pos, min_partner_pos) - PAD, l_end = pmax(pos, max_partner_pos) + PAD)]
cat("leads per eGene x tissue:\n"); print(L[, .N, by = .(tissue, gene_base)][, .(genes = .N, one = sum(N == 1), two = sum(N == 2), three_plus = sum(N >= 3),
  median = median(N), max = max(N))])
P <- fread(cmd = paste("zcat", file.path(O, "pairs.tsv.gz")))[egene == TRUE]
## does the rank-1 clump lead equal the GTEx eGenes top variant used by the main analysis?
top <- unique(P[, .(tissue, gene_base, eqtl_top)])
chk <- merge(top, L[rank == 1, .(tissue, gene_base, lead)], by = c("tissue", "gene_base"), all.x = TRUE)
cat("rank-1 lead == GTEx top variant:", chk[lead == eqtl_top, .N], "of", nrow(chk), "; missing:", chk[is.na(lead), .N], "\n")
X <- merge(P[, .(tissue, promoter, type, k, gene_base, chr, signal_snp, signal_pos, eqtl_top, p_start, p_end)],
           L[, .(tissue, gene_base, lead, rank, lead_p, n_absorbed, n_partners, l_start, l_end)], by = c("tissue", "gene_base"), allow.cartesian = TRUE)
X[, `:=`(int_start = pmax(p_start, l_start), int_end = pmin(p_end, l_end), uni_start = pmin(p_start, l_start), uni_end = pmax(p_end, l_end))]
X[, `:=`(overlap_bp = pmax(0L, int_end - int_start + 1L), lead_dist = abs(as.integer(sub("^chr[0-9]+_([0-9]+)_.*$", "\\1", lead)) - signal_pos))]
setorder(X, tissue, promoter, type, k, rank)
fwrite(X, file.path(SO, "pairs.tsv.gz"), sep = "\t")
W <- X[, .(start = min(uni_start), end = max(uni_end)), by = .(tissue, gene_base, chr)]
fwrite(W, file.path(SO, "eqtl_windows.tsv"), sep = "\t")
cat("signals:", uniqueN(X[, .(tissue, promoter, type, k)]), "of", nrow(P), "eGene signals; signal x lead pairs:", nrow(X),
    "; with positional overlap:", X[overlap_bp > 0, .N], "\n")
cat("signals whose windows overlap no lead (positionally):", X[, .(any = any(overlap_bp > 0)), by = .(tissue, promoter, type, k)][any == FALSE, .N], "\n")
