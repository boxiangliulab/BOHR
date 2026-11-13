## Figure 4 panel e, LD-window coloc, step 3 (test, 2026-10-12): windows of the GTEx eGene top variants absent from the
## MAF >= 0.05 panel. Their r2 >= 0.1 partners were computed on
## NSCC from the unfiltered v9 VCF (944-donor freeze, --maf 0.005, tissue donors; job rarelead.pbs). Three windows per
## rare lead, all = [min(lead, partners) - 50 kb, max(lead, partners) + 50 kb]:
##   A    : partners of any frequency (MAF >= 0.005) -- the author's "top variant LD" with the full VCF;
##   Apan : partners restricted to the panel (MAF >= 0.05) -- the same partner universe as every other lead;
##   C    : no partner, lead +/- 50 kb (my first suggestion).
## Validity proxy per puQTL signal x eQTL pair: panel SNPs inside the intersection of the two windows (>= 5). The real
## rule counts SNPs with statistics in both datasets; panel SNPs are an upper bound.
suppressPackageStartupMessages(library(data.table))
O <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/results/coloc_eqtl_ldwin"
PAD <- 50000
L <- fread(file.path(O, "lead_windows.tsv"))
R <- L[absent == TRUE, .(tissue, snp, chr, pos = as.integer(sub("^chr[0-9]+_([0-9]+)_.*$", "\\1", snp)))]
fs <- list.files(file.path(O, "rarelead"), "\\.ld$", full.names = TRUE)
D <- rbindlist(lapply(fs, function(f) { x <- fread(f); if (!nrow(x)) return(NULL)
  x[, tissue := sub("\\.chr[0-9]+\\.ld$", "", basename(f))][, .(tissue, snp = SNP_A, bp_b = BP_B, snp_b = SNP_B, r2 = R2)] }))
PAN <- fread(cmd = paste("zcat", file.path(O, "panel/panel_snps.tsv.gz")), header = FALSE, col.names = c("chrn", "id", "bp"))
D[, in_panel := snp_b %in% PAN$id]
found <- unique(D[, .(tissue, snp)])            # every typed lead is reported with itself (r2 = 1)
R[, typed := paste(tissue, snp) %in% found[, paste(tissue, snp)]]
cat("rare leads:", nrow(R), "; typed in the full VCF (MAF >= 0.005):", R[typed == TRUE, .N], "\n")
P <- D[snp_b != snp]
S <- function(x) x[, .(n = .N, lo = min(bp_b), hi = max(bp_b)), by = .(tissue, snp)]
A  <- S(P); Ap <- S(P[in_panel == TRUE])
R <- merge(R, A[, .(tissue, snp, nA = n, loA = lo, hiA = hi)], by = c("tissue", "snp"), all.x = TRUE)
R <- merge(R, Ap[, .(tissue, snp, nAp = n, loAp = lo, hiAp = hi)], by = c("tissue", "snp"), all.x = TRUE)
R[is.na(nA), `:=`(nA = 0L, loA = pos, hiA = pos)]; R[is.na(nAp), `:=`(nAp = 0L, loAp = pos, hiAp = pos)]
R[, `:=`(sA = pmin(pos, loA) - PAD, eA = pmax(pos, hiA) + PAD, sAp = pmin(pos, loAp) - PAD, eAp = pmax(pos, hiAp) + PAD,
         sC = pos - PAD, eC = pos + PAD)]
cat("\npartners (typed leads): any MAF / panel only\n"); print(R[typed == TRUE, .(no_partner_A = mean(nA == 0), no_partner_Apan = mean(nAp == 0),
  med_nA = median(nA), med_nAp = median(nAp))])
q <- function(x) round(quantile(x / 1e3, c(.1, .25, .5, .75, .9, 1)))
cat("width kb  A:\n"); print(q(R[typed == TRUE, eA - sA])); cat("width kb  Apan:\n"); print(q(R[typed == TRUE, eAp - sAp]))
cat("width kb  present leads (reference):\n"); print(q(L[absent == FALSE, win_end - win_start]))
fwrite(R, file.path(O, "rarelead_windows.tsv"), sep = "\t")

## validity proxy on the affected signal x eQTL pairs
W <- fread(cmd = paste("zcat", file.path(O, "pair_windows.tsv.gz")))
X <- merge(W[e_absent == TRUE], R[, .(tissue, eqtl_lead = snp, typed, sA, eA, sAp, eAp, sC, eC)], by = c("tissue", "eqtl_lead"))
POS <- split(sort(PAN$bp), PAN$chrn[order(PAN$bp)])   # sorted panel positions per chromosome
npan <- function(chr, s, e) { k <- sub("chr", "", chr); out <- integer(length(s))
  for (c in unique(k)) { i <- which(k == c); p <- POS[[c]]
    out[i] <- ifelse(s[i] > e[i], 0L, findInterval(e[i], p) - findInterval(s[i] - 1, p)) }
  out }
for (v in c("A", "Ap", "C")) {
  s <- pmax(X$p_start, X[[paste0("s", v)]]); e <- pmin(X$p_end, X[[paste0("e", v)]])
  X[, (paste0("n_int_", v)) := npan(sub("_.*$", "", eqtl_lead), s, e)]
}
X[, gene_tissue := paste(tissue, gene_base)]
cat("\nsignal x eQTL pairs with a rare eQTL lead:", nrow(X), "(eGene:", X[egene == TRUE, .N], ")\n")
cat("valid (>= 5 panel SNPs in the intersection), eGene pairs:\n")
print(X[egene == TRUE, .(n = .N, valid_A = sum(n_int_A >= 5), valid_Apan = sum(n_int_Ap >= 5), valid_C = sum(n_int_C >= 5))])
cat("eGene gene x tissue with no valid signal pair:\n")
print(X[egene == TRUE, .(A = any(n_int_A >= 5), Ap = any(n_int_Ap >= 5), C = any(n_int_C >= 5)), by = gene_tissue][,
  .(n = .N, none_A = sum(!A), none_Apan = sum(!Ap), none_C = sum(!C))])
Wp <- W[e_absent == FALSE & egene == TRUE]
Wp[, n_int := npan(sub("_.*$", "", eqtl_lead), pmax(p_start, e_start), pmin(p_end, e_end))]
cat("reference, eGene pairs with a panel eQTL lead:\n")
print(Wp[, .(n = .N, valid = sum(n_int >= 5), share = round(mean(n_int >= 5), 3))])
fwrite(X, file.path(O, "rarelead_pairs.tsv"), sep = "\t")
