## Mirror rule v3 coloc, step 4: complete pairs for every new row (step 1 rows + the rows whose windows / leads came from the
## NSCC step 3), with the window rules of coloc_eqtl_ldwin/02_windows.R (farthest r2 >= 0.1 partner +/- 50 kb; lead +/- 50 kb
## without a partner; panel leads from the panel bfile, absent leads option A') and the clumped-lead windows of sens_clump/s3.
## Output: results/coloc_eqtl_v3m/{pairs_top_new_all.tsv.gz, pairs_clump_new_all.tsv.gz, eqtl_windows_{top,clump}_new.tsv}
suppressPackageStartupMessages(library(data.table))
R2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"; O <- file.path(R2, "results/coloc_eqtl_v3m"); LW <- file.path(R2, "results/coloc_eqtl_ldwin")
NS <- file.path(O, "nscc"); PAD <- 50000L
pos_of <- function(s) as.integer(sub("^chr[0-9]+_([0-9]+)_.*$", "\\1", s))
A <- readRDS(file.path(O, "analysis_set.rds")); N <- A[carried == FALSE]
## windows: existing (lead_windows panel leads + rarelead A') and new (NSCC step 3)
L <- fread(file.path(LW, "lead_windows.tsv")); RW <- fread(file.path(LW, "rarelead_windows.tsv"))
W0 <- rbind(L[absent == FALSE, .(tissue, snp, n_partners, w_start = win_start, w_end = win_end, rare = FALSE)],
            RW[, .(tissue, snp, n_partners = nAp, w_start = sAp, w_end = eAp, rare = TRUE)])
nw <- rbindlist(lapply(Sys.glob(file.path(NS, "out/*.win.tsv")), function(f) { b <- basename(f); x <- fread(f, colClasses = list(character = "snp"))
  if (!nrow(x)) return(NULL); x[, `:=`(tissue = sub("\\.chr.*$", "", b), rare = grepl("\\.absent\\.", b))] }))
nw[, pos := pos_of(snp)]; nw[n_partners == 0, `:=`(min_partner_pos = pos, max_partner_pos = pos)]
W1 <- nw[, .(tissue, snp, n_partners, w_start = pmin(pos, min_partner_pos) - PAD, w_end = pmax(pos, max_partner_pos) + PAD, rare)]
need <- fread(file.path(O, "new_inputs/ld_snps.tsv"))
cat("new windows:", nrow(W1), "for", uniqueN(need[, .(tissue, snp)]), "needed (tissue x SNP); missing:", nrow(fsetdiff(unique(need[, .(tissue, snp)]), W1[, .(tissue, snp)])), "\n")
W <- rbind(W0, W1[!paste(tissue, snp) %in% W0[, paste(tissue, snp)]]); stopifnot(!anyDuplicated(W[, .(tissue, snp)]))
## GTEx top variant per gene
EL <- rbind(fread(file.path(LW, "eqtl_leads.tsv"))[, .(tissue, gene_base, eqtl_top = eqtl_lead)], fread(file.path(O, "new_inputs/egene_top.tsv"))[, .(tissue, gene_base, eqtl_top = eqtl_lead)])
EL <- unique(EL, by = c("tissue", "gene_base"))
P <- merge(N[, .(tissue, promoter, type, k, gene_base, chr, signal_snp, signal_pos, egene)], EL, by = c("tissue", "gene_base"), all.x = TRUE)
P <- merge(P, W[, .(tissue, signal_snp = snp, p_partners = n_partners, p_start = w_start, p_end = w_end)], by = c("tissue", "signal_snp"), all.x = TRUE)
P <- merge(P, W[, .(tissue, eqtl_top = snp, e_rare = rare, e_partners = n_partners, e_start = w_start, e_end = w_end)], by = c("tissue", "eqtl_top"), all.x = TRUE)
stopifnot(nrow(P) == nrow(N), !anyNA(P$eqtl_top), !anyNA(P$p_start), !anyNA(P$e_start))
P[, `:=`(int_start = pmax(p_start, e_start), int_end = pmin(p_end, e_end), uni_start = pmin(p_start, e_start), uni_end = pmax(p_end, e_end))]
P[, `:=`(eqtl_top_pos = pos_of(eqtl_top), overlap_bp = pmax(0L, int_end - int_start + 1L))]; P[, lead_dist := abs(eqtl_top_pos - signal_pos)]
## step-1 consistency: rows built in step 1 must agree
P1 <- fread(cmd = paste("zcat", file.path(O, "pairs_top_new.tsv.gz")))
ck <- merge(P1[, .(tissue, promoter, type, k, a = uni_start, b = uni_end)], P[, .(tissue, promoter, type, k, c = uni_start, d = uni_end)], by = c("tissue", "promoter", "type", "k"))
cat("top pairs: step-1 rows", nrow(P1), "identical unions", ck[a == c & b == d, .N], "\n")
fwrite(P, file.path(O, "pairs_top_new_all.tsv.gz"), sep = "\t")
## clumped leads: existing (sens_clump pairs) + new (NSCC step 3 clump_out)
CL0 <- unique(fread(cmd = paste("zcat", file.path(LW, "sens_clump/pairs.tsv.gz")))[, .(tissue, gene_base, lead, rank, lead_p, n_absorbed, n_partners, l_start, l_end)])
cg <- fread(file.path(O, "new_inputs/clump_genes.tsv"))
CL1 <- rbindlist(lapply(unique(cg$tissue), function(t) { l <- fread(file.path(NS, "clump_out", t, "leads.tsv"))[, tissue := t]
  q <- fread(file.path(NS, "clump_out", t, "lead_partners.tsv"))[, tissue := t]; merge(l[!is.na(lead)], q[, .(tissue, lead, n_partners, min_partner_pos, max_partner_pos)], by = c("tissue", "lead"), all.x = TRUE) }))
stopifnot(!anyNA(CL1$n_partners))
CL1[, pos := pos_of(lead)]; CL1[n_partners == 0, `:=`(min_partner_pos = pos, max_partner_pos = pos)]
CL1[, `:=`(l_start = pmin(pos, min_partner_pos) - PAD, l_end = pmax(pos, max_partner_pos) + PAD)]
cat("new clumped genes:", uniqueN(CL1[, .(tissue, gene_base)]), "of", nrow(cg), "; leads", nrow(CL1), "\n")
CL <- rbind(CL0, CL1[, .(tissue, gene_base, lead, rank, lead_p, n_absorbed, n_partners, l_start, l_end)])
X <- merge(P[egene == TRUE, .(tissue, promoter, type, k, gene_base, chr, signal_snp, signal_pos, eqtl_top, p_start, p_end)], CL, by = c("tissue", "gene_base"), allow.cartesian = TRUE)
X[, `:=`(int_start = pmax(p_start, l_start), int_end = pmin(p_end, l_end), uni_start = pmin(p_start, l_start), uni_end = pmax(p_end, l_end))]
X[, `:=`(overlap_bp = pmax(0L, int_end - int_start + 1L), lead_dist = abs(pos_of(lead) - signal_pos))]
setorder(X, tissue, promoter, type, k, rank)
cat("eGene new rows:", P[egene == TRUE, .N], "; with clumped pairs", uniqueN(X[, .(tissue, promoter, type, k)]), "; pairs", nrow(X), "\n")
fwrite(X, file.path(O, "pairs_clump_new_all.tsv.gz"), sep = "\t")
## eQTL extraction windows for genes not (fully) covered by the existing subsets
EW <- list(top = fread(file.path(LW, "eqtl_windows.tsv")), clump = fread(file.path(LW, "sens_clump/eqtl_windows.tsv")))
sp <- list(top = P[, .(start = min(uni_start), end = max(uni_end)), by = .(tissue, gene_base, chr)],
           clump = X[, .(start = min(uni_start), end = max(uni_end)), by = .(tissue, gene_base, chr)])
for (s in names(sp)) { cv <- merge(sp[[s]], EW[[s]], by = c("tissue", "gene_base", "chr"), all.x = TRUE, suffixes = c("", "_have"))
  ex <- cv[is.na(start_have) | start < start_have | end > end_have, .(tissue, gene_base, chr, start, end)]
  cat("eQTL extraction (", s, "): tissue x gene to extract", nrow(ex), "\n"); fwrite(ex, file.path(O, paste0("eqtl_windows_", s, "_new.tsv")), sep = "\t") }
