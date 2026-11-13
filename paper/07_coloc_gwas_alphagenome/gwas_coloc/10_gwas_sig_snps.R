#!/usr/bin/env Rscript
# LD-window GWAS coloc (author 2026-10-12), step 10: the significant SNPs of one GWAS file,
# typed in GTEx v9, as input for recursive clumping on NSCC (step 11).
# - Significant = P < 1e-5, the loosest tier. Recursive clumping in P order gives the same leads below any stricter
#   threshold (a SNP's lead status depends only on more significant SNPs), so one run serves 5e-8 / 1e-7 / 1e-6 / 1e-5.
# - Typed = present in the unfiltered v9 bfile (bfile_full_maf0.005; option A'), matched by chr:pos (the GWAS files carry
#   no alleles). A position with several GTEx variants takes the one in the MAF >= 0.05 panel, then a SNV, then the
#   smallest ID.
# - Ties in P (incl. P underflow to 0) are broken by |BETA/SE| (larger first), then by ID, in step 11.
# Output: gwas_sig/<trait>.tsv.gz (gwas_id, chr, pos, Pvalue, absz, variant_id, in_panel, n_ids_at_pos) and
#         gwas_sig/<trait>.counts.tsv (per threshold: significant, typed, in panel).
# Usage: Rscript 10_gwas_sig_snps.R <trait_index>    (idempotent)
source("/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/scripts/coloc_final_lib.R")
D  <- file.path(NDIR, "coloc_ldunion_20261012")
tr <- TRAITS[as.integer(commandArgs(TRUE)[1])]
OF <- file.path(D, "gwas_sig", sub("\\.txt$", ".tsv.gz", tr))
if (file.exists(OF)) { cat("exists, skip\n"); quit(save = "no") }
dir.create(dirname(OF), showWarnings = FALSE)
P_MAX <- 1e-5
gw <- read_gwas(tr)[Pvalue < P_MAX]
gw[, c("chr", "pos") := tstrsplit(gwas_id, ":", fixed = TRUE)]; gw[, pos := as.integer(pos)]
gw <- gw[chr %in% as.character(1:22)]
gw[, absz := abs(BETA / SE)]
rd_bim <- function(dir) rbindlist(lapply(unique(gw$chr), function(cn)
  fread(file.path(D, "ref_copy_nscc", dir, paste0("chr", cn, ".bim")), header = FALSE, select = c(2, 4, 5, 6),
        col.names = c("variant_id", "pos", "a1", "a2"))[pos %in% gw[chr == cn, pos]][, chr := cn]))
B <- rd_bim("bim_full_maf0.005")
pan <- rbindlist(lapply(unique(gw$chr), function(cn)
  fread(file.path(D, "ref_copy_nscc", "bim_panel", paste0("chr", cn, ".bim")), header = FALSE, select = c(2, 4),
        col.names = c("variant_id", "pos"))[pos %in% gw[chr == cn, pos]]))$variant_id
B[, `:=`(in_panel = variant_id %chin% pan, snv = nchar(a1) == 1L & nchar(a2) == 1L)]
B[, n_ids_at_pos := .N, by = .(chr, pos)]
setorder(B, chr, pos, -in_panel, -snv, variant_id)
B <- B[!duplicated(B[, .(chr, pos)])]
X <- merge(gw[, .(gwas_id, chr, pos, Pvalue, absz)], B[, .(chr, pos, variant_id, in_panel, n_ids_at_pos)], by = c("chr", "pos"), all.x = TRUE)
X[is.na(in_panel), `:=`(in_panel = FALSE, n_ids_at_pos = 0L)]
setorder(X, Pvalue, -absz, gwas_id)
fwrite(X, OF, sep = "\t", quote = FALSE)
C <- rbindlist(lapply(c(5e-8, 1e-7, 1e-6, 1e-5), function(t) X[Pvalue < t, .(trait = tr, threshold = t, n_sig = .N,
  n_typed = sum(!is.na(variant_id)), n_panel = sum(in_panel), n_lowfreq = sum(!is.na(variant_id) & !in_panel), n_untyped = sum(is.na(variant_id)))]))
fwrite(C, sub("\\.tsv\\.gz$", ".counts.tsv", OF), sep = "\t", quote = FALSE)
print(C)
