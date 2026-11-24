## Figure 4 panel e, LD-window coloc, sensitivity "eQTL leads by recursive clumping" (2026-10-12),
## step s2 (NSCC, one tissue per call): for every eGene x tissue of the analysis set,
##   1. its GTEx significant eQTL variants (s1) typed in the unfiltered v9 bfile (bfile_full_maf0.005), tissue donors;
##   2. r2 among them with PLINK 1.9 --r2 (the estimator of every LD window; PLINK --clump uses a different, haplotype-EM
##      r2, so the clumping is done here on --r2 pairs), +-1 Mb, r2 >= 0.1;
##   3. recursive clumping (the rule used for the GWAS leads): sort by P, the most significant unassigned variant
##      is a lead and takes every unassigned variant with r2 >= 0.1; repeat until all are assigned;
##   4. LD window of every lead as in the main analysis (option A'): r2 >= 0.1 partners restricted to the MAF >= 0.05
##      panel, +-1 Mb, tissue donors; window = farthest partner +/- 50 kb (computed locally in s3).
## Output (OUTDIR/<tissue>/): leads.tsv (gene_base, chr, lead, lead_p, rank, n_sig, n_typed, n_absorbed),
##   lead_partners.tsv (lead, pos, n_partners, min_partner_pos, max_partner_pos).
## Usage: Rscript s2_clump_nscc.R <tissue> <sig_dir> <outdir> <n_cores>
suppressPackageStartupMessages({ library(data.table); library(parallel) })
a <- commandArgs(TRUE); TIS <- a[1]; SIG <- a[2]; OUTD <- a[3]; NC <- as.integer(a[4])
PL <- "/home/users/nus/e1327989/bin/plink"
REF <- "/home/users/nus/e1327989/YuhaoDong/GTEx/LD_ref_v9_plink19"
FULL <- file.path(REF, "bfile_full_maf0.005"); PAN <- file.path(REF, "bfile"); DON <- file.path(REF, "donors", paste0(TIS, ".txt"))
W <- file.path(OUTD, TIS); dir.create(file.path(W, "tmp"), recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(W, "DONE"))) quit(save = "no")
S <- fread(cmd = paste("zcat", file.path(SIG, paste0(TIS, ".tsv.gz"))))
run <- function(args, log) { st <- system2(PL, args, stdout = log, stderr = log); if (st != 0) stop("plink failed: ", log) }
LEADS <- list(); PART <- list()
for (cs in unique(S$chr)) {
  Sc <- S[chr == cs]; cn <- sub("chr", "", cs)
  ## 1. small bfile of the chromosome's significant variants, tissue donors
  u <- file.path(W, "tmp", paste0(cs, ".sig.txt")); writeLines(unique(Sc$variant_id), u)
  sb <- file.path(W, "tmp", paste0(cs, ".small"))
  run(c("--bfile", file.path(FULL, cs), "--keep-allele-order", "--keep", DON, "--extract", u, "--make-bed", "--threads", "1",
        "--memory", "3000", "--out", sb), paste0(sb, ".stdout"))
  typed <- fread(paste0(sb, ".bim"), header = FALSE)$V2
  ## 2-3. per gene: --r2 among its typed significant variants, then recursive clumping
  one_gene <- function(g) {
    x <- Sc[gene_base == g][order(pval_nominal, variant_id)]; n_sig <- nrow(x)
    x <- x[variant_id %in% typed]
    if (!nrow(x)) return(data.table(gene_base = g, chr = cs, lead = NA_character_, lead_p = NA_real_, rank = NA_integer_,
                                    n_sig = n_sig, n_typed = 0L, n_absorbed = NA_integer_))
    gl <- file.path(W, "tmp", paste0(g, ".txt")); writeLines(x$variant_id, gl)
    go <- file.path(W, "tmp", g)
    run(c("--bfile", sb, "--keep-allele-order", "--extract", gl, "--r2", "--ld-window", "99999", "--ld-window-kb", "1000",
          "--ld-window-r2", "0.1", "--threads", "1", "--memory", "2000", "--out", go), paste0(go, ".stdout"))
    ld <- if (file.exists(paste0(go, ".ld"))) fread(paste0(go, ".ld"), select = c("SNP_A", "SNP_B", "R2")) else data.table(SNP_A = character(), SNP_B = character(), R2 = numeric())
    ld[, `:=`(SNP_A = as.character(SNP_A), SNP_B = as.character(SNP_B))]          # an empty table reads as logical
    ld <- ld[SNP_A != SNP_B & R2 >= 0.1]
    adj <- rbind(ld[, .(a = SNP_A, b = SNP_B)], ld[, .(a = SNP_B, b = SNP_A)]); setkey(adj, a)
    assigned <- setNames(rep(FALSE, nrow(x)), x$variant_id); out <- list(); r <- 0L
    for (i in seq_len(nrow(x))) { s <- x$variant_id[i]; if (assigned[[s]]) next
      nb <- adj[.(s), b, nomatch = 0L]; nb <- nb[!assigned[nb]]
      assigned[c(s, nb)] <- TRUE; r <- r + 1L
      out[[r]] <- data.table(gene_base = g, chr = cs, lead = s, lead_p = x$pval_nominal[i], rank = r, n_sig = n_sig, n_typed = nrow(x),
                             n_absorbed = length(nb)) }
    unlink(paste0(go, c(".ld", ".log", ".stdout", ".nosex"))); unlink(gl)
    rbindlist(out) }
  res <- mclapply(unique(Sc$gene_base), one_gene, mc.cores = NC)
  bad <- vapply(res, function(z) inherits(z, "try-error") || is.null(z), TRUE)
  if (any(bad)) stop(TIS, " ", cs, ": ", sum(bad), " genes failed, e.g. ", as.character(res[[which(bad)[1]]]))
  Lc <- rbindlist(res)
  LEADS[[cs]] <- Lc
  ## 4. partners of every lead among the panel variants (tissue donors, +-1 Mb, r2 >= 0.1)
  ld_list <- unique(na.omit(Lc$lead))
  if (length(ld_list)) {
    pl <- file.path(W, "tmp", paste0(cs, ".leads.txt")); writeLines(ld_list, pl)
    pan <- fread(file.path(PAN, paste0(cs, ".bim")), header = FALSE, select = 2)$V2
    ex <- file.path(W, "tmp", paste0(cs, ".extract.txt")); writeLines(unique(c(pan, ld_list)), ex)
    po <- file.path(W, "tmp", paste0(cs, ".part"))
    run(c("--bfile", file.path(FULL, cs), "--keep-allele-order", "--keep", DON, "--extract", ex, "--ld-snp-list", pl, "--r2",
          "--ld-window", "99999", "--ld-window-kb", "1000", "--ld-window-r2", "0.1", "--threads", "4", "--memory", "8000", "--out", po),
        paste0(po, ".stdout"))
    P <- fread(paste0(po, ".ld"), select = c("SNP_A", "BP_A", "SNP_B", "BP_B", "R2"))
    P <- P[SNP_A %in% ld_list]
    P[, keep := SNP_B %in% pan & SNP_B != SNP_A]                      # partners: panel variants only, not the lead itself
    Q <- P[, .(pos = BP_A[1], n_partners = sum(keep), min_partner_pos = if (any(keep)) min(BP_B[keep]) else NA_integer_,
               max_partner_pos = if (any(keep)) max(BP_B[keep]) else NA_integer_), by = .(lead = SNP_A)]
    miss <- setdiff(ld_list, Q$lead)                                  # a lead always pairs with itself; none should be missing
    if (length(miss)) Q <- rbind(Q, data.table(lead = miss, pos = NA_integer_, n_partners = NA_integer_, min_partner_pos = NA_integer_, max_partner_pos = NA_integer_))
    PART[[cs]] <- Q
    unlink(paste0(po, ".ld")); unlink(c(ex, pl))
  }
  unlink(paste0(sb, c(".bed", ".bim", ".fam"))); unlink(u)
  cat(TIS, cs, "genes", uniqueN(Sc$gene_base), "leads", sum(!is.na(Lc$lead)), "\n")
}
L <- rbindlist(LEADS); Q <- rbindlist(PART)
fwrite(L, file.path(W, "leads.tsv"), sep = "\t"); fwrite(Q, file.path(W, "lead_partners.tsv"), sep = "\t")
cat(TIS, "eGenes", uniqueN(L$gene_base), "; leads", sum(!is.na(L$lead)), "; genes with >= 2 leads", L[!is.na(lead), uniqueN(gene_base[rank >= 2])],
    "; significant variants not typed", sum(L[!duplicated(gene_base), n_sig - n_typed]), "\n")
writeLines(format(Sys.time()), file.path(W, "DONE"))
