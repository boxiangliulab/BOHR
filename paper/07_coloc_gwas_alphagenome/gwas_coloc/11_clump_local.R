## LD-window GWAS coloc, step 11 (one GWAS file per call; paths from the environment): GWAS leads by recursive clumping and
## their LD windows. Mirrors Figure 4's eQTL version (QTL_figure/rerun_v2/code/coloc_eqtl_ldwin/sens_clump/s2_clump_nscc.R).
##   1. The file's significant SNPs (step 10: P < 1e-5, typed in the unfiltered v9 bfile bfile_full_maf0.005); the donors of
##      CLUMP_DONORS (all 943 donors, the ld_v3 standard, or the 799 European donors for the European GWAS files).
##   2. r2 among them with PLINK 1.9 --r2 (the estimator of every LD window; --clump would use a haplotype-EM r2),
##      +-1 Mb, r2 >= 0.1.
##   3. Recursive clumping: sort by P (ties: larger |BETA/SE|, then ID); the most significant unassigned
##      SNP is a lead and takes every unassigned significant SNP with r2 >= 0.1; repeat until all are assigned.
##      Leads below a stricter threshold are exactly the leads of clumping at that threshold.
##   4. LD partners of every lead among the MAF >= 0.05 panel variants (option A'), the same donors, +-1 Mb, r2 >= 0.1.
##      Window (step 12, local) = [min(lead, partners) - 50 kb, max(lead, partners) + 50 kb].
## Output (OUTDIR/<trait>/): leads.tsv (chr, lead, gwas_id, lead_p, lead_absz, rank, n_sig_chr, n_absorbed, in_panel),
##   lead_partners.tsv (lead, pos, n_partners, min_partner_pos, max_partner_pos), DONE.
## Usage: Rscript 11_clump_local.R <trait_basename_without_.txt> <sig_dir> <outdir> <n_threads>   (CLUMP_PLINK, CLUMP_REF, CLUMP_DONORS)
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); TR <- a[1]; SIG <- a[2]; OUTD <- a[3]; NT <- as.integer(a[4])
PL  <- Sys.getenv("CLUMP_PLINK", "/home/users/nus/e1327989/bin/plink")   # the NSCC script of the earlier run, with its three paths taken from the environment
REF <- Sys.getenv("CLUMP_REF", "/home/users/nus/e1327989/YuhaoDong/GTEx/LD_ref_v9_plink19")
FULL <- file.path(REF, "bfile_full_maf0.005"); PAN <- file.path(REF, "bfile"); DON <- Sys.getenv("CLUMP_DONORS", file.path(REF, "donors", "ALL.txt"))
W <- file.path(OUTD, TR); dir.create(file.path(W, "tmp"), recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(W, "DONE"))) quit(save = "no")
S <- fread(cmd = paste("zcat", file.path(SIG, paste0(TR, ".tsv.gz"))), colClasses = list(character = c("chr", "gwas_id", "variant_id")))
S <- S[!is.na(variant_id) & variant_id != ""]
run <- function(args, log) { st <- system2(PL, args, stdout = log, stderr = log); if (st != 0) stop("plink failed: ", log) }
LEADS <- list(); PART <- list()
for (cn in unique(S$chr)) {
  cs <- paste0("chr", cn)
  Sc <- S[chr == cn][order(Pvalue, -absz, variant_id)]
  ## 1-2. r2 among the chromosome's significant typed SNPs (donors of CLUMP_DONORS)
  u <- file.path(W, "tmp", paste0(cs, ".sig.txt")); writeLines(unique(Sc$variant_id), u)
  go <- file.path(W, "tmp", paste0(cs, ".sig"))
  run(c("--bfile", file.path(FULL, cs), "--keep-allele-order", "--keep", DON, "--extract", u, "--r2", "--ld-window", "99999",
        "--ld-window-kb", "1000", "--ld-window-r2", "0.1", "--threads", NT, "--memory", "16000", "--out", go), paste0(go, ".stdout"))
  ld <- if (file.exists(paste0(go, ".ld"))) fread(paste0(go, ".ld"), select = c("SNP_A", "SNP_B", "R2")) else data.table(SNP_A = character(), SNP_B = character(), R2 = numeric())
  ld[, `:=`(SNP_A = as.character(SNP_A), SNP_B = as.character(SNP_B))]
  ld <- ld[SNP_A != SNP_B & R2 >= 0.1]
  ## a SNP absent from the bfile cannot be clumped (step 10 already restricts to typed SNPs; this is a guard)
  bimc <- fread(file.path(FULL, paste0(cs, ".bim")), header = FALSE, select = 2)$V2
  Sc[, in_bfile := variant_id %chin% bimc]
  x <- Sc[in_bfile == TRUE]                                            # already in clumping order
  ## 3. recursive clumping on integer indices (x row = P rank)
  ia <- match(c(ld$SNP_A, ld$SNP_B), x$variant_id); ib <- match(c(ld$SNP_B, ld$SNP_A), x$variant_id)
  ok <- !is.na(ia) & !is.na(ib); nbl <- split(ib[ok], factor(ia[ok], levels = seq_len(nrow(x))))
  assigned <- rep(FALSE, nrow(x)); lead_i <- integer(0); n_abs <- integer(0)
  for (i in seq_len(nrow(x))) { if (assigned[i]) next
    nb <- nbl[[i]]; nb <- nb[!assigned[nb]]
    assigned[c(i, nb)] <- TRUE; lead_i <- c(lead_i, i); n_abs <- c(n_abs, length(nb)) }
  stopifnot(all(assigned))
  Lc <- x[lead_i, .(chr = cn, lead = variant_id, gwas_id, lead_p = Pvalue, lead_absz = absz, rank = seq_along(lead_i),
                    n_sig_chr = nrow(Sc), n_absorbed = n_abs, in_panel)]
  LEADS[[cs]] <- Lc
  ## 4. partners of every lead among the panel variants (same donors, +-1 Mb, r2 >= 0.1)
  if (nrow(Lc)) {
    pl <- file.path(W, "tmp", paste0(cs, ".leads.txt")); writeLines(Lc$lead, pl)
    pan <- fread(file.path(PAN, paste0(cs, ".bim")), header = FALSE, select = 2)$V2
    ex <- file.path(W, "tmp", paste0(cs, ".extract.txt")); writeLines(unique(c(pan, Lc$lead)), ex)
    po <- file.path(W, "tmp", paste0(cs, ".part"))
    run(c("--bfile", file.path(FULL, cs), "--keep-allele-order", "--keep", DON, "--extract", ex, "--ld-snp-list", pl, "--r2",
          "--ld-window", "99999", "--ld-window-kb", "1000", "--ld-window-r2", "0.1", "--threads", NT, "--memory", "16000", "--out", po),
        paste0(po, ".stdout"))
    P <- fread(paste0(po, ".ld"), select = c("SNP_A", "BP_A", "SNP_B", "BP_B", "R2"))
    P[, `:=`(SNP_A = as.character(SNP_A), SNP_B = as.character(SNP_B))]
    P <- P[SNP_A %chin% Lc$lead]
    P[, keep := SNP_B %chin% pan & SNP_B != SNP_A]                     # partners: panel variants only, not the lead itself
    Q <- P[, .(pos = BP_A[1], n_partners = sum(keep), min_partner_pos = if (any(keep)) min(BP_B[keep]) else NA_integer_,
               max_partner_pos = if (any(keep)) max(BP_B[keep]) else NA_integer_), by = .(lead = SNP_A)]
    miss <- setdiff(Lc$lead, Q$lead)                                   # a lead always pairs with itself; none should be missing
    if (length(miss)) Q <- rbind(Q, data.table(lead = miss, pos = NA_integer_, n_partners = NA_integer_, min_partner_pos = NA_integer_, max_partner_pos = NA_integer_))
    PART[[cs]] <- Q
    unlink(c(paste0(po, c(".ld", ".nosex")), ex, pl))
  }
  unlink(c(paste0(go, c(".ld", ".nosex")), u))
  cat(TR, cs, "significant typed", nrow(Sc), "; in bfile", nrow(x), "; leads", nrow(Lc), "; leads P<5e-8", Lc[lead_p < 5e-8, .N], "\n")
}
L <- rbindlist(LEADS); Q <- rbindlist(PART)
fwrite(L, file.path(W, "leads.tsv"), sep = "\t"); fwrite(Q, file.path(W, "lead_partners.tsv"), sep = "\t")
cat(TR, "leads", nrow(L), "; P<5e-8", L[lead_p < 5e-8, .N], "; P<1e-7", L[lead_p < 1e-7, .N], "; P<1e-6", L[lead_p < 1e-6, .N],
    "; leads not in the panel (A')", L[in_panel == FALSE, .N], "; leads with missing partner rows", Q[is.na(pos), .N], "\n")
writeLines(format(Sys.time()), file.path(W, "DONE"))
