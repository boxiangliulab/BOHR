## make_site_lists.R: per tissue and chromosome, the variants BOHR does not test, so that QTLtools tests the same variant set.
## BOHR tests a variant for a promoter only if each of the three genotype classes has at least 5 donors among the promoter's
## analysis samples (MIN_GENO = 5; run_final2.R line 63); QTLtools tested every variant of the MAF >= 0.05 VCF. To give QTLtools
## the same variant set, this script lists, per tissue and chromosome, the variants that FAIL that rule among the tissue's donors
## (the freeze object's U), for QTLtools --exclude-sites. Genotypes are read and missing calls imputed exactly as in bohr2_data.R
## (gt_dosage; modal dosage among the tissue's donors, ties to the smaller dosage). The rule is applied on all donors of the
## tissue, so it equals BOHR's for a promoter with no missing value and is slightly looser for promoters with missing values.
## Output: sites/<Tissue>.chr<N>.exclude.txt (variant IDs) and sites/summary.tsv (counts only). No genotype is written.
## Usage: Rscript make_site_lists.R <ncores>
suppressPackageStartupMessages({ library(VariantAnnotation); library(data.table); library(parallel) })
NC <- as.integer(commandArgs(TRUE)[1]); if (is.na(NC)) NC <- 11
W <- "/projects/group/bliu/e1327989/GTEx/nofloor_full"; V2 <- file.path(W, "GTEx/QTL_figure/rerun_v4")
VCF <- file.path(W, "GTEx/Rawdata/v9_vcf_MAF_0.05"); OUT <- file.path(W, "consult/fig2e_20261119/matched_variants/sites")
MIN_GENO <- 5L; BLOCK <- 2e6L
LEN <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622,
         133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468)   # GRCh38 chr1-22
tis <- sub("[.]rds$", "", list.files(file.path(V2, "work/final/freeze"), "[.]rds$"))
US <- lapply(setNames(tis, tis), function(t) readRDS(file.path(V2, "work/final/freeze", paste0(t, ".rds")))$U)
all_s <- sort(unique(unlist(US))); Mm <- sapply(US, function(u) as.numeric(all_s %in% u)); rownames(Mm) <- all_s   # donors x tissues
gt_dosage <- function(gt) { D <- matrix(NA_integer_, nrow(gt), ncol(gt), dimnames = dimnames(gt))
  D[gt == "0|0" | gt == "0/0"] <- 0L; D[gt == "0|1" | gt == "1|0" | gt == "0/1" | gt == "1/0"] <- 1L; D[gt == "1|1" | gt == "1/1"] <- 2L; D }
one_chr <- function(k) { chr <- paste0("chr", k); tf <- TabixFile(file.path(VCF, paste0(chr, ".vcf.gz")))
  ex <- lapply(tis, function(t) character(0)); names(ex) <- tis; n_all <- 0L; n_keep <- setNames(numeric(length(tis)), tis)
  for (from in seq(1L, LEN[k], by = BLOCK)) { to <- min(from + BLOCK - 1L, LEN[k])
    v <- readVcf(tf, "hg38", ScanVcfParam(which = GRanges(chr, IRanges(from, to)), samples = all_s, fixed = NA, info = NA, geno = "GT"))
    if (nrow(v) == 0) next
    st <- start(rowRanges(v)); own <- st >= from & st <= to                 # a record belongs to the block that holds its start
    if (!any(own)) next
    gt <- geno(v)$GT[own, all_s, drop = FALSE]; ids <- rownames(gt); D <- gt_dosage(gt); rm(gt, v)
    C0 <- (D == 0L); C0[is.na(C0)] <- FALSE; C1 <- (D == 1L); C1[is.na(C1)] <- FALSE; C2 <- (D == 2L); C2[is.na(C2)] <- FALSE; NA_ <- is.na(D)
    c0 <- C0 %*% Mm; c1 <- C1 %*% Mm; c2 <- C2 %*% Mm; nm <- NA_ %*% Mm       # variants x tissues
    md <- ifelse(c0 >= c1 & c0 >= c2, 0L, ifelse(c1 >= c2, 1L, 2L))            # modal class, ties to the smaller dosage
    c0 <- c0 + nm * (md == 0L); c1 <- c1 + nm * (md == 1L); c2 <- c2 + nm * (md == 2L)
    keep <- pmin(c0, c1, c2) >= MIN_GENO
    n_all <- n_all + length(ids); n_keep <- n_keep + colSums(keep)
    for (t in tis) ex[[t]] <- c(ex[[t]], ids[!keep[, t]]) }
  for (t in tis) writeLines(unique(ex[[t]]), file.path(OUT, paste0(t, ".", chr, ".exclude.txt")))
  data.table(chr = chr, tissue = tis, variants = n_all, kept = as.numeric(n_keep[tis])) }
R <- rbindlist(mclapply(1:22, one_chr, mc.cores = NC, mc.preschedule = FALSE))
S <- R[, .(variants = sum(variants), kept = sum(kept)), by = tissue][, kept_share := round(kept / variants, 4)]
S[, donors := lengths(US)[tissue]]; setorder(S, donors); fwrite(S, file.path(OUT, "summary.tsv"), sep = "\t"); fwrite(R, file.path(OUT, "summary_by_chr.tsv"), sep = "\t")
options(width = 160, datatable.print.nrows = 60); print(S); cat("DONE\n")
