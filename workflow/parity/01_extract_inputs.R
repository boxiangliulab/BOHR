## Real-data parity, step A (original R + data.table): rebuild the inputs of the stored run
## results/final/freeze/<TIS>/real/<CHR>.rds exactly as run_final2.R built them, as plain matrices for bohr_cis().
## Output stays in the scratchpad (never in the package repository).
## Usage: Rscript extract_inputs.R <tissue> <chr> <out.rds>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); TIS <- a[1]; CHR <- a[2]; OUT <- a[3]
V2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"
CODE <- "/projects/group/bliu/e1327989/GTEx/TATA_Enhancer/_tssQTL_repo/Figure4_puQTL_landscape_eQTL_coloc/rerun_v2/code"
VCF <- file.path("/home/e1327989/project/GTEx/Rawdata/v9_vcf_MAF_0.05", paste0(CHR, ".vcf.gz"))
## gt_dosage() of bohr2_data.R (definition only; the file itself loads VariantAnnotation)
for (e in parse(file.path(CODE, "bohr2_data.R"), keep.source = FALSE))
  if (is.call(e) && identical(e[[1]], as.name("<-")) && identical(as.character(e[[2]]), "gt_dosage")) eval(e)
td <- readRDS(file.path(V2, "work", "final", "freeze", paste0(TIS, ".rds"))); K <- td$K          # run_final2.R:18
A <- td$act; P <- td$peer                                                                     # arm real (run_final2.R:19)
prs <- colnames(td$act)[td$info[colnames(td$act), "seqnames"] == CHR]                          # run_final2.R:24
prs <- prs[order(td$info[prs, "start"])]
## covariates = cbind(PC1-5, PEER 1..K) (cov_matrix, bohr2_data.R:131-132); bohr_cis() standardises them
Z <- cbind(td$pcs, P[, seq_len(K), drop = FALSE])
## genotypes: every record of the chromosome for the analysis samples td$U; dosage as read_geno() (missing kept as NA:
## bohr_cis() fills them with the per-variant mode over the same samples, like impute_modal() in read_geno())
hdr <- names(fread(cmd = paste("zcat", VCF, "| grep -v '^##' | head -1")))
stopifnot(all(td$U %in% hdr))
v <- fread(cmd = paste("zcat", VCF, "| grep -v '^##'"), select = c("#CHROM", "POS", "ID", "REF", td$U), colClasses = "character")
gt <- as.matrix(v[, td$U, with = FALSE]); rownames(gt) <- v$ID
D <- gt_dosage(gt)
stopifnot(!anyDuplicated(v$ID))
variants <- data.frame(id = v$ID, chr = v[["#CHROM"]], pos = as.numeric(v$POS), end = as.numeric(v$POS) + nchar(v$REF) - 1, stringsAsFactors = FALSE)
promoters <- data.frame(id = prs, chr = CHR, pos = td$info[prs, "start"], gene = td$info[prs, "geneId"], stringsAsFactors = FALSE)
## canonical member of mirror pairs: pooled mean over all freeze samples of all 50 tissues (mirror_pairs.R:22-24, 40-41)
M <- fread(cmd = paste("zcat", file.path(V2, "results/mirror/mirror_pairs.tsv.gz")))[tissue == TIS]
major <- c(setNames(M$pooled_mean_1, M$promoter_1), setNames(M$pooled_mean_2, M$promoter_2))
saveRDS(list(phenotype = t(A[, prs, drop = FALSE]), genotype = D, variants = variants, promoters = promoters, covariates = Z,
             major_score = major, tissue = TIS, chr = CHR, K = K), OUT)
cat(sprintf("%s %s: %d promoters, %d samples, %d covariates (K = %d), %d variants (%d indel/MNP), %d missing genotypes, %d mirror pairs in tissue\n",
            TIS, CHR, length(prs), length(td$U), ncol(Z), K, nrow(D), sum(variants$end > variants$pos), sum(is.na(D)), nrow(M)))
