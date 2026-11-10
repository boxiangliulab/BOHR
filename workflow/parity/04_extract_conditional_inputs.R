## Real-data parity of bohr_conditional(), step A (gtexfig R: data.table + VariantAnnotation): the inputs of the stored
## conditional runs for the significant units of one tissue (all chromosomes), with the stored cis units as the `cis` input.
## Genotypes: the VCF records overlapping each unit promoter's window [TSS - 5e5, TSS + 5e5 - 1], read with readVcf() as the
## original read_geno() did (missing calls kept; bohr_conditional() fills them per variant over the tissue's samples).
## Output (scratch only, never the repository): inputs as plain matrices + the stored units, threshold and tested pairs.
## Usage: Rscript 04_extract_conditional_inputs.R <tissue> <out.rds> [units.tsv: columns tissue, promoter, type]
suppressPackageStartupMessages({ library(data.table); library(VariantAnnotation) })
a <- commandArgs(TRUE); TIS <- a[1]; OUT <- a[2]; SEL <- if (length(a) > 2) a[3] else ""
V2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"
CODE <- "/projects/group/bliu/e1327989/GTEx/TATA_Enhancer/_tssQTL_repo/Figure4_puQTL_landscape_eQTL_coloc/rerun_v2/code"
VCFD <- "/home/e1327989/project/GTEx/Rawdata/v9_vcf_MAF_0.05"
for (e in parse(file.path(CODE, "bohr2_data.R"), keep.source = FALSE))                 # gt_dosage() only
  if (is.call(e) && identical(e[[1]], as.name("<-")) && identical(as.character(e[[2]]), "gt_dosage")) eval(e)
td <- readRDS(file.path(V2, "work/final/freeze", paste0(TIS, ".rds"))); K <- td$K
U <- as.data.table(readRDS(file.path(V2, "results/final_nogate_v3m/freeze", paste0(TIS, ".real.units.rds"))))
thr <- U[q <= 0.05, max(p)]                                                               # the tissue's discovery threshold
S <- U[q <= 0.05]
if (nzchar(SEL)) { w <- fread(SEL)[tissue == TIS]; S <- S[paste(promoter, type) %in% paste(w$promoter, w$type)] }
prs <- unique(S$promoter)
info <- td$info[prs, , drop = FALSE]
G <- list(); VV <- list()
for (ch in unique(as.character(info$seqnames))) {
  pos <- info$start[as.character(info$seqnames) == ch]
  gr <- reduce(GRanges(ch, IRanges(pmax(1L, pos - 500000L), pos + 500000L - 1L)))
  v <- readVcf(TabixFile(file.path(VCFD, paste0(ch, ".vcf.gz"))), "hg38", ScanVcfParam(which = gr, samples = td$U, fixed = NA, info = NA, geno = "GT"))
  gt <- geno(v)$GT[, td$U, drop = FALSE]; rr <- rowRanges(v)
  keep <- !duplicated(rownames(gt))
  G[[ch]] <- gt_dosage(gt[keep, , drop = FALSE])
  VV[[ch]] <- data.frame(id = rownames(gt)[keep], chr = ch, pos = start(rr)[keep], end = end(rr)[keep], stringsAsFactors = FALSE)
}
D <- do.call(rbind, G); variants <- do.call(rbind, VV); stopifnot(!anyDuplicated(variants$id), identical(rownames(D), variants$id))
## tested pairs of the stored main run (the original driver checks its tested set against them)
snp <- rbindlist(lapply(unique(as.character(info$seqnames)), function(ch) {
  x <- as.data.table(readRDS(file.path(V2, "results/final/freeze", TIS, "real", paste0(ch, ".rds")))$snp)
  x[promoter %in% prs, .(promoter, variant = snp)] }))
units <- S[, .(promoter, gene, component = type, family, n, n_mid, n0, n1, p, leader, leader_log_p = leader_logp, estimable, boundary_exact,
               interior_exact, mirror_role, mirror_partner = fifelse(mirror_partner == "", NA_character_, mirror_partner), dedup_drop, q,
               significant = TRUE)]
saveRDS(list(phenotype = t(td$act[, prs, drop = FALSE]), genotype = D, variants = variants,
             promoters = data.frame(id = prs, chr = as.character(info$seqnames), pos = info$start, gene = info$geneId, stringsAsFactors = FALSE),
             covariates = cbind(td$pcs, td$peer[, seq_len(K), drop = FALSE]), units = as.data.frame(units), threshold = thr,
             snp = as.data.frame(snp), tissue = TIS, K = K), OUT)
cat(sprintf("%s: %d significant units (%s) on %d promoters, threshold %.4g; %d samples, %d covariates; %d variants on %d chromosomes, %d missing calls\n",
            TIS, nrow(units), paste(names(table(units$component)), table(units$component), collapse = " "), length(prs), thr, ncol(t(td$act)),
            5 + K, nrow(D), length(G), sum(is.na(D))))
