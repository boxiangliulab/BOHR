#!/usr/bin/env Rscript
# 2026-10-09: two European files without GRCh38 coordinates, prepared to the pipeline format
# (SNPID_hg19 SNPID_hg38 BETA SE Pvalue A1FREQ R2 ALTfreq; the QC of the other converters: drop missing values, SE > 0).
#  * Rheumatoid arthritis, Okada et al. 2014, European meta-analysis (RA_GWASmeta_European_v2.txt.gz, the authors' site;
#    14,361 cases / 43,923 controls): GRCh37 -> GRCh38 with the UCSC chain (rtracklayer::liftOver, unique mappings only);
#    BETA = log(OR), SE = (log(CIup) - log(CIlow)) / (2 x 1.96), as convert_gwas_catalog(use_log_or = TRUE) does.
#  * Adiponectin, Dastani et al. 2012 (ADIPOGen), European discovery meta-analysis (GWAS Catalog GCST001465, 29,347):
#    the file carries rsID and chromosome only, so GRCh38 positions come from the rsID -> (hm_chrom, hm_pos) map of the
#    GWAS Catalog harmonised MAGIC fasting-glucose file (GCST90002232, GRCh38), required to agree on chromosome.
suppressPackageStartupMessages({ library(data.table); library(rtracklayer); library(GenomicRanges) })
W <- "/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/gwas_eur_20261009"; RAW <- file.path(W, "raw")
wr <- function(d, name) { of <- file.path(W, "prepared", name)
  fwrite(d[, .(SNPID_hg19, SNPID_hg38, BETA = signif(BETA, 6), SE = signif(SE, 6), Pvalue, A1FREQ, R2 = NA, ALTfreq = NA)], of, sep = " ", na = "NA", quote = FALSE)
  cat(name, nrow(d), "variants\n") }
# ---- RA (Okada 2014, European) ----
ra <- fread(cmd = paste("zcat", file.path(RAW, "RA_GWASmeta_European_v2.txt.gz")))
setnames(ra, c("SNPID", "chr", "pos", "A1", "A2", "OR", "lo", "hi", "P"))
n0 <- nrow(ra)
ra <- ra[!is.na(OR) & !is.na(lo) & !is.na(hi) & !is.na(P) & OR > 0 & lo > 0 & hi > 0]
# The authors' European file gives OR and its CI to 2 decimals and P to 2 significant digits. An SE from the rounded CI
# carries up to ~10% error, so SE is taken from P and BETA (SE = |BETA| / z_P), which keeps each variant's Z as precise
# as its P; where that is undefined (OR rounded to 1.00, or P >= 0.99) the CI-based SE is used (SE_source).
ra[, `:=`(BETA = log(OR), SE_ci = (log(hi) - log(lo)) / (2 * 1.96), z = qnorm(P / 2, lower.tail = FALSE))]
ra[, SE := fifelse(BETA != 0 & P < 0.99 & is.finite(z) & z > 0, abs(BETA) / z, SE_ci)]
ra[, SE_source := fifelse(BETA != 0 & P < 0.99 & is.finite(z) & z > 0, "P", "CI")]
ra <- ra[is.finite(SE) & SE > 0]
print(ra[, .N, by = SE_source]); print(ra[P < 5e-8, summary(SE / SE_ci)])
ch <- import.chain(file.path(RAW, "hg19ToHg38.over.chain"))
gr <- GRanges(paste0("chr", fifelse(as.character(ra$chr) == "23", "X", as.character(ra$chr))), IRanges(ra$pos, width = 1)); gr$i <- seq_len(nrow(ra))   # chr 23 = X
lo <- liftOver(gr, ch); n1 <- lengths(lo); ok <- which(n1 == 1); l1 <- unlist(lo[ok])
ra <- ra[l1$i][, `:=`(SNPID_hg19 = paste0(fifelse(as.character(chr) == "23", "X", as.character(chr)), ":", pos), SNPID_hg38 = paste0(sub("^chr", "", as.character(seqnames(l1))), ":", start(l1)))]
ra <- ra[!grepl("_|^Y|^M", SNPID_hg38)]
cat(sprintf("RA: %d rows; %d after QC; %d lifted uniquely\n", n0, length(gr), nrow(ra)))
ra[, `:=`(Pvalue = P, A1FREQ = NA)]
wr(ra, "Rheumatoid_Arthritis_Okada_2014_Nature_hg38.txt")
# ---- Adiponectin (Dastani 2012, ADIPOGen European) ----
ad <- rbindlist(lapply(Sys.glob(file.path(RAW, "ADIPOGen_EUR_part*.txt")), fread))
n0 <- nrow(ad); ad <- ad[!is.na(beta) & !is.na(se) & !is.na(pvalue) & se > 0]
mp <- fread(cmd = paste("zcat", file.path(RAW, "MAGIC_FG_EUR_GCST90002232.h.tsv.gz"), "| cut -f2-4"), col.names = c("rsid", "c38", "p38"))
mp <- unique(mp[grepl("^rs", rsid) & !is.na(p38)]); mp <- mp[!rsid %in% mp[duplicated(rsid), rsid]]   # unique rsid -> position
ad <- merge(ad, mp, by.x = "marker", by.y = "rsid"); ad <- ad[as.character(chr) == as.character(c38)]
ad[, `:=`(SNPID_hg19 = NA_character_, SNPID_hg38 = paste0(c38, ":", p38), BETA = beta, SE = se, Pvalue = pvalue, A1FREQ = maf)]
cat(sprintf("Adiponectin: %d rows; %d placed on GRCh38 by rsID (chromosome agrees)\n", n0, nrow(ad)))
wr(ad, "Adiponectin_Dastani_2012_PLoS_Genetics_hg38.txt")
