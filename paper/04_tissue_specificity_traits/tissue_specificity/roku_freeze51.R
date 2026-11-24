## Tissue specificity (Fig. 3A,B; Supplementary Tables 3g,h): one mean per promoter and tissue, ROKU calls and the 0.10
## effect filter, with the tissue means taken over the GTEx v11 analysis-freeze samples (SMAFRZE = RNASEQ) of the
## genotyped donors. The paper uses the run freeze56 (56 tissues). The steps are those of the earlier run on all samples:
## BOHR::tissue_means() on the no-floor relative activity of the genotyped donors -> TCC::ROKU(upper.limit = 0.25,
## sort = TRUE) -> the 0.10 effect filter (delta_filter() below, the code of the earlier roku_effect_filter.R).
## Checks first:
##   C0  ROKU on the stored means gives the stored NSCC ROKU result (outlier matrix identical)
##   C1  tissue_means() with the earlier run's samples gives the stored means
## Runs: all56 (the earlier calls, recomputed), all51, freeze56 (the paper's calls) and freeze51. The 51-tissue runs
## leave out the three stomach LCM sub-tissues and the two cell lines. Usage: Rscript roku_freeze51.R <out_dir> <ncores>
suppressMessages({ library(data.table); library(tidyr); library(dplyr) })
a <- commandArgs(TRUE); OUT <- a[1]; NC <- as.integer(a[2])
W <- "/projects/group/bliu/e1327989/GTEx/nofloor_full"
ACT <- file.path(W, "GTEx/v11/filter_5_recal_att_filter.RData")
VCF <- file.path(W, "GTEx/Rawdata/v9_vcf_MAF_0.05/chr22.vcf.gz")
OLD <- file.path(W, "GTEx/fit")                       # the earlier run's means, ROKU calls and 0.10 filter (all samples of the genotyped donors)
DS  <- "/home/e1327989/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"
SRC <- "/projects/group/bliu/e1327989/GTEx/manuscript/BOHR"
DROP5 <- c("Stomach_Mixed_Cell", "Stomach_Mucosa", "Stomach_Muscularis", "Cells_EBV-transformed_lymphocytes", "Cells_Cultured_fibroblasts")
ROKU <- TCC::ROKU   # as run, this line sourced a copy of R/ROKU.R of TCC 1.52.0 (see README)
suppressMessages(pkgload::load_all(SRC, quiet = TRUE))
say <- function(...) cat(sprintf(...), "\n", sep = "")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ROKU with the input built as the earlier ROKU.R built it (read.table, pivot_wider in order of appearance)
roku_from_file <- function(f) {
  t <- read.table(f, sep = "\t", header = TRUE)
  tt <- t %>% pivot_wider(id_cols = promoter, names_from = tissue, values_from = overall_mean)
  df <- as.data.frame(tt); row.names(df) <- df$promoter; d <- df[, -1]
  ROKU(d, upper.limit = 0.25, sort = TRUE)
}
## the 0.10 filter, the code of the earlier roku_effect_filter.R (DELTA = 0.10), returning the same object structure
delta_filter <- function(r, f, DELTA = 0.10) {
  t <- fread(f, select = c("tissue", "promoter", "overall_mean"))
  Wd <- dcast(t, promoter ~ tissue, value.var = "overall_mean")
  X <- as.matrix(Wd[, -1]); rownames(X) <- Wd$promoter
  O <- r$outlier; X <- X[rownames(O), colnames(O)]
  K <- O; K[!is.na(K)] <- 0
  Dm <- matrix(NA_real_, nrow(O), ncol(O), dimnames = dimnames(O))
  for (i in seq_len(nrow(O))) {
    x <- X[i, ]; o <- O[i, ]
    out <- which(!is.na(o) & o != 0 & !is.na(x)); if (!length(out)) next
    nor <- which(!is.na(o) & o == 0 & !is.na(x)); if (length(nor) < 2) next
    d <- x[out] - median(x[nor]); Dm[i, out] <- d
    K[i, out] <- ifelse(abs(d) >= DELTA, o[out], 0)
  }
  res <- r; res$true_outliers <- K; res$delta_matrix <- Dm; res$p_values <- NULL
  res$filter <- sprintf("keep ROKU call if |overall_mean - median(non-outlier tissues)| >= %.2f (2026-10-07; replaces the z-test)", DELTA)
  res
}
cnt <- function(M) { M[is.na(M)] <- 0; c(up_calls = sum(M == 1), down_calls = sum(M == -1), promoters_any = sum(rowSums(M != 0) > 0),
  promoters_up = sum(rowSums(M == 1) > 0), promoters_down = sum(rowSums(M == -1) > 0), promoters_both = sum(rowSums(M == 1) > 0 & rowSums(M == -1) > 0)) }

## C0: ROKU on the stored means = the stored NSCC result (TCC::ROKU, earlier run)
r_old <- readRDS(file.path(OLD, "v11_tissue_specific_ROKU_results.rds"))
r0 <- roku_from_file(file.path(OLD, "v11_raw_fit_fit_results.txt"))
o_new <- r0$outlier[rownames(r_old$outlier), colnames(r_old$outlier)]
## The calls are identical, but modH differs by <= 9e-16 (R 4.5.3 here, R 4.4.3 on NSCC), which reorders
## near-tied ranks (1,810 ties in the stored rank) and so the sorted row order. Only the calls are used downstream, so C0 checks
## the calls (values and NA pattern, aligned by name) and reports the rest.
c0 <- identical(is.na(o_new), is.na(r_old$outlier)) && isTRUE(all.equal(unname(o_new), unname(r_old$outlier), tolerance = 0))
say("C0 ROKU calls (TCC::ROKU) on the stored means vs the stored NSCC result: %s (%d promoters x %d tissues); max |modH diff| %.2g, max |Tbw diff| %.2g; rank and row order %s",
    if (c0) "IDENTICAL" else "DIFFERENT", nrow(o_new), ncol(o_new), max(abs(r0$modH[names(r_old$modH)] - r_old$modH)),
    max(abs(r0$Tbw[names(r_old$Tbw)] - r_old$Tbw)), if (identical(rownames(r0$outlier), rownames(r_old$outlier))) "same" else "differ (near-ties)")
if (!c0) stop("C0 failed: the local ROKU does not reproduce the NSCC calls")
k_old <- readRDS(file.path(OLD, "roku_20261007/v11_tissue_specific_ROKU_results_delta010.rds"))
say("stored 0.10 filter: %s", paste(names(cnt(k_old$true_outliers)), cnt(k_old$true_outliers), sep = "=", collapse = " "))

## sample sets
e <- new.env(); load(ACT, envir = e); x <- e[["recal_att_f"]]; rm(e); invisible(gc())
pc <- grep("^PromoterID_", names(x))
con <- pipe(paste("zcat", VCF, "| grep -m 1 '^#CHROM'")); donors <- strsplit(readLines(con), "\t")[[1]][-(1:9)]; close(con)
att <- fread(DS, select = c("SAMPID", "SMAFRZE"), quote = "")
fz <- att[SMAFRZE == "RNASEQ", SAMPID]
say("activity rows %d; row names that are GTEx sample IDs: %d; freeze rows: %d", nrow(x), sum(rownames(x) %in% att$SAMPID), sum(rownames(x) %in% fz))
stopifnot(sum(rownames(x) %in% att$SAMPID) == nrow(x))
means <- function(keep, tag) {
  ph <- t(as.matrix(x[keep, pc])); tis <- x$SMTSD[keep]; colnames(ph) <- paste0("s", seq_len(ncol(ph)))
  say("[%s] samples of %d genotyped donors: %d; promoters %d; tissues %d", tag, length(donors), sum(keep), nrow(ph), length(unique(tis)))
  tm <- as.data.table(tissue_means(ph, tissue = tis, ncores = NC))
  tm[!is.na(mean), .(tissue, promoter = feature, model_type = family, status = fifelse(method == "fit", "Success", "Skip_fit"), reason,
                     mu, p0, p1, overall_mean = mean, sample_size = n, n0, n1, n_mid, converged)]
}
## C1: the earlier run's sample set reproduces the stored means
m_all <- means(x$sample %in% donors, "all")
st <- fread(file.path(OLD, "v11_raw_fit_fit_results.txt"), select = c("tissue", "promoter", "overall_mean", "sample_size"))
cmp <- merge(st, m_all[, .(tissue, promoter, om = overall_mean, n = sample_size)], by = c("tissue", "promoter"), all = TRUE)
say("C1 tissue_means() with the paper's samples vs the stored means: rows %d vs %d, unmatched %d, max |diff| %.3g, sample sizes equal: %s",
    nrow(st), nrow(m_all), cmp[is.na(om) | is.na(overall_mean), .N], cmp[, max(abs(om - overall_mean), na.rm = TRUE)], cmp[, all(n == sample_size, na.rm = TRUE)])
m_fz <- means(x$sample %in% donors & rownames(x) %in% fz, "freeze")

run <- function(m, tag) {
  d <- file.path(OUT, tag); dir.create(file.path(d, "roku_20261007"), recursive = TRUE, showWarnings = FALSE)
  f <- file.path(d, "v11_raw_fit_fit_results.txt"); fwrite(m, f, sep = "\t")
  r <- roku_from_file(f); saveRDS(r, file.path(d, "v11_tissue_specific_ROKU_results.rds"))
  k <- delta_filter(r, f); saveRDS(k, file.path(d, "roku_20261007", "v11_tissue_specific_ROKU_results_delta010.rds"))
  S <- rbind(roku_raw = cnt(r$outlier), delta010 = cnt(k$true_outliers))
  fwrite(data.table(set = rownames(S), S), file.path(d, "roku_20261007", "summary_delta_filter.tsv"), sep = "\t")
  as.data.table(c(list(run = tag, tissues = ncol(r$outlier), promoters_assessable = nrow(r$outlier)), as.list(cnt(k$true_outliers))))
}
SUM <- rbindlist(list(
  run(m_all, "all56"),
  run(m_all[!tissue %in% DROP5], "all51"),
  run(m_fz, "freeze56"),
  run(m_fz[!tissue %in% DROP5], "freeze51")))
fwrite(SUM, file.path(OUT, "summary_runs.tsv"), sep = "\t"); print(SUM)
## the earlier run recomputed must equal the stored filter output
k_all <- readRDS(file.path(OUT, "all56/roku_20261007/v11_tissue_specific_ROKU_results_delta010.rds"))
say("all56 0.10 calls identical to the stored ones: %s",
    identical(k_all$true_outliers[rownames(k_old$true_outliers), colnames(k_old$true_outliers)], k_old$true_outliers))
## per-tissue sample counts of the freeze run
ns <- data.table(tissue = x$SMTSD, s = x$sample %in% donors, f = x$sample %in% donors & rownames(x) %in% fz)[, .(all = sum(s), freeze = sum(f)), by = tissue]
fwrite(ns[order(tissue)], file.path(OUT, "samples_per_tissue.tsv"), sep = "\t")
say("samples (genotyped donors): all %d, freeze %d; freeze in the 51 tissues %d", ns[, sum(all)], ns[, sum(freeze)], ns[!tissue %in% DROP5, sum(freeze)])
say("DONE")
