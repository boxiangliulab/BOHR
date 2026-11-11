## Real-data parity of tissue_means() (stage 3b-1; R with BOHR installed, or BOHR_LOAD_ALL=<package dir> for the working tree).
## The paper's tissue means (fit/v11_raw_fit_fit_results.txt, one row per promoter x tissue with >= 10 values) came from
## raw_fit() of code_orig/gamlss_functions_optimized.R:324-444 (gamlss, intercept only) on the v11 relative activity
## (v11/filter_5_recal_att_filter.RData), over the samples of the 944 donors with genotypes (the donors of the VCF; this gives
## the stored sample sizes in every tissue, and no sample of the CML cell line). This script runs tissue_means() on the same
## activity and samples and compares: the pairs and
## their sample sizes, which pairs get the sample mean, the means, the families, and the ROKU effect filter (|tissue mean -
## median of the tissues without a ROKU call| >= 0.10, roku_20261007/roku_effect_filter.R:25-31) with ROKU's stored calls held
## fixed. Prints aggregate numbers only (no sample or donor IDs). Writes the package table to <out_dir>.
## Usage: Rscript 11_tissue_means_parity.R <out_dir> <ncores> [lib] [tissue ...]
a <- commandArgs(TRUE); OUT <- a[1]; NC <- as.integer(a[2]); LIB <- if (length(a) > 2) a[3] else ""; TIS <- if (length(a) > 3) a[-(1:3)] else NULL
if (nzchar(Sys.getenv("BOHR_LOAD_ALL"))) pkgload::load_all(Sys.getenv("BOHR_LOAD_ALL"), quiet = TRUE) else if (nzchar(LIB)) library(BOHR, lib.loc = LIB) else library(BOHR)
options(width = 200)
V11 <- "/home/e1327989/project/GTEx/v11/filter_5_recal_att_filter.RData"
FIT <- "/projects/group/bliu/e1327989/GTEx/fit"
VCF <- "/home/e1327989/project/GTEx/Rawdata/v9_vcf_MAF_0.05/chr22.vcf.gz"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
e <- new.env(); load(V11, envir = e); x <- e[[ls(e)[1]]]; rm(e)
pc <- grep("^PromoterID_", names(x))
con <- pipe(paste("zcat", VCF, "| grep -m 1 '^#CHROM'")); donors <- strsplit(readLines(con), "\t")[[1]][-(1:9)]; close(con)
keep <- x$sample %in% donors & (if (is.null(TIS)) TRUE else x$SMTSD %in% TIS)
cat(sprintf("samples of the %d genotyped donors: %d of %d\n", length(donors), sum(keep), nrow(x)))
ph <- t(as.matrix(x[keep, pc])); tis <- x$SMTSD[keep]; colnames(ph) <- paste0("s", seq_len(ncol(ph)))     # anonymous sample IDs
rm(x); invisible(gc())
cat(sprintf("activity: %d promoters x %d samples in %d tissues\n", nrow(ph), ncol(ph), length(unique(tis))))
t0 <- proc.time()
tm <- tissue_means(ph, tissue = tis, ncores = NC)
el <- proc.time() - t0
cat(sprintf("tissue_means(): %.0f s wall, %.0f s CPU on %d cores\n", el[["elapsed"]], sum(el[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE), NC))
saveRDS(tm, file.path(OUT, "tissue_means_package.rds"))

## ---- stored gamlss table ----
S <- utils::read.delim(file.path(FIT, "v11_raw_fit_fit_results.txt"), stringsAsFactors = FALSE)
if (!is.null(TIS)) S <- S[S$tissue %in% TIS, ]
P <- tm[!is.na(tm$mean), ]
ks <- paste(S$tissue, S$promoter); kp <- paste(P$tissue, P$feature)
cat(sprintf("\npairs with a mean: stored %d (%d tissues), package %d (%d tissues); in both %d; stored only %d; package only %d (of which in tissues absent from the stored table: %d)\n",
            nrow(S), length(unique(S$tissue)), nrow(P), length(unique(P$tissue)), sum(ks %in% kp), sum(!ks %in% kp), sum(!kp %in% ks), sum(!kp %in% ks & !P$tissue %in% S$tissue)))
M <- merge(S, P, by.x = c("tissue", "promoter"), by.y = c("tissue", "feature"), suffixes = c(".st", ".pk"))
cat(sprintf("sample sizes equal: %d of %d\n", sum(M$sample_size == M$n), nrow(M)))
## ---- which pairs use the sample mean ----
M$st_rule <- ifelse(M$status == "Success", "fit (gamlss)", ifelse(M$reason.st == "all_same_value", "all_same_value",
                    ifelse(M$reason.st == "insufficient_diversity", "insufficient_diversity", "gamlss error")))
M$pk_rule <- ifelse(M$method == "fit", "fit", M$reason.pk)
cat("\nstored rule (rows) x package rule (columns):\n"); print(table(stored = M$st_rule, package = M$pk_rule))
sm_st <- M$status != "Success"; sm_pk <- M$method == "sample_mean"
cat(sprintf("sample mean in both %d; fitted in both %d; sample mean only in the stored table %d; only in the package %d\n",
            sum(sm_st & sm_pk), sum(!sm_st & !sm_pk), sum(sm_st & !sm_pk), sum(!sm_st & sm_pk)))
d <- abs(M$overall_mean - M$mean)
q <- function(v) { na <- sum(!is.finite(v)); v <- v[is.finite(v)]
  paste0(if (length(v)) sprintf("max %.3g; 99.9%% %.3g; 99%% %.3g; median %.3g", max(v), stats::quantile(v, 0.999), stats::quantile(v, 0.99), stats::median(v)) else "none",
         if (na) sprintf(" (%d not finite)", na) else "") }
cat(sprintf("|delta mean|, sample mean in both: %s\n", q(d[sm_st & sm_pk])))
fb <- !sm_st & !sm_pk
cat(sprintf("|delta mean|, fitted in both, gamlss converged (%d): %s\n", sum(fb & M$converged.st %in% TRUE), q(d[fb & M$converged.st %in% TRUE])))
cat(sprintf("|delta mean|, fitted in both, gamlss not converged (%d): %s\n", sum(fb & !M$converged.st %in% TRUE), q(d[fb & !M$converged.st %in% TRUE])))
cat(sprintf("|delta mu| (beta mean), fitted in both, gamlss converged: %s\n", q(abs(M$mu.st - M$mu.pk)[fb & M$converged.st %in% TRUE])))
cat(sprintf("|delta p0| and |delta p1|, fitted in both: max %.3g and %.3g\n", max(abs(M$p0.st - M$p0.pk)[fb]), max(abs(M$p1.st - M$p1.pk)[fb])))
cat(sprintf("|delta mean| where the rules differ: %s\n", q(d[sm_st != sm_pk])))
cat(sprintf("families equal among pairs fitted in both: %d of %d\n", sum(M$model_type[fb] == M$family[fb]), sum(fb)))
w <- M$st_rule != "fit (gamlss)" | M$pk_rule != "fit"
cat("\nrule differences by the number of values strictly between 0 and 1 (n_mid; 0-1 = few_mid) and gamlss convergence:\n")
print(stats::aggregate(list(pairs = d[w & sm_st != sm_pk]), list(stored = M$st_rule[w & sm_st != sm_pk], package = M$pk_rule[w & sm_st != sm_pk],
      n_mid = pmin(M$n_mid[w & sm_st != sm_pk], 4), gamlss_converged = M$converged.st[w & sm_st != sm_pk]), length))
print(stats::aggregate(list(max_abs_delta_mean = d[w & sm_st != sm_pk]), list(stored = M$st_rule[w & sm_st != sm_pk], package = M$pk_rule[w & sm_st != sm_pk]),
      function(v) signif(max(v, na.rm = TRUE), 3)))
ff <- which(M$pk_rule == "fit_failed")
if (length(ff)) { cat("package fit failures (the sample mean; stored rule, family, n, n_mid):\n")
  print(data.frame(stored = M$st_rule[ff], family = M$model_type[ff], n = M$n[ff], n_mid = M$n_mid[ff]), row.names = FALSE) }
big <- which(fb & d > 1e-3)
if (length(big)) { cat("fitted pairs with |delta mean| > 1e-3 (family, gamlss converged, n, n_mid):\n")
  print(utils::head(data.frame(family = M$family[big], gamlss_converged = M$converged.st[big], n = M$n[big], n_mid = M$n_mid[big],
                               mean_stored = signif(M$overall_mean[big], 5), mean_package = signif(M$mean[big], 5))[order(-d[big]), ], 10), row.names = FALSE) }
saveRDS(M, file.path(OUT, "tissue_means_vs_stored.rds"))

## ---- ROKU effect filter with ROKU's stored calls held fixed ----
r <- readRDS(file.path(FIT, "roku_20261007", "v11_tissue_specific_ROKU_results_delta010.rds"))
O <- r$outlier
wide <- function(tab, val, feat) { W <- matrix(NA_real_, nrow(O), ncol(O), dimnames = dimnames(O))
  i <- match(tab[[feat]], rownames(O)); j <- match(tab$tissue, colnames(O)); ok <- !is.na(i) & !is.na(j); W[cbind(i[ok], j[ok])] <- tab[[val]][ok]; W }
XS <- wide(S, "overall_mean", "promoter"); XP <- wide(tm, "mean", "feature")
delta <- function(X) { D <- matrix(NA_real_, nrow(O), ncol(O), dimnames = dimnames(O))
  for (i in seq_len(nrow(O))) { x <- X[i, ]; o <- O[i, ]
    out <- which(!is.na(o) & o != 0 & !is.na(x)); if (!length(out)) next
    nor <- which(!is.na(o) & o == 0 & !is.na(x)); if (length(nor) < 2) next
    D[i, out] <- x[out] - stats::median(x[nor]) }
  D }
DS <- delta(XS); DP <- delta(XP)
cat(sprintf("\nROKU (stored calls held fixed; %d promoters x %d tissues): calls with an effect: stored means %d, package means %d\n",
            nrow(O), ncol(O), sum(!is.na(DS)), sum(!is.na(DP))))
cat(sprintf("recomputed effects equal the stored delta_matrix: max |diff| %.3g\n", max(abs(DS - r$delta_matrix), na.rm = TRUE)))
both <- !is.na(DS) & !is.na(DP)
kS <- abs(DS) >= 0.10; kP <- abs(DP) >= 0.10
cat(sprintf("calls kept at |effect| >= 0.10: stored %d, package %d; kept by one only: %d (stored only %d, package only %d); max |delta effect| %.3g\n",
            sum(kS[both]), sum(kP[both]), sum(kS[both] != kP[both]), sum(kS[both] & !kP[both]), sum(!kS[both] & kP[both]), max(abs(DS - DP)[both])))
if (any(kS[both] != kP[both])) { w <- which(both & kS != kP, arr.ind = TRUE)
  print(data.frame(effect_stored = signif(DS[w], 4), effect_package = signif(DP[w], 4))) }
