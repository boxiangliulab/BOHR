## 184 (2026-11-11) -- P values for major vs minor TSS marks, on exactly the
## units of Fig. 5A / Source Data Fig5a (script 134: every tested (tissue, promoter) of a two-promoter gene with a TSS element, once;
## per unit = share of the tissue's matched ENCODE samples with the mark called, z > 1.64; major/minor = pooled usage, mirror_pairs).
## Marks: DNase, H3K27ac, H3K4me3; plus "H3K27ac only" (k27_balance: share of samples with H3K27ac and not H3K4me3 among samples with either).
## Tests, each two-sided, major vs minor:
##   (1) Welch t-test of the per-unit fractions (t, df, P; P also from log-scale pt so underflow is visible);
##   (2) gene-clustered bootstrap (B = 2000, seed 20261003 as 122a/134): genes resampled with replacement over the pooled major + minor
##       table (a gene's major and minor units travel together); statistic = mean(major) - mean(minor);
##       P = max(2 * min(share of b <= 0, share of b >= 0), 1/B) -> floor 5e-4 (convention of 71/110, Methods §9).
## Does not change any figure or Source Data. Output: out/PROMOTER_LEVEL_v4/184_fig5a_major_minor_pvalues.tsv
F5S <- Sys.getenv("F5S", "/home/e1327989/f5_stage_20261106_v6")
src <- readLines(file.path(F5S, "code/134_major_minor_pooled_chromatin.R"))
stop_at <- grep("^bt <- function", src) - 1
eval(parse(text = src[1:stop_at]))                     # 134 up to the merge with roles (MK, with 122a's B, SEED, gboot)
stopifnot(B == 2000, SEED == 20261003)
res <- rbindlist(lapply(c("DNase", "H3K27ac", "H3K4me3", "k27_balance"), function(m) {
  d <- MK[is.finite(get(m)), .(gene0, role, v = get(m))]
  x <- d[role == "major"]$v; y <- d[role == "minor"]$v
  tt <- t.test(x, y, alternative = "two.sided", var.equal = FALSE)
  t <- unname(tt$statistic); df <- unname(tt$parameter)
  log10p <- (log(2) + pt(-abs(t), df, log.p = TRUE)) / log(10)
  f <- function(z) mean(z$v[z$role == "major"]) - mean(z$v[z$role == "minor"])
  b <- gboot(d, f)
  data.table(mark = fifelse(m == "k27_balance", "H3K27ac only", m), n_major = length(x), n_minor = length(y),
             genes = uniqueN(d$gene0), mean_major = mean(x), mean_minor = mean(y), diff = mean(x) - mean(y),
             welch_t = t, welch_df = df, welch_P = tt$p.value, welch_P_underflow = tt$p.value == 0, welch_log10P = log10p,
             boot_B = length(b), boot_diff_lo = unname(quantile(b, 0.025)), boot_diff_hi = unname(quantile(b, 0.975)),
             boot_n_le0 = sum(b <= 0), boot_n_ge0 = sum(b >= 0), boot_P = max(2 * min(mean(b <= 0), mean(b >= 0)), 1 / length(b)))
}))
fwrite(res, file.path(OUTD, "184_fig5a_major_minor_pvalues.tsv"), sep = "\t")
options(width = 250, digits = 6); print(res)
cat("check vs Source Data Fig5a means:\n"); print(fread(file.path(OUTD, "134_major_minor_pooled_chromatin.tsv"))[kind == "mark", .(role, what, units, est)])
