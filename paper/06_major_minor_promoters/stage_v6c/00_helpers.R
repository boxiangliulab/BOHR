## Shared helpers for the SNP-level "regulatory route" analyses (2026-11-03).
## Gene-clustered bootstrap for ratio-of-sums statistics: genes are resampled with replacement, all signals of a
## resampled gene come along (signals of one gene across tissues and types are not independent).
suppressPackageStartupMessages(library(data.table))
S   <- Sys.getenv("F5S", "/home/e1327989/f5_stage_20261107_v6c")
IN  <- file.path(S, "in")
OUTROOT <- file.path(S, "out")

## d: data.table(gene, group, num, den). Returns per-gene matrices (genes x groups).
gene_mats <- function(d, groups) {
  genes <- sort(unique(d$gene))
  a <- d[group %in% groups, .(num = sum(num), den = sum(den)), by = .(gene, group)]
  NUM <- matrix(0, length(genes), length(groups), dimnames = list(genes, groups)); DEN <- NUM
  ii <- cbind(match(a$gene, genes), match(a$group, groups))
  NUM[ii] <- a$num; DEN[ii] <- a$den
  list(NUM = NUM, DEN = DEN, genes = genes)
}
boot_weights <- function(G, B = 2000L, seed = 1L) {
  set.seed(seed); W <- matrix(0, B, G)
  for (b in seq_len(B)) W[b, ] <- tabulate(sample.int(G, G, replace = TRUE), G)
  W
}
## ratio-of-sums per group with percentile CIs; pairwise differences with two-sided bootstrap P
ratio_table <- function(d, groups, B = 2000L, seed = 1L, label = "") {
  m <- gene_mats(d, groups); W <- boot_weights(length(m$genes), B, seed)
  est <- colSums(m$NUM) / colSums(m$DEN)
  bs  <- (W %*% m$NUM) / (W %*% m$DEN)
  one <- data.table(label = label, group = groups, n_signals = as.numeric(colSums(m$DEN)), estimate = est,
                    lo = apply(bs, 2, quantile, 0.025, na.rm = TRUE), hi = apply(bs, 2, quantile, 0.975, na.rm = TRUE))
  pr <- if (length(groups) > 1) rbindlist(combn(groups, 2, simplify = FALSE, FUN = function(g) {
    dd <- bs[, g[1]] - bs[, g[2]]; e <- est[g[1]] - est[g[2]]
    p <- 2 * min(mean(dd <= 0, na.rm = TRUE), mean(dd >= 0, na.rm = TRUE))
    data.table(label = label, contrast = paste(g[1], "-", g[2]), difference = e,
               lo = quantile(dd, 0.025, na.rm = TRUE), hi = quantile(dd, 0.975, na.rm = TRUE), p_boot = min(1, max(p, 1 / B)))  # capped at 1 (ties at 0 counted in both tails); 1/B = resolution floor
  })) else NULL
  list(groups = one, contrasts = pr, boot = bs, est = est)
}
## statistic = sum(num_a)/sum(den_a) / (sum(num_b)/sum(den_b)) etc. can be built from $boot of ratio_table.
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
pct <- function(x, d = 1) paste0(formatC(100 * x, format = "f", digits = d), "%")
nlog10p <- function(lnp) -lnp / log(10)
writeT <- function(x, f) fwrite(x, f, sep = "\t", quote = FALSE)
## fread of .gz without R.utils (not installed in gtexfig): stream through zcat
fr <- function(f, ...) if (grepl("\\.gz$", f)) fread(cmd = sprintf("zcat '%s'", f), ...) else fread(f, ...)
