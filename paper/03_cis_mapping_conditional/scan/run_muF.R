## Amendment 7 record: per-SNP mu p-values before AND after the small-sample F correction, kept side by side so the
## correction can be revised later without touching the main run. For every promoter with interior values the design is
## rebuilt exactly as in run_final2.R (prepare_y -> orth_design) to get the interior-sample rank (p = rank + 1, phi not
## counted); no genotypes are needed because W (mu_stat) is read from the main run.
##   mu_chi_logp = main run's chi2(1) log p (unchanged);  mu_F_logp = log P(F(1, n_mid - p) > (n_mid - p) expm1(W / n_mid))
## Usage: Rscript run_muF.R <ncores> [def/tissue/arm ...]   (default: every complete arm under results/final)
## Output: results/final_muF/<def>/<tissue>/<arm>/<chr>.rds = list(snp = (promoter, snp, mu_chi_logp, mu_F_logp),
##         prom = (promoter, n_mid, rank_mid, ncol_design, p_mean), meta)
suppressPackageStartupMessages(library(parallel))
CODE <- Sys.getenv("BOHR2_CODE", "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/code")
source(file.path(CODE, "bohr2_data.R"))
a <- commandArgs(TRUE); NC <- as.integer(a[1])
FR <- file.path(V2, "results", "final"); OUTR <- file.path(V2, "results", "final_muF")
arms <- if (length(a) > 1) a[-1] else { d <- Sys.glob(file.path(FR, "*", "*", "*")); sub(paste0(FR, "/"), "", d) }
arms <- arms[vapply(arms, function(x) all(file.exists(file.path(FR, x, sprintf("chr%d.rds", 1:22)))), TRUE)]
one_arm <- function(x) { parts <- strsplit(x, "/")[[1]]; def <- parts[1]; tis <- parts[2]; arm <- parts[3]
  td <- readRDS(file.path(V2, "work", "final", def, paste0(tis, ".rds")))
  if (arm == "real") { A <- td$act; P <- td$peer } else {
    m <- fread(file.path(V2, "perm_maps", "final", def, paste0(tis, ".", arm, ".perm_map.tsv")))
    stopifnot(setequal(m$slot, td$U), setequal(m$donor, td$U), !anyDuplicated(m$donor))
    d <- setNames(m$donor, m$slot)[td$U]; A <- td$act[d, , drop = FALSE]; P <- td$peer[d, , drop = FALSE]; rownames(A) <- rownames(P) <- td$U }
  C_all <- cov_matrix(td, P, td$K); n_def <- 0L
  for (chr in sprintf("chr%d", 1:22)) {
    out <- file.path(OUTR, def, tis, arm, paste0(chr, ".rds")); if (file.exists(out)) next
    M <- readRDS(file.path(FR, def, tis, arm, paste0(chr, ".rds")))
    pr <- M$prom[is.na(skip) & !is.na(family)]
    PR <- rbindlist(lapply(pr$promoter, function(p) { z <- prepare_y(A[, p], ncol(C_all) - 1L); if (!is.na(z$skip)) return(NULL)
      y <- unname(z$y); Cs <- orth_design(C_all[names(z$y), , drop = FALSE]); mid <- y > 0 & y < 1
      rk <- if (sum(mid) > 0) qr(Cs[mid, , drop = FALSE])$rank else NA_integer_
      data.table(promoter = p, n_mid = sum(mid), rank_mid = rk, ncol_design = ncol(Cs), p_mean = rk + 1L) }))
    if (nrow(PR)) n_def <- n_def + PR[rank_mid < ncol_design, .N]
    S <- M$snp[, .(promoter, snp, W = mu_stat, mu_chi_logp = mu_lrt_logp)]
    if (nrow(PR)) S <- merge(S, PR[, .(promoter, n_mid, p_mean)], by = "promoter", all.x = TRUE) else S[, `:=`(n_mid = NA_integer_, p_mean = NA_integer_)]
    S[, mu_F_logp := fifelse(is.finite(W) & n_mid > p_mean + 1L, pf((n_mid - p_mean) * expm1(pmax(W, 0) / n_mid), 1, n_mid - p_mean, lower.tail = FALSE, log.p = TRUE), NA_real_)]
    dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
    tmp <- paste0(out, ".tmp"); saveRDS(list(snp = S[, .(promoter, snp, mu_chi_logp, mu_F_logp)], prom = PR,
      meta = list(def = def, tissue = tis, arm = arm, chr = chr, formula = "F = (n_mid - p) * expm1(W / n_mid) ~ F(1, n_mid - p); p = interior rank + 1",
                  date = format(Sys.time()))), tmp); file.rename(tmp, out) }
  data.table(arm = x, rank_deficient_promoters = n_def) }
R <- rbindlist(mclapply(arms, function(x) tryCatch(one_arm(x), error = function(e) data.table(arm = x, rank_deficient_promoters = NA_integer_, err = conditionMessage(e))),
                        mc.cores = NC, mc.preschedule = FALSE), fill = TRUE)
nerr <- if ("err" %in% names(R)) sum(!is.na(R$err)) else 0L
print(R[order(-rank_deficient_promoters)][1:10]); cat("arms", nrow(R), " errors", nerr, " rank-deficient promoter-chromosome cases", sum(R$rank_deficient_promoters, na.rm = TRUE), "\n")
