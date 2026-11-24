## build_v4.R -- inputs of the final promoter rules (v4): no SD <= 0.05 skip; the filter SD >= 0.005 stays.
## Run after promoter_sets_rule.R (NF_SD_SKIP=none) and prep_freeze_nofloor.R into the lab staging root rerun_v4.
##  1. checks every v4 freeze object against the v2 one of the same tissue: same donors U, K, PCs, PEER, and the same
##     activity on every v2 promoter (the v4 set only adds promoters);
##  2. the delta to scan: promoters new in the v4 set, and v2-set promoters that the earlier prepare_y skipped for
##     SD <= 0.05 on the scan samples (the set rule looks at all samples of the tissue, the scan at one per donor);
##  3. writes the delta freeze objects (the v4 object with only the delta columns) to the root rerun_v4_delta, with the
##     permutation maps, so that run_final2.R / run_firth_addon.R scan exactly the delta, unchanged;
##  4. core/v4/delta_promoters.tsv (tissue, promoter, reason, v4 skip on the scan samples) and delta_summary.tsv.
## Usage: Rscript build_v4.R   (lab compute node)
suppressPackageStartupMessages(library(data.table))
W <- "/projects/group/bliu/e1327989/GTEx/nofloor_full"; G <- file.path(W, "GTEx/QTL_figure")
R2 <- file.path(G, "rerun_v2"); R4 <- file.path(G, "rerun_v4"); RD <- file.path(G, "rerun_v4_delta")
source(file.path(W, "code/bohr2_data.R"))                       # the earlier prepare_y (SD <= 0.05 skip) and cov_matrix
OUT <- file.path(W, "core/v4"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(RD, "work/final/freeze"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RD, "perm_maps/final/freeze"), recursive = TRUE, showWarnings = FALSE)
TIS <- sub("[.]rds$", "", list.files(file.path(R4, "work/final/freeze"), "[.]rds$"))
stopifnot(length(TIS) == 50)
rows <- list(); dl <- list()
for (t in TIS) {
  a <- readRDS(file.path(R2, "work/final/freeze", paste0(t, ".rds"))); b <- readRDS(file.path(R4, "work/final/freeze", paste0(t, ".rds")))
  stopifnot(identical(a$U, b$U), identical(a$row_id, b$row_id), a$K == b$K, isTRUE(all.equal(a$pcs, b$pcs)), isTRUE(all.equal(a$peer, b$peer)))
  stopifnot(all(colnames(a$act) %in% colnames(b$act)))
  d <- b$act[, colnames(a$act), drop = FALSE] - a$act; same_act <- all(is.na(d) == (is.na(a$act) | is.na(b$act[, colnames(a$act), drop = FALSE]))) && max(abs(d), na.rm = TRUE) == 0
  stopifnot(same_act)
  p_cov <- ncol(cov_matrix(b, b$peer, b$K)) - 1L
  why <- function(y) { y <- y[!is.na(y)]; n <- length(y)
    if (n < 20 || n < p_cov) return("n"); p0 <- mean(y == 0); p1 <- mean(y == 1)
    if (p0 >= 0.95 || p1 >= 0.95) return("share") ; if (sd(y) < 0.005) return("sd<0.005"); if (sd(y) <= 0.05) return("sd<=0.05"); "tested" }
  old_why <- vapply(colnames(a$act), function(p) why(b$act[, p]), "")
  new <- setdiff(colnames(b$act), colnames(a$act)); sd_old <- names(old_why)[old_why == "sd<=0.05"]
  delta <- c(new, sd_old)
  v4_why <- vapply(delta, function(p) { w <- why(b$act[, p]); if (w == "sd<=0.05") "tested" else w }, "")
  td <- b; td$act <- b$act[, delta, drop = FALSE]; td$info <- b$info[delta, , drop = FALSE]; td$def <- "freeze (v4 delta)"
  saveRDS(td, file.path(RD, "work/final/freeze", paste0(t, ".rds")))
  for (arm in c("perm0", "perm1", "perm2")) { f <- paste0(t, ".", arm, ".perm_map.tsv")
    file.copy(file.path(R2, "perm_maps/final/freeze", f), file.path(RD, "perm_maps/final/freeze", f), overwrite = TRUE) }
  dl[[t]] <- data.table(tissue = t, promoter = delta, reason = c(rep("new in the v4 set", length(new)), rep("v2 set, SD <= 0.05 on the scan samples", length(sd_old))), v4_scan = v4_why)
  rows[[t]] <- data.table(tissue = t, n = length(b$U), K = b$K, v2_promoters = ncol(a$act), v4_promoters = ncol(b$act), new = length(new),
                          v2_sd_skipped = length(sd_old), delta = length(delta), delta_tested_v4 = sum(v4_why == "tested"),
                          v2_tested_paper = sum(old_why == "tested"))
  cat(t, ": v2", ncol(a$act), " v4", ncol(b$act), " delta", length(delta), "(", sum(v4_why == "tested"), "tested )\n")
}
D <- rbindlist(dl); S <- rbindlist(rows)
fwrite(D, file.path(OUT, "delta_promoters.tsv"), sep = "\t"); fwrite(S, file.path(OUT, "delta_summary.tsv"), sep = "\t")
fwrite(S[, .(tissue, n, K)], file.path(RD, "work/final/tissues.tsv"), sep = "\t")
cat("\nall 50 tissues: v2 promoters", sum(S$v2_promoters), "| v4", sum(S$v4_promoters), "| delta", sum(S$delta), "(new", sum(S$new), ", v2 SD-skipped on scan samples", sum(S$v2_sd_skipped), ")",
    "| delta tested under v4", sum(S$delta_tested_v4), "\n")
print(D[, .N, by = .(reason, v4_scan)])
