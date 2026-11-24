## mirror_abs_labels.R: ONE-TIME step that decides the boundary label of every exact mirror pair by absolute
## activity and SAVES it: only the mirror pairs need this judgement; it runs once and saves the labels, and an unclear
## case is recorded as p0. The unit step (units_v3b_from_labels.R) and every downstream
## script read the saved table; nothing refits during routine runs.
## Scope: REAL arm only; ALL primary, estimable boundary units on boundary_exact pairs (v3m: p0 on the silent TSS S, O = its
##   partner; p0(S) is the same per-SNP test as p1(O)), significant or not (two OLS per unit, genotypes read once
##   per tissue x chromosome). Permuted arms are not fitted: no permuted count depends on the label.
## Fit (same as mirror_v3_abs_rule_check.R, which it reproduces): OLS of each TSS's absolute activity A (proActiv log2 normalised
##   count, v11/tabs_att_filter_In.rds $abs, all freeze donors, 0-read donors = 0) on the unit's leader dosage, with the discovery
##   covariates (cov_matrix() of the engine: intercept, scaled genotype PCs, PEER 1..K) and, for significant units, the unit's
##   other backward-kept conditional signals (results/conditional_final_v3m); if that design is rank-deficient the conditioning
##   SNPs are dropped (flag cond_dropped). Non-significant units have no conditional run: leader + covariates only.
##   regulated = |z| >= 1.96.
## Label: S only -> p0(S); O only -> p1(O); both or neither (or not estimable) -> p0(S).
## Option B (one rule with Figure 3): every boundary event of an exact pair is TESTED as p0(S) and BH is
##   unchanged; the rule label is applied to SIGNIFICANT units only (q <= 0.05). Non-significant units keep p0(S); their fit
##   outcome is stored for reference only (column outcome). label_applied says which.
## Output: results/mirror_abs_labels/labels_real.tsv.gz (one row per judged unit), log.txt. Env COND_FINAL / LABEL_DIR choose the conditional
## table and the output folder (final run: COND_FINAL=results/conditional_final_v3c LABEL_DIR=results/mirror_abs_labels_v3c).
suppressPackageStartupMessages({ library(data.table); library(parallel) })
CODE <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/code"; source(file.path(CODE, "bohr2_data.R"))
OUT <- file.path(V2, Sys.getenv("LABEL_DIR", "results/mirror_abs_labels")); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
COND_FINAL <- Sys.getenv("COND_FINAL", "results/conditional_final_v3m")   # final run: results/conditional_final_v3c
NC <- as.integer(Sys.getenv("NC", "16")); ZC <- 1.96; t0 <- proc.time()[3]
z <- function(f) fread(cmd = paste("zcat", f))
U <- rbindlist(lapply(Sys.glob(file.path(V2, "results/final_nogate_v3m/freeze/*.real.units.rds")), function(f) as.data.table(readRDS(f))[
  type == "p0" & boundary_exact %in% TRUE & dedup_drop == FALSE & estimable == TRUE,
  .(tissue, S = promoter, O = mirror_partner, gene, leader, p, q, sig = q <= 0.05)]))
## the partner's p1 row (the same per-SNP test; its own +-500 kb window) for the record
P1 <- rbindlist(lapply(Sys.glob(file.path(V2, "results/final_nogate_v3m/freeze/*.real.units.rds")), function(f) as.data.table(readRDS(f))[
  type == "p1" & boundary_exact %in% TRUE, .(tissue, O = promoter, p_p1_O = p, leader_p1_O = leader, estimable_p1_O = estimable)]))
U <- merge(U, P1, by = c("tissue", "O"), all.x = TRUE)
stopifnot(!anyNA(U$O), !anyDuplicated(U[, .(tissue, S)]))
TIS <- Sys.getenv("TIS", ""); if (nzchar(TIS)) { U <- U[tissue %in% strsplit(TIS, " ")[[1]]]; OUT <- file.path(OUT, "test"); dir.create(OUT, showWarnings = FALSE) }   # smoke test only
cat("units to judge (real arm, primary estimable p0 on boundary_exact pairs):", nrow(U), "; significant:", U[sig == TRUE, .N], "\n")
S3 <- z(file.path(V2, COND_FINAL, "freeze_real_signals.tsv.gz"))[type == "p0" & backward_keep %in% TRUE, .(tissue, S = promoter, snp)]
abs_all <- readRDS(file.path(GTEX, "v11/tabs_att_filter_In.rds"))$abs
abs_all <- abs_all[, intersect(colnames(abs_all), unique(c(U$S, U$O))), drop = FALSE]; invisible(gc())
stopifnot(all(c(U$S, U$O) %in% colnames(abs_all)))
pos_of <- function(id) as.integer(vapply(strsplit(id, "_", fixed = TRUE), `[`, "", 2))
## genotypes of a set of SNP ids on one chromosome, one tabix query (same dosage + modal imputation as read_geno(), per variant)
read_ids <- function(chr, ids, samples) {
  gr <- reduce(GRanges(chr, IRanges(pos_of(ids), width = 1L)))
  v <- readVcf(TabixFile(file.path(PATHS$vcf, paste0(chr, ".vcf.gz"))), "hg38",
               ScanVcfParam(which = gr, samples = samples, fixed = NA, info = NA, geno = "GT"))
  gt <- geno(v)$GT[, samples, drop = FALSE]; gt <- gt[!duplicated(rownames(gt)) & rownames(gt) %in% ids, , drop = FALSE]
  D <- impute_modal(gt_dosage(gt)); t(D) }
fit <- function(X, y) { f <- lm.fit(X, y); p <- ncol(X); if (f$rank < p) return(c(b = NA, se = NA, z = NA, sdy = sd(y)))
  s2 <- sum(f$residuals^2) / (length(y) - p); V <- chol2inv(qr.R(f$qr)); j <- which(f$qr$pivot == p)
  b <- unname(f$coefficients[p]); se <- sqrt(V[j, j] * s2); c(b = b, se = se, z = b / se, sdy = sd(y)) }
one_tissue <- function(t) {
  td <- readRDS(file.path(V2, "work/final/freeze", paste0(t, ".rds"))); C <- cov_matrix(td, td$peer, td$K)
  u <- U[tissue == t]; A <- abs_all[td$row_id, , drop = FALSE]; rownames(A) <- td$U
  cs <- S3[tissue == t]; u[, chr := as.character(td$info[S, "seqnames"])]
  need <- unique(rbind(u[, .(chr, snp = leader)], merge(u[sig == TRUE, .(S, chr)], cs, by = "S")[, .(chr, snp)]))
  sp <- split(need$snp, need$chr); G <- do.call(cbind, unname(Map(function(ch, ids) read_ids(ch, ids, td$U), names(sp), sp)))
  res <- vector("list", nrow(u))
  for (i in seq_len(nrow(u))) { r <- u[i]
    if (!r$leader %in% colnames(G)) { res[[i]] <- data.table(tissue = t, S = r$S, status = "leader genotype missing"); next }
    cond <- if (isTRUE(r$sig)) intersect(setdiff(cs[S == r$S, snp], r$leader), colnames(G)) else character(0)
    dose <- as.numeric(G[, r$leader]); Xc <- cbind(C, G[, cond, drop = FALSE], dose); dropped <- FALSE
    if (qr(Xc)$rank < ncol(Xc)) { Xc <- cbind(C, dose); dropped <- length(cond) > 0 }
    fS <- fit(Xc, A[, r$S]); fO <- fit(Xc, A[, r$O])
    res[[i]] <- data.table(tissue = t, S = r$S, status = "ok", n_cond = if (dropped) 0L else length(cond), cond_dropped = dropped,
      cond_snps = if (dropped || !length(cond)) "" else paste(cond, collapse = ";"), n_donors = length(dose),
      z_S = fS[["z"]], z_O = fO[["z"]], b_S = fS[["b"]], b_O = fO[["b"]], se_S = fS[["se"]], se_O = fO[["se"]]) }
  rbindlist(res, fill = TRUE) }
R <- rbindlist(mclapply(sort(unique(U$tissue)), function(t) tryCatch(one_tissue(t), error = function(e) { message("ERROR ", t, ": ", conditionMessage(e)); NULL }),
                        mc.cores = NC, mc.preschedule = FALSE), fill = TRUE)
X <- merge(U, R, by = c("tissue", "S"), all.x = TRUE)
X[, regS := abs(z_S) >= ZC]; X[, regO := abs(z_O) >= ZC]
X[, outcome := fcase(is.na(z_S) | is.na(z_O), "not estimable", regS & !regO, "S only", regO & !regS, "O only", regS & regO, "both", default = "neither")]
X[, label_applied := sig %in% TRUE]
X[, `:=`(label_promoter = fifelse(label_applied & outcome == "O only", O, S), label_type = fifelse(label_applied & outcome == "O only", "p1", "p0"))]
X[, label_unit := paste0(label_promoter, ":", label_type)][, v3m_unit := paste0(S, ":p0")][, changed := label_unit != v3m_unit]
X[, `:=`(mirror_label_rule = "abs-activity, significant units only (author 2026-11-06, option B)", arm = "real")]
setcolorder(X, c("tissue", "arm", "S", "O", "gene", "leader", "sig", "p", "q", "z_S", "z_O", "outcome", "label_unit", "v3m_unit", "changed"))
fwrite(X, file.path(OUT, "labels_real.tsv.gz"), sep = "\t")
sink(file.path(OUT, "log.txt"), split = TRUE)
cat("mirror_abs_labels.R,", format(Sys.time()), "; cores", NC, "; conditioning signals from", COND_FINAL, "; output", OUT, "\n")
cat("units judged:", nrow(X), "; fits ok:", X[status == "ok", .N], "; not ok:", X[is.na(status) | status != "ok", .N], "\n")
cat("\noutcome by significance:\n"); print(dcast(X[, .N, by = .(sig, outcome)], outcome ~ sig, value.var = "N", fill = 0))
cat("\nlabel differences vs v3m (significant units, O only -> p1(O)):", X[changed == TRUE, .N], "; non-significant units with outcome O only (reference, label stays p0(S)):", X[sig == FALSE & outcome == "O only", .N], "\n")
print(X[changed == TRUE & sig == TRUE, .(tissue, v3m_unit, label_unit, z_S, z_O, n_cond, cond_snps, p, p_p1_O)])
cat("\nconditioning (significant units): with conditioning SNPs", X[sig == TRUE & n_cond > 0, .N], "; dropped for rank", X[cond_dropped == TRUE, .N], "\n")
## reproduction check against the earlier check table (same rule, same conditioning, significant units)
if (file.exists(file.path(V2, "results/mirror_v3/abs_rule_check/abs_rule_units.tsv"))) {   # only when that table exists (it does not in the final run)
OLD <- fread(file.path(V2, "results/mirror_v3/abs_rule_check/abs_rule_units.tsv"))[, .(tissue, S, oz_S = z_S, oz_O = z_O)]
M <- merge(X[sig == TRUE, .(tissue, S, z_S, z_O)], OLD, by = c("tissue", "S"))
cat("\nreproduces mirror_v3_abs_rule_check.R: matched", nrow(M), "of", X[sig == TRUE, .N], "significant; max |dz|",
    signif(M[, max(abs(c(z_S - oz_S, z_O - oz_O)), na.rm = TRUE)], 3), "\n") } else cat("\nno abs-rule check table: comparison skipped\n")
cat("elapsed:", round((proc.time()[3] - t0) / 60, 1), "min\n")
sink()
cat("DONE\n")
