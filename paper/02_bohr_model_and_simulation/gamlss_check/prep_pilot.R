## Stage 0: read the shared inputs once, cache one compact object per pilot tissue, and check the
## sample set / design before any fitting. Output: work/tissue/<T>.rds, results/stage0_prep_checks.tsv
source("/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/code/bohr2_data.R")
setDTthreads(8)
t0 <- Sys.time()
pcs <- as.matrix(read.table(PATHS$pcs, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)[, paste0("PC", 1:N_PC)])
info <- fread(PATHS$info, data.table = FALSE); rownames(info) <- info[[1]]
proms_all <- unique(unlist(lapply(PILOT$tissue, function(t) unlist(lapply(CHRS, function(ch) {
  f <- file.path(PATHS$proms, sprintf("%s.%s_promoters.txt", t, ch)); if (file.exists(f)) readLines(f) })))))
hdr <- names(fread(PATHS$act, nrows = 0))
act_all <- fread(PATHS$act, select = c("sample", "SMTSD", intersect(hdr, proms_all)), data.table = FALSE)
act_all <- act_all[act_all$SMTSD %in% PILOT$tissue, , drop = FALSE]
pc <- grep("^PromoterID_", names(act_all), value = TRUE)
nonnum <- 0L
for (p in pc) if (!is.numeric(act_all[[p]])) { x <- suppressWarnings(as.numeric(act_all[[p]]))
  nonnum <- nonnum + sum(is.na(x) & !is.na(act_all[[p]])); act_all[[p]] <- x }
cat("activity read:", nrow(act_all), "rows,", length(pc), "promoter columns; non-numeric cells coerced:", nonnum,
    "; elapsed", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
vcf_s <- samples(scanVcfHeader(file.path(PATHS$vcf, "chr1.vcf.gz")))

rows <- list()
for (i in seq_len(nrow(PILOT))) {
  tis <- PILOT$tissue[i]; k_now <- PILOT$k_now[i]
  td <- build_tissue(tis, k_now, act_all, pcs, info)
  m0 <- fread(map_file(tis, "perm0"))
  stopifnot(setequal(m0$slot, td$U))                        # U == published perm slot set
  n <- length(td$U); A <- td$act
  Ks <- sapply(c("C1", "C2", "C3"), peer_k, k_now = k_now, n = n)
  C <- cov_matrix(td, td$peer, k_now)
  fam <- table(factor(vapply(colnames(A), function(p) {
    r <- prepare_y(A[, p], ncol(C) - 1L); if (is.na(r$skip)) r$family else paste0("skip_", r$skip) }, ""),
    levels = c("BEINF", "BEZI", "BEOI", "BE", "skip_n", "skip_prop_or_sd")))
  psd <- apply(td$peer[, seq_len(k_now), drop = FALSE], 2, sd)
  rows[[i]] <- data.table(tissue = tis, n_U = n, U_order_equals_perm0 = identical(m0$slot, td$U),
    U_in_vcf = all(td$U %in% vcf_s), n_dup_rows = td$n_dup_rows,
    n_prom_listed = ncol(A) + length(td$missing_promoter_cols), n_prom_found = ncol(A),
    n_near1 = sum(A > 1 - 5e-7 & A < 1, na.rm = TRUE), n_exact0 = sum(A == 0, na.rm = TRUE),
    n_exact1 = sum(A == 1, na.rm = TRUE),
    k_now = k_now, K_C1 = Ks[["C1"]], K_C2 = Ks[["C2"]], K_C3 = Ks[["C3"]],
    design_cols = ncol(C), design_rank = qr(C, tol = 1e-7)$rank, design_kappa = kappa(C, exact = TRUE),
    peer_sd_min = min(psd), peer_sd_max = max(psd),
    fam_BEINF = fam[["BEINF"]], fam_BEZI = fam[["BEZI"]], fam_BEOI = fam[["BEOI"]], fam_BE = fam[["BE"]],
    skip_n = fam[["skip_n"]], skip_prop_or_sd = fam[["skip_prop_or_sd"]])
  saveRDS(td, file.path(PATHS$work, "tissue", paste0(tis, ".rds")))
  print(rows[[i]])
}
out <- rbindlist(rows)
fwrite(out, file.path(V2, "results", "stage0_prep_checks.tsv"), sep = "\t")
cat("done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
