## make_version_inputs.R: promoter rules of three versions; nofloor_reads5_E80 (r5e80) is the rule set of the final run,
## and its promoters_excluded table defines the units of the final run (build_r5e80_units.R).
## Three versions, all with the endpoint skip at 80% (a promoter is skipped in a tissue when its share of 0 or of 1 among
## the tissue's analysed samples is >= 0.80; the earlier rule was 0.95):
##   paper_E80           floor data (QTL_figure/rerun_v2)
##   nofloor_E80         no-floor data (nofloor_full/GTEx/QTL_figure/rerun_v4)
##   nofloor_reads5_E80  no-floor data, plus genes with median first-intron reads < 5 in the tissue dropped (whole gene)
## BH again within tissue on the kept tested units. Outputs per version, in figures_v4/inputs/<version>/:
##   promoters_excluded.tsv.gz      tissue, promoter, gene, reason, share0, share1, med_D
##   units_keep.tsv.gz              significant before and after (their signals, leads and coloc exist): tissue, promoter, type
##   units_new_topup.tsv.gz         significant only after the re-BH: need the conditional run and coloc (top-up)
##   units_lost.tsv.gz              significant before, kept, not significant after the re-BH
##   signals_keep.tsv.gz            delivery v7 signals of units_keep: tissue, promoter, type, k
##   fig2e/                         Fig. 2E tables rebuilt with the logic of fig2e_deliver_freeze_v3m_pairs.R on the kept
##                                  promoters (bohr/, qtltools/<mode>/, counts.tsv, counts_display.tsv)
##   summary.tsv
## Aggregates and identifiers only. Usage: Rscript make_version_inputs.R
suppressPackageStartupMessages({ library(data.table); library(parallel) })
G <- "/projects/group/bliu/e1327989/GTEx"; W <- file.path(G, "nofloor_full"); C <- file.path(W, "consult/version_choice_20261118/main")
OUTR <- file.path(W, "figures_v4/inputs"); dir.create(OUTR, FALSE, TRUE)
fr <- function(f, ...) fread(cmd = paste("zcat", shQuote(f)), ...)
TREE <- c(paper = file.path(G, "QTL_figure/rerun_v2"), v4 = file.path(W, "GTEx/QTL_figure/rerun_v4"))
VERS <- list(paper_E80 = list(tree = "paper", gene_cut = 0), nofloor_E80 = list(tree = "v4", gene_cut = 0),
             nofloor_reads5_E80 = list(tree = "v4", gene_cut = 5))
ECUT <- 0.80; ARMS <- c("real", "perm0", "perm1", "perm2"); MODES <- c("rintfirst", "raw", "normal")
tis <- sub("[.]real[.]tsv[.]gz$", "", list.files(file.path(TREE[["v4"]], "results/fig2e_freeze_v3m/bohr"), "[.]real[.]tsv[.]gz$"))
DS <- fr(file.path(C, "out/promoter_depth_v4.tsv.gz"))[, .(tissue, promoter, med_D)]
info <- unique(rbind(fread(file.path(W, "GTEx/v11/Promoter_info_filter_5.txt"))[, .(promoter = V1, gene = geneId)],
                     fread(file.path(G, "v11/Promoter_info_filter_5.txt"))[, .(promoter = V1, gene = geneId)]), by = "promoter")
GD <- merge(DS, info, by = "promoter")[, .(med_D = max(med_D, na.rm = TRUE)), by = .(tissue, gene)]
simes <- function(p) { p <- sort(p); k <- length(p); min(1, min(k * p / seq_len(k))) }

## promoter shares per tree and tissue (all promoters of the freeze act)
shares <- function(tr, t) { a <- readRDS(file.path(TREE[[tr]], "work/final/freeze", paste0(t, ".rds")))$act
  data.table(tissue = t, promoter = colnames(a), nfin = colSums(!is.na(a)), share0 = colSums(a == 0, na.rm = TRUE) / colSums(!is.na(a)),
             share1 = colSums(a == 1, na.rm = TRUE) / colSums(!is.na(a))) }
SH <- list(paper = rbindlist(mclapply(tis, shares, tr = "paper", mc.cores = 10)), v4 = rbindlist(mclapply(tis, shares, tr = "v4", mc.cores = 10)))
excl_of <- function(v) { s <- copy(SH[[v$tree]]); s <- merge(s, info, by = "promoter", all.x = TRUE); s <- merge(s, GD, by = c("tissue", "gene"), all.x = TRUE)
  s[, e := pmax(share0, share1) >= ECUT][, g := v$gene_cut > 0 & (is.na(med_D) | med_D < v$gene_cut)]
  s[e | g, .(tissue, promoter, gene, reason = fifelse(e & g, "endpoint+gene_reads", fifelse(e, "endpoint", "gene_reads")),
             share0 = round(share0, 4), share1 = round(share1, 4), med_D)] }

## Fig. 2E on the kept promoters (logic of fig2e_deliver_freeze_v3m_pairs.R)
fig2e_one <- function(tr, t, arm, ex, OUT) {
  U <- as.data.table(readRDS(file.path(TREE[[tr]], "results/final_nogate_v3m/freeze", paste(t, arm, "units.rds", sep = "."))))[dedup_drop == FALSE & estimable == TRUE & is.finite(p)]
  U <- U[!promoter %in% ex]
  U[, unit := fifelse(boundary_exact %in% TRUE & interior_exact %in% TRUE, mirror_pair_id, promoter)]
  U[, promoters_in_unit := paste(sort(unique(promoter)), collapse = ";"), by = unit]
  U[, promoter_real := promoter][, promoter := unit]
  B <- U[, .(gene = gene[1], promoters_in_unit = promoters_in_unit[1], p_simes = simes(p),
             types = paste(sort(paste0(type, if (uniqueN(promoter_real) > 1) paste0("@", promoter_real) else "")), collapse = ","),
             p_mu = if (any(type == "mu")) min(p[type == "mu"]) else NA_real_, p_p0 = if (any(type == "p0")) min(p[type == "p0"]) else NA_real_,
             p_p1 = if (any(type == "p1")) min(p[type == "p1"]) else NA_real_), by = promoter]
  B[, q := p.adjust(p_simes, "BH")][, significant := q <= 0.05]
  m <- nrow(B); R <- B[significant == TRUE, .N]; lvl <- 0.05 * R / m
  L <- U[promoter %in% B[significant == TRUE]$promoter][, rej := p.adjust(p, "BH") <= lvl, by = promoter][,
         .(layer_class = { s <- type[rej]; if (all(s == "mu")) "mu only" else if (!any(s == "mu")) "p0/p1 only" else "both" }), by = promoter]
  B <- merge(B, L, by = "promoter", all.x = TRUE); res <- list()
  for (mode in MODES) {
    d <- fr(file.path(TREE[[tr]], "results/fig2e_freeze_v3m/qtltools", mode, paste(t, arm, "tsv.gz", sep = ".")), select = c("promoter", "p", "promoters_in_unit"))
    d <- d[promoter %in% B$promoter][, q := p.adjust(p, "BH")][, significant := q <= 0.05]
    B[, (paste0("qtltools_", mode, "_significant")) := promoter %in% d[significant == TRUE]$promoter]
    res[[mode]] <- data.table(tissue = t, arm = arm, qtltools_mode = mode, promoters = m, bohr = R,
      bohr_mu_only = L[layer_class == "mu only", .N], bohr_both = L[layer_class == "both", .N], bohr_p0p1_only = L[layer_class == "p0/p1 only", .N],
      qtltools = d[significant == TRUE, .N], qtltools_tested = nrow(d))
    fwrite(d, file.path(OUT, "qtltools", mode, paste(t, arm, "tsv.gz", sep = ".")), sep = "\t") }
  fwrite(B[order(p_simes)], file.path(OUT, "bohr", paste(t, arm, "tsv.gz", sep = ".")), sep = "\t")
  rbindlist(res) }

SUM <- list()
for (vn in names(VERS)) {
  v <- VERS[[vn]]; tr <- v$tree; OD <- file.path(OUTR, vn); dir.create(OD, FALSE, TRUE)
  EX <- excl_of(v); fwrite(EX, file.path(OD, "promoters_excluded.tsv.gz"), sep = "\t")
  ## units v3c (the final units of each tree) and re-BH
  U <- rbindlist(mclapply(tis, function(t) { x <- as.data.table(readRDS(file.path(TREE[[tr]], "results/final_nogate_v3c/freeze", paste0(t, ".real.units.rds"))))
    x[, tissue := t]; x[, .(tissue, promoter, type, chr = sub("_.*$", "", leader), p, q, estimable, dedup_drop)] }, mc.cores = 10))
  U[, tested := estimable %in% TRUE & dedup_drop %in% FALSE & !is.na(p)][, sig_old := tested & !is.na(q) & q <= 0.05]
  U[, kept := !paste(tissue, promoter) %in% paste(EX$tissue, EX$promoter)]
  U[tested & kept, qn := p.adjust(p, "BH"), by = tissue][, sig_new := tested & kept & !is.na(qn) & qn <= 0.05]
  fwrite(U[sig_old & sig_new, .(tissue, promoter, type)], file.path(OD, "units_keep.tsv.gz"), sep = "\t")
  fwrite(U[!sig_old & sig_new, .(tissue, promoter, type, chr, p, q_before = q, q_after = qn)], file.path(OD, "units_new_topup.tsv.gz"), sep = "\t")
  fwrite(U[sig_old & kept & !sig_new, .(tissue, promoter, type, q_before = q, q_after = qn)], file.path(OD, "units_lost.tsv.gz"), sep = "\t")
  S <- fr(file.path(TREE[[tr]], "results/fig5_delivery_v7/signals.tsv.gz"), select = c("tissue", "promoter", "type", "k"))
  K <- U[sig_old & sig_new, .(tissue, promoter, type)]
  SK <- merge(S, K, by = c("tissue", "promoter", "type")); fwrite(SK, file.path(OD, "signals_keep.tsv.gz"), sep = "\t")
  SUM[[vn]] <- U[, .(version = vn, sig_before = sum(sig_old), removed_by_filter = sum(sig_old & !kept), kept_and_sig = sum(sig_old & sig_new),
                     lost_after_rebh = sum(sig_old & kept & !sig_new), new_after_rebh = sum(!sig_old & sig_new), sig_after = sum(sig_new)), by = type]
  SUM[[paste0(vn, "_signals")]] <- data.table(version = vn, type = "signals", sig_before = nrow(S), removed_by_filter = nrow(S) - nrow(SK), kept_and_sig = nrow(SK),
                                              lost_after_rebh = NA, new_after_rebh = NA, sig_after = NA)
  ## Fig. 2E
  F2 <- file.path(OD, "fig2e"); for (d in c("bohr", paste0("qtltools/", MODES))) dir.create(file.path(F2, d), FALSE, TRUE)
  exl <- split(EX$promoter, EX$tissue)
  Rt <- rbindlist(mclapply(tis, function(t) rbindlist(lapply(ARMS, function(a) fig2e_one(tr, t, a, exl[[t]], F2))), mc.cores = 10))
  fwrite(Rt, file.path(F2, "counts.tsv"), sep = "\t")
  L2 <- rbind(unique(Rt[qtltools_mode == "rintfirst", .(tissue, arm, method = "BOHR", n = bohr, promoters)]),
              Rt[, .(tissue, arm, method = paste0("QTLtools ", qtltools_mode), n = qtltools, promoters)])
  D <- L2[, .(promoters = promoters[arm == "real"], real = n[arm == "real"], n_perm = sum(arm != "real"), perm_mean = round(mean(n[arm != "real"]), 2),
              perm_min = min(n[arm != "real"]), perm_max = max(n[arm != "real"])), by = .(tissue, method)]
  D <- merge(D, unique(Rt[arm == "real" & qtltools_mode == "rintfirst", .(tissue, bohr_mu_only, bohr_both, bohr_p0p1_only)]), by = "tissue")
  D[method != "BOHR", c("bohr_mu_only", "bohr_both", "bohr_p0p1_only") := NA_integer_]
  fwrite(D, file.path(F2, "counts_display.tsv"), sep = "\t")
  cat(vn, ": Fig. 2E real totals\n"); print(D[, .(real = sum(real), perm_mean = sum(perm_mean)), by = method]) }
SM <- rbindlist(SUM, fill = TRUE); fwrite(SM, file.path(OUTR, "summary.tsv"), sep = "\t"); options(width = 200); print(SM)
## check: Fig. 2E totals equal the filter-only computation (fig2e_variants.R: E80 paper / v4, G5E80 v4)
FV <- fread(file.path(C, "out/fig2e_variants/fig2e_variants.tsv"))
for (x in list(c("paper_E80", "E80", "paper"), c("nofloor_E80", "E80", "v4"), c("nofloor_reads5_E80", "G5E80", "v4"))) {
  D <- fread(file.path(OUTR, x[1], "fig2e/counts_display.tsv"))
  a <- D[, .(real = sum(real)), by = method][order(method)]; b <- FV[variant == x[2] & version == x[3], .(method = c(BOHR = "BOHR", raw = "QTLtools raw", normal = "QTLtools normal", rint = "QTLtools rintfirst")[method], real)][order(method)]
  cat(x[1], "matches fig2e_variants:", isTRUE(all.equal(a$real, b$real)), "\n"); print(merge(a, b, by = "method", suffixes = c("_rebuilt", "_filter"))) }
