## v6c 2026-11-07: >= 3 promoters use the pooled major (see the v6c block below).
## 101_K_overlap_multi_location.R -- shared leads, genes with >= 3 promoters and lead location. Frame in plain terms: in the
## analysed groups the switching promoter is always the MINOR promoter and the other one is the MAJOR promoter
## (p0 on the minor, p1 on the major = the minor switching, mu: minor by convention); the signals where the major switches
## (p0 on the major, p1 on the minor) stay out, as before.
##  K1 overlap: drop signals whose lead is shared with a signal of the other kind (p0/p1 vs mu) in the same gene x tissue
##     (identical lead or lead inside the other signal's r2 >= 0.8 proxy set); also keep only gene x tissue with one kind only.
##  K2 genes with >= 3 promoters: minor side vs the rest of the gene (p0 on a minor promoter: minor = own, rest = all others;
##     p1 on the major: minor side = all other promoters together, rest = the major; mu: own if own is minor, else all others).
##  K3 location of the lead: distance to the minor and to the major TSS (partner TSS, not the nearest promoter), by class and
##     inside the same promoter x tissue pairs (F2).
## Read-outs per set: minor only, both together (same direction), opposite (minor up, major down), major only (minor unchanged),
## gene output moves with the minor, mean gene-output effect. Gene bootstrap B = 2000.
source(file.path(Sys.getenv("F5S", "/home/e1327989/f5_stage_20261107_v6c"), "code/00_helpers.R"))
OUT <- file.path(OUTROOT, "K_overlap_multi_location"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE); B <- 2000L; ZC <- 1.96
X <- fr(file.path(OUTROOT, "B1v2_multipromoter/B1v2_signal_patterns.tsv.gz"))
a1 <- fr(file.path(OUTROOT, "A1_location/A1_signal_location.tsv"), select = c("sid", "own_minor"))
X <- merge(X, a1, by = "sid", all.x = TRUE)
sig <- fr(file.path(IN, "delivery/signals.tsv.gz"), select = c("tissue", "promoter", "type", "k", "snp"))
sig[, sid := paste(tissue, promoter, type, k, sep = "|")]
X <- merge(X, sig[, .(sid, lead = snp)], by = "sid")
## v5 (Methods §9): code by |z| regardless of direction (z oriented to the allele raising the minor promoter's usage)
patS <- function(s, n) fifelse(abs(s) >= ZC & abs(n) < ZC, "minor_only", fifelse(abs(s) < ZC & abs(n) >= ZC, "major_only",
  fifelse(abs(s) >= ZC & abs(n) >= ZC, fifelse(sign(s) == sign(n), "together", "opposite"), "neither")))
## two-promoter frame (as 74b): keep minor-switching groups
A2 <- X[two == TRUE & route %in% c("onoff_S_minor", "mu"), .(sid, gene, tissue, type, lead, minor_is_own = S_is_own, partner, promoter,
        zM = zS, zJ = zN, zT = zTS, dT, kind = fifelse(type == "mu", "mu", "p0p1"), frame = "two")]
## >= 3 promoters: own_mean (lead genotype groups) is kept only to compare the old rule (own_mean < 0.5) with the major promoter set below
pxl <- fr(file.path(IN, "delivery/genotype_pheno_proxies.tsv.gz"), select = c("tissue", "promoter", "type", "k", "signal_snp", "proxy_snp", "n0", "n1", "n2", "m0", "m1", "m2"))[signal_snp == proxy_snp]
pxl[, sid := paste(tissue, promoter, type, k, sep = "|")]; pxl <- unique(pxl[, .(own_mean = (n0 * m0 + n1 * m1 + n2 * m2) / (n0 + n1 + n2)), by = sid], by = "sid")
X <- merge(X, pxl, by = "sid", all.x = TRUE)
## v6c (2026-11-07; the largest promoter is the major one): >= 3 tested promoters use the same pooled definition as the
## two-promoter genes (Figure 4's mirror canonical). Pooled usage of a promoter = its mean usage (freeze td$act) over all freeze
## samples of the tissues where it is tested; this reproduces mirror_pairs pooled_mean exactly (checked below,
## stopifnot). Major = the gene's tested promoter with the highest pooled usage; every other tested promoter is minor.
## Replaces "own mean usage in that tissue < 0.5" (v5/v6), which called every promoter minor when none had usage >= 0.5.
TTf <- fr(file.path(IN, "freeze/tissues.tsv"))[def == "freeze"]; stopifnot(nrow(TTf) == 50)
PS <- rbindlist(lapply(TTf$tissue, function(t) { a <- readRDS(file.path(IN, "freeze", paste0(t, ".rds")))$act
  data.table(pid = colnames(a), s = colSums(a, na.rm = TRUE), n = colSums(!is.na(a))) }))[, .(pooled_usage = sum(s) / sum(n)), by = pid]
mp0 <- fr(file.path(IN, "mirror/mirror_pairs.tsv.gz"), select = c("promoter_1", "pooled_mean_1"))
ck <- merge(unique(mp0, by = "promoter_1"), PS, by.x = "promoter_1", by.y = "pid"); stopifnot(max(abs(ck$pooled_mean_1 - ck$pooled_usage)) < 1e-9)
f5g <- fread(file.path(IN, "v11/Promoter_info_filter_5.txt")); setnames(f5g, 1, "pid")
g3 <- unique(X[two == FALSE & n_prom_gene >= 3, gene])
PU <- merge(f5g[geneId %in% g3, .(gene = geneId, pid)], PS, by = "pid")
MJ <- PU[order(gene, -pooled_usage)][, .(major_p = pid[1], major_pooled = pooled_usage[1], second_pooled = pooled_usage[2], n_tested = .N), by = gene]
cat("genes with >= 3 tested promoters:", nrow(MJ), "; pooled usage of the major: median", round(median(MJ$major_pooled), 3), "; major < 0.5 in", MJ[major_pooled < 0.5, .N], "genes\n")
fwrite(merge(PU, MJ[, .(gene, major_p)], by = "gene")[, is_major := pid == major_p][], file.path(OUT, "K2_pooled_major_three_plus.tsv.gz"), sep = "\t")
X <- merge(X, MJ[, .(gene, major_p)], by = "gene", all.x = TRUE)
X[two == FALSE, own_minor_v6 := own_mean < 0.5]
X[two == FALSE, own_minor := fifelse(is.na(major_p), NA, promoter != major_p)]
cat("own_minor, old (tissue mean < 0.5) vs new (pooled major), >= 3 tested promoters:\n"); print(X[two == FALSE & n_prom_gene >= 3, .N, by = .(type, old_minor = own_minor_v6, new_minor = own_minor)][order(type, old_minor, new_minor)])
T3 <- X[two == FALSE & !is.na(own_minor) & n_prom_gene >= 3]   # v5: genes with >= 3 TESTED promoters only (signals in one-promoter genes dropped)
T3[, minor_is_own := fifelse(type == "p0", TRUE, fifelse(type == "p1", FALSE, own_minor))]
T3 <- T3[(type == "p0" & own_minor == TRUE) | (type == "p1" & own_minor == FALSE) | type == "mu"]
T3[, sM := fifelse(minor_is_own, dir_c, -dir_c)]
T3[, `:=`(zM = fifelse(minor_is_own, z_t, z_O) * sM, zJ = fifelse(minor_is_own, z_O, z_t) * sM, zT = z_T * sM, dT = b_T / sd_T * sM)]
print(T3[, .N, by = .(type, own_minor)])
A3 <- T3[, .(sid, gene, tissue, type, lead, minor_is_own, partner, promoter, zM, zJ, zT, dT, kind = fifelse(type == "mu", "mu", "p0p1"), frame = "three_plus")]
A <- rbind(A2, A3)[is.finite(zM) & is.finite(zJ)]
A[, `:=`(pat = patS(zM, zJ), out_with = as.numeric(zT >= ZC))]
cat("signals: two", A[frame == "two", .N], "; three_plus", A[frame == "three_plus", .N], "\n"); print(A[, .N, by = .(frame, type)])
## K1 shared leads (identical or r2 >= 0.8), same gene x tissue, between p0/p1 and mu
px <- fr(file.path(IN, "delivery/genotype_pheno_proxies.tsv.gz"), select = c("tissue", "promoter", "type", "k", "proxy_snp", "r2"))
px[, sid := paste(tissue, promoter, type, k, sep = "|")]
SP <- unique(rbind(merge(px[r2 >= 0.8, .(sid, snp = proxy_snp)], A[, .(sid, gene, tissue, kind)], by = "sid"),
                   A[, .(sid, snp = lead, gene, tissue, kind)]))
L <- A[, .(sid, gene, tissue, kind, snp = lead)]
hit <- merge(L, SP, by = c("gene", "tissue", "snp"), suffixes = c("", ".o"), allow.cartesian = TRUE)[kind != kind.o]
sh_ids <- unique(c(hit$sid, hit$sid.o))
A[, shared := sid %in% sh_ids]
A[, gt := paste(gene, tissue)]; A[, mixed_unit := uniqueN(kind) > 1, by = .(gt, frame)]
cat("shared leads:", A[, sum(shared)], "of", nrow(A), "\n"); print(A[, .(n = .N, shared = mean(shared), mixed_unit = mean(mixed_unit)), by = .(frame, type)])
fwrite(A, file.path(OUT, "K_signals.tsv.gz"), sep = "\t")
tabs <- list()
add <- function(d, groups, lab, seed) { groups <- intersect(groups, unique(d$group)); if (length(groups) < 1 || nrow(d) < 5) return(invisible())
  r <- ratio_table(d, groups, B, seed, lab)
  tabs[[length(tabs) + 1]] <<- rbind(r$groups[, .(label, group, n_signals, estimate, lo, hi, contrast = NA_character_, p_boot = NA_real_)],
    if (!is.null(r$contrasts)) r$contrasts[, .(label, group = NA_character_, n_signals = NA_real_, estimate = difference, lo, hi, contrast, p_boot)]) }
sets <- list("two-promoter genes, all" = A[frame == "two"],
             "two-promoter genes, shared leads removed" = A[frame == "two" & shared == FALSE],
             "two-promoter genes, gene x tissue with one kind only" = A[frame == "two" & mixed_unit == FALSE],
             ">= 3 promoters, all" = A[frame == "three_plus"],
             ">= 3 promoters, shared leads removed" = A[frame == "three_plus" & shared == FALSE],
             "all genes, shared leads removed" = A[shared == FALSE])
mets <- c(together = "both together (same direction)", major_only = "major only (minor unchanged)", opposite = "opposite (minor up, major down)", minor_only = "minor only")
for (sn in names(sets)) { d <- sets[[sn]]; s0 <- match(sn, names(sets)) * 100L
  for (m in names(mets)) { add(d[, .(gene, group = kind, num = as.numeric(pat == m), den = 1)], c("p0p1", "mu"), paste(sn, "|", mets[[m]]), s0 + 1L)
                           add(d[, .(gene, group = type, num = as.numeric(pat == m), den = 1)], c("p0", "p1", "mu"), paste(sn, "|", mets[[m]], "| by type"), s0 + 2L) }
  add(d[, .(gene, group = kind, num = out_with, den = 1)], c("p0p1", "mu"), paste(sn, "| gene output moves with the minor"), s0 + 3L)
  add(d[, .(gene, group = type, num = out_with, den = 1)], c("p0", "p1", "mu"), paste(sn, "| gene output moves with the minor | by type"), s0 + 4L)
  add(d[is.finite(dT), .(gene, group = kind, num = dT, den = 1)], c("p0p1", "mu"), paste(sn, "| mean gene-output effect (SD)"), s0 + 5L) }
## K3 location (two-promoter frame): distance of the lead to the minor and the major TSS
info <- fread(file.path(IN, "v11/Promoter_info_filter_5.txt")); setnames(info, 1, "pid")
D <- A[frame == "two"]; D[, pos := as.integer(sub("^chr[^_]+_([0-9]+)_.*$", "\\1", lead))]
D[, `:=`(minor_p = fifelse(minor_is_own, promoter, partner), major_p = fifelse(minor_is_own, partner, promoter))]
D <- merge(D, info[, .(minor_p = pid, tss_minor = start)], by = "minor_p"); D <- merge(D, info[, .(major_p = pid, tss_major = start)], by = "major_p")
D[, `:=`(d_minor = abs(pos - tss_minor), d_major = abs(pos - tss_major))]
D[, `:=`(near_minor1k = as.numeric(d_minor <= 1000), near_major1k = as.numeric(d_major <= 1000), nearer_major = as.numeric(d_major < d_minor),
         at_minor_core = as.numeric(d_minor <= 300), at_major_core = as.numeric(d_major <= 300))]
fwrite(D[, .(sid, type, kind, pat, d_minor, d_major, dT)], file.path(OUT, "K3_lead_location.tsv.gz"), sep = "\t")
locm <- c(at_minor_core = "lead within 300 bp of the minor TSS", near_minor1k = "lead within 1 kb of the minor TSS",
          at_major_core = "lead within 300 bp of the major TSS", near_major1k = "lead within 1 kb of the major TSS", nearer_major = "lead nearer the major than the minor TSS")
for (m in names(locm)) { add(D[, .(gene, group = type, num = get(m), den = 1)], c("p0", "p1", "mu"), paste("K3 location |", locm[[m]], "| by type"), 700L)
                         add(D[, .(gene, group = kind, num = get(m), den = 1)], c("p0p1", "mu"), paste("K3 location |", locm[[m]]), 701L)
                         add(D[pat == "major_only" | kind == "p0p1", .(gene, group = fifelse(kind == "mu", "mu_major_only", "p0p1"), num = get(m), den = 1)], c("p0p1", "mu_major_only"),
                             paste("K3 location |", locm[[m]], "| mu signals with major only"), 702L) }
TAB <- rbindlist(tabs); writeT(TAB, file.path(OUT, "K_summary.tsv"))
## K3b F2 pairs (same promoter x tissue): paired location
Q <- fr(file.path(OUTROOT, "I_continuous_effects/I1_F2_pairs.tsv.gz"), select = c("sid_b", "sid_mu", "gene"))
Q <- merge(Q, D[, .(sid_b = sid, b_minor = d_minor, b_major = d_major)], by = "sid_b"); Q <- merge(Q, D[, .(sid_mu = sid, m_minor = d_minor, m_major = d_major)], by = "sid_mu")
pf <- function(v_b, v_m, lab, seed) { ug <- unique(Q$gene); f <- function(d) c(mean(v_b(d)), mean(v_m(d))); set.seed(seed); est <- f(Q)
  bt <- t(replicate(B, { g <- sample(ug, length(ug), TRUE); f(Q[data.table(gene = g), on = "gene", allow.cartesian = TRUE]) }))
  data.table(label = lab, p0p1 = est[1], mu = est[2], diff = est[1] - est[2], lo = quantile(bt[, 1] - bt[, 2], 0.025), hi = quantile(bt[, 1] - bt[, 2], 0.975),
             p_boot = min(1, max(2 * min(mean(bt[, 1] - bt[, 2] <= 0), mean(bt[, 1] - bt[, 2] >= 0)), 1 / B)), n_pairs = nrow(Q)) }
F2L <- rbind(pf(function(d) d$b_minor <= 300, function(d) d$m_minor <= 300, "F2 pairs: lead within 300 bp of the minor TSS", 801L),
             pf(function(d) d$b_minor <= 1000, function(d) d$m_minor <= 1000, "F2 pairs: lead within 1 kb of the minor TSS", 802L),
             pf(function(d) d$b_major <= 1000, function(d) d$m_major <= 1000, "F2 pairs: lead within 1 kb of the major TSS", 803L),
             pf(function(d) d$b_major < d$b_minor, function(d) d$m_major < d$m_minor, "F2 pairs: lead nearer the major TSS", 804L),
             pf(function(d) log10(d$b_minor + 10), function(d) log10(d$m_minor + 10), "F2 pairs: mean log10 distance to the minor TSS", 805L))
writeT(F2L, file.path(OUT, "K3_F2_location.tsv")); print(F2L)
fmtrow <- function(x) if (is.na(x$contrast)) sprintf("  - %s: %s [%s, %s] (n = %s)", x$group, if (grepl("mean gene-output", x$label)) sprintf("%.3f", x$estimate) else pct(x$estimate),
  if (grepl("mean gene-output", x$label)) sprintf("%.3f", x$lo) else pct(x$lo), if (grepl("mean gene-output", x$label)) sprintf("%.3f", x$hi) else pct(x$hi), format(round(x$n_signals, 1))) else
  sprintf("  - %s: %+.1f pp [%+.1f, %+.1f], P = %s", x$contrast, 100 * x$estimate, 100 * x$lo, 100 * x$hi, formatC(x$p_boot, format = "g", digits = 2))
Lr <- c("# K -- overlap removal, genes with >= 3 promoters, lead location (minor / major frame)", "",
        sprintf("Signals: two-promoter %d (shared leads %d); >= 3 promoters %d (shared %d).", A[frame == "two", .N], A[frame == "two", sum(shared)], A[frame == "three_plus", .N], A[frame == "three_plus", sum(shared)]), "")
for (lab in unique(TAB$label)) { Lr <- c(Lr, paste0("- ", lab)); x <- TAB[label == lab]; for (i in seq_len(nrow(x))) Lr <- c(Lr, fmtrow(x[i])) }
Lr <- c(Lr, "", "F2 pairs (same promoter x tissue):", F2L[, sprintf("- %s: p0/p1 %.3f vs mu %.3f; diff %+.3f [%+.3f, %+.3f], P = %s (n = %d pairs)", label, p0p1, mu, diff, lo, hi, formatC(p_boot, format = "g", digits = 2), n_pairs)])
writeLines(Lr, file.path(OUT, "K_report.md")); cat(Lr, sep = "\n")
