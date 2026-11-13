## 71_A1_location.R -- SNP-level regulatory routes, A1: where the candidate variants of each signal sit relative to
## ALL promoters of the gene (not only the promoter the signal was called on).
## Signals: all signals of the staged delivery. Candidate variants: lead + r2 >= 0.8 proxies (PLINK standard),
## each signal's variants weighted 1/n so every signal counts once.
## Other promoters: the gene's promoters in the same tissue's freeze activity matrix (= the promoters the mirror code
## and n_prom_gene use). Major/minor: pooled canonical/partner for two-promoter genes (mirror code); tissue mean usage
## rank for >= 3-promoter genes.
## Categories at threshold T (exclusive): own = within T of the signal's TSS and not closer to another promoter;
## other = within T of another promoter of the gene and closer to it; gene body = inside the GENCODE v26 gene span;
## outside = the rest. T = 1 kb (main), 300 bp and 5 kb (sensitivity).
## Two-promoter genes, aligned frame: in a two-promoter gene "minor = 0" is the same event as "major = 1", so an on/off
## (p0/p1) signal is described by the promoter that switches (S): S = own for p0, S = the other promoter for p1.
## mu signals are described by minor / major.
## Outputs: out/A1_location/*.tsv and A1_report.md
source(file.path(Sys.getenv("F5S", "/home/e1327989/f5_stage_20261106_v6"), "code/00_helpers.R"))
OUT <- file.path(OUTROOT, "A1_location"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setDTthreads(2)
B <- 2000L
cat("start", format(Sys.time()), "\n")

sig <- fr(file.path(IN, "delivery/signals_partner.tsv.gz"))
s2  <- fr(file.path(IN, "delivery/signals.tsv.gz"), select = c("tissue", "promoter", "type", "k", "logp_at_selection", "maf", "backward_keep"))
sig <- merge(sig, s2, by = c("tissue", "promoter", "type", "k"), all.x = TRUE)
stopifnot(nrow(sig) == 14823)   # delivery v6 (conditional refit 2026-11-06); 14,825 in v4, 14,808 in v3
fl <- fr(file.path(IN, "mirror/promoter_flags.tsv.gz"))[, .(tissue, promoter, gene, n_prom_gene, mirror_role)]
info <- fread(file.path(IN, "v11/Promoter_info_filter_5.txt")); setnames(info, 1, "pid")
info <- info[, .(promoter = pid, chr = seqnames, tss = start, pstrand = strand)]
TT <- fread(file.path(IN, "freeze/tissues.tsv"))[def == "freeze"]
use <- rbindlist(lapply(TT$tissue, function(t) { td <- readRDS(file.path(IN, "freeze", paste0(t, ".rds")))
  data.table(tissue = t, promoter = colnames(td$act), mean_use = colMeans(td$act, na.rm = TRUE)) }))
P <- merge(fl, info, by = "promoter")
P <- merge(P, use, by = c("tissue", "promoter"), all.x = TRUE)
P[, rank_in_gene := frank(-mean_use, ties.method = "first"), by = .(tissue, gene)]
cat("tissue x promoter rows", nrow(P), "\n")

sig <- merge(sig, P[, .(tissue, promoter, gene, n_prom_gene, own_rank = rank_in_gene, own_use = mean_use)], by = c("tissue", "promoter"), all.x = TRUE)
stopifnot(!anyNA(sig$gene))
sig[, gene_base := sub("\\..*$", "", gene)]
sig[, sid := paste(tissue, promoter, type, k, sep = "|")]

## gene spans (GENCODE v26; version-stripped match)
gtf <- fread(cmd = sprintf("zcat %s | awk -F'\\t' '$3==\"gene\"' | cut -f1,4,5,9", file.path(IN, "gtf/gencode.v26.annotation.gtf.gz")),
             header = FALSE, sep = "\t", col.names = c("gchr", "gstart", "gend", "attr"))
gtf[, gene_base := sub('^gene_id "([^".]+).*$', "\\1", attr)]; gtf <- gtf[!duplicated(gene_base), .(gene_base, gstart, gend)]
sig <- merge(sig, gtf, by = "gene_base", all.x = TRUE)
cat("signals with a GENCODE v26 gene span:", sum(!is.na(sig$gstart)), "of", nrow(sig), "\n")

## candidate variants
gp <- fr(file.path(IN, "delivery/genotype_pheno_proxies.tsv.gz"), select = c("tissue", "promoter", "type", "k", "signal_snp", "proxy_snp", "r2"))[r2 >= 0.8]
gp[, pos := as.integer(sub("^chr[^_]+_([0-9]+)_.*$", "\\1", proxy_snp))]
gp[, sid := paste(tissue, promoter, type, k, sep = "|")]
stopifnot(all(sig$sid %in% gp$sid))
gp[, w := 1 / .N, by = sid]
V <- merge(gp[, .(sid, proxy_snp, pos, w, is_lead = proxy_snp == signal_snp)],
           sig[, .(sid, tissue, gene, promoter, tss_own = tss, gstart, gend)], by = "sid")
V[, vid := .I]
OP <- P[, .(tissue, gene, oprom = promoter, otss = tss, orank = rank_in_gene)]
X <- merge(V[, .(vid, tissue, gene, promoter, pos)], OP, by = c("tissue", "gene"), allow.cartesian = TRUE)[oprom != promoter]
X[, d := abs(pos - otss)]
X <- X[X[, .I[which.min(d)], by = vid]$V1]
V <- merge(V, X[, .(vid, d_oth = d, near_oprom = oprom, near_orank = orank)], by = "vid", all.x = TRUE)
V[, d_own := abs(pos - tss_own)]
V[, in_gene := !is.na(gstart) & pos >= gstart & pos <= gend]
catT <- function(T) fifelse(V$d_own <= T & (is.na(V$d_oth) | V$d_own <= V$d_oth), "own",
                    fifelse(!is.na(V$d_oth) & V$d_oth <= T & V$d_oth < V$d_own, "other", fifelse(V$in_gene, "gene_body", "outside")))
for (T in c(300L, 1000L, 5000L)) V[, (paste0("cat_", T)) := catT(T)]
cat("variant rows", nrow(V), "\n")
fwrite(V[, .(sid, proxy_snp, pos, w, is_lead, d_own, d_oth, near_oprom, near_orank, in_gene, cat_300, cat_1000, cat_5000)], file.path(OUT, "A1_variant_location.tsv.gz"), sep = "\t")

## per-signal weighted fractions
agg <- V[, .(n_var = .N,
  f_own = sum(w * (cat_1000 == "own")), f_other = sum(w * (cat_1000 == "other")),
  f_gene = sum(w * (cat_1000 == "gene_body")), f_out = sum(w * (cat_1000 == "outside")),
  f_own300 = sum(w * (cat_300 == "own")), f_other300 = sum(w * (cat_300 == "other")),
  f_own5k = sum(w * (cat_5000 == "own")), f_other5k = sum(w * (cat_5000 == "other")),
  lead_cat = cat_1000[is_lead][1], lead_d_own = d_own[is_lead][1], lead_d_oth = d_oth[is_lead][1],
  f_near_major = sum(w * (cat_1000 == "other" & near_orank == 1))), by = sid]
S1 <- merge(sig, agg, by = "sid")
S1[, spacing := abs(tss - partner_tss)]
## author 2026-11-04: major / minor = higher / lower mean usage POOLED over all tissues (= the mirror table's canonical / partner, as Figure 4).
S1[, own_minor := fifelse(n_prom_gene == 2, mirror_role == "partner", NA)]
S1[, nlp := pmin(nlog10p(logp_at_selection), 50)]
S1 <- merge(S1, TT[, .(tissue, n_donors = n)], by = "tissue")
## aligned frame (two-promoter genes)
S1[n_prom_gene == 2 & type %in% c("p0", "p1"), `:=`(
  route = fifelse((type == "p0") == own_minor, "onoff_S_minor", "onoff_S_major"),
  f_S = fifelse(type == "p0", f_own, f_other), f_N = fifelse(type == "p0", f_other, f_own),
  f_S300 = fifelse(type == "p0", f_own300, f_other300), f_N300 = fifelse(type == "p0", f_other300, f_own300),
  f_S5k = fifelse(type == "p0", f_own5k, f_other5k), f_N5k = fifelse(type == "p0", f_other5k, f_own5k))]
S1[n_prom_gene == 2 & type == "mu", `:=`(route = "mu",
  f_S = fifelse(own_minor, f_own, f_other), f_N = fifelse(own_minor, f_other, f_own),
  f_S300 = fifelse(own_minor, f_own300, f_other300), f_N300 = fifelse(own_minor, f_other300, f_own300),
  f_S5k = fifelse(own_minor, f_own5k, f_other5k), f_N5k = fifelse(own_minor, f_other5k, f_own5k))]
## for mu, "S" = the minor promoter and "N" = the major promoter
writeT(S1[, .(sid, tissue, promoter, type, k, gene, n_prom_gene, mirror_role, exact_mirror, own_minor, route, spacing, n_var,
             f_own, f_other, f_gene, f_out, f_own300, f_other300, f_own5k, f_other5k, f_S, f_N, f_S300, f_N300, f_S5k, f_N5k,
             f_near_major, lead_cat, lead_d_own, lead_d_oth, maf, nlp, n_donors)], file.path(OUT, "A1_signal_location.tsv"))

## ---- 1. class level, all signals
D <- rbindlist(lapply(c("f_own", "f_other", "f_gene", "f_out"), function(v) S1[, .(gene, group = type, num = get(v), den = 1, var = v)]))
res1 <- rbindlist(lapply(split(D, D$var), function(d) { r <- ratio_table(d, c("mu", "p0", "p1"), B, 11L, d$var[1]); r$groups }))
con1 <- rbindlist(lapply(split(D, D$var), function(d) { r <- ratio_table(d, c("mu", "p0", "p1"), B, 11L, d$var[1]); r$contrasts }))
writeT(res1, file.path(OUT, "A1_class_fractions.tsv")); writeT(con1, file.path(OUT, "A1_class_contrasts.tsv"))
## by class x n_prom_gene (2 vs >= 3)
S1[, gene_kind := fifelse(n_prom_gene == 2, "two_promoter", "three_plus")]
res1b <- rbindlist(lapply(c("two_promoter", "three_plus"), function(gk) rbindlist(lapply(c("f_own", "f_other"), function(v) {
  r <- ratio_table(S1[gene_kind == gk, .(gene, group = type, num = get(v), den = 1)], c("mu", "p0", "p1"), B, 12L, paste(gk, v))
  rbind(r$groups[, .(label, group, n_signals, estimate, lo, hi, contrast = NA_character_, p_boot = NA_real_)],
        r$contrasts[, .(label, group = NA_character_, n_signals = NA_real_, estimate = difference, lo, hi, contrast, p_boot)]) }))))
writeT(res1b, file.path(OUT, "A1_class_by_genekind.tsv"))

## ---- 2. aligned frame, two-promoter genes
A <- S1[!is.na(route)]
alg <- rbindlist(lapply(list(c("f_S", "f_N", "1kb"), c("f_S300", "f_N300", "300bp"), c("f_S5k", "f_N5k", "5kb")), function(v) {
  dS <- A[, .(gene, group = route, num = get(v[1]), den = 1)]; dN <- A[, .(gene, group = route, num = get(v[2]), den = 1)]
  dsh <- A[get(v[1]) + get(v[2]) > 0, .(gene, group = route, num = get(v[2]), den = get(v[1]) + get(v[2]))]
  g <- c("onoff_S_minor", "onoff_S_major", "mu")
  r1 <- ratio_table(dS, g, B, 21L, paste0("at S (switching; mu: minor) ", v[3]))
  r2 <- ratio_table(dN, g, B, 21L, paste0("at N (non-switching; mu: major) ", v[3]))
  r3 <- ratio_table(dsh, g, B, 21L, paste0("share at N among promoter-proximal ", v[3]))
  rbind(rbindlist(lapply(list(r1, r2, r3), function(r) r$groups[, .(label, group, n_signals, estimate, lo, hi, contrast = NA_character_, p_boot = NA_real_)])),
        rbindlist(lapply(list(r1, r2, r3), function(r) r$contrasts[, .(label, group = NA_character_, n_signals = NA_real_, estimate = difference, lo, hi, contrast, p_boot)])))
}))
writeT(alg, file.path(OUT, "A1_aligned_two_promoter.tsv"))
## exact-mirror vs not (two-promoter genes)
alg_ex <- rbindlist(lapply(c(TRUE, FALSE), function(ex) {
  dsh <- A[exact_mirror == ex & f_S + f_N > 0, .(gene, group = route, num = f_N, den = f_S + f_N)]
  r <- ratio_table(dsh, c("onoff_S_minor", "mu"), B, 22L, paste("share at N, exact_mirror =", ex))
  rbind(r$groups[, .(label, group, n_signals, estimate, lo, hi, contrast = NA_character_, p_boot = NA_real_)],
        r$contrasts[, .(label, group = NA_character_, n_signals = NA_real_, estimate = difference, lo, hi, contrast, p_boot)]) }))
writeT(alg_ex, file.path(OUT, "A1_aligned_by_exactness.tsv"))

## ---- 3. strata of TSS spacing (two-promoter genes)
A[, sp_bin := cut(spacing, c(-1, 2000, 10000, 50000, Inf), labels = c("<=2kb", "2-10kb", "10-50kb", ">50kb"))]
strat <- rbindlist(lapply(levels(A$sp_bin), function(b) {
  d <- A[sp_bin == b & route %in% c("onoff_S_minor", "mu")]
  rS <- ratio_table(d[, .(gene, group = route, num = f_S, den = 1)], c("onoff_S_minor", "mu"), 1000L, 31L, paste("at S", b))
  rN <- ratio_table(d[, .(gene, group = route, num = f_N, den = 1)], c("onoff_S_minor", "mu"), 1000L, 31L, paste("at N", b))
  rsh <- ratio_table(d[f_S + f_N > 0, .(gene, group = route, num = f_N, den = f_S + f_N)], c("onoff_S_minor", "mu"), 1000L, 31L, paste("share at N", b))
  rbind(rbindlist(lapply(list(rS, rN, rsh), function(r) r$groups[, .(label, group, n_signals, estimate, lo, hi, contrast = NA_character_, p_boot = NA_real_)])),
        rbindlist(lapply(list(rS, rN, rsh), function(r) r$contrasts[, .(label, group = NA_character_, n_signals = NA_real_, estimate = difference, lo, hi, contrast, p_boot)])))
}))
writeT(strat, file.path(OUT, "A1_aligned_by_spacing.tsv"))

## ---- 4. quasi-binomial model for the share at N (major) among promoter-proximal candidates
M <- A[route %in% c("onoff_S_minor", "mu") & f_S + f_N > 0]
M[, `:=`(y = f_N / (f_S + f_N), wt = f_S + f_N, onoff = as.integer(route != "mu"), lsp = log10(spacing + 1))]
fitq <- function(d) coef(suppressWarnings(glm(y ~ onoff + lsp + maf + nlp + log(n_donors), family = quasibinomial(), weights = wt, data = d)))["onoff"]
est <- fitq(M); ug <- unique(M$gene); set.seed(41)
bt <- replicate(500, { g <- sample(ug, length(ug), TRUE); d <- M[data.table(gene = g), on = "gene", allow.cartesian = TRUE]; tryCatch(fitq(d), error = function(e) NA) })
glm_res <- data.table(model = "share at major promoter ~ on/off + log10 spacing + MAF + -log10 P + log n donors (quasibinomial; two-promoter genes; on/off with minor switching vs mu)",
                      n_signals = nrow(M), coef_onoff_logodds = est, odds_ratio = exp(est), or_lo = exp(quantile(bt, 0.025, na.rm = TRUE)), or_hi = exp(quantile(bt, 0.975, na.rm = TRUE)),
                      p_boot = max(2 * min(mean(bt <= 0, na.rm = TRUE), mean(bt >= 0, na.rm = TRUE)), 1 / 500))
writeT(glm_res, file.path(OUT, "A1_glm_share_major.tsv"))

## ---- 5. >= 3-promoter genes: own vs other by class; nearest other = major?
T3 <- S1[n_prom_gene >= 3]
r3 <- rbindlist(lapply(c("f_own", "f_other", "f_near_major"), function(v) {
  r <- ratio_table(T3[, .(gene, group = type, num = get(v), den = 1)], c("mu", "p0", "p1"), B, 51L, paste("three_plus", v))
  rbind(r$groups[, .(label, group, n_signals, estimate, lo, hi, contrast = NA_character_, p_boot = NA_real_)],
        r$contrasts[, .(label, group = NA_character_, n_signals = NA_real_, estimate = difference, lo, hi, contrast, p_boot)]) }))
writeT(r3, file.path(OUT, "A1_three_plus_promoter_genes.tsv"))

## ---- 6. within-unit pairs (panel-e set, 'within' pairs): mu and boundary signal of the same promoter x tissue
pr <- fr(file.path(IN, "delivery/crosstype_conditional_pairs_all.tsv.gz"), select = c("pair_id", "pair_kind", "tissue", "mu_promoter", "b_promoter", "b_type", "mu_k", "b_k"))[pair_kind %in% c("within", "exact_cross")]
pt <- fr(file.path(IN, "delivery/signals_pertype.tsv.gz"), select = c("tissue", "promoter", "type", "k", "mu_F_logp", "nu_firth_logp", "tau_firth_logp"))
pr[, sid_mu := paste(tissue, mu_promoter, "mu", mu_k, sep = "|")][, sid_b := paste(tissue, b_promoter, b_type, b_k, sep = "|")]
pt[, sid := paste(tissue, promoter, type, k, sep = "|")]
pr <- merge(pr, pt[, .(sid_mu = sid, cross_at_mu_p0 = nu_firth_logp, cross_at_mu_p1 = tau_firth_logp)], by = "sid_mu", all.x = TRUE)
pr <- merge(pr, pt[, .(sid_b = sid, cross_at_b = mu_F_logp)], by = "sid_b", all.x = TRUE)
pr[, cross_at_mu := fifelse(b_type == "p0", cross_at_mu_p0, cross_at_mu_p1)]
pr[, type_specific := exp(cross_at_mu) > 0.05 & exp(cross_at_b) > 0.05]
pm <- S1[, .(sid, gene, n_prom_gene, own_minor, f_own, f_other, f_S, f_N)]
pr <- merge(pr, pm[, .(sid_mu = sid, gene, n_prom_gene, mu_own = f_own, mu_other = f_other, mu_S = f_S, mu_N = f_N)], by = "sid_mu")
pr <- merge(pr, pm[, .(sid_b = sid, b_own = f_own, b_other = f_other, b_S = f_S, b_N = f_N)], by = "sid_b")
writeT(pr, file.path(OUT, "A1_within_unit_pairs.tsv"))
pair_stat <- function(d, a, b, lab, seed) {   # paired mean difference, gene bootstrap
  d <- d[!is.na(get(a)) & !is.na(get(b))]; if (nrow(d) < 10) return(NULL)
  dd <- d[, .(gene, x = get(a) - get(b))]; e <- mean(dd$x); ug <- unique(dd$gene); set.seed(seed)
  bt <- replicate(2000, { g <- sample(ug, length(ug), TRUE); mean(dd[data.table(gene = g), on = "gene", allow.cartesian = TRUE]$x) })
  data.table(comparison = lab, n_pairs = nrow(d), mean_a = mean(d[[a]]), mean_b = mean(d[[b]]), diff = e,
             lo = quantile(bt, 0.025), hi = quantile(bt, 0.975), p_boot = max(2 * min(mean(bt <= 0), mean(bt >= 0)), 1 / 2000))
}
ps <- rbind(
  pair_stat(pr, "b_own", "mu_own", "all within pairs: boundary vs mu, at own promoter", 61L),
  pair_stat(pr, "b_other", "mu_other", "all within pairs: boundary vs mu, at another promoter", 62L),
  pair_stat(pr[n_prom_gene == 2], "b_S", "mu_S", "two-promoter pairs: boundary at S vs mu at minor", 63L),
  pair_stat(pr[n_prom_gene == 2], "b_N", "mu_N", "two-promoter pairs: boundary at N vs mu at major", 64L),
  pair_stat(pr[type_specific == TRUE], "b_own", "mu_own", "type-specific pairs: at own promoter", 65L),
  pair_stat(pr[type_specific == TRUE], "b_other", "mu_other", "type-specific pairs: at another promoter", 66L))
writeT(ps, file.path(OUT, "A1_within_unit_pair_tests.tsv"))

## ---- 7. de-duplicated sensitivity: one signal per gene and route group (strongest), aligned shares
Ad <- A[order(-nlp)][!duplicated(paste(gene, route))]
rd <- ratio_table(Ad[f_S + f_N > 0, .(gene, group = route, num = f_N, den = f_S + f_N)], c("onoff_S_minor", "mu"), B, 71L, "dedup one per gene: share at N")
rdS <- ratio_table(Ad[, .(gene, group = route, num = f_S, den = 1)], c("onoff_S_minor", "mu"), B, 72L, "dedup one per gene: at S")
rdN <- ratio_table(Ad[, .(gene, group = route, num = f_N, den = 1)], c("onoff_S_minor", "mu"), B, 73L, "dedup one per gene: at N")
writeT(rbind(rbindlist(lapply(list(rdS, rdN, rd), function(r) r$groups[, .(label, group, n_signals, estimate, lo, hi, contrast = NA_character_, p_boot = NA_real_)])),
             rbindlist(lapply(list(rdS, rdN, rd), function(r) r$contrasts[, .(label, group = NA_character_, n_signals = NA_real_, estimate = difference, lo, hi, contrast, p_boot)]))),
       file.path(OUT, "A1_dedup_one_per_gene.tsv"))

## ---- 8. raw class x role table (the preview of 2026-11-03, now with candidate variants and CIs)
S1[, role := fifelse(n_prom_gene == 2, fifelse(own_minor, "minor", "major"), "three_plus")]
rr <- rbindlist(lapply(c("f_own", "f_other"), function(v) {
  d <- S1[, .(gene, group = paste(type, role, sep = "_"), num = get(v), den = 1)]
  g <- sort(unique(d$group)); ratio_table(d, g, 1000L, 81L, v)$groups }))
writeT(rr, file.path(OUT, "A1_class_by_role.tsv"))

## ---- report
g1 <- function(lab, grp) { x <- alg[label == lab & group == grp]; sprintf("%s [%s, %s] (n = %d)", pct(x$estimate), pct(x$lo), pct(x$hi), as.integer(x$n_signals)) }
c1 <- function(lab, ctr) { x <- alg[label == lab & contrast == ctr]; sprintf("%+.1f pp [%+.1f, %+.1f], P = %s", 100 * x$estimate, 100 * x$lo, 100 * x$hi, format(x$p_boot, digits = 2)) }
L <- c("# A1 -- variant location relative to all promoters of the gene", "",
  sprintf("Signals: %d; candidate variant rows: %d (lead + r2 >= 0.8 proxies, weighted 1/n per signal). Gene-clustered bootstrap, B = %d.", nrow(S1), nrow(V), B), "",
  "## Two-promoter genes, aligned frame (1 kb, exclusive)",
  paste0("- on/off, minor switching: at S ", g1("at S (switching; mu: minor) 1kb", "onoff_S_minor"), "; at N ", g1("at N (non-switching; mu: major) 1kb", "onoff_S_minor")),
  paste0("- mu: at minor ", g1("at S (switching; mu: minor) 1kb", "mu"), "; at major ", g1("at N (non-switching; mu: major) 1kb", "mu")),
  paste0("- share at the major/non-switching promoter among promoter-proximal candidates: on/off ", g1("share at N among promoter-proximal 1kb", "onoff_S_minor"),
         "; mu ", g1("share at N among promoter-proximal 1kb", "mu"), "; difference ", c1("share at N among promoter-proximal 1kb", "onoff_S_minor - mu")),
  sprintf("- quasi-binomial OR (on/off vs mu) for the share at major, adjusted for spacing, MAF, -log10 P, n donors: %.2f [%.2f, %.2f], P = %s",
          glm_res$odds_ratio, glm_res$or_lo, glm_res$or_hi, format(glm_res$p_boot, digits = 2)), "",
  "See A1_*.tsv for class-level fractions, >= 3-promoter genes, spacing strata, exactness split, within-unit pairs and de-duplication.")
writeLines(L, file.path(OUT, "A1_report.md"))
cat(L, sep = "\n"); cat("\nend", format(Sys.time()), "\n")
