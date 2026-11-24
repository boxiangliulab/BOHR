## 72_A2_eqtl_direction.R -- SNP-level regulatory routes, A2: gene-level consequence of each signal.
## (a) Direction: does the allele that raises the usage of the signal's promoter raise or lower total gene expression?
##     eQTL = GTEx v11 all-pairs (rows cut by the Figure 4e coloc step: genes of the signals, variants within the signal
##     windows; eqtl_subset/<tissue>.tsv.gz). Lead used when present, otherwise the best r2 >= 0.8 proxy; the proxy's slope
##     is put on the lead's ALT allele by the sign of r (PLINK r, ALT dosages).
##     Usage direction of a signal = sign of its own-type estimate: mu: sign(mu_est) (interior mean, logit); p0: -sign(nu_est)
##     (more zeros = lower usage); p1: sign(tau_est) (more ones = higher usage). Checked against the empirical genotype-group means.
## (b) Aligned frame for two-promoter genes: the direction of the SWITCHING promoter's usage (p0: own; p1: the other, whose
##     usage is 1 - own) and, for mu, of the minor promoter's usage.
## (c) Splice-site QC (GENCODE v26): candidate variants inside the splice region of the first intron of the signal's own
##     promoter or of another promoter of the gene (donor: last 3 exonic + first 8 intronic bp; acceptor: last 20 intronic +
##     first 3 exonic bp). proActiv quantifies a promoter from its first-intron junction reads.
## Outputs: out/A2_eqtl_sqtl/*.tsv, A2_report.md
source(file.path(Sys.getenv("F5S", "/home/e1327989/f5_stage_20261106_v6"), "code/00_helpers.R"))
OUT <- file.path(OUTROOT, "A2_eqtl_sqtl"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setDTthreads(2); B <- 2000L
cat("start", format(Sys.time()), "\n")

sig <- fr(file.path(IN, "delivery/signals.tsv.gz"), select = c("tissue", "promoter", "type", "k", "snp", "mu_est", "nu_est", "tau_est", "logp_at_selection", "maf"))
sp  <- fr(file.path(IN, "delivery/signals_partner.tsv.gz"), select = c("tissue", "promoter", "type", "k", "mirror_role", "exact_mirror", "tss", "partner_tss"))
sig <- merge(sig, sp, by = c("tissue", "promoter", "type", "k"))
fl  <- fr(file.path(IN, "mirror/promoter_flags.tsv.gz"))[, .(tissue, promoter, gene, n_prom_gene)]
sig <- merge(sig, fl, by = c("tissue", "promoter"))
stopifnot(nrow(sig) > 0)   # the per-signal table of delivery v6 must not be empty
sig[, gene_base := sub("\\..*$", "", gene)]
sig[, sid := paste(tissue, promoter, type, k, sep = "|")]
sig[, dir := fifelse(type == "mu", sign(mu_est), fifelse(type == "p0", -sign(nu_est), sign(tau_est)))]
## empirical check from the genotype-group means at the lead
gpm <- fr(file.path(IN, "delivery/genotype_pheno.tsv.gz"))
gpm[, sid := paste(tissue, promoter, type, k, sep = "|")]
gpm[, emp := sign(((n1 * m1 + 2 * n2 * m2) / (n0 + n1 + n2) - (n1 + 2 * n2) / (n0 + n1 + n2) * (n0 * m0 + n1 * m1 + n2 * m2) / (n0 + n1 + n2)))]
sig <- merge(sig, gpm[, .(sid, emp)], by = "sid", all.x = TRUE)
agree <- sig[, .(n = .N, n_dir_na = sum(is.na(dir) | dir == 0), agree_with_group_means = mean(dir == emp, na.rm = TRUE)), by = type]
print(agree); writeT(agree, file.path(OUT, "A2_direction_check.tsv"))
## major / minor = higher / lower mean usage POOLED over all tissues (= the mirror table's canonical / partner, as Figure 4).
sig[, own_minor := fifelse(n_prom_gene == 2, mirror_role == "partner", NA)]
sig[n_prom_gene == 2, route := fifelse(type == "mu", "mu", fifelse((type == "p0") == own_minor, "onoff_S_minor", "onoff_S_major"))]
sig[n_prom_gene == 2, dirS := fifelse(type == "mu", fifelse(own_minor, dir, -dir), fifelse(type == "p0", dir, -dir))]
sig[, role := fifelse(n_prom_gene == 2, fifelse(own_minor, "minor", "major"), "three_plus")]

## ---- eQTL lookup
gp <- fr(file.path(IN, "delivery/genotype_pheno_proxies.tsv.gz"), select = c("tissue", "promoter", "type", "k", "signal_snp", "proxy_snp", "r", "r2"))[r2 >= 0.8]
gp[, sid := paste(tissue, promoter, type, k, sep = "|")]
gp <- merge(gp, sig[, .(sid, gene_base)], by = "sid")
EQ <- rbindlist(lapply(sort(unique(sig$tissue)), function(t) {
  f <- file.path(IN, "coloc/eqtl_subset", paste0(t, ".tsv.gz")); if (!file.exists(f)) { cat("no eQTL subset for", t, "\n"); return(NULL) }
  e <- fr(f); g <- gp[tissue == t]
  m <- merge(g, e, by.x = c("gene_base", "proxy_snp"), by.y = c("gene_base", "variant_id"))
  m[order(-(proxy_snp == signal_snp), -r2)][!duplicated(sid)] }))
cat("signals with an eQTL row at the lead or a proxy:", nrow(EQ), "of", nrow(sig), "\n")
EQ[, `:=`(z_eq = sign(r) * slope / slope_se, eq_at_lead = proxy_snp == signal_snp)]
sig <- merge(sig, EQ[, .(sid, eq_snp = proxy_snp, eq_r2 = r2, eq_at_lead, z_eq, eq_p = pval_nominal)], by = "sid", all.x = TRUE)
co <- if (file.exists(file.path(IN, "coloc/coloc_eqtl_signals.tsv.gz"))) fr(file.path(IN, "coloc/coloc_eqtl_signals.tsv.gz"), select = c("tissue", "promoter", "type", "k", "call", "PP.H4", "egene")) else   # without the eQTL coloc table: an empty table with the same columns
  data.table(tissue = character(), promoter = character(), type = character(), k = integer(), call = character(), PP.H4 = numeric(), egene = logical())
co[, sid := paste(tissue, promoter, type, k, sep = "|")]
sig <- merge(sig, co[, .(sid, coloc_call = call, PP.H4, egene)], by = "sid", all.x = TRUE)
sig[, `:=`(conc_own = sign(z_eq) == dir, conc_S = sign(z_eq) == dirS)]
writeT(sig[, .(sid, tissue, promoter, type, k, gene, n_prom_gene, role, route, dir, emp, dirS, eq_snp, eq_r2, eq_at_lead, z_eq, eq_p, coloc_call, PP.H4, egene, conc_own, conc_S, maf)],
       file.path(OUT, "A2_signal_eqtl.tsv"))

## ---- summaries
tabs <- list()
add <- function(d, groups, lab, seed) { r <- ratio_table(d, groups, B, seed, lab)
  tabs[[length(tabs) + 1]] <<- rbind(r$groups[, .(label, group, n_signals, estimate, lo, hi, contrast = NA_character_, p_boot = NA_real_)],
                                      r$contrasts[, .(label, group = NA_character_, n_signals = NA_real_, estimate = difference, lo, hi, contrast, p_boot)]) }
cls <- c("mu", "p0", "p1")
x <- sig[!is.na(z_eq)]
add(x[, .(gene, group = type, num = as.numeric(abs(z_eq) >= 1.96), den = 1)], cls, "eQTL nominal (|z| >= 1.96) at lead/proxy", 1L)
add(x[, .(gene, group = type, num = as.numeric(eq_p < 1e-4), den = 1)], cls, "eQTL P < 1e-4", 2L)
for (thr in c(1.96, 3.89)) {
  y <- x[abs(z_eq) >= thr & !is.na(dir) & dir != 0]
  add(y[, .(gene, group = type, num = as.numeric(conc_own), den = 1)], cls, sprintf("concordant with OWN usage | |z_eq| >= %.2f", thr), 10L + thr)
  add(y[, .(gene, group = paste(type, role, sep = "_"), num = as.numeric(conc_own), den = 1)], sort(unique(y[, paste(type, role, sep = "_")])),
      sprintf("concordant with OWN usage by role | |z_eq| >= %.2f", thr), 20L + thr)
  z <- y[!is.na(route)]
  add(z[, .(gene, group = route, num = as.numeric(conc_S), den = 1)], c("onoff_S_minor", "onoff_S_major", "mu"),
      sprintf("two-promoter aligned: concordant with SWITCHING (mu: minor) promoter usage | |z_eq| >= %.2f", thr), 30L + thr)
}
if (nrow(co) > 0) {   # without the eQTL coloc table the coloc summaries are skipped
y <- x[PP.H4 >= 0.8 & !is.na(dir) & dir != 0]
add(y[, .(gene, group = type, num = as.numeric(conc_own), den = 1)], cls, "concordant with OWN usage | coloc H4 >= 0.8", 41L)
z <- y[!is.na(route)]
add(z[, .(gene, group = route, num = as.numeric(conc_S), den = 1)], c("onoff_S_minor", "onoff_S_major", "mu"), "two-promoter aligned: concordant with switching/minor | coloc H4 >= 0.8", 42L)
## coloc calls by class
cc <- sig[!is.na(coloc_call), .N, by = .(type, coloc_call)][, share := N / sum(N), by = type]
writeT(dcast(cc, coloc_call ~ type, value.var = "share"), file.path(OUT, "A2_coloc_calls_by_class.tsv"))
} else cat("no eQTL coloc table: coloc summaries skipped\n")
TAB <- rbindlist(tabs); writeT(TAB, file.path(OUT, "A2_eqtl_direction_summary.tsv"))

## ---- (c) splice-site QC
info <- fread(file.path(IN, "v11/Promoter_info_filter_5.txt")); setnames(info, 1, "pid")
ex <- fread(cmd = sprintf("zcat %s | awk -F'\\t' '$3==\"exon\"' | cut -f1,4,5,7,9", file.path(IN, "gtf/gencode.v26.annotation.gtf.gz")),
            header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "strand", "attr"))
ex[, tx := sub('^.*transcript_id "([^"]+)".*$', "\\1", attr)][, exn := as.integer(sub('^.*exon_number ([0-9]+);.*$', "\\1", attr))][, attr := NULL]
ex[, txb := sub("\\..*$", "", tx)]
e12 <- ex[exn %in% 1:2]
w1 <- dcast(e12, txb + chr + strand ~ exn, value.var = c("start", "end"))
setnames(w1, c("start_1", "end_1", "start_2", "end_2"), c("s1", "e1", "s2", "e2"))
w1 <- w1[!is.na(s2)]
## donor / acceptor regions of the first intron
w1[strand == "+", `:=`(don_lo = e1 - 2L, don_hi = e1 + 8L, acc_lo = s2 - 20L, acc_hi = s2 + 2L)]
w1[strand == "-", `:=`(don_lo = s1 - 8L, don_hi = s1 + 2L, acc_lo = e2 - 2L, acc_hi = e2 + 20L)]
ptx <- info[, .(txb = sub("\\..*$", "", unlist(strsplit(txId, "|", fixed = TRUE)))), by = .(promoter = pid, geneId)]
ptx <- merge(ptx, w1, by = "txb")
cat("promoters with a first intron from GENCODE v26:", uniqueN(ptx$promoter), "of", nrow(info), "\n")
gpv <- gp[, .(sid, tissue, gene_base, proxy_snp, chr = sub("_.*$", "", proxy_snp), pos = as.integer(sub("^chr[^_]+_([0-9]+)_.*$", "\\1", proxy_snp)))]
gpv[, w := 1 / .N, by = sid]
gpv <- merge(gpv, sig[, .(sid, own = promoter, type, gene)], by = "sid")
fl2 <- fl[, .(tissue, gene, promoter)]
reg <- rbind(ptx[, .(promoter, chr, lo = don_lo, hi = don_hi, site = "donor")], ptx[, .(promoter, chr, lo = acc_lo, hi = acc_hi, site = "acceptor")])
reg <- unique(reg)
## variants x promoters of the same gene in the tissue
vp <- merge(gpv, fl2, by = c("tissue", "gene"), allow.cartesian = TRUE)
vp <- merge(vp, reg, by = c("promoter", "chr"), allow.cartesian = TRUE)[pos >= lo & pos <= hi]
vp[, which_prom := fifelse(promoter == own, "own", "other")]
hit <- unique(vp[, .(sid, proxy_snp, which_prom, site)])
ss <- merge(gpv[, .(sid, proxy_snp, w, type, gene)], dcast(hit[, .(sid, proxy_snp, k = paste(which_prom, site, sep = "_"))], sid + proxy_snp ~ k, fun.aggregate = length),
            by = c("sid", "proxy_snp"), all.x = TRUE)
for (v in setdiff(names(ss), c("sid", "proxy_snp", "w", "type", "gene"))) set(ss, which(is.na(ss[[v]])), v, 0L)
ss[, any_own := 0L]; ss[, any_other := 0L]
for (v in grep("^own_", names(ss), value = TRUE)) ss[, any_own := pmax(any_own, as.integer(get(v) > 0))]
for (v in grep("^other_", names(ss), value = TRUE)) ss[, any_other := pmax(any_other, as.integer(get(v) > 0))]
sss <- ss[, .(f_splice_own = sum(w * any_own), f_splice_other = sum(w * any_other), any_cand_own = max(any_own), any_cand_other = max(any_other)), by = .(sid, type, gene)]
writeT(sss, file.path(OUT, "A2_splice_region_by_signal.tsv"))
spl <- rbindlist(lapply(c("f_splice_own", "f_splice_other", "any_cand_own", "any_cand_other"), function(v) ratio_table(sss[, .(gene, group = type, num = as.numeric(get(v)), den = 1)], cls, B, 51L, v)$groups))
writeT(spl, file.path(OUT, "A2_splice_region_summary.tsv"))

## ---- report
gv <- function(lab, grp) { x <- TAB[label == lab & group == grp]; if (!nrow(x)) return("NA"); sprintf("%s [%s, %s] (n = %d)", pct(x$estimate), pct(x$lo), pct(x$hi), as.integer(x$n_signals)) }
L <- c("# A2 -- gene-level consequence (eQTL direction) and splice-site QC", "",
  sprintf("Signals with an eQTL row at the lead or an r2 >= 0.8 proxy: %d of %d.", sum(!is.na(sig$z_eq)), nrow(sig)),
  sprintf("Own-type estimate sign agrees with the genotype-group mean slope: %s.", paste(sprintf("%s %s", agree$type, pct(agree$agree_with_group_means)), collapse = "; ")), "",
  "## Concordance of the eQTL with the signal's OWN promoter usage (|z_eq| >= 1.96)",
  paste0("- ", cls, ": ", sapply(cls, function(g) gv("concordant with OWN usage | |z_eq| >= 1.96", g))), "",
  "## Two-promoter genes, aligned: concordance with the SWITCHING promoter (mu: minor) (|z_eq| >= 1.96)",
  paste0("- ", c("onoff_S_minor", "onoff_S_major", "mu"), ": ", sapply(c("onoff_S_minor", "onoff_S_major", "mu"), function(g) gv("two-promoter aligned: concordant with SWITCHING (mu: minor) promoter usage | |z_eq| >= 1.96", g))), "",
  "## Splice-region QC (fraction of signals with any candidate variant in a first-intron splice region)",
  paste0("- own promoter: ", paste(sprintf("%s %s", spl[label == "any_cand_own"]$group, pct(spl[label == "any_cand_own"]$estimate, 2)), collapse = "; ")),
  paste0("- another promoter: ", paste(sprintf("%s %s", spl[label == "any_cand_other"]$group, pct(spl[label == "any_cand_other"]$estimate, 2)), collapse = "; ")))
writeLines(L, file.path(OUT, "A2_report.md")); cat(L, sep = "\n"); cat("\nend", format(Sys.time()), "\n")
