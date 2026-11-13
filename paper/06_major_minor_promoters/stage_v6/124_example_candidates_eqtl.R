## 124_example_candidates_eqtl.R -- example loci for Fig. 5 h-j, re-chosen on v4 (author 2026-11-04): the gene-output column must
## show its own test, and each example must say whether the lead SNP is a GTEx eQTL of the gene (nominal P and slope of the
## lead, or the best r2 >= 0.8 proxy, from the GTEx v11 eQTL subset) -- not a coloc call. Fixed rules (two-promoter genes,
## shared leads removed, lead MAF >= 0.2, >= 10 donors per genotype group):
##   h  p0/p1, minor promoter only (pat minor_only); gene output changes (|z_gene| >= 4, covariate-adjusted, from B1v2);
##      GTEx eQTL nominal P < 1e-5, same direction as gene output;
##   i  μ, major promoter only; |z_gene| >= 4; GTEx eQTL P < 1e-5, same direction;
##   j  μ, both promoters in opposite directions (pat opposite) with |z_minor|, |z_major| >= 3; gene output flat (|z_gene| < 1);
##      GTEx eQTL P > 0.05 (not an eQTL of the gene).
## Ranked by min(|z| of the defining effects); top 5 distinct genes each; contact sheets with all read-outs printed.
suppressPackageStartupMessages({ library(data.table); library(VariantAnnotation); library(ggplot2); library(patchwork); library(scales) })
source(file.path(Sys.getenv("F5S", "/home/e1327989/f5_stage_20261106_v6"), "code/00_helpers.R"))
OUT <- file.path(OUTROOT, "FIG5_example_candidates_v4"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
K <- fr(file.path(OUTROOT, "K_overlap_multi_location/K_signals.tsv.gz"))[frame == "two" & shared == FALSE]
a1 <- fr(file.path(OUTROOT, "A1_location/A1_signal_location.tsv"), select = c("sid", "maf", "nlp"))
B1 <- fr(file.path(OUTROOT, "B1v2_multipromoter/B1v2_signal_patterns.tsv.gz"), select = c("sid", "z_T", "b_T", "se_T"))
K <- merge(merge(K, a1, by = "sid"), B1, by = "sid")
K[, gene_base := sub("\\..*$", "", gene)]
## GTEx eQTL of the lead (or best r2 >= 0.8 proxy) for this gene in this tissue
PX <- fr(file.path(IN, "delivery/proxies_pertype.tsv.gz"), select = c("tissue", "promoter", "type", "k", "signal_snp", "proxy_snp", "r2"))[r2 >= 0.8]
PX[, sid := paste(tissue, promoter, type, k, sep = "|")]
eq_of <- function(d) { rbindlist(lapply(split(d, d$tissue), function(x) { f <- file.path(IN, "coloc/eqtl_subset", paste0(x$tissue[1], ".tsv.gz"))
  if (!file.exists(f)) return(NULL); E <- fr(f)[gene_base %in% x$gene_base]
  cand <- rbind(x[, .(sid, gene_base, snp = lead, r2 = 1)], merge(x[, .(sid, gene_base)], PX[, .(sid, snp = proxy_snp, r2)], by = "sid"))
  m <- merge(cand, E, by.x = c("gene_base", "snp"), by.y = c("gene_base", "variant_id"))
  m[order(sid, -r2, pval_nominal)][, .SD[1], by = sid][, .(sid, eqtl_snp = snp, eqtl_r2 = r2, eqtl_p = pval_nominal, eqtl_slope = slope)] })) }
CND <- rbind(K[kind == "p0p1" & pat == "minor_only" & maf >= 0.2 & abs(z_T) >= 4][, set := "h_p0p1_minor_switch"][, score := pmin(abs(zM), abs(z_T))],
             K[kind == "mu" & pat == "major_only" & maf >= 0.2 & abs(z_T) >= 4][, set := "i_mu_major_only"][, score := pmin(abs(zJ), abs(z_T))],
             K[kind == "mu" & pat == "opposite" & maf >= 0.2 & abs(z_T) < 1 & abs(zM) >= 3 & abs(zJ) >= 3][, set := "j_mu_redistribution"][, score := pmin(abs(zM), abs(zJ))])
CND <- merge(CND, eq_of(CND), by = "sid", all.x = TRUE)
CND[, eq_ok := fifelse(set == "j_mu_redistribution", is.na(eqtl_p) | eqtl_p > 0.05, !is.na(eqtl_p) & eqtl_p < 1e-5 & sign(eqtl_slope) == sign(b_T))]
cat("candidates passing:\n"); print(CND[, .(all = .N, eq_ok = sum(eq_ok, na.rm = TRUE), eqtl_missing = sum(is.na(eqtl_p))), by = set])
SEL <- CND[eq_ok == TRUE][order(set, -score)][, .SD[!duplicated(gene)][1:min(5, .N)], by = set]
info <- fread(file.path(IN, "v11/Promoter_info_filter_5.txt")); setnames(info, 1, "pid")
x <- readRDS(file.path(IN, "v11/tabs_att_filter_In.rds")); abs <- x$abs; rm(x); invisible(gc())
theme_set(theme_classic(base_size = 6.5, base_family = "Arial") + theme(axis.text = element_text(size = 5.6, colour = "black"), plot.title = element_text(size = 6, hjust = 0.5),
  plot.subtitle = element_text(size = 5.4, colour = "grey25"), legend.position = "none", plot.margin = margin(1, 2, 1, 2)))
UST <- c("unused" = "grey85", "partly" = "#9AA8B8", "exclusive" = "#3F4F63"); rows <- list(); RD <- list()
for (i in seq_len(nrow(SEL))) { r <- SEL[i]; td <- readRDS(file.path(IN, "freeze", paste0(r$tissue, ".rds")))
  ch <- sub("_.*$", "", r$lead); pos <- as.integer(sub("^chr[^_]+_([0-9]+)_.*$", "\\1", r$lead))
  v <- readVcf(TabixFile(file.path(IN, "vcf", paste0(ch, ".vcf.gz"))), "hg38", ScanVcfParam(which = GRanges(ch, IRanges(pos, pos)), samples = td$U, fixed = NA, info = NA, geno = "GT"))
  gt <- geno(v)$GT[, td$U, drop = FALSE]; gt <- gt[rownames(gt) == r$lead, , drop = FALSE][1, ]
  dose <- fifelse(gt %in% c("0|0", "0/0"), 0L, fifelse(gt %in% c("0|1", "1|0", "0/1", "1/0"), 1L, fifelse(gt %in% c("1|1", "1/1"), 2L, NA_integer_)))
  minor_p <- if (r$minor_is_own) r$promoter else r$partner; major_p <- if (r$minor_is_own) r$partner else r$promoter
  pr <- union(c(minor_p, major_p), intersect(info[geneId == r$gene & internalPromoter == FALSE, pid], colnames(abs)))
  M <- as.matrix(abs[td$row_id, pr, drop = FALSE]); ok <- !is.na(dose)
  W <- data.table(dose = dose[ok], minor = M[ok, minor_p], major = M[ok, major_p], tot = rowSums(M[ok, , drop = FALSE]), gene_out = log2(rowSums(2^M[ok, , drop = FALSE] - 1) + 1))
  W[, usage := fifelse(tot > 0, minor / tot, NA_real_)]
  if (min(table(factor(W$dose, levels = 0:2))) < 10) { rm(td); next }
  ref <- strsplit(r$lead, "_")[[1]][3]; alt <- strsplit(r$lead, "_")[[1]][4]; gl <- c(paste0(ref, "/", ref), paste0(ref, "/", alt), paste0(alt, "/", alt))
  W[, glab := factor(gl[dose + 1], levels = gl)]
  colx <- c(h_p0p1_minor_switch = "#726CB1", i_mu_major_only = "#8E4618", j_mu_redistribution = "#D6874C")[[r$set]]
  act <- function(v, ttl) ggplot(W, aes(glab, get(v))) + geom_jitter(width = 0.2, height = 0, size = 0.3, alpha = 0.4, colour = colx, stroke = 0) +
    geom_boxplot(outlier.shape = NA, width = 0.48, linewidth = 0.25, fill = NA) + labs(x = NULL, y = NULL, title = ttl)
  S <- W[!is.na(usage), .(st = fifelse(usage == 0, "unused", fifelse(usage == 1, "exclusive", "partly")), glab)][, .N, by = .(glab, st)][, f := N / sum(N), by = glab]
  S[, st := factor(st, levels = rev(names(UST)))]
  p3 <- ggplot(S, aes(glab, f, fill = st)) + geom_col(width = 0.6) + scale_fill_manual(values = UST) + scale_y_continuous(labels = percent, breaks = c(0, 0.5, 1)) + labs(x = NULL, y = NULL, title = "minor usage")
  G <- W[, .(m = mean(gene_out), lo = mean(gene_out) - 1.96 * sd(gene_out) / sqrt(.N), hi = mean(gene_out) + 1.96 * sd(gene_out) / sqrt(.N)), by = glab]
  p4 <- ggplot(W, aes(glab, gene_out)) + geom_violin(fill = colx, alpha = 0.18, colour = NA, scale = "width", width = 0.75) + geom_line(data = G, aes(glab, m, group = 1), linewidth = 0.35) +
    geom_pointrange(data = G, aes(glab, m, ymin = lo, ymax = hi), size = 0.1, linewidth = 0.35) + labs(x = NULL, y = NULL, title = "gene output")
  pg <- 2 * pnorm(-abs(r$z_T))
  sub <- sprintf("%s %d · %s · %s · MAF %.2f | z minor %+.1f, major %+.1f, gene %+.1f (P = %.1e) | GTEx eQTL: %s", r$set, i, gsub("_", " ", r$tissue), r$lead, r$maf, r$zM, r$zJ, r$z_T, pg,
                 if (is.na(r$eqtl_p)) "lead not in GTEx subset" else sprintf("P = %.1e, slope %+.2f%s", r$eqtl_p, r$eqtl_slope, if (r$eqtl_r2 < 1) sprintf(" (proxy r2 %.2f)", r$eqtl_r2) else ""))
  rows[[length(rows) + 1]] <- wrap_elements(full = (act("minor", "minor activity") | act("major", "major activity") | p3 | p4) + plot_annotation(subtitle = sub))
  names(rows)[length(rows)] <- r$set
  RD[[length(RD) + 1]] <- r[, .(set, sid, gene, tissue, lead, maf, zM, zJ, z_T, p_gene = pg, eqtl_snp, eqtl_r2, eqtl_p, eqtl_slope)]
  rm(td); invisible(gc()) }
R <- rbindlist(RD); writeT(R, file.path(OUT, "candidates_v4_eqtl.tsv")); print(R)
for (s in unique(names(rows))) { ii <- which(names(rows) == s)
  ggsave(file.path(OUT, paste0("candidates_", s, ".png")), wrap_plots(rows[ii], ncol = 1), width = 183, height = 34 * length(ii), units = "mm", dpi = 200, bg = "white") }
cat("done\n")
