## 104_M_genotype_usage.R -- the second dimension (2026-11-03 night): fix ONE promoter and compare the usage
## distribution between genotype groups. No new statistics on the QTL itself; this shows what the three model parameters mean
## in donors.
##  M1 all signals, from the delivered tables: per genotype group g (n_g, mean usage m_g from genotype_pheno; interior
##     count and mean from signals_pe_moments) the exclusive count = n_g m_g - n_mid_g m_mid_g and the unused count =
##     n_g - n_mid_g - exclusive (checked for being whole numbers). Genotypes oriented so the usage-raising allele (dir_c) is
##     counted; group means over signals per class with gene-bootstrap CIs; groups with < 10 donors left out.
##  M2 example loci, chosen by a fixed rule (two-promoter genes, MAF >= 0.25, lead not shared with the other kind unless the
##     shared example, strongest -log10 P): p0 minor only; p1 minor only; p1 both together; mu minor only; mu major only;
##     one lead that is both a mu and a p0/p1 lead. Donor level: absolute activity (proActiv log2 normalised count) of every
##     promoter of the gene by genotype, and the target promoter's usage state and interior usage by genotype.
suppressPackageStartupMessages({ library(data.table); library(VariantAnnotation); library(Matrix); library(ggplot2); library(patchwork); library(ragg); library(scales) })
source(file.path(Sys.getenv("F5S", "/home/e1327989/f5_stage_20261106_v6"), "code/00_helpers.R"))
OUT <- file.path(OUTROOT, "M_genotype_usage"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE); B <- 2000L
CL <- c(p0 = "#4D9E9B", p1 = "#726CB1", mu = "#D6874C"); LABT <- c(p0 = "p0", p1 = "p1", mu = "μ")
th <- theme_classic(base_size = 8.5) + theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3), strip.background = element_blank(),
  strip.text = element_text(face = "bold", size = 8), plot.title = element_text(face = "bold", size = 9), plot.title.position = "plot",
  plot.subtitle = element_text(size = 7, colour = "grey30"), legend.title = element_blank(), legend.key.size = unit(3, "mm"))
save <- function(p, name, w, h) { ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 220, device = ragg::agg_png, bg = "white")
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf); cat("wrote", name, "\n") }

## ---------------- M1
gp <- fr(file.path(IN, "delivery/genotype_pheno.tsv.gz"))
pm <- fr(file.path(IN, "delivery/signals_pe_moments.tsv.gz"), select = c("tissue", "promoter", "type", "k", "n0_mid", "n1_mid", "n2_mid", "m0_mid", "m1_mid", "m2_mid"))
M <- merge(gp, pm, by = c("tissue", "promoter", "type", "k"))
M[, sid := paste(tissue, promoter, type, k, sep = "|")]
P <- fr(file.path(OUTROOT, "B1v2_multipromoter/B1v2_signal_patterns.tsv.gz"), select = c("sid", "gene", "two", "dir_c"))
M <- merge(M, P, by = "sid")
L <- rbindlist(lapply(0:2, function(g) M[, .(sid, gene, type, two, dir_c, g = g, n = get(paste0("n", g)), m = get(paste0("m", g)),
                                              nmid = get(paste0("n", g, "_mid")), mmid = get(paste0("m", g, "_mid")))]))
L[is.na(nmid), nmid := 0]
L[, n_one := n * m - nmid * fifelse(nmid > 0, mmid, 0)][, n_zero := n - nmid - n_one]
chk <- L[n > 0, .(max_dev_one = max(abs(n_one - round(n_one)), na.rm = TRUE), min_one = min(n_one, na.rm = TRUE), min_zero = min(n_zero, na.rm = TRUE))]
cat("whole-number check:\n"); print(chk)
L[, `:=`(n_one = pmax(round(n_one), 0), n_zero = pmax(round(n_zero), 0))]
L[, `:=`(f_zero = n_zero / n, f_one = n_one / n, f_mid = nmid / n, m_mid = fifelse(nmid >= 5, mmid, NA_real_))]
L[, go := fifelse(dir_c >= 0, g, 2L - g)]          # oriented: go = number of usage-raising alleles
L <- L[n >= 10]
fwrite(L, file.path(OUT, "M1_signal_genotype_states.tsv.gz"), sep = "\t")
gb <- function(d, v, seed) { ug <- unique(d$gene); x <- d[is.finite(get(v)), .(gene, y = get(v))]; set.seed(seed)
  est <- mean(x$y); bt <- replicate(B, { s <- sample(ug, length(ug), TRUE); mean(x[data.table(gene = s), on = "gene", allow.cartesian = TRUE, nomatch = 0]$y) })
  c(est = est, lo = unname(quantile(bt, 0.025, na.rm = TRUE)), hi = unname(quantile(bt, 0.975, na.rm = TRUE)), n = nrow(x)) }
vars <- c(f_zero = "not used\n(usage = 0)", f_mid = "partly used\n(0 < usage < 1)", f_one = "used exclusively\n(usage = 1)", m_mid = "usage when\npartly used")
S1 <- rbindlist(lapply(c("all", "two"), function(st) rbindlist(lapply(c("p0", "p1", "mu"), function(ty) rbindlist(lapply(0:2, function(gg) rbindlist(lapply(names(vars), function(v) {
  d <- L[type == ty & go == gg & (st == "all" | two == TRUE)]; b <- gb(d, v, 11L)
  data.table(set = st, type = ty, go = gg, var = v, est = b[["est"]], lo = b[["lo"]], hi = b[["hi"]], n = b[["n"]]) }))))))))
writeT(S1, file.path(OUT, "M1_class_genotype_means.tsv"))
## per signal change from the low to the high homozygote (or het when a homozygote group is too small)
W <- dcast(L, sid + gene + type + two ~ go, value.var = c("f_zero", "f_one", "m_mid"))
W[, `:=`(d_zero = f_zero_2 - f_zero_0, d_one = f_one_2 - f_one_0, d_mid = m_mid_2 - m_mid_0)]
W[is.na(d_zero), d_zero := f_zero_1 - f_zero_0]; W[is.na(d_one), d_one := f_one_1 - f_one_0]; W[is.na(d_mid), d_mid := m_mid_1 - m_mid_0]
fwrite(W, file.path(OUT, "M1_signal_changes.tsv.gz"), sep = "\t")
S2 <- W[, .(n = .N, zero_down_5pp = mean(d_zero <= -0.05, na.rm = TRUE), one_up_5pp = mean(d_one >= 0.05, na.rm = TRUE),
            mid_up_005 = mean(d_mid >= 0.05, na.rm = TRUE), med_d_zero = median(d_zero, na.rm = TRUE), med_d_one = median(d_one, na.rm = TRUE),
            med_d_mid = median(d_mid, na.rm = TRUE)), by = type][order(type)]
writeT(S2, file.path(OUT, "M1_change_shares.tsv")); print(S2)

## figure M1
S1p <- S1[set == "all"]; S1p[, `:=`(var = factor(vars[var], levels = vars), type = factor(LABT[type], levels = LABT),
                                     gl = factor(c("0", "1", "2")[go + 1], levels = c("0", "1", "2")))]
pM1 <- ggplot(S1p, aes(x = gl, y = est, colour = type, group = type)) + geom_line(linewidth = 0.6) + geom_point(size = 1.4) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, linewidth = 0.35) + facet_wrap(~var, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = setNames(CL, LABT)) + scale_y_continuous(labels = percent) +
  labs(x = "copies of the allele that raises the promoter's usage", y = "share of donors / usage",
       title = "One promoter, three genotype groups: how the usage distribution changes",
       subtitle = sprintf("All %s signals. Each line averages the signals of one class; each signal is read on its own promoter.", format(uniqueN(L$sid), big.mark = ","))) +
  th + theme(legend.position = "top")
save(pM1, "figM1_genotype_usage_states", 190, 80)

## ---------------- M2 example loci
K <- fr(file.path(OUTROOT, "K_overlap_multi_location/K_signals.tsv.gz"))[frame == "two"]
a1 <- fr(file.path(OUTROOT, "A1_location/A1_signal_location.tsv"), select = c("sid", "maf", "nlp", "n_donors"))
K <- merge(K, a1, by = "sid")
pick <- function(d, lab) { x <- d[maf >= 0.25][order(-nlp)][1]; x[, ex := lab]; x }
EX <- rbind(pick(K[type == "p0" & pat == "minor_only" & shared == FALSE], "p0: minor promoter only"),
            pick(K[type == "p1" & pat == "minor_only" & shared == FALSE], "p1: minor promoter only"),
            pick(K[type == "p1" & pat == "together" & shared == FALSE], "p1: both promoters together"),
            pick(K[type == "mu" & pat == "minor_only" & shared == FALSE], "μ: minor promoter only"),
            pick(K[type == "mu" & pat == "major_only" & shared == FALSE], "μ: major promoter only"))
## shared lead: a mu signal whose lead is identical to a p0/p1 lead in the same gene x tissue
sh <- merge(K[kind == "mu" & shared == TRUE & maf >= 0.25, .(sid_mu = sid, gene, tissue, lead, nlp)], K[kind == "p0p1", .(sid_b = sid, gene, tissue, lead, type_b = type)], by = c("gene", "tissue", "lead"))
if (nrow(sh)) { s1 <- sh[order(-nlp)][1]; x <- K[sid == s1$sid_mu]; x[, ex := sprintf("shared lead: μ and %s", s1$type_b)]; EX <- rbind(EX, x) }
fwrite(EX[, .(ex, sid, gene, tissue, type, lead, pat, maf, nlp, minor_is_own, promoter, partner)], file.path(OUT, "M2_examples.tsv"), sep = "\t"); print(EX[, .(ex, sid, lead, pat, maf, nlp)])

DF <- file.path(OUT, "M2_donor_values.tsv.gz")
if (!file.exists(DF)) {
e <- new.env(); load(file.path(IN, "v11/Promoter_info_all.RData"), envir = e); ia <- as.data.table(e$info_all, keep.rownames = "pid"); rm(e)
x <- readRDS(file.path(IN, "v11/tabs_att_filter_In.rds")); abs <- x$abs; rm(x); invisible(gc())
ia <- ia[internalPromoter == FALSE & pid %in% colnames(abs)]
gsym <- tryCatch(fread(file.path(IN, "v11/Promoter_info_filter_5.txt"))[, .(gene = geneId, sym = if ("geneName" %in% names(.SD)) geneName else geneId)], error = function(e) NULL)
donor <- list()
for (i in seq_len(nrow(EX))) { r <- EX[i]
  td <- readRDS(file.path(IN, "freeze", paste0(r$tissue, ".rds")))
  ch <- sub("_.*$", "", r$lead); pos <- as.integer(sub("^chr[^_]+_([0-9]+)_.*$", "\\1", r$lead))
  v <- readVcf(TabixFile(file.path(IN, "vcf", paste0(ch, ".vcf.gz"))), "hg38", ScanVcfParam(which = GRanges(ch, IRanges(pos, pos)), samples = td$U, fixed = NA, info = NA, geno = "GT"))
  gt <- geno(v)$GT[, td$U, drop = FALSE]; gt <- gt[rownames(gt) == r$lead, , drop = FALSE][1, ]
  dose <- fifelse(gt %in% c("0|0", "0/0"), 0L, fifelse(gt %in% c("0|1", "1|0", "0/1", "1/0"), 1L, fifelse(gt %in% c("1|1", "1/1"), 2L, NA_integer_)))
  pr <- intersect(unique(c(r$promoter, r$partner, ia[geneId == r$gene, pid])), colnames(abs)); cat(r$ex, r$sid, "promoters:", length(pr), "\n")
  if (!(r$promoter %in% pr)) { cat("  target not in matrix, skipped\n"); next }
  A <- as.matrix(abs[td$row_id, pr, drop = FALSE]); tot <- rowSums(A)
  minor_p <- if (r$minor_is_own) r$promoter else r$partner; major_p <- if (r$minor_is_own) r$partner else r$promoter
  rolev <- setNames(fifelse(pr == minor_p, "minor promoter", fifelse(pr == major_p, "major promoter", "other promoter (low)")), pr)
  u <- ifelse(tot > 0, A[, r$promoter] / tot, NA_real_)
  ref <- strsplit(r$lead, "_")[[1]][3]; alt <- strsplit(r$lead, "_")[[1]][4]
  glab <- c(paste0(ref, "/", ref), paste0(ref, "/", alt), paste0(alt, "/", alt))
  donor[[i]] <- rbind(data.table(ex = r$ex, sid = r$sid, what = "activity", prom = rep(pr, each = nrow(A)), role = rep(rolev, each = nrow(A)), dose = rep(dose, length(pr)), y = as.vector(A)),
                      data.table(ex = r$ex, sid = r$sid, what = "usage", prom = r$promoter, role = rolev[[r$promoter]], dose = dose, y = u), fill = TRUE)
  donor[[i]][, `:=`(glab = factor(glab[dose + 1], levels = glab), target = r$promoter, target_role = rolev[[r$promoter]], gene = r$gene, tissue = r$tissue, lead = r$lead, cls = r$type)]
  rm(td); invisible(gc()) }
DD <- rbindlist(donor, fill = TRUE)[!is.na(dose)]
fwrite(DD, DF, sep = "\t")
}
DD <- fr(DF)
## one figure per example: left activity of every promoter by genotype; right usage states + interior usage of the target
for (exn in unique(DD$ex)) { d <- DD[ex == exn]; cl <- d$cls[1]; colx <- CL[[cl]]
  d[, glab := factor(as.character(glab), levels = unique(as.character(glab)[order(dose)]))]
  da <- d[what == "activity"]; da[, role := factor(role, levels = c("minor promoter", "major promoter", "other promoter (low)"))]
  da <- da[role != "other promoter (low)" | prom %in% da[, .(m = mean(y)), by = prom][m > 0.5, prom]]
  pa <- ggplot(da, aes(x = glab, y = y)) + geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.3, fill = NA, colour = "grey35") +
    geom_jitter(width = 0.18, height = 0, size = 0.35, alpha = 0.35, colour = colx) + facet_wrap(~ role + prom, nrow = 1, labeller = label_wrap_gen(multi_line = TRUE)) +
    labs(x = NULL, y = "activity (log2 count)", title = "a  Absolute activity, every promoter") + th
  du <- d[what == "usage" & !is.na(y)]
  du[, state := factor(fifelse(y <= 0, "not used", fifelse(y >= 1, "used exclusively", "partly used")), levels = c("not used", "partly used", "used exclusively"))]
  st <- du[, .N, by = .(glab, state)][, f := N / sum(N), by = glab]; nn <- du[, .N, by = glab]
  pb <- ggplot(st, aes(x = glab, y = f, fill = state)) + geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
    geom_text(data = nn, aes(x = glab, y = 1.04, label = paste0("n = ", N)), inherit.aes = FALSE, size = 2.2) +
    scale_fill_manual(values = c("not used" = "grey85", "partly used" = alpha(colx, 0.55), "used exclusively" = colx)) +
    scale_y_continuous(labels = percent, limits = c(0, 1.08), breaks = c(0, 0.5, 1)) + labs(x = NULL, y = "donors", title = "b  Usage state of the target") + scale_x_discrete(drop = FALSE) + th + theme(legend.position = "top", legend.text = element_text(size = 6.5))
  pc <- ggplot(du[state == "partly used"], aes(x = glab, y = y)) + geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.3, fill = NA, colour = "grey35") +
    geom_jitter(width = 0.18, height = 0, size = 0.4, alpha = 0.4, colour = colx) + scale_y_continuous(limits = c(0, 1)) +
    scale_x_discrete(drop = FALSE) + labs(x = NULL, y = "usage", title = "c  Usage when partly used") + th
  p <- (pa | pb | pc) + plot_layout(widths = c(max(2, uniqueN(da$prom)) * 0.9, 1, 1)) +
    plot_annotation(title = sprintf("%s  ·  %s  ·  %s", exn, d$gene[1], gsub("_", " ", d$tissue[1])),
                    subtitle = sprintf("lead %s; target promoter %s (%s)", d$lead[1], d$target[1], d$target_role[1]),
                    theme = theme(plot.title = element_text(face = "bold", size = 9.5), plot.subtitle = element_text(size = 7.5, colour = "grey30")))
  save(p, paste0("figM2_example_", gsub("[^A-Za-z0-9]+", "_", sub("μ", "mu", exn))), 185, 72) }
cat("done\n")
