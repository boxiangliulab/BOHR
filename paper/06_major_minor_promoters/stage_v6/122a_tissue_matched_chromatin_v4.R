## 122a -- copy of Figure5_rerun_v3/scripts/28 rebuilt on the v3m unit set (mirror rule v3); inputs staged; outputs to out/PROMOTER_LEVEL_v4.
#!/usr/bin/env Rscript
# 28_tissue_matched_tss_chromatin.R -- TSS chromatin in the QTL's OWN tissue (PI request 2026-10-03).
# Unit level: (GTEx tissue, promoter). For every tested (tissue, promoter) of Figure 4's v3 freeze (estimable, not
# dedup-dropped), the state of the cCRE covering the promoter's TSS (25a) is read in the ENCODE v4 core-collection
# biosamples MATCHED to that GTEx tissue (mapping fixed below before any result was looked at):
#   exact       = the same tissue or cell type;
#   approximate = a related tissue of the same organ, or the main cell type (non-cortex brain -> dorsolateral
#                 prefrontal cortex; GE junction -> esophagus muscularis mucosa; whole blood -> CD14+ monocyte; skin ->
#                 keratinocyte);
#   unmatched   = adipose, bladder, fibroblasts, kidney, minor salivary gland, ovary, pituitary.
# State per unit, over its matched biosamples (same biosample for all marks; Z > 1.64): promoter-like when the share of
# accessible biosamples with high H3K4me3 >= the share with high H3K27ac and low H3K4me3 (> 0); enhancer-like when the
# latter is larger; accessible-unmarked; or no accessible element (no cCRE, or never accessible in the matched samples).
# Groups: (tissue, promoter) units with a significant mu / p0 / p1 QTL in THAT tissue, and tested units with none.
# Primary = exact + approximate; sensitivity = exact only. Gene-clustered bootstrap.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) }); setDTthreads(8)
F5S <- Sys.getenv("F5S", "/home/e1327989/f5_stage_20261106_v6"); WD <- file.path(F5S, "in/fig5_v3data"); FZ <- file.path(F5S, "in/freeze_v3c")
OUTD <- file.path(F5S, "out/PROMOTER_LEVEL_v4"); dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
rz <- function(f) fread(cmd = sprintf("zcat '%s'", f)); B <- 2000; SEED <- 20261003; TH <- 1.64
MAP <- rbindlist(list(
  data.table(gtex = "Adrenal_Gland", term = "adrenal gland", match = "exact"),
  data.table(gtex = "Artery_Aorta", term = c("ascending aorta", "thoracic aorta"), match = "exact"),
  data.table(gtex = "Artery_Coronary", term = "coronary artery", match = "exact"),
  data.table(gtex = "Artery_Tibial", term = "tibial artery", match = "exact"),
  data.table(gtex = c("Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Anterior_cingulate_cortex_BA24"), term = "dorsolateral prefrontal cortex", match = "exact"),
  data.table(gtex = c("Brain_Amygdala", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Hippocampus", "Brain_Hypothalamus",
                      "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra"),
             term = "dorsolateral prefrontal cortex", match = "approximate"),
  data.table(gtex = "Breast_Mammary_Tissue", term = "breast epithelium", match = "exact"),
  data.table(gtex = "Esophagus_Mucosa", term = "esophagus squamous epithelium", match = "exact"),
  data.table(gtex = "Esophagus_Muscularis", term = "esophagus muscularis mucosa", match = "exact"),
  data.table(gtex = "Esophagus_Gastroesophageal_Junction", term = "esophagus muscularis mucosa", match = "approximate"),
  data.table(gtex = "Heart_Left_Ventricle", term = c("heart left ventricle", "left ventricle myocardium inferior"), match = "exact"),
  data.table(gtex = "Heart_Atrial_Appendage", term = "right atrium auricular region", match = "exact"),
  data.table(gtex = "Colon_Sigmoid", term = "sigmoid colon", match = "exact"),
  data.table(gtex = "Colon_Transverse", term = c("transverse colon", "colonic mucosa"), match = "exact"),
  data.table(gtex = "Liver", term = "right lobe of liver", match = "exact"),
  data.table(gtex = "Lung", term = c("left lung", "lower lobe of left lung", "upper lobe of left lung"), match = "exact"),
  data.table(gtex = "Muscle_Skeletal", term = "gastrocnemius medialis", match = "exact"),
  data.table(gtex = "Nerve_Tibial", term = "tibial nerve", match = "exact"),
  data.table(gtex = "Pancreas", term = c("body of pancreas", "pancreas"), match = "exact"),
  data.table(gtex = "Prostate", term = "prostate gland", match = "exact"),
  data.table(gtex = "Small_Intestine_Terminal_Ileum", term = "Peyer's patch", match = "exact"),
  data.table(gtex = "Spleen", term = "spleen", match = "exact"),
  data.table(gtex = "Stomach", term = "stomach", match = "exact"),
  data.table(gtex = "Testis", term = "testis", match = "exact"),
  data.table(gtex = "Thyroid", term = "thyroid gland", match = "exact"),
  data.table(gtex = "Uterus", term = "uterus", match = "exact"),
  data.table(gtex = "Vagina", term = "vagina", match = "exact"),
  data.table(gtex = "Cells_EBV-transformed_lymphocytes", term = "GM12878", match = "exact"),
  data.table(gtex = "Whole_Blood", term = "CD14-positive monocyte", match = "approximate"),
  data.table(gtex = c("Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg"), term = "keratinocyte", match = "approximate")))
Z <- lapply(c(k4 = "h3k4me3", k27 = "h3k27ac", dn = "dnase"), function(m) list(z = rz(file.path(WD, sprintf("zscores_%s_alltested.tsv.gz", m))),
                                                                             c = fread(file.path(WD, sprintf("zscores_%s_alltested_columns.tsv", m)))))
CM <- Z$k4$c[, .(biosample_value, term = sub(",.*$", "", celltype))]
common <- Reduce(intersect, lapply(Z, function(x) x$c$biosample_value))
MAPb <- merge(MAP, CM[biosample_value %in% common], by = "term", allow.cartesian = TRUE)
mat <- function(k) { cols <- Z[[k]]$c[match(common, biosample_value)]$col; m <- as.matrix(Z[[k]]$z[, ..cols]); rownames(m) <- Z[[k]]$z$acc; colnames(m) <- common; m }
A <- mat("k4"); K <- mat("k27")[rownames(A), ]; D <- mat("dn")[rownames(A), ]
OPEN <- D > TH; PROM <- OPEN & A > TH; ENH <- OPEN & K > TH & A <= TH
## tested (tissue, promoter) units and their QTL classes in that tissue
U <- rbindlist(lapply(list.files(FZ, pattern = "real.units.rds$", full.names = TRUE), function(f)
  as.data.table(readRDS(f))[estimable == TRUE & dedup_drop == FALSE, .(tissue, promoter, gene, type, q)]))
TP <- U[, .(gene0 = sub("\\..*$", "", gene[1]), mu = any(type == "mu" & !is.na(q) & q <= .05), p0 = any(type == "p0" & !is.na(q) & q <= .05),
            p1 = any(type == "p1" & !is.na(q) & q <= .05)), by = .(tissue, promoter)]
EL <- fread(file.path(WD, "alltested_tss_elements.tsv"))[, .(promoter, acc)]
TP <- merge(TP, EL, by = "promoter", all.x = TRUE)
state_for <- function(tp, match_set) {
  mp <- MAPb[match %in% match_set]
  out <- rbindlist(lapply(split(tp, tp$tissue), function(x) { bs <- mp[gtex == x$tissue[1]]$biosample_value
    if (!length(bs)) return(NULL)
    has <- !is.na(x$acc) & x$acc %in% rownames(A); a <- x$acc[has]
    no <- integer(nrow(x)); pr <- en <- rep(NA_real_, nrow(x))
    if (length(a)) { o <- rowSums(OPEN[a, bs, drop = FALSE]); no[has] <- o
      pr[has] <- rowSums(PROM[a, bs, drop = FALSE]) / pmax(o, 1); en[has] <- rowSums(ENH[a, bs, drop = FALSE]) / pmax(o, 1) }
    x[, `:=`(n_bios = length(bs), n_open = no, prom = pr, enh = en)] }))
  out[, state := fcase(is.na(prom) | n_open == 0, "no accessible element", prom >= enh & prom > 0, "promoter-like",
                       enh > prom, "enhancer-like", default = "accessible, unmarked")][] }
GL <- c("tested, no QTL", "mu (size)", "p0 (freq)", "p1 (freq)")
grp <- function(x) rbind(x[!mu & !p0 & !p1][, group := GL[1]], x[mu == TRUE][, group := GL[2]], x[p0 == TRUE][, group := GL[3]], x[p1 == TRUE][, group := GL[4]])
gboot <- function(d, f) { gi <- split(seq_len(nrow(d)), d$gene0); g <- names(gi); set.seed(SEED)
  replicate(B, f(d[unlist(gi[sample(g, length(g), TRUE)], use.names = FALSE)])) }
summ <- function(G, tag) rbindlist(lapply(GL, function(g) { d <- G[group == g]
  rbindlist(lapply(c("enhancer-like", "promoter-like", "no accessible element"), function(st) { f <- function(x) 100 * mean(x$state == st); b <- gboot(d, f)
    data.table(set = tag, group = g, state = st, units = nrow(d), genes = uniqueN(d$gene0), tissues = uniqueN(d$tissue),
               pct = f(d), lo = quantile(b, .025), hi = quantile(b, .975)) })) }))
G1 <- grp(state_for(TP, c("exact", "approximate"))); G2 <- grp(state_for(TP, "exact"))
RES <- rbind(summ(G1, "matched: exact + approximate"), summ(G2, "matched: exact only"))
for (j in c("pct", "lo", "hi")) set(RES, j = j, value = round(RES[[j]], 2))
fwrite(RES, file.path(OUTD, "28v4_tissue_matched_tss_chromatin.tsv"), sep = "\t")
fwrite(MAPb[, .N, by = .(gtex, match)][order(gtex)], file.path(OUTD, "28v4_tissue_map_used.tsv"), sep = "\t")
options(width = 220)
cat(sprintf("mapped GTEx tissues: exact %d, approximate %d; ENCODE biosamples used: %d\n", MAP[match == "exact", uniqueN(gtex)],
            MAP[match == "approximate", uniqueN(gtex)], uniqueN(MAPb$biosample_value)))
print(dcast(RES[, v := sprintf("%.1f [%.1f, %.1f]", pct, lo, hi)], set + group + units ~ state, value.var = "v"))

## figure
source(file.path(F5S, "code/theme_nature.R"))
BASE <- 7; ann5 <- pt2mm(5)
theme_set(theme_nature(BASE) + theme(plot.title = element_text(size = BASE, face = "bold"), plot.subtitle = element_text(size = 5.5, colour = "grey35"),
                                     strip.text = element_text(size = BASE - 1)))
CCOL <- c("tested, no QTL" = "grey55", "mu (size)" = "#4E79A7", "p0 (freq)" = "#F28E2B", "p1 (freq)" = "#59A14F")
pd <- RES[set == "matched: exact + approximate"][, `:=`(group = factor(group, GL), state = factor(state, c("promoter-like", "enhancer-like", "no accessible element")))]
p <- ggplot(pd, aes(group, pct, fill = group)) +
  geom_col(width = 0.68, colour = "grey25", linewidth = 0.2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, linewidth = 0.3, colour = "grey25") +
  geom_text(aes(y = hi, label = sprintf("%.1f%%", pct)), vjust = -0.5, size = ann5, colour = "grey20") +
  facet_wrap(~state, scales = "free_y", nrow = 1) + scale_fill_manual(values = CCOL, guide = "none") +
  scale_x_discrete(labels = c("tested,\nno QTL", "mu\n(size)", "p0\n(freq)", "p1\n(freq)")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), limits = c(0, NA)) +
  labs(x = NULL, y = "QTL units (%)", title = "TSS chromatin in the QTL's own tissue (ENCODE samples matched to the GTEx tissue)",
       subtitle = sprintf("units = (tissue, promoter); %d GTEx tissues matched (exact + approximate); gene-clustered 95%% CI",
                          uniqueN(G1$tissue)))
ggsave(file.path(OUTD, "28v4_tissue_matched_tss_chromatin.png"), p, width = 180, height = 68, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(OUTD, "28v4_tissue_matched_tss_chromatin.pdf"), p, width = 180, height = 68, units = "mm", device = cairo_pdf)
cat("written figures/bcd_exploration/tissue_matched_tss_chromatin_v3.{png,pdf}\n")
