library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(scales)

setwd("/data/projects/11003054/YuhaoDong/GTEx/version_comparison")

# ============================================================
# 1. Read CSV data
# ============================================================
csv_data <- read.csv("GTEx_promoter.csv", check.names = FALSE)
colnames(csv_data)[1] <- "Filter"

# Build long-format data
gencode_row <- csv_data[1, -1]
count_data <- csv_data[-1, ]  # remove gencode row

# Version ordering (critical: define once, use everywhere)
ver_levels <- c("v7", "v8", "v10", "v11")
ver_gencode <- c("v19", "v26", "v39", "v47")
ver_labels <- paste0(ver_levels, " (GENCODE ", ver_gencode, ")")

# Long format for filter stages
df_long <- count_data %>%
  pivot_longer(cols = -Filter, names_to = "version", values_to = "count") %>%
  mutate(
    count = as.numeric(count),
    version = factor(version, levels = ver_levels),
    Filter = factor(Filter, levels = count_data$Filter)
  )

# ============================================================
# 2. Load RData for deeper analysis
# ============================================================
all_info <- list()
for (v in ver_levels) {
  load(paste0("Promoter_info_all_", v, ".RData"))
  info_all$version <- v
  all_info[[v]] <- info_all
}

# Summary statistics
summary_stats <- do.call(rbind, lapply(ver_levels, function(v) {
  d <- all_info[[v]]
  n_tx <- sapply(d$txId, length)
  data.frame(
    version = v,
    total_promoters = nrow(d),
    unique_genes = length(unique(d$geneId)),
    external_promoters = sum(!d$internalPromoter),
    internal_promoters = sum(d$internalPromoter),
    mean_tx_per_promoter = mean(n_tx),
    mean_promoters_per_gene = nrow(d) / length(unique(d$geneId))
  )
}))
summary_stats$version <- factor(summary_stats$version, levels = ver_levels)

# ============================================================
# Color palette - v11 highlighted in red, others in gray gradient
# ============================================================
version_colors <- c("v7" = "#B0BEC5", "v8" = "#90A4AE",
                     "v10" = "#78909C", "v11" = "#E53935")

# ============================================================
# Plot 1: Filter pipeline (grouped bar chart)
# ============================================================
p1 <- ggplot(df_long, aes(x = Filter, y = count, fill = version)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = format(count, big.mark = ",")),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 2.8, fontface = "bold") +
  scale_fill_manual(values = version_colors, labels = ver_labels) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Promoter Counts Across GTEx Versions at Each Filtering Stage",
       subtitle = "GTEx v11 (GENCODE v47) consistently yields the highest promoter counts",
       x = NULL, y = "Number of Promoters", fill = "GTEx Version") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "#E53935", size = 11),
    axis.text.x = element_text(angle = 15, hjust = 1, size = 10),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

ggsave("Fig1_filter_pipeline.pdf", p1, width = 14, height = 7)
ggsave("Fig1_filter_pipeline.png", p1, width = 14, height = 7, dpi = 300)

# ============================================================
# Plot 2: Retention rate through filters
# ============================================================
retention <- df_long %>%
  group_by(version) %>%
  mutate(retention_pct = count / count[Filter == "all"] * 100) %>%
  ungroup()

p2 <- ggplot(retention, aes(x = Filter, y = retention_pct,
                             color = version, group = version)) +
  geom_line(aes(linewidth = version), show.legend = TRUE) +
  geom_point(aes(size = version), show.legend = TRUE) +
  scale_color_manual(values = version_colors, labels = ver_labels) +
  scale_linewidth_manual(values = c("v7" = 0.8, "v8" = 0.8, "v10" = 0.8, "v11" = 1.5), guide = "none") +
  scale_size_manual(values = c("v7" = 2, "v8" = 2, "v10" = 2, "v11" = 3.5), guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 105)) +
  geom_text(data = retention %>% filter(Filter == "Filter Single"),
            aes(label = paste0(round(retention_pct, 1), "%")),
            vjust = -1.2, size = 3.5, fontface = "bold", show.legend = FALSE) +
  labs(title = "Promoter Retention Rate Through Filtering Pipeline",
       subtitle = "v10 and v11 show markedly higher retention; v11 yields the most promoters in absolute count",
       x = NULL, y = "Retention Rate (% of All Promoters)",
       color = "GTEx Version") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "#E53935", size = 11),
    axis.text.x = element_text(angle = 15, hjust = 1, size = 10),
    legend.position = "top"
  )

ggsave("Fig2_retention_rate.pdf", p2, width = 12, height = 7)
ggsave("Fig2_retention_rate.png", p2, width = 12, height = 7, dpi = 300)

# ============================================================
# Plot 3: Gene coverage & promoters per gene (triple panel)
# ============================================================
p3a <- ggplot(summary_stats, aes(x = version, y = unique_genes, fill = version)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = format(unique_genes, big.mark = ",")),
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = version_colors) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Gene Coverage",
       subtitle = "Unique genes with identified promoters",
       x = NULL, y = "Number of Genes") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

p3b <- ggplot(summary_stats, aes(x = version, y = mean_promoters_per_gene, fill = version)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = round(mean_promoters_per_gene, 2)),
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = version_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Promoter Density",
       subtitle = "Mean promoters per gene",
       x = NULL, y = "Promoters / Gene") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

p3c <- ggplot(summary_stats, aes(x = version, y = mean_tx_per_promoter, fill = version)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = round(mean_tx_per_promoter, 2)),
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = version_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Transcript Support",
       subtitle = "Mean transcripts per promoter",
       x = NULL, y = "Transcripts / Promoter") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

p3 <- p3a + p3b + p3c +
  plot_annotation(
    title = "Promoter Annotation Quality Metrics Across GTEx Versions",
    subtitle = "v11 achieves the best gene coverage and transcript support",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "#E53935", size = 12)
    )
  )

ggsave("Fig3_quality_metrics.pdf", p3, width = 14, height = 6)
ggsave("Fig3_quality_metrics.png", p3, width = 14, height = 6, dpi = 300)

# ============================================================
# Plot 4: Internal vs External promoter composition (stacked)
# ============================================================
comp_data <- summary_stats %>%
  select(version, external_promoters, internal_promoters) %>%
  pivot_longer(cols = c(external_promoters, internal_promoters),
               names_to = "type", values_to = "count") %>%
  mutate(
    type = ifelse(type == "external_promoters", "First Exon (External)", "Internal"),
    version = factor(version, levels = ver_levels)
  )

fill_colors <- c("First Exon (External)" = "#42A5F5", "Internal" = "#FFA726")

p4 <- ggplot(comp_data, aes(x = version, y = count, fill = type)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = format(count, big.mark = ",")),
            position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold") +
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Promoter Composition: External vs Internal",
       subtitle = "v11 identifies substantially more external (first exon) promoters",
       x = "GTEx Version", y = "Number of Promoters", fill = "Promoter Type") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "#E53935", size = 11),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

ggsave("Fig4_composition.pdf", p4, width = 10, height = 7)
ggsave("Fig4_composition.png", p4, width = 10, height = 7, dpi = 300)

# ============================================================
# Plot 5: Promoter position distribution per version
# ============================================================
pos_data <- do.call(rbind, lapply(ver_levels, function(v) {
  d <- all_info[[v]]
  pos_tab <- as.data.frame(table(pmin(d$promoterPosition, 6)))
  colnames(pos_tab) <- c("position", "count")
  pos_tab$version <- v
  pos_tab$position <- factor(pos_tab$position, levels = 1:6,
                              labels = c("1st", "2nd", "3rd", "4th", "5th", "6th+"))
  pos_tab
}))
pos_data$version <- factor(pos_data$version, levels = ver_levels)

p5 <- ggplot(pos_data, aes(x = position, y = count, fill = version)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = version_colors, labels = ver_labels) +
  scale_y_continuous(labels = comma) +
  labs(title = "Promoter Position Distribution Within Genes",
       subtitle = "v11 captures more alternative promoters at every position rank",
       x = "Promoter Position (Ranked by Gene Structure)",
       y = "Number of Promoters", fill = "GTEx Version") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "#E53935", size = 11),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

ggsave("Fig5_position_distribution.pdf", p5, width = 12, height = 7)
ggsave("Fig5_position_distribution.png", p5, width = 12, height = 7, dpi = 300)

# ============================================================
# Plot 6: Fold-change improvement relative to v7
# ============================================================
fold_data <- df_long %>%
  group_by(Filter) %>%
  mutate(fold_vs_v7 = count / count[version == "v7"]) %>%
  ungroup() %>%
  filter(version != "v7") %>%
  mutate(version = factor(version, levels = c("v8", "v10", "v11")))

fold_colors <- c("v8" = "#90A4AE", "v10" = "#78909C", "v11" = "#E53935")
fold_labels <- c("v8" = "v8 vs v7", "v10" = "v10 vs v7", "v11" = "v11 vs v7")

p6 <- ggplot(fold_data, aes(x = Filter, y = fold_vs_v7, fill = version)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_text(aes(label = paste0(round(fold_vs_v7, 2), "x")),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 3, fontface = "bold") +
  scale_fill_manual(values = fold_colors, labels = fold_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Fold Improvement Over GTEx v7 Baseline",
       subtitle = "v11 shows the largest improvement across all filtering stages",
       x = NULL, y = "Fold Change (relative to v7)", fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "#E53935", size = 11),
    axis.text.x = element_text(angle = 15, hjust = 1, size = 10),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

ggsave("Fig6_fold_improvement.pdf", p6, width = 12, height = 7)
ggsave("Fig6_fold_improvement.png", p6, width = 12, height = 7, dpi = 300)

# ============================================================
# Plot 8: Single- vs Multi-promoter gene proportion
# ============================================================
gene_class <- do.call(rbind, lapply(ver_levels, function(v) {
  d <- all_info[[v]]
  gene_counts <- table(d$geneId)
  n_single <- sum(gene_counts == 1)
  n_multi  <- sum(gene_counts > 1)
  data.frame(
    version = v,
    type = c("Single-promoter\n(not usable for AP analysis)",
             "Multi-promoter\n(usable for AP analysis)"),
    count = c(n_single, n_multi)
  )
}))
gene_class$version <- factor(gene_class$version, levels = ver_levels)
gene_class$type <- factor(gene_class$type,
  levels = c("Single-promoter\n(not usable for AP analysis)",
             "Multi-promoter\n(usable for AP analysis)"))

# Compute percentages
gene_class <- gene_class %>%
  group_by(version) %>%
  mutate(total = sum(count), pct = count / total * 100) %>%
  ungroup()

# --- Panel A: stacked bar (absolute counts) with percentage labels ---
ap_fill <- c("Single-promoter\n(not usable for AP analysis)" = "#BDBDBD",
             "Multi-promoter\n(usable for AP analysis)" = "#E53935")

p8a <- ggplot(gene_class, aes(x = version, y = count, fill = type)) +
  geom_bar(stat = "identity", width = 0.55) +
  geom_text(aes(label = paste0(format(count, big.mark = ","),
                               "\n(", round(pct, 1), "%)")),
            position = position_stack(vjust = 0.5),
            size = 3.5, fontface = "bold", color = "white") +
  scale_fill_manual(values = ap_fill) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.03))) +
  labs(title = "Gene Classification for Alternative Promoter (AP) Analysis",
       subtitle = "Single-promoter genes must be excluded; v11 provides the most multi-promoter genes",
       x = "GTEx Version", y = "Number of Genes", fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "#E53935", size = 11),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

# --- Panel B: multi-promoter gene count only (highlight v11 advantage) ---
multi_only <- gene_class %>% filter(grepl("Multi", type))

p8b <- ggplot(multi_only, aes(x = version, y = count, fill = version)) +
  geom_bar(stat = "identity", width = 0.55) +
  geom_text(aes(label = format(count, big.mark = ",")),
            vjust = -0.5, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = version_colors) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Genes Usable for AP Analysis",
       subtitle = "Genes with >= 2 promoters",
       x = "GTEx Version", y = "Number of Genes") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

p8 <- p8a + p8b + plot_layout(widths = c(1.6, 1)) +
  plot_annotation(
    title = "Single- vs Multi-Promoter Gene Distribution Across GTEx Versions",
    subtitle = "Although v11 includes more single-promoter genes from expanded annotation, it yields 24,710 multi-promoter genes -- the most across all versions",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "#E53935", size = 11)
    )
  )

ggsave("Fig8_single_vs_multi_promoter.pdf", p8, width = 14, height = 7)
ggsave("Fig8_single_vs_multi_promoter.png", p8, width = 14, height = 7, dpi = 300)

# ============================================================
# Plot 9: Source of promoter changes (new genes vs more promoters/gene)
# ============================================================
strip_ver <- function(x) sub("\\..*", "", x)

# Add base gene ID to each version
for (v in ver_levels) {
  all_info[[v]]$geneId_base <- strip_ver(all_info[[v]]$geneId)
}

# ---- Panel A: consecutive version decomposition (waterfall-style) ----
transitions <- list(c("v7","v8"), c("v8","v10"), c("v10","v11"))
decomp_data <- do.call(rbind, lapply(transitions, function(pair) {
  v_old <- pair[1]; v_new <- pair[2]
  genes_old <- unique(all_info[[v_old]]$geneId_base)
  genes_new <- unique(all_info[[v_new]]$geneId_base)
  shared <- intersect(genes_old, genes_new)
  only_new <- setdiff(genes_new, genes_old)
  only_old <- setdiff(genes_old, genes_new)

  pc_old <- table(all_info[[v_old]]$geneId_base)
  pc_new <- table(all_info[[v_new]]$geneId_base)
  diff_shared <- as.numeric(pc_new[shared]) - as.numeric(pc_old[shared])

  prom_new_genes   <- sum(all_info[[v_new]]$geneId_base %in% only_new)
  prom_lost_genes  <- sum(all_info[[v_old]]$geneId_base %in% only_old)
  gained_in_shared <- sum(diff_shared[diff_shared > 0])
  lost_in_shared   <- sum(abs(diff_shared[diff_shared < 0]))

  label <- paste0(v_old, " \u2192 ", v_new)
  data.frame(
    transition = label,
    source = c("New genes\n(added)", "Old genes\n(removed)",
               "Shared genes\n(gained promoters)", "Shared genes\n(lost promoters)"),
    promoters = c(prom_new_genes, -prom_lost_genes, gained_in_shared, -lost_in_shared),
    n_genes = c(length(only_new), length(only_old),
                sum(diff_shared > 0), sum(diff_shared < 0))
  )
}))

decomp_data$transition <- factor(decomp_data$transition,
  levels = c("v7 \u2192 v8", "v8 \u2192 v10", "v10 \u2192 v11"))
decomp_data$source <- factor(decomp_data$source,
  levels = c("New genes\n(added)", "Shared genes\n(gained promoters)",
             "Shared genes\n(lost promoters)", "Old genes\n(removed)"))

source_colors <- c(
  "New genes\n(added)" = "#E53935",
  "Shared genes\n(gained promoters)" = "#FF8A65",
  "Shared genes\n(lost promoters)" = "#90A4AE",
  "Old genes\n(removed)" = "#546E7A"
)

p9a <- ggplot(decomp_data, aes(x = source, y = promoters, fill = source)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste0(ifelse(promoters >= 0, "+", ""),
                               format(promoters, big.mark = ","),
                               "\n(", format(n_genes, big.mark = ","), " genes)")),
            vjust = ifelse(decomp_data$promoters >= 0, -0.2, 1.2),
            size = 3, fontface = "bold") +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  facet_wrap(~ transition, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = source_colors) +
  scale_y_continuous(labels = comma) +
  labs(title = "Sources of Promoter Changes Between Consecutive Versions",
       subtitle = "Where do new promoters come from: newly annotated genes, or existing genes gaining more promoters?",
       x = NULL, y = "Promoter Count Change") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10.5),
    legend.position = "none",
    axis.text.x = element_text(size = 8),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

# ---- Panel B: percentage stacked bar (positive contributions only) ----
pos_data <- decomp_data %>%
  filter(promoters > 0) %>%
  mutate(source_simple = ifelse(grepl("New genes", source),
                                "New genes", "Existing genes\n(more promoters)")) %>%
  group_by(transition, source_simple) %>%
  summarise(promoters = sum(promoters), .groups = "drop") %>%
  group_by(transition) %>%
  mutate(pct = promoters / sum(promoters) * 100) %>%
  ungroup()

pos_data$source_simple <- factor(pos_data$source_simple,
  levels = c("New genes", "Existing genes\n(more promoters)"))

pct_colors <- c("New genes" = "#E53935",
                "Existing genes\n(more promoters)" = "#FF8A65")

p9b <- ggplot(pos_data, aes(x = transition, y = pct, fill = source_simple)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(label = paste0(round(pct, 1), "%")),
            position = position_stack(vjust = 0.5),
            size = 4.5, fontface = "bold", color = "white") +
  scale_fill_manual(values = pct_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Relative Contribution to Promoter Gains",
       subtitle = "Proportion from new genes vs existing genes gaining promoters",
       x = NULL, y = "% of Gained Promoters", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

# ---- Panel C: cumulative vs v7 baseline ----
cumul_data <- do.call(rbind, lapply(c("v8", "v10", "v11"), function(v) {
  genes_v7 <- unique(all_info[["v7"]]$geneId_base)
  genes_v  <- unique(all_info[[v]]$geneId_base)
  shared <- intersect(genes_v7, genes_v)
  only_new <- setdiff(genes_v, genes_v7)

  pc_v7 <- table(all_info[["v7"]]$geneId_base)
  pc_v  <- table(all_info[[v]]$geneId_base)
  diff_s <- as.numeric(pc_v[shared]) - as.numeric(pc_v7[shared])

  prom_new <- sum(all_info[[v]]$geneId_base %in% only_new)
  gained_s <- sum(diff_s[diff_s > 0])

  data.frame(
    version = v,
    source = c("New genes", "Existing genes\n(more promoters)"),
    promoters = c(prom_new, gained_s)
  )
}))
cumul_data$version <- factor(cumul_data$version, levels = c("v8", "v10", "v11"))
cumul_data$source <- factor(cumul_data$source,
  levels = c("New genes", "Existing genes\n(more promoters)"))

p9c <- ggplot(cumul_data, aes(x = version, y = promoters, fill = source)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  geom_text(aes(label = format(promoters, big.mark = ",")),
            position = position_dodge(width = 0.6),
            vjust = -0.4, size = 3.8, fontface = "bold") +
  scale_fill_manual(values = pct_colors) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Cumulative Promoter Gains vs v7",
       subtitle = "Total new promoters from each source compared to v7 baseline",
       x = "GTEx Version", y = "Additional Promoters vs v7", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

p9 <- (p9a / (p9b | p9c)) +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(
    title = "Decomposition of Promoter Increases Across GTEx Versions",
    subtitle = "Both new gene annotations AND deeper promoter discovery in existing genes contribute to the growth",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(color = "#E53935", size = 12)
    )
  )

ggsave("Fig9_promoter_source_decomposition.pdf", p9, width = 16, height = 12)
ggsave("Fig9_promoter_source_decomposition.png", p9, width = 16, height = 12, dpi = 300)

# ============================================================
# Plot 7: Combined summary figure
# ============================================================
p_combined <- (p1 / p2) | (p3a / p3b / p3c)
p_combined <- p_combined +
  plot_annotation(
    title = "GTEx Version Comparison: Promoter Identification Performance",
    subtitle = "GTEx v11 with GENCODE v47 delivers superior promoter identification across all metrics",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(color = "#E53935", size = 12, hjust = 0.5)
    )
  ) +
  plot_layout(widths = c(2, 1))

ggsave("Fig7_combined_summary.pdf", p_combined, width = 20, height = 14)
ggsave("Fig7_combined_summary.png", p_combined, width = 20, height = 14, dpi = 300)

cat("\n=== All figures saved successfully ===\n")
