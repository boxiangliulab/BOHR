### theme_nature.R — portable Nature-standard ggplot2 theme + export helper.
### Drop into any project: source("theme_nature.R").
### Requires: ggplot2; for save_pub() also svglite, ragg (cairo_pdf is base grDevices).
suppressMessages({library(ggplot2)})

## Project palette — red = primary method, blue = baseline. Override per project.
PALETTE <- c(primary = "#C0392B", baseline = "#2E5A88")
NATURE_FONT <- "Arial"   # set explicitly; "sans" resolves to DejaVu, not Arial

## Nature-standard theme. base_size in pt (5-7 typical; 7 = this project default).
theme_nature <- function(base_size = 7, base_family = NATURE_FONT) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line   = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks  = element_line(linewidth = 0.3, colour = "black"),
      axis.title  = element_text(size = base_size),
      axis.text   = element_text(size = base_size - 1, colour = "black"),
      legend.title = element_blank(),
      legend.text  = element_text(size = base_size - 1),
      legend.key.size = unit(3.2, "mm"),
      legend.position = "bottom",
      strip.text  = element_text(size = base_size - 0.5, face = "bold"),
      strip.background = element_rect(fill = "grey92", colour = NA),
      plot.title  = element_text(size = base_size + 1, face = "bold"),
      plot.tag    = element_text(size = 9, face = "bold", family = base_family),
      panel.grid  = element_blank())
}

## geom_text() size is in mm. Use pt2mm(5) to place a 5 pt annotation.
pt2mm <- function(pt) pt / 2.845

## Submission export bundle: vector SVG + Arial-embedded PDF + 600-dpi TIFF + PNG preview.
save_pub <- function(plot, filename, width_mm = 183, height_mm = 120,
                     dpi = 600, family = NATURE_FONT) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = family); print(plot); dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi, compression = "lzw"); print(plot); dev.off()
  ggsave(paste0(filename, ".png"), plot, width = w, height = h, dpi = 200)
}
