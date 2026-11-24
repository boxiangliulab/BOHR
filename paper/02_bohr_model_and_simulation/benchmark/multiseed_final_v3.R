#!/usr/bin/env Rscript
### Supplementary per-seed and pooled permuted-null QQ with the FINAL Fig 2 statistics. BOHR = Simes (or Fisher, see COMBINE) over the
### SNP's existing tests (mu F-corrected LRT, p0 / p1 Firth), the rule of Fig 2a-c and Fig 2e; comparator = the version's
### QTLtools mode (vsRaw: raw; vsRINTfirst: RINT-first). Five independent sets of 1,000 pre-drawn permutations per unit:
###   seed 1 = base 20261001 (the figure), 2 = 20262025, 3 = 20282025, 4 = 20292025, 5 = 20302025.
### Same calibration units as the figure (n_mid >= 5; raw RNA-seq exact mirror pairs counted once); DATASET env MS_DATASET.
### Framework display names: "Method 1" = direct activity simulation, "Method 2" = raw RNA-seq
### simulation; the data keep the keys M1 / M2.
### The seed figures are NOT in the supplement; they are kept as a reserve. The tables feed Supplementary Table 2b.
### Outputs per comparator version (supplementary_v3_<DS>/seed_reserve_S9/<suf>/):
###   S9a_perseed_calibration_direct_activity.*  per-seed QQ, direct activity simulation (rows: BOHR, comparator; columns: seeds)
###   S9b_perseed_calibration_raw_rnaseq.*       per-seed QQ, raw RNA-seq simulation
###   S9c_pooled_seed_calibration.*              five seeds pooled: BOHR Simes, its three component tests, the comparator
### Tables (data_v3_<DS>/seed_replicates/final_stats/): lambda_tail_by_seed.csv, qq_points_per_seed.csv, qq_points_pooled.csv.
### For the paper these replace supplementary_v3_<DS>/multiseed/QQ_multiseed_* (mu test and RINT-first), which stay as
### diagnostics; multiseed_qq_v3.R is not changed. PLOT_ONLY=1 re-plots from the saved tables.
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
ROOT <- "/projects/group/bliu/e1327989/GTEx/new_simulation/figures_bohr2"; pd <- readLines(file.path(ROOT, "scripts", "prepare_data_v3.R"))
DS <- Sys.getenv("MS_DATASET", "nothr_realjunc"); PLOT_ONLY <- Sys.getenv("PLOT_ONLY", "0") == "1"
## COMBINE (the final runs: fisher_sidak): the promoter-level BOHR P is taken from prepare_data_v3.R with the same switch; tables and
## reserve plots go to data_v3_<DS><suffix>/ and supplementary_v3_<DS><suffix>/ (suffix "" for simes, the default, which reproduces the earlier outputs).
CMB <- Sys.getenv("COMBINE", "simes"); stopifnot(CMB %in% c("simes", "fisher", "fisher_sidak")); MSFX <- c(simes = "", fisher = "_fisher", fisher_sidak = "_fishersidak")[[CMB]]
TAGS <- c("Seed 1 (figure)" = "", "Seed 2" = "_seed2", "Seed 3" = "_seed3", "Seed 4" = "_seed4", "Seed 5" = "_seed5")
TESTS <- c("BOHR_simes", "BOHR", "BOHR_p0", "BOHR_p1", "raw", "RINTfirst")
TLAB <- c(BOHR_simes = c(simes = "BOHR (Simes)", fisher = "BOHR (Fisher)", fisher_sidak = "BOHR (Fisher, Šidák parts)")[[CMB]], BOHR = "BOHR μ test", BOHR_p0 = "BOHR p0 test", BOHR_p1 = "BOHR p1 test",
          raw = "QTLtools raw", RINTfirst = "QTLtools RINT-first")
FW <- c(M1 = "Direct activity simulation", M2 = "Raw RNA-seq simulation"); FWS <- c(M1 = "Direct activity", M2 = "Raw RNA-seq")
TD <- file.path(ROOT, sprintf("data_v3_%s%s", DS, MSFX), "seed_replicates", "final_stats"); dir.create(TD, recursive = TRUE, showWarnings = FALSE)
lam <- function(p) median(qchisq(p, 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(.5, 1)
qq_of <- function(p) { p <- sort(p[!is.na(p)]); n <- length(p)                    # tail kept, body thinned
  d <- data.table(rank = seq_len(n), n = n, observed = -log10(pmax(p, 1e-300)), expected = -log10(seq_len(n) / (n + 1)))
  d[, `:=`(band_lo = -log10(qbeta(.975, rank, n - rank + 1)), band_hi = -log10(qbeta(.025, rank, n - rank + 1)))]
  d[(rank <= 2000) | (rank %% pmax(1L, as.integer(n / 3000)) == 0L)] }

if (!PLOT_ONLY) {
  build_nul <- function(tag) { Sys.setenv(SEEDTAG = tag, DATASET = DS); e <- new.env()
    eval(parse(text = pd[1:(grep("^MODES <- ", pd) - 1)]), envir = e)             # NUL incl. BOHR_simes and the QTLtools modes
    d <- e$OUT                                                                    # prepare_data_v3.R's dir.create (data_v3_<DS><tag><suffix>); seed data live in seed_replicates/
    if (nzchar(tag) && dir.exists(d) && length(list.files(d, all.files = TRUE, no.. = TRUE)) == 0) unlink(d, recursive = TRUE)
    e$NUL[, c("m", "N", "id", "perm", TESTS), with = FALSE] }
  NN <- rbindlist(lapply(names(TAGS), function(s) build_nul(TAGS[[s]])[, seed := s])); Sys.setenv(SEEDTAG = "")
  cat("units and tests per framework x seed:\n")
  print(NN[, .(units = uniqueN(paste(N, id)), tests = .N, simes_na = sum(is.na(BOHR_simes)), raw_na = sum(is.na(raw)), rint_na = sum(is.na(RINTfirst))), by = .(m, seed)])
  TB <- rbindlist(lapply(TESTS, function(t) {
    f <- function(by) NN[!is.na(get(t)), .(test = t, tests = .N, lambda = lam(get(t)), r1e3 = mean(get(t) < 1e-3) / 1e-3,
                                           r1e4 = mean(get(t) < 1e-4) / 1e-4, r1e5 = mean(get(t) < 1e-5) / 1e-5), by = by]
    rbind(f(c("m", "seed")), f("m")[, seed := "Pooled"]) }))
  fwrite(TB, file.path(TD, "lambda_tail_by_seed.csv"))
  QS <- rbindlist(lapply(c("M1", "M2"), function(mm) rbindlist(lapply(names(TAGS), function(s) { X <- NN[m == mm & seed == s]
    rbindlist(lapply(c("BOHR_simes", "raw", "RINTfirst"), function(t) qq_of(X[[t]])[, `:=`(m = mm, seed = s, test = t)])) }))))
  QP <- rbindlist(lapply(c("M1", "M2"), function(mm) { X <- NN[m == mm]; rbindlist(lapply(TESTS, function(t) qq_of(X[[t]])[, `:=`(m = mm, test = t)])) }))
  fwrite(QS, file.path(TD, "qq_points_per_seed.csv")); fwrite(QP, file.path(TD, "qq_points_pooled.csv"))
} else {
  TB <- fread(file.path(TD, "lambda_tail_by_seed.csv")); QS <- fread(file.path(TD, "qq_points_per_seed.csv")); QP <- fread(file.path(TD, "qq_points_pooled.csv"))
}
op <- options(width = 200); cat("\nlambda and observed / expected, per seed and pooled:\n")
print(TB[, .(m, seed, test, tests, lambda = round(lambda, 3), r1e3 = round(r1e3, 2), r1e4 = round(r1e4, 2), r1e5 = round(r1e5, 2))][order(m, test, seed)]); options(op)

FAM <- "Arial"; CAP <- 8
## Arial plain text only -- no plotmath; every label is one Unicode string (italic / subscript set by hand in
## Illustrator); SVG without textLength (fix_text_size = FALSE); PNG/TIFF via cairo, hinting off (FONTCONFIG_FILE=scripts/fonts_nohint.conf).
update_geom_defaults("text", list(family = FAM)); update_geom_defaults("label", list(family = FAM))
theme_set(theme_classic(base_size = 7, base_family = FAM) +
  theme(axis.line = element_line(linewidth = .3), axis.ticks = element_line(linewidth = .3), axis.text = element_text(size = 6, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 6), legend.position = "bottom",
        strip.text = element_text(size = 6.5, face = "bold"), strip.background = element_rect(fill = "grey92", colour = NA),
        plot.title = element_text(size = 8, face = "bold"), panel.grid = element_blank()))
save3 <- function(p, f, w_mm, h_mm) { w <- w_mm / 25.4; h <- h_mm / 25.4
  svglite::svglite(paste0(f, ".svg"), width = w, height = h, fix_text_size = FALSE); print(p); dev.off()
  grDevices::cairo_pdf(paste0(f, ".pdf"), width = w, height = h, family = FAM); print(p); dev.off()
  grDevices::png(paste0(f, ".png"), width = w, height = h, units = "in", res = 300, type = "cairo", family = FAM); print(p); dev.off() }
COLT <- c(BOHR_simes = "#C0392B", BOHR = "#C0392B", BOHR_p0 = "#C0392B", BOHR_p1 = "#C0392B", raw = "#4C78A8", RINTfirst = "#1F4E79")
XL <- "−log10(expected P)"; YL <- "−log10(observed P)"
qqbase <- function(D) ggplot(D, aes(expected, observed)) + geom_ribbon(aes(ymin = band_lo, ymax = band_hi), fill = "grey85", alpha = .6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45", linewidth = .3) +
  geom_point(aes(colour = test), size = .3, alpha = .55, show.legend = FALSE) +                     # all points drawn alike
  geom_point(data = D[capped == TRUE], aes(colour = test), shape = 17, size = .9, show.legend = FALSE) +
  scale_colour_manual(values = COLT) + coord_cartesian(ylim = c(0, CAP)) + labs(x = XL, y = YL)
for (v in c("vsRaw", "vsRINTfirst")) { cp <- c(vsRaw = "raw", vsRINTfirst = "RINTfirst")[[v]]
  OUT <- file.path(ROOT, sprintf("supplementary_v3_%s%s", DS, MSFX), paste0("seed_reserve_S9", Sys.getenv("OUT_TAG", "")), v); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)   # reserve, not the supplement
  ## per seed, one file per framework
  for (mm in c("M1", "M2")) { K <- c("BOHR_simes", cp)
    D <- QS[m == mm & test %in% K][, `:=`(test = factor(test, levels = K), seed = factor(seed, levels = names(TAGS)))][, capped := observed > CAP][, observed := pmin(observed, CAP)]
    L <- TB[m == mm & seed != "Pooled" & test %in% K][, `:=`(test = factor(test, levels = K), seed = factor(seed, levels = names(TAGS)))][, txt := sprintf("λ = %.2f", lambda)]
    p <- qqbase(D) + geom_text(data = L, aes(-Inf, Inf, label = txt), hjust = -.12, vjust = 1.6, size = 1.9, colour = "grey20") +
      facet_grid(test ~ seed, labeller = labeller(test = TLAB)) + labs(title = sprintf("%s: permuted-null QQ, each of five permutation seeds", FW[[mm]]))
    save3(p, file.path(OUT, sprintf("S9%s_perseed_calibration_%s", c(M1 = "a", M2 = "b")[[mm]], c(M1 = "direct_activity", M2 = "raw_rnaseq")[[mm]])), 183, 88) }
  ## five seeds pooled: BOHR Simes, its component tests, the comparator
  K <- c("BOHR_simes", "BOHR", "BOHR_p0", "BOHR_p1", cp)
  D <- QP[test %in% K][, `:=`(test = factor(test, levels = K), fw = factor(FWS[m], levels = FWS))][, capped := observed > CAP][, observed := pmin(observed, CAP)]
  L <- TB[seed == "Pooled" & test %in% K][, `:=`(test = factor(test, levels = K), fw = factor(FWS[m], levels = FWS))][
    , `:=`(txt = sprintf("λ = %.2f", lambda), txt2 = sprintf("obs/exp at 1e-3, 1e-4, 1e-5:\n%.2f, %.2f, %.2f", r1e3, r1e4, r1e5))]
  p <- qqbase(D) + geom_text(data = L, aes(-Inf, Inf, label = txt), hjust = -.12, vjust = 1.6, size = 1.9, colour = "grey20") +
    geom_text(data = L, aes(-Inf, Inf, label = txt2), hjust = -.07, vjust = 2.2, size = 1.6, lineheight = .9, family = FAM, colour = "grey30") +
    facet_grid(fw ~ test, labeller = labeller(test = TLAB)) + labs(title = "Permuted-null QQ, five permutation seeds pooled")
  save3(p, file.path(OUT, "S9c_pooled_seed_calibration"), 183, 92)
  cat("saved S9a-c ->", OUT, "\n") }
cat("MULTISEED_FINAL_DONE\n")
