#!/usr/bin/env Rscript
### Source Data for Fig. 2B: detection power from the existing final runs, binned
### by simulated effect size at a fixed N = 300.
### Detection power by simulated effect size at ONE sample size, from the same final runs, tests and permutation-calibrated
### thresholds as Fig. 2A: BOHR = promoter-level Simes(μ F-LRT, p0 Firth, p1 Firth); QTLtools modes = their single test; a call
### is p <= the alpha = 0.05 quantile of the statistic's own null at that N (thresholds.csv). The setup and helpers (sources,
### mirror de-duplication, n_mid >= 5 filter, Simes) are taken verbatim from prepare_data_v3.R, so the units are Fig. 2A's.
### Effect size = the gene's simulated effect_size (after the generator's scale factor), x 2 for the dosage 0 -> 2 scale (the
### same scale as Fig. 2D's x axis). It is the magnitude of the change in every component the effect type
### moves (μ, p0 and/or p1; equal size, generator line 201 / v12 line 209). Bins of 0.05 from 0 to 0.6, then >= 0.6.
### Built-in check: pooled over bins, power reproduces power.csv (Fig. 2A) at this N for every framework and model.
### Output: data_v3_<DATASET>/power_by_effect_N<N>/ {per_promoter.csv, power_by_effect_bin.csv, CHECKS.txt}. Reads only.
suppressMessages(library(data.table))
Sys.setenv(DATASET = Sys.getenv("DATASET", "nothr_realjunc"), SEEDTAG = "")
FIGROOT <- "/projects/group/bliu/e1327989/GTEx/new_simulation/figures_bohr2"
pd <- readLines(file.path(FIGROOT, "scripts", "prepare_data_v3.R"))
i_nul <- grep("^## ================= nulls", pd); i_s3 <- grep("^simes3 <- ", pd)
eval(parse(text = pd[1:(i_nul - 1)])); eval(parse(text = pd[i_s3:(i_s3 + 2)]))      # defines DATASET, NM, SRC, FfromLog, tag, mirror_map, dedup_ids, relabel, simes3
NN <- as.integer(Sys.getenv("EFFECT_N", "300")); stopifnot(NN %in% Ns)
DATA <- file.path(FIGROOT, paste0("data_v3_", DATASET)); OUTD <- file.path(DATA, sprintf("power_by_effect_N%d", NN)); dir.create(OUTD, showWarnings = FALSE)
TH <- fread(file.path(DATA, "thresholds.csv"))
lab <- c(middle_value = "mu", zero_proportion = "p0", one_proportion = "p1", mixed_effect = "mixed")
PWD <- rbindlist(lapply(SRC, function(src) {
  x <- readRDS(file.path(NM, src$pow, sprintf("N%d.rds", NN))); b3 <- readRDS(file.path(NM, src$b3, "power", sprintf("N%d.rds", NN)))$res[, perm := NULL]
  h <- x$header[is.na(skip) & effect_size > 0 & effect_type %in% names(lab)]
  if (src$m == "M2") { mp <- mirror_map(src, NN); h <- relabel(h, mp); h <- h[(n - n0 - n1) >= MIN_NMID]; h <- h[id %in% dedup_ids(h$id, mp)] }
  h[, set := lab[effect_type]]
  r <- merge(merge(x$res, h[, .(id, set, effect_size, nmid = n - n0 - n1)], by = "id"), b3, by = "id")
  r <- r[nmid >= MIN_NMID]
  data.table(m = src$m, id = tag(src$b, r$id), set = r$set, effect_0to2 = 2 * r$effect_size,
             BOHR = FfromLog(r$mu_lrt_logp, r$nmid), BOHR_p0 = exp(r$nu_firth_logp), BOHR_p1 = exp(r$tau_firth_logp),
             RINTfirst = r$rintfirst, raw = r$raw, normal = r$normal) }))
PWD[, BOHR_simes := simes3(BOHR, BOHR_p0, BOHR_p1)]
MODELS <- c(BOHR = "BOHR_simes", RINTfirst = "RINTfirst", raw = "raw", normal = "normal")
L <- rbindlist(lapply(names(MODELS), function(md) PWD[, .(m, N = NN, id, set, effect_0to2, model = md, stat = MODELS[[md]], p = get(MODELS[[md]]))]))
L <- merge(L, TH[, .(m, N, stat, t05)], by = c("m", "N", "stat")); L[, call := !is.na(p) & p <= t05]
BR <- c(seq(0, 0.6, 0.05), Inf); BL <- c(sprintf("%.2f-%.2f", head(BR, -2), BR[2:13]), "0.60 and above")
L[, bin := factor(BL[findInterval(effect_0to2, BR)], levels = BL)]
fwrite(L[order(m, model, id), .(framework = c(M1 = "Simulated activity", M2 = "Simulated RNA-seq")[m], N, promoter = id, effect_type = set,
                                 simulated_effect_0to2 = effect_0to2, effect_bin = bin, model, test = stat, p, threshold_alpha05 = t05, called = call)],
       file.path(OUTD, "per_promoter.csv"))
wilson <- function(k, n, z = qnorm(.975)) { p <- k / n; c0 <- p + z^2 / (2 * n); s <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)); d <- 1 + z^2 / n
  list((c0 - s) / d, (c0 + s) / d) }
S <- L[, .(promoters = .N, called = sum(call), mean_effect_0to2 = mean(effect_0to2)), by = .(m, model, bin)][order(m, model, bin)]
S[, power := called / promoters][, c("ci_lo", "ci_hi") := wilson(called, promoters)]
fwrite(S[, .(framework = c(M1 = "Simulated activity", M2 = "Simulated RNA-seq")[m], m, N = NN, model, effect_bin = bin, promoters, called,
             mean_effect_0to2 = round(mean_effect_0to2, 4), power_pct = round(100 * power, 2), ci95_lo_pct = round(100 * ci_lo, 2), ci95_hi_pct = round(100 * ci_hi, 2))],
       file.path(OUTD, "power_by_effect_bin.csv"))
pw <- fread(file.path(DATA, "power.csv"))[set == "all" & N == NN]
ck <- merge(L[, .(scan = mean(call), n = .N), by = .(m, model)], pw[, .(m, model, power_cal05, promoters)], by = c("m", "model"))
ok <- ck[, all(abs(scan - power_cal05) < 1e-12 & n == promoters)] && nrow(ck) == 8
chk <- c(sprintf("prepare_power_by_effect.R, DATASET %s, N = %d, %s", DATASET, NN, format(Sys.time(), "%Y-%m-%d %H:%M")),
         sprintf("[%s] pooled over bins = power.csv (Fig. 2A) for %d framework x model cells; promoters %s", if (ok) "PASS" else "FAIL", nrow(ck),
                 paste(ck[model == "BOHR", sprintf("%s %d", m, n)], collapse = ", ")),
         sprintf("[%s] every promoter in exactly one bin", if (L[, !anyNA(bin)]) "PASS" else "FAIL"))
writeLines(chk, file.path(OUTD, "CHECKS.txt")); cat(chk, sep = "\n"); stopifnot(ok)
print(dcast(S[model %in% c("BOHR", "raw")], m + bin ~ model, value.var = c("power", "promoters"))[, .(m, bin, n = promoters_BOHR, BOHR = round(100 * power_BOHR, 1), raw = round(100 * power_raw, 1))], nrows = 40)
