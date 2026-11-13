## Ungated units v2 (author 2026-10-01): mu units need n_mid >= 80 (estimability rule; perm0 tails were 2-5x at 1e-4 for
## n_mid < 80 and all six perm0 discoveries were such mu units). p0/p1 get NO count rule: their perm0 tails are flat
## across minor-class counts (Firth). Primary set = dedup_drop == FALSE & estimable; q = BH on it (q_v1 kept).
## Input: results/final_nogate/<def>/*.units.rds (v1, untouched: the running NSCC conditional reads those).
## Output: results/final_nogate_v2/<def>/<tissue>.<arm>.units.rds
suppressPackageStartupMessages(library(data.table))
V2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"; NMID_MIN <- 80L
for (def in c("freeze", "pilot")) { IN <- file.path(V2, "results/final_nogate", def); if (!dir.exists(IN)) next
  OUT <- file.path(V2, "results/final_nogate_v2", def); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  for (f in list.files(IN, pattern = "units[.]rds$", full.names = TRUE)) { U <- as.data.table(readRDS(f))
    U[, q_v1 := q]; U[, estimable := !(type == "mu" & n_mid < NMID_MIN)]
    U[, q := NA_real_]; U[dedup_drop == FALSE & estimable == TRUE, q := p.adjust(p, "BH")]
    U[, minor_count := fifelse(type == "mu", n_mid, fifelse(type == "p0", pmin(n0, n - n0), pmin(n1, n - n1)))]
    tmp <- file.path(OUT, paste0(basename(f), ".tmp")); saveRDS(U, tmp); file.rename(tmp, file.path(OUT, basename(f))) } }
S <- rbindlist(lapply(Sys.glob(file.path(V2, "results/final_nogate_v2/freeze/*.units.rds")), function(f) readRDS(f)[, .(arm = arm[1], sig = sum(q <= .05, na.rm = TRUE))]))
print(S[, .(files = .N, significant = sum(sig)), by = arm])
