## Ungated units v3: estimability rules applied BEFORE BH, as information-content rules:
##   mu    : n_mid >= 80                                        (v2 rule, perm0 calibration)
##   p0/p1 : minor class >= 10%, i.e. 0.10 <= boundary share <= 0.90 (p0: f0; p1: f1)   (NEW; symmetric)
## Primary set = dedup_drop == FALSE & estimable; q = BH within tissue on it. q_v1 (no rule) and q_v2 (mu rule only) kept.
## The high display thresholds (e.g. mu mid >= 0.70, p0/p1 >= 0.35) are NOT applied here: they are a credibility filter
## applied AFTER BH.
## Input : results/final_nogate/<def>/*.units.rds (v1, untouched)
## Output: results/final_nogate_v3/<def>/<tissue>.<arm>.units.rds
suppressPackageStartupMessages(library(data.table))
V2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"; NMID_MIN <- 80L; SHARE_MIN <- 0.10
for (def in c("freeze", "pilot")) { IN <- file.path(V2, "results/final_nogate", def); if (!dir.exists(IN)) next
  OUT <- file.path(V2, "results/final_nogate_v3", def); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  for (f in list.files(IN, pattern = "units[.]rds$", full.names = TRUE)) { U <- as.data.table(readRDS(f))
    U[, q_v1 := q]
    U[, estimable_mu_rule := !(type == "mu" & n_mid < NMID_MIN)]
    U[, boundary_share := fifelse(type == "p0", f0, fifelse(type == "p1", f1, NA_real_))]
    U[, estimable := estimable_mu_rule & (type == "mu" | (boundary_share >= SHARE_MIN & boundary_share <= 1 - SHARE_MIN))]
    U[, q_v2 := NA_real_]; U[dedup_drop == FALSE & estimable_mu_rule == TRUE, q_v2 := p.adjust(p, "BH")]
    U[, q := NA_real_]; U[dedup_drop == FALSE & estimable == TRUE, q := p.adjust(p, "BH")]
    U[, minor_count := fifelse(type == "mu", n_mid, fifelse(type == "p0", pmin(n0, n - n0), pmin(n1, n - n1)))]
    tmp <- file.path(OUT, paste0(basename(f), ".tmp")); saveRDS(U, tmp); file.rename(tmp, file.path(OUT, basename(f))) } }
S <- rbindlist(lapply(Sys.glob(file.path(V2, "results/final_nogate_v3/freeze/*.units.rds")), function(f) readRDS(f)[, .(arm = arm[1], tissue = tissue[1],
  sig_v3 = sum(q <= .05, na.rm = TRUE), sig_v2 = sum(q_v2 <= .05, na.rm = TRUE), thr_v3 = max(p[q <= .05], na.rm = TRUE), thr_v2 = max(p[q_v2 <= .05], na.rm = TRUE))]))
print(S[, .(files = .N, sig_v2 = sum(sig_v2), sig_v3 = sum(sig_v3), thr_looser = sum(thr_v3 > thr_v2, na.rm = TRUE), thr_stricter = sum(thr_v3 < thr_v2, na.rm = TRUE), thr_same = sum(thr_v3 == thr_v2, na.rm = TRUE)), by = arm])
fwrite(S, file.path(V2, "results/final_nogate_v3/freeze_thresholds_v2_v3.tsv"), sep = "\t")
