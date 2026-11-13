## fig2e_deliver_freeze_v3m_pairs.R (2026-11-04): Fig 2e on the mirror rule v3 units (results/final_nogate_v3m; author 2026-11-04
## final: exact switches are p0 on the silent TSS), with EXACT MIRROR PAIRS COUNTED ONCE (author 2026-11-04: count a mirror pair once).
## Why: for an exact pair (boundary_exact & interior_exact) y_B = 1 - y_A in every donor, so QTLtools' linear test is identical for
## the two promoters; under v3 both promoters carry primary units, and promoter-level counting would count one switch twice
## (QTLtools 2,258 such pairs vs BOHR 892; code/fig2e_v3m_pair_collapse.R). Comparison unit = the exact pair, else the promoter:
##   BOHR unit p = Simes over the primary units of the unit's promoters; QTLtools unit p = the smallest p of its promoters
##   (identical for exact pairs); BH over units; layer classes (Benjamini-Bogomolov) within the unit.
## Everything else as fig2e_deliver_freeze.R. Output: results/fig2e_freeze_v3m/ (promoter-level reference in by_promoter_reference/).
## Figure 2e delivery on the FREEZE definition (author 2026-10-02: one definition, the real data; plan B added perm1 +
## perm2 for all 50 tissues). Same decisions and output layout as fig2e_deliver.R (version 2, pilot), which stays as is:
##   * BOHR promoter p = Simes over the promoter's PRIMARY units (dedup_drop == FALSE & estimable), BH over promoters;
##     layer classes from Benjamini-Bogomolov stage 2 (BH within the selected promoter at level q * R / m).
##   * units v3 (results/final_nogate_v3m/freeze): mu n_mid >= 80; p0/p1 minor class >= 10%; both before BH.
##   * QTLtools as run (permutation pass, adj_beta_pval = V19) in three modes, restricted to the BOHR-tested promoters
##     before BH. Inputs: QTLtools_linear_compare/freeze_all/results/{rint,raw,normal}/<T>/prm.<arm>.chrN.txt (BOHR's
##     donors, promoters, covariates and permutation maps).
## Arms: real, perm0, perm1, perm2. Tissues: env FIG2E_TIS ("ALL" = every tissue complete in all arms and modes;
## default the three tissues of the current panel).
## Output: results/fig2e_freeze/ (counts.tsv, counts_display.tsv, bohr/<T>.<arm>.tsv.gz, qtltools/<mode>/<T>.<arm>.tsv.gz)
suppressPackageStartupMessages(library(data.table))
V2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"; QA <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/QTLtools_linear_compare/freeze_all/results"
OUT <- file.path(V2, "results/fig2e_freeze_v3m"); MODES <- c(rintfirst = "rint", raw = "raw", normal = "normal")
for (d in c("bohr", paste0("qtltools/", names(MODES)))) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
simes <- function(p) { p <- sort(p); k <- length(p); min(1, min(k * p / seq_len(k))) }
ARMS <- c("real", "perm0", "perm1", "perm2")
qfiles <- function(tis, arm, mode) Sys.glob(file.path(QA, MODES[[mode]], tis, paste0("prm.", arm, ".chr*.txt")))
ufile <- function(tis, arm) file.path(V2, "results/final_nogate_v3m/freeze", paste(tis, arm, "units.rds", sep = "."))
ts <- Sys.getenv("FIG2E_TIS", "Brain_Hippocampus Adipose_Visceral_Omentum Skin_Sun_Exposed_Lower_leg")
TIS <- if (ts == "ALL") sort(basename(Sys.glob(file.path(QA, "rint", "*")))) else strsplit(ts, " ")[[1]]
ok <- vapply(TIS, function(t) all(file.exists(ufile(t, ARMS))) && all(vapply(names(MODES), function(md) all(vapply(ARMS, function(a) length(qfiles(t, a, md)) == 22, TRUE)), TRUE)), TRUE)
if (any(!ok)) cat("skipped (incomplete arm or mode):", TIS[!ok], "\n"); TIS <- TIS[ok]
res <- list()
for (tis in TIS) for (arm in ARMS) {
  U <- as.data.table(readRDS(ufile(tis, arm)))[dedup_drop == FALSE & estimable == TRUE & is.finite(p)]
  U[, unit := fifelse(boundary_exact %in% TRUE & interior_exact %in% TRUE, mirror_pair_id, promoter)]
  U[, promoters_in_unit := paste(sort(unique(promoter)), collapse = ";"), by = unit]
  U[, promoter_real := promoter][, promoter := unit]                                   # the comparison unit replaces the promoter below
  B <- U[, .(gene = gene[1], promoters_in_unit = promoters_in_unit[1], p_simes = simes(p), types = paste(sort(paste0(type, if (uniqueN(promoter_real) > 1) paste0("@", promoter_real) else "")), collapse = ","),
             p_mu = if (any(type == "mu")) min(p[type == "mu"]) else NA_real_, p_p0 = if (any(type == "p0")) min(p[type == "p0"]) else NA_real_,
             p_p1 = if (any(type == "p1")) min(p[type == "p1"]) else NA_real_), by = promoter]
  B[, q := p.adjust(p_simes, "BH")][, significant := q <= 0.05]
  m <- nrow(B); R <- B[significant == TRUE, .N]; lvl <- 0.05 * R / m
  L <- U[promoter %in% B[significant == TRUE]$promoter][, rej := p.adjust(p, "BH") <= lvl, by = promoter][,
         .(layer_class = { s <- type[rej]; if (all(s == "mu")) "mu only" else if (!any(s == "mu")) "p0/p1 only" else "both" }), by = promoter]
  B <- merge(B, L, by = "promoter", all.x = TRUE)
  for (mode in names(MODES)) {
    d <- rbindlist(lapply(qfiles(tis, arm, mode), function(f) if (file.size(f) > 0) fread(f, header = FALSE) else NULL))
    d <- d[is.finite(suppressWarnings(as.numeric(V19))), .(promoter = V1, p = as.numeric(V19))]; n_bed <- nrow(d)
    um <- unique(U[, .(promoter_real, unit = promoter)])                              # QTLtools promoter -> comparison unit
    d <- merge(d, um, by.x = "promoter", by.y = "promoter_real")[, .(p = min(p), promoters_in_unit = paste(sort(promoter), collapse = ";")), by = .(promoter = unit)]
    d <- d[promoter %in% B$promoter][, q := p.adjust(p, "BH")][, significant := q <= 0.05]
    B[, (paste0("qtltools_", mode, "_significant")) := promoter %in% d[significant == TRUE]$promoter]
    res[[length(res) + 1]] <- data.table(tissue = tis, arm = arm, qtltools_mode = mode, promoters = m, bohr = R,
      bohr_mu_only = L[layer_class == "mu only", .N], bohr_both = L[layer_class == "both", .N], bohr_p0p1_only = L[layer_class == "p0/p1 only", .N],
      qtltools = d[significant == TRUE, .N], qtltools_tested = nrow(d), qtltools_bed_promoters = n_bed)
    fwrite(d, file.path(OUT, "qtltools", mode, paste(tis, arm, "tsv.gz", sep = ".")), sep = "\t") }
  fwrite(B[order(p_simes)], file.path(OUT, "bohr", paste(tis, arm, "tsv.gz", sep = ".")), sep = "\t") }
Rt <- rbindlist(res); fwrite(Rt, file.path(OUT, "counts.tsv"), sep = "\t")
L2 <- rbind(unique(Rt[qtltools_mode == "rintfirst", .(tissue, arm, method = "BOHR", n = bohr, promoters)]),
            Rt[, .(tissue, arm, method = paste0("QTLtools ", qtltools_mode), n = qtltools, promoters)])
D <- L2[, .(promoters = promoters[arm == "real"], real = n[arm == "real"], n_perm = sum(arm != "real"), perm_mean = round(mean(n[arm != "real"]), 2),
            perm_min = min(n[arm != "real"]), perm_max = max(n[arm != "real"])), by = .(tissue, method)]
D <- merge(D, unique(Rt[arm == "real" & qtltools_mode == "rintfirst", .(tissue, bohr_mu_only, bohr_both, bohr_p0p1_only)]), by = "tissue")
D[method != "BOHR", c("bohr_mu_only", "bohr_both", "bohr_p0p1_only") := NA_integer_]
fwrite(D, file.path(OUT, "counts_display.tsv"), sep = "\t"); cat("tissues:", length(TIS), "\n"); print(D[tissue %in% head(TIS, 3)])
