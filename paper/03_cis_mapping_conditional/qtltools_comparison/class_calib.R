## class_calib.R: BOHR's and QTLtools' calls per unit, by the part of BOHR's model that is significant (Fig. 2E).
## Fig. 2E uses UNIT_MODE=merge_all: the unit map, the real-arm calls of scheme B (B_bh), column near and the no-fill QTLtools calls.
## BOHR's discoveries two ways, each at BH 5% and at a permutation-matched FDR of 5%, per tissue:
##   S (the earlier Fig. 2E rule): unit P = Simes over the unit's component tests; a unit is a discovery by its unit P; its class = the
##       components that pass BH inside the unit at the unit-level level (Benjamini-Bogomolov: 0.05 R / m for BH, the unit P
##       threshold for the matched set);
##   D (direct): every component test (mu, p0, p1 of every unit) is a hypothesis of its own; a component is significant by
##       its own P; a unit is a discovery if any of its components is; class = which components.
##   BH: thresholds from BH 5% on the real data (S over unit P values, D over component P values).
##   matched: the largest set whose permutation-estimated FDR is <= 5% (mean of the permuted arms that are in).
## QTLtools (unit P = the smallest P of the unit's promoters): BH 5% and matched FDR 5% over units; modes raw, rint, normal;
##   kinds: pairs (missing values mean-filled; real + 3 permuted arms) and nofill (donors without a value left out; the arms
##   that are in for all 50 tissues).
## Units and component tests as fig2e_deliver_freeze_v3m_pairs.R (final_nogate_v3m units: dedup_drop == FALSE & estimable).
## Output (folder O): K_units_<set>.tsv.gz (flags per unit), K_summary_<set>.tsv (long: BOHR discoveries by
## class, how many QTLtools also finds, the units only QTLtools finds by BOHR's smallest component P), K_schemes_<set>.tsv.
## UNIT_MODE: a pair of promoters whose 0/1
## pattern is an exact mirror but whose interior values are not ("partial pair") shares one boundary test in BOHR (kept at one
## promoter, dropped at the other as a duplicate), while QTLtools tests both promoters and sees the same switch twice.
##   base  = as the earlier Fig. 2E rule: only pairs exact in both parts are one unit; a partial pair stays two units;
##   merge = a partial pair is one unit as well (BOHR: all its primary tests; QTLtools: the smaller P of its promoters);
##   drop  = the promoters of partial pairs are left out for both methods;
##   droppairs = every promoter of an exact or partial pair is left out (single promoters only).
## Output files of a mode other than base carry the suffix _<mode>.
## A unit only QTLtools finds is placed in a class by BOHR's components (column near): "both" when the smallest mu P and the
## smallest 0/1 P of the unit are both under 0.05, else the class of the smallest component P.
## Estimability (the rule is BOHR's own later control, minor class
## >= 10% for p0 / p1 and >= 80 interior samples for mu, and QTLtools has no counterpart). A suffix on UNIT_MODE:
##   (none)   = only BOHR's estimable component tests count;
##   _all     = every component test BOHR computed counts (finite P, not a dropped duplicate), estimable or not;
##   _dropne  = units with a computed test that is not estimable are left out for both methods.
## The units are those of Fig. 2E in every case (at least one estimable test), the promoters QTLtools was run on.
## Usage: Rscript class_calib.R <v4|paper> <ncores>
suppressPackageStartupMessages({ library(data.table); library(parallel) })
a <- commandArgs(TRUE); ver <- a[1]; NC <- as.integer(a[2]); if (is.na(NC)) NC <- 4
G <- "/projects/group/bliu/e1327989/GTEx"; W <- file.path(G, "nofloor_full"); K <- file.path(W, "consult/fig2e_20261119"); MV <- file.path(K, "matched_variants"); O <- file.path(W, "consult/why_fold_20261120/out")
V2 <- c(v4 = file.path(W, "GTEx/QTL_figure/rerun_v4"), paper = file.path(G, "QTL_figure/rerun_v2"))[[ver]]; set_name <- if (ver == "v4") "r5e80" else "paper"
ARMS <- c("real", "perm0", "perm1", "perm2"); MODES <- c("raw", "rint", "normal"); KD <- c(pairs = "pairs", nofill = "pairs_nofill")
UMF <- Sys.getenv("UNIT_MODE", "base"); EM <- if (grepl("_all$", UMF)) "all" else if (grepl("_dropne$", UMF)) "dropne" else "est"; UM <- sub("_(all|dropne)$", "", UMF)
stopifnot(UM %in% c("base", "merge", "drop", "droppairs")); SFX <- if (UMF == "base") "" else paste0("_", UMF)
fr <- function(f, ...) fread(cmd = paste("zcat", shQuote(f)), ...)
X <- fr(file.path(K, "main/out/units_paper_r5e80.tsv.gz"), select = c("set", "tissue", "unit", "med_D"))[set == set_name]; tis <- sort(unique(X$tissue))
qf <- function(kd, md, t, arm) file.path(MV, KD[[kd]], ver, md, paste(t, arm, "promoters.tsv.gz", sep = "."))
PA <- lapply(KD, function(d) ARMS[-1][vapply(ARMS[-1], function(arm) all(file.exists(file.path(MV, d, ver, "raw", paste(tis, arm, "promoters.tsv.gz", sep = ".")))), TRUE)])   # permuted arms that are in, per kind
cat(ver, ": permuted arms of QTLtools in for all tissues: pairs", PA$pairs, "| nofill", PA$nofill, "\n")
simes <- function(p) { p <- sort(p); k <- length(p); min(1, min(k * p / seq_len(k))) }
thr_bh <- function(p) { q <- p.adjust(p, "BH"); if (any(q <= 0.05)) max(p[q <= 0.05]) else -Inf }
thr_m <- function(pr, pp) { if (!length(pp)) return(NA_real_); s <- sort(pr); V <- Reduce(`+`, lapply(pp, function(q) findInterval(s, sort(q)))) / length(pp); ok <- which(V / seq_along(s) <= 0.05); if (length(ok)) s[max(ok)] else -Inf }
cls <- function(ty) if (!length(ty)) NA_character_ else if (all(ty == "mu")) "mu only" else if (!any(ty == "mu")) "p0/p1 only" else "both"
one <- function(t) { keep <- X[tissue == t, unit]
  CA <- lapply(setNames(ARMS, ARMS), function(arm) { U <- as.data.table(readRDS(file.path(V2, "results/final_nogate_v3m/freeze", paste(t, arm, "units.rds", sep = "."))))[dedup_drop == FALSE & is.finite(p)]; if (EM == "est") U <- U[estimable == TRUE]
    U[, unit0 := fifelse(boundary_exact %in% TRUE & interior_exact %in% TRUE, mirror_pair_id, promoter)]; U <- U[unit0 %in% keep]            # the units of Fig. 2E
    U[, partial := boundary_exact %in% TRUE & !(interior_exact %in% TRUE)]                                                                  # exact in the 0/1 pattern only
    if (UM == "drop") U <- U[partial == FALSE]; if (UM == "droppairs") U <- U[!(boundary_exact %in% TRUE | interior_exact %in% TRUE)]
    U[, unit := if (UM == "merge") fifelse(partial, mirror_pair_id, unit0) else unit0]; U[, .(unit, unit0, promoter, type, p, gene, partial, estimable)] })
  CA <- lapply(CA, function(d) d[unit %in% d[estimable == TRUE, unit]])                    # the units of Fig. 2E: at least one estimable test
  n_ne <- uniqueN(CA$real[estimable == FALSE, unit]); if (EM == "dropne") { ne <- unique(CA$real[estimable == FALSE, unit]); CA <- lapply(CA, function(d) d[!unit %in% ne]) }
  SA <- lapply(CA, function(d) d[, .(p = simes(p)), by = unit]); m <- nrow(SA$real)
  b <- fr(file.path(V2, "results/fig2e_freeze_v3m/bohr", paste0(t, ".real.tsv.gz")), select = c("promoter", "p_simes", "promoters_in_unit")); setnames(b, c("unit", "p_stored", "piu")); b <- b[unit %in% keep]
  if (UMF == "base") { ck <- merge(SA$real, b, by = "unit"); stopifnot(nrow(ck) == m, nrow(b) == m, max(abs(ck$p - ck$p_stored) / pmax(ck$p_stored, 1e-300)) < 1e-6) }       # the unit P of Fig. 2E, reproduced
  um0 <- b[, .(promoter = unlist(strsplit(piu, ";", fixed = TRUE))), by = .(unit0 = unit)]; um <- unique(merge(um0, unique(CA$real[, .(unit0, unit)]), by = "unit0")[, .(unit, promoter)])   # promoters QTLtools tests for the unit
  ui <- CA$real[, .(gene = gene[1], partial = any(partial), unit0s = uniqueN(unit0)), by = unit]; ui <- merge(ui, merge(unique(CA$real[, .(unit, unit0)]), X[tissue == t, .(unit0 = unit, med_D)], by = "unit0")[, .(med_D = max(med_D)), by = unit], by = "unit")
  ## thresholds
  tS <- c(bh = thr_bh(SA$real$p), m = thr_m(SA$real$p, lapply(SA[-1], `[[`, "p"))); R <- sum(SA$real$p <= tS["bh"]); lvlS <- c(bh = 0.05 * R / m, m = unname(tS["m"]))
  umin <- function(d) d[, .(p = min(p)), by = unit]$p                                # under D a unit is a discovery when its smallest component P passes the threshold
  tD <- c(bh = thr_bh(CA$real$p), m = thr_m(umin(CA$real), lapply(CA[-1], umin)))   # matched at the unit level, as for QTLtools (matching over component tests left the units at about 6%)
  Cr <- CA$real; Ct <- Cr[, .(p = min(p)), by = .(unit, type)]; Ct[, padj := p.adjust(p, "BH"), by = unit]      # S: the class as in the stored unit table (one P per component type, BH inside the unit)
  u <- SA$real[, .(unit, pS = p)]
  for (h in c("bh", "m")) {
    cs <- Ct[unit %in% u[pS <= tS[h], unit] & padj <= lvlS[h], .(cl = cls(type)), by = unit]; u[, paste0("S_", h) := cs$cl[match(unit, cs$unit)]]
    cd <- Cr[p <= tD[h], .(cl = cls(type)), by = unit]; u[, paste0("D_", h) := cd$cl[match(unit, cd$unit)]] }
  Cr[, qc := p.adjust(p, "BH")]                                              # D: BH over all component tests of the tissue
  nr <- Cr[, .(near = if (type[which.min(p)] == "mu") "mu only" else "p0/p1 only", qD = min(qc), pD = min(p), p_mu = if (any(type == "mu")) min(p[type == "mu"]) else NA_real_,
               p_b = if (any(type != "mu")) min(p[type != "mu"]) else NA_real_, n_comp = .N), by = unit]
  nr[!is.na(p_mu) & !is.na(p_b) & p_mu < 0.05 & p_b < 0.05, near := "both"]      # rule: both components have P < 0.05
  u <- merge(u, nr, by = "unit"); u <- merge(u, ui, by = "unit")
  ## permuted discoveries of BOHR under each scheme (units), for the estimated FDR at the unit level
  fp <- c(S_bh = mean(vapply(SA[-1], function(d) sum(d$p <= tS["bh"]), 0)), S_m = mean(vapply(SA[-1], function(d) sum(d$p <= tS["m"]), 0)),
          D_bh = mean(vapply(CA[-1], function(d) uniqueN(d[p <= tD["bh"], unit]), 0)), D_m = mean(vapply(CA[-1], function(d) uniqueN(d[p <= tD["m"], unit]), 0)))
  ex <- function(d) if (EM != "all") 0 else d[, .(pe = min(c(Inf, p[estimable == TRUE])), pn = min(c(Inf, p[estimable == FALSE]))), by = unit][pn <= tD["bh"] & pe > tD["bh"], .N]   # units found through a not estimable test only
  fp <- c(fp, NE_real = ex(CA$real), NE_perm = mean(vapply(CA[-1], ex, 0)), NE_units = n_ne, tests = nrow(CA$real), tests_ne = CA$real[estimable == FALSE, .N])
  u[, ne_only := unit %in% CA$real[, .(pe = min(c(Inf, p[estimable == TRUE])), pn = min(c(Inf, p[estimable == FALSE]))), by = unit][pn <= tD["bh"] & pe > tD["bh"], unit]]
  ## U and J. U: a unit is a discovery when a component passes on its own (D) or
  ## its Simes P passes (S), each at its own BH threshold. J: one BH family per tissue of every component test and, for units with two or more tests, the unit's
  ## Simes P; a unit is a discovery when any of its tests in the family passes. A unit found through the combined test only takes the class of column near.
  u[, U_bh := fifelse(!is.na(D_bh), D_bh, fifelse(!is.na(S_bh), near, NA_character_))]
  multi <- function(d) d[, .N, by = unit][N >= 2, unit]; tJ <- thr_bh(c(Cr$p, SA$real[unit %in% multi(Cr), p]))
  cj <- Cr[p <= tJ, .(cl = cls(type)), by = unit]; u[, J_bh := cj$cl[match(unit, cj$unit)]]; u[is.na(J_bh) & n_comp >= 2 & pS <= tJ, J_bh := near]
  ## B: D, plus the units
  ## in which no component passes on its own, the Simes P passes, and the mu part and the 0/1 part both have P < 0.05 (near == "both"); these go to the class "both".
  u[, B_bh := fifelse(!is.na(D_bh), D_bh, fifelse(!is.na(S_bh) & near == "both", "both", NA_character_))]
  pu <- function(d, sm) { a <- merge(sm[, .(unit, ps = p)], d[, .(pm = min(p), k = .N, pmu = min(c(Inf, p[type == "mu"])), pb = min(c(Inf, p[type != "mu"]))), by = unit], by = "unit")
    c(U = a[ps <= tS["bh"] | pm <= tD["bh"], .N], J = a[pm <= tJ | (k >= 2 & ps <= tJ), .N], B = a[pm <= tD["bh"] | (ps <= tS["bh"] & pmu < 0.05 & pb < 0.05), .N]) }
  pj <- vapply(ARMS[-1], function(arm) pu(CA[[arm]], SA[[arm]]), c(U = 0, J = 0, B = 0)); fp <- c(fp, U_bh = mean(pj["U", ]), J_bh = mean(pj["J", ]), B_bh = mean(pj["B", ]))
  ## QTLtools
  qfp <- list()
  for (kd in names(KD)) for (md in MODES) { rd1 <- function(arm) merge(fr(qf(kd, md, t, arm), select = c("promoter", "p")), um, by = "promoter")[!is.na(p), .(p = min(p)), by = unit]
    qr <- rd1("real"); pp <- lapply(PA[[kd]], function(arm) rd1(arm)$p); tq <- c(bh = thr_bh(qr$p), m = thr_m(qr$p, pp))
    for (h in c("bh", "m")) { u[, paste("Q", kd, md, h, sep = "_") := if (is.na(tq[h])) NA else qr$p[match(unit, qr$unit)] <= tq[h]]; set(u, which(is.na(u[[paste("Q", kd, md, h, sep = "_")]]) & !is.na(tq[h])), paste("Q", kd, md, h, sep = "_"), FALSE)
      qfp[[paste(kd, md, h, sep = "_")]] <- if (is.na(tq[h]) || !length(pp)) NA_real_ else mean(vapply(pp, function(q) sum(q <= tq[h]), 0)) } }
  list(u = u[, `:=`(set = set_name, tissue = t)], map = um[, `:=`(set = set_name, tissue = t)], fp = data.table(set = set_name, tissue = t, what = c(names(fp), names(qfp)), fp = c(fp, unlist(qfp))),
       thr = data.table(set = set_name, tissue = t, m_units = m, m_comp = nrow(Cr), tS_bh = tS["bh"], tS_m = tS["m"], tD_bh = tD["bh"], tD_m = tD["m"], tJ_bh = tJ, m_simes = length(multi(Cr)))) }
RS <- mclapply(tis, one, mc.cores = NC); bad <- vapply(RS, inherits, TRUE, "try-error"); if (any(bad)) { print(RS[bad][[1]]); stop("a tissue failed") }
U <- rbindlist(lapply(RS, `[[`, "u"), use.names = TRUE, fill = TRUE); FP <- rbindlist(lapply(RS, `[[`, "fp")); TH <- rbindlist(lapply(RS, `[[`, "thr"))
fwrite(U, file.path(O, paste0("K_units_", set_name, SFX, ".tsv.gz")), sep = "\t"); fwrite(TH, file.path(O, paste0("K_thresholds_", set_name, SFX, ".tsv")), sep = "\t")
fwrite(rbindlist(lapply(RS, `[[`, "map")), file.path(O, paste0("K_map_", set_name, SFX, ".tsv.gz")), sep = "\t")
xx <- FP[what %in% c("NE_real", "NE_perm", "NE_units", "tests", "tests_ne"), .(v = sum(fp)), by = what]; cat("estimability mode", EM, ":", paste(xx$what, round(xx$v, 1), collapse = " | "), "\n"); fwrite(xx[, set := set_name], file.path(O, paste0("K_ne_", set_name, SFX, ".tsv")), sep = "\t")
cat("unit mode", UM, ": units", nrow(U), "| units of a partial pair", sum(U$partial), "| merged from two", sum(U$unit0s > 1), "\n")
options(width = 250)
## 1. the two schemes: discoveries, classes, estimated FDR at the unit level
sc <- rbindlist(lapply(c("S_bh", "D_bh", "S_m", "D_m", "U_bh", "J_bh", "B_bh"), function(k) { x <- U[[k]]; data.table(set = set_name, scheme = k, units = sum(!is.na(x)), mu_only = sum(x == "mu only", na.rm = TRUE), p0p1_only = sum(x == "p0/p1 only", na.rm = TRUE),
  both = sum(x == "both", na.rm = TRUE), permuted_units = round(FP[what == k, sum(fp)], 1)) }))[, est_fdr := round(permuted_units / units, 4)]
cat("\n== BOHR under the two schemes (S = Simes per unit, D = each component on its own), BH 5% and matched FDR 5%\n"); print(sc); fwrite(sc, file.path(O, paste0("K_schemes_", set_name, SFX, ".tsv")), sep = "\t")
if (UMF == "base") { XS <- fr(file.path(K, "main/out/units_paper_r5e80.tsv.gz"), select = c("set", "tissue", "unit", "sig_BOHR", "class"))[set == set_name]; ck <- merge(U[, .(tissue, unit, S_bh)], XS, by = c("tissue", "unit"))
  cat("\ncheck against the stored unit table: same discovery status", ck[, sum(!is.na(S_bh) == sig_BOHR)], "of", nrow(ck), "| same class among discoveries", ck[sig_BOHR == TRUE, sum(S_bh == class, na.rm = TRUE)], "of", ck[sig_BOHR == TRUE, .N], "\n") }
cat("\n== same units? S against D at BH 5% (rows S, columns D; 'none' = not a discovery)\n"); print(dcast(U[, .N, by = .(S = fifelse(is.na(S_bh), "none", S_bh), D = fifelse(is.na(D_bh), "none", D_bh))], S ~ D, value.var = "N", fill = 0))
## 2. by class against QTLtools
L <- list()
for (sch in c("S", "D", "U", "J", "B")) for (h in c("bh", "m")) for (kd in names(KD)) for (md in MODES) { qc <- paste("Q", kd, md, h, sep = "_"); if (all(is.na(U[[qc]])) || !paste0(sch, "_", h) %in% names(U)) next
  bc <- U[[paste0(sch, "_", h)]]; q <- U[[qc]]
  for (k in c("mu only", "p0/p1 only", "both")) L[[length(L) + 1]] <- data.table(set = set_name, scheme = sch, thr = h, kind = kd, mode = md, class = k, bohr = sum(bc == k, na.rm = TRUE), also = sum(bc == k & q, na.rm = TRUE),
                                                                                  only = sum(is.na(bc) & q & U$near == k)) }
Z <- rbindlist(L); Z[, `:=`(fold_also = round(bohr / also, 3), fold_all = round(bohr / (also + only), 3))]
tot <- Z[, .(class = "all", bohr = sum(bohr), also = sum(also), only = sum(only)), by = .(set, scheme, thr, kind, mode)][, `:=`(fold_also = round(bohr / also, 3), fold_all = round(bohr / (also + only), 3))]
Z <- rbind(Z, tot)[order(scheme, thr, kind, match(mode, MODES), match(class, c("mu only", "p0/p1 only", "both", "all")))]
qf2 <- FP[grepl("^(pairs|nofill)_", what), .(fp = sum(fp)), by = what]; Z[, q_fp := qf2$fp[match(paste(kind, mode, thr, sep = "_"), qf2$what)]]
cat("\n== by class (also = QTLtools finds it too; only = QTLtools only, by BOHR's smallest component P; fold_also = bohr / also; fold_all = bohr / (also + only))\n"); print(Z[, !"q_fp"])
fwrite(Z, file.path(O, paste0("K_summary_", set_name, SFX, ".tsv")), sep = "\t"); cat("DONE\n")
