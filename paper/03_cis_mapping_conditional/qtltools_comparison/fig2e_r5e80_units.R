## fig2e_r5e80_units.R: the units of the Fig. 2E comparison under the final promoter rules, with tables of the stored calls.
## The steps that follow read only the unit list out/units_paper_r5e80.tsv.gz (set, tissue, unit, med_D); the other tables describe the calls.
## Unit-level tables of the stored Fig. 2E calls for the earlier run ("paper", all units) and the final rules ("r5e80": v4 units with gene median reads >= 5 and
## neither share of 0 nor share of 1 >= 80%; BH again within tissue, method and arm; no refit), as in
## the earlier filter-only computation (variant G5E80). Then:
##   A  totals (they reproduce the filter-only totals);
##   B  paper -> r5e80 split of every method's discoveries: units tested in both, in the paper only, in r5e80 only;
##   C  BOHR discoveries by component class and the share QTLtools also finds; QTLtools-only discoveries;
##   D  permuted discoveries by the unit's largest share at one value and by tissue;
##   E  discoveries / permuted discoveries with exact Poisson limits on the permuted count;
##   F  discoveries at a matched permutation FDR (step-up on mean permuted count / real count, per tissue);
##   H  per-tissue table for the sample-size display (four methods).
## Aggregates only. Usage: Rscript fig2e_r5e80_units.R <ncores>
suppressPackageStartupMessages({ library(data.table); library(parallel) })
NC <- as.integer(commandArgs(TRUE)[1]); if (is.na(NC)) NC <- 12
G <- "/projects/group/bliu/e1327989/GTEx"; W <- file.path(G, "nofloor_full")
K <- file.path(W, "consult/fig2e_20261119/main"); OUT <- file.path(K, "out"); dir.create(OUT, FALSE, TRUE)
C0 <- file.path(W, "consult/version_choice_20261118/main")
fr <- function(f, ...) fread(cmd = paste("zcat", shQuote(f)), ...)
ROOT <- c(paper = file.path(G, "QTL_figure/rerun_v2/results/fig2e_freeze_v3m"), v4 = file.path(W, "GTEx/QTL_figure/rerun_v4/results/fig2e_freeze_v3m"))
FRZ <- c(paper = file.path(G, "QTL_figure/rerun_v2/work/final/freeze"), v4 = file.path(W, "GTEx/QTL_figure/rerun_v4/work/final/freeze"))
METH <- c(BOHR = "bohr", raw = "qtltools/raw", normal = "qtltools/normal", rint = "qtltools/rintfirst")
ARMS <- c("real", "perm0", "perm1", "perm2")
tis <- sub("[.]real[.]tsv[.]gz$", "", list.files(file.path(ROOT[["v4"]], "bohr"), "[.]real[.]tsv[.]gz$"))
DS <- fr(file.path(C0, "out/promoter_depth_v4.tsv.gz"))[, .(tissue, promoter, med_D)]
info <- unique(rbind(fread(file.path(W, "GTEx/v11/Promoter_info_filter_5.txt"))[, .(promoter = V1, gene = geneId)],
                     fread(file.path(G, "v11/Promoter_info_filter_5.txt"))[, .(promoter = V1, gene = geneId)]), by = "promoter")
GD <- merge(DS, info, by = "promoter")[, .(med_D = max(med_D, na.rm = TRUE)), by = .(tissue, gene)]

one <- function(t) {
  P <- list(); U <- list()
  for (ver in names(ROOT)) {
    a <- readRDS(file.path(FRZ[[ver]], paste0(t, ".rds")))$act
    sh <- data.table(p1st = colnames(a), nfin = colSums(!is.na(a)), n0 = colSums(a == 0, na.rm = TRUE), n1 = colSums(a == 1, na.rm = TRUE)); ns <- nrow(a); rm(a)
    b <- fr(file.path(ROOT[[ver]], "bohr", paste0(t, ".real.tsv.gz")), select = c("promoter", "gene", "promoters_in_unit", "p_mu", "p_p0", "p_p1"))
    u <- merge(b[, .(unit = promoter, gene, p1st = sub(";.*$", "", promoters_in_unit), mirror = grepl(";", promoters_in_unit), p_mu, p_p0, p_p1)], sh, by = "p1st", all.x = TRUE)
    u <- merge(u, GD[tissue == t, .(gene, med_D)], by = "gene", all.x = TRUE)
    U[[ver]] <- u[, .(version = ver, tissue = t, n_samples = ns, unit, gene, mirror, p_mu, p_p0, p_p1, med_D, s0 = n0 / nfin, s1 = n1 / nfin)]
    for (m in names(METH)) for (arm in ARMS) {
      x <- fr(file.path(ROOT[[ver]], METH[[m]], paste0(t, ".", arm, ".tsv.gz")), select = c("promoter", if (m == "BOHR") "p_simes" else "p"))
      setnames(x, c("unit", "p")); P[[paste(ver, m, arm)]] <- x[, `:=`(version = ver, tissue = t, method = m, arm = arm)] } }
  list(P = rbindlist(P), U = rbindlist(U)) }
R <- mclapply(tis, one, mc.cores = NC)
bad <- which(vapply(R, inherits, TRUE, "try-error")); if (length(bad)) { print(R[bad]); stop("failed: ", paste(tis[bad], collapse = ",")) }
P <- rbindlist(lapply(R, `[[`, "P")); U <- rbindlist(lapply(R, `[[`, "U")); rm(R); gc()
P <- merge(P, U[, .(version, tissue, unit, med_D, s0, s1)], by = c("version", "tissue", "unit"), all.x = TRUE); P <- P[!is.na(p)]

## the two sets: paper as run; r5e80 = v4 filtered (gene median reads >= 5, neither share of 0 nor of 1 >= 80%)
P[, keep := (version == "paper") | (!is.na(med_D) & med_D >= 5 & pmax(s0, s1) < 0.80)]
D <- P[keep == TRUE][, set := fifelse(version == "paper", "paper", "r5e80")]
V4 <- P[version == "v4"][, set := "v4"]                                  # v4 as run, for section D only
D <- rbind(D, V4)
D[, q := p.adjust(p, "BH"), by = .(set, tissue, method, arm)]
thr <- D[arm == "real", .(thr = if (any(q <= 0.05)) max(p[q <= 0.05]) else -Inf), by = .(set, tissue, method)]
D <- merge(D, thr, by = c("set", "tissue", "method"))
D[, `:=`(sig = q <= 0.05, hit = p <= thr, smax = pmax(s0, s1), bshare = s0 + s1)]
options(width = 250, datatable.print.nrows = 300)
show <- function(nm, x) { fwrite(x, file.path(OUT, paste0(nm, ".tsv")), sep = "\t"); cat("\n==", nm, "\n"); print(x) }

## A totals
bt <- D[, .(sig = sum(sig), hit = sum(hit), tests = .N), by = .(set, tissue, method, arm)]
A <- bt[, .(tests = sum(tests[arm == "real"]), real = sum(sig[arm == "real"]), perm_bh = sum(sig[arm != "real"]) / 3, perm_events = sum(sig[arm != "real"]),
            fp = sum(hit[arm != "real"]) / 3), by = .(set, method)][, rFDR := round(fp / real, 4)]
A[, BOHR_over := round(real[method == "BOHR"] / real, 3), by = set]
show("A_totals", A[order(set, method)])

## wide real-arm table (one row per set x tissue x unit)
Wd <- dcast(D[arm == "real"], set + tissue + unit ~ method, value.var = c("p", "sig"))
Wd <- merge(Wd, U[, .(version, tissue, unit, gene, mirror, n_samples, p_mu, p_p0, p_p1, med_D, s0, s1)][, set_v := version], by = c("tissue", "unit"), allow.cartesian = TRUE)
Wd <- Wd[(set == "paper" & set_v == "paper") | (set %in% c("r5e80", "v4") & set_v == "v4")][, c("set_v", "version") := NULL]
for (m in names(METH)) { s <- paste0("sig_", m); Wd[is.na(get(s)), (s) := FALSE] }
## BOHR component class under the set's own BH (fig2e_deliver_freeze_v3m_pairs.R: within a discovery, components rejected by BH at level 0.05 R / m;
## here from the stored per-type minimum P values, so mirror pairs with two tests of one type are approximated)
lv <- Wd[, .(lvl = 0.05 * sum(sig_BOHR) / .N), by = .(set, tissue)]
Wd <- merge(Wd, lv, by = c("set", "tissue"))
cls <- function(pm, p0, p1, lvl) { p <- c(mu = pm, p0 = p0, p1 = p1); p <- p[!is.na(p)]; if (!length(p)) return(NA_character_)
  r <- names(p)[p.adjust(p, "BH") <= lvl]; if (!length(r)) return("none") ; if (all(r == "mu")) "mu only" else if (!any(r == "mu")) "p0/p1 only" else "both" }
Wd[sig_BOHR == TRUE, class := mapply(cls, p_mu, p_p0, p_p1, lvl)]
fwrite(Wd[set != "v4"], file.path(OUT, "units_paper_r5e80.tsv.gz"), sep = "\t")

## B split paper -> r5e80
X <- Wd[set != "v4"]
kp <- X[set == "paper", paste(tissue, unit)]; kr <- X[set == "r5e80", paste(tissue, unit)]; both <- intersect(kp, kr)
X[, part := fifelse(paste(tissue, unit) %in% both, "tested in both", fifelse(set == "paper", "paper only", "r5e80 only"))]
B <- X[, .(tests = .N, BOHR = sum(sig_BOHR), raw = sum(sig_raw), rint = sum(sig_rint), normal = sum(sig_normal),
           boundary_share_mean = round(mean(s0 + s1, na.rm = TRUE), 3), gene_reads_median = as.numeric(median(med_D, na.rm = TRUE)),
           raw_found_of_BOHR = round(sum(sig_BOHR & sig_raw) / sum(sig_BOHR), 3)), by = .(set, part)][order(set, part)]
B[, `:=`(BOHR_pct = round(100 * BOHR / tests, 2), raw_pct = round(100 * raw / tests, 2), BOHR_over_raw = round(BOHR / raw, 3), BOHR_over_rint = round(BOHR / rint, 3))]
show("B_split", B)
## why units leave or enter
pv4 <- unique(P[version == "v4" & arm == "real" & method == "BOHR", .(tissue, unit, med_D, smax = pmax(s0, s1))])
po <- merge(X[set == "paper" & part == "paper only", .(tissue, unit)], pv4, by = c("tissue", "unit"), all.x = TRUE)
show("B_paper_only_reason", po[, .(units = .N), by = .(reason = fifelse(is.na(smax) & is.na(med_D), "not a v4 unit",
       fifelse(!is.na(med_D) & med_D < 5, "gene reads < 5", fifelse(!is.na(smax) & smax >= 0.8, "share at one value >= 80%", "other"))))])
## same units, joint state paper -> r5e80 (BOHR / raw)
st <- function(b, r) fifelse(b & r, "both", fifelse(b, "BOHR only", fifelse(r, "raw only", "neither")))
S <- merge(X[set == "paper" & part == "tested in both", .(tissue, unit, paper = st(sig_BOHR, sig_raw), b_paper = s0 + s1)],
           X[set == "r5e80" & part == "tested in both", .(tissue, unit, r5e80 = st(sig_BOHR, sig_raw), b_r5e80 = s0 + s1)], by = c("tissue", "unit"))
show("B_joint_state", dcast(S[, .N, by = .(paper, r5e80)], paper ~ r5e80, value.var = "N", fill = 0))
show("B_boundary_share_same_units", S[, .(units = .N, boundary_share_paper = round(mean(b_paper, na.rm = TRUE), 3), boundary_share_r5e80 = round(mean(b_r5e80, na.rm = TRUE), 3))])
## r5e80 by gene reads and by boundary share
X[, depth := cut(med_D, c(-Inf, 5, 10, 20, 50, 100, Inf), right = FALSE, labels = c("<5", "5-10", "10-20", "20-50", "50-100", ">=100"))]
X[, bbin := cut(s0 + s1, c(-Inf, 0.05, 0.1, 0.3, 0.5, 0.8, Inf), right = FALSE, labels = c("<5%", "5-10%", "10-30%", "30-50%", "50-80%", ">=80%"))]
sm <- function(d, by) d[, .(tests = .N, BOHR = sum(sig_BOHR), raw = sum(sig_raw), rint = sum(sig_rint), normal = sum(sig_normal),
                            BOHR_over_raw = round(sum(sig_BOHR) / sum(sig_raw), 3), BOHR_over_rint = round(sum(sig_BOHR) / sum(sig_rint), 3),
                            raw_found_of_BOHR = round(sum(sig_BOHR & sig_raw) / sum(sig_BOHR), 3)), by = by][order(get(by[1]), get(by[2]))]
show("B_by_gene_reads", sm(X, c("set", "depth"))); show("B_by_boundary_share", sm(X, c("set", "bbin")))
show("B_by_part_boundary", sm(X, c("set", "part", "bbin"))[order(set, part, bbin)])

## C component classes and overlap
Cc <- X[sig_BOHR == TRUE, .(BOHR = .N, raw_also = sum(sig_raw), raw_share = round(mean(sig_raw), 3), rint_also = sum(sig_rint), rint_share = round(mean(sig_rint), 3),
                            normal_also = sum(sig_normal), normal_share = round(mean(sig_normal), 3)), by = .(set, class)][order(set, class)]
show("C_class_overlap", Cc)
show("C_overlap", X[, .(BOHR = sum(sig_BOHR), raw = sum(sig_raw), both = sum(sig_BOHR & sig_raw), BOHR_only = sum(sig_BOHR & !sig_raw), raw_only = sum(sig_raw & !sig_BOHR),
                        union = sum(sig_BOHR | sig_raw), rint = sum(sig_rint), BOHR_and_rint = sum(sig_BOHR & sig_rint), BOHR_not_rint = sum(sig_BOHR & !sig_rint),
                        rint_not_BOHR = sum(sig_rint & !sig_BOHR)), by = set])
## raw-only discoveries: where is BOHR?
bq <- merge(X[sig_raw & !sig_BOHR, .(set, tissue, unit, p_BOHR, s0, s1, med_D)], D[arm == "real" & method == "BOHR", .(set, tissue, unit, q)], by = c("set", "tissue", "unit"), all.x = TRUE)
show("C_raw_only_where_is_BOHR", bq[, .(raw_only = .N, BOHR_untested = sum(is.na(q)), BOHR_q_le_0.10 = sum(q <= 0.10, na.rm = TRUE), BOHR_q_le_0.20 = sum(q <= 0.20, na.rm = TRUE),
                                         BOHR_q_gt_0.5 = sum(q > 0.5, na.rm = TRUE), boundary_share_mean = round(mean(s0 + s1, na.rm = TRUE), 3)), by = set])

## D permuted discoveries by the largest share at one value, and by tissue
D[, sbin := cut(smax, c(-Inf, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, Inf), right = FALSE, labels = c("<30%", "30-50%", "50-60%", "60-70%", "70-80%", "80-90%", ">=90%"))]
Dd <- D[, .(tests = sum(arm == "real"), real = sum(sig[arm == "real"]), perm_bh = round(sum(sig[arm != "real"]) / 3, 2), fp = round(sum(hit[arm != "real"]) / 3, 1)), by = .(set, method, sbin)]
Dd[, rFDR := round(fp / real, 4)]
show("D_perm_by_share", dcast(Dd, set + sbin ~ method, value.var = c("tests", "real", "perm_bh", "rFDR"))[order(set, sbin)][, c("tests_normal", "tests_raw", "tests_rint") := NULL])
ns <- unique(U[, .(version, tissue, n_samples)])[, set := fifelse(version == "paper", "paper", "r5e80")]
Dt <- merge(bt[set != "v4", .(real = sum(sig[arm == "real"]), tests = sum(tests[arm == "real"]), perm_bh = round(sum(sig[arm != "real"]) / 3, 2), fp = round(sum(hit[arm != "real"]) / 3, 1)),
               by = .(set, tissue, method)], ns[, .(set, tissue, n_samples)], by = c("set", "tissue"))
fwrite(Dt[order(set, n_samples, tissue, method)], file.path(OUT, "H_per_tissue_four_methods.tsv"), sep = "\t")
Dt[, nbin := cut(n_samples, c(0, 250, 450, Inf), labels = c("< 250", "250-449", ">= 450"))]
show("D_perm_by_tissue_size", Dt[, .(tissues = uniqueN(tissue), tests = sum(tests), real = sum(real), perm_bh = sum(perm_bh), fp = sum(fp)), by = .(set, method, nbin)][, rFDR := round(fp / real, 4)][order(set, method, nbin)])

## E discoveries / permuted discoveries, with exact Poisson limits on the permuted total over the three arms
E <- A[set != "v4", .(set, method, real, perm_events, perm_mean = round(perm_events / 3, 2))]
E[, `:=`(perm_lo = round(qgamma(0.025, perm_events) / 3, 2), perm_hi = round(qgamma(0.975, perm_events + 1) / 3, 2))]
E[, `:=`(ratio = round(real / (perm_events / 3)), ratio_lo = round(real / perm_hi), ratio_hi = round(real / pmax(perm_lo, 1e-9)), perm_pct_of_real = round(100 * (perm_events / 3) / real, 3))]
show("E_ratio", E[order(set, method)])

## F matched permutation FDR
mfdr <- function(pr, pp, a = 0.05) { pr <- sort(pr[!is.na(pr)]); if (!length(pr)) return(0L)
  V <- Reduce(`+`, lapply(pp, function(q) findInterval(pr, sort(q[!is.na(q)])))) / length(pp)
  ok <- which(V / seq_along(pr) <= a); if (length(ok)) max(ok) else 0L }
Fm <- D[set != "v4", { pr <- p[arm == "real"]; pp <- split(p[arm != "real"], arm[arm != "real"])
  .(bh = sum(sig[arm == "real"]), m03 = mfdr(pr, pp, 0.03), m05 = mfdr(pr, pp, 0.05), m06 = mfdr(pr, pp, 0.06)) }, by = .(set, tissue, method)]
fwrite(Fm, file.path(OUT, "F_matched_by_tissue.tsv"), sep = "\t")
Ft <- Fm[, .(bh = sum(bh), m03 = sum(m03), m05 = sum(m05), m06 = sum(m06)), by = .(set, method)]
Ft[, `:=`(BOHR_over_bh = round(bh[method == "BOHR"] / bh, 3), BOHR_over_m03 = round(m03[method == "BOHR"] / m03, 3), BOHR_over_m05 = round(m05[method == "BOHR"] / m05, 3),
          BOHR_over_m06 = round(m06[method == "BOHR"] / m06, 3)), by = set]
show("F_matched", Ft[order(set, method)])
Fw <- dcast(Fm, set + tissue ~ method, value.var = c("bh", "m05"))
show("F_matched_tissue_spread", Fw[, .(tissues = .N, BOHR_gt_raw_bh = sum(bh_BOHR > bh_raw), BOHR_gt_raw_m05 = sum(m05_BOHR > m05_raw), BOHR_gt_rint_bh = sum(bh_BOHR > bh_rint),
                                        BOHR_gt_rint_m05 = sum(m05_BOHR > m05_rint), median_ratio_raw_m05 = round(median(m05_BOHR / m05_raw), 3),
                                        median_ratio_rint_m05 = round(median(m05_BOHR / m05_rint), 3)), by = set])
cat("\nDONE\n")
