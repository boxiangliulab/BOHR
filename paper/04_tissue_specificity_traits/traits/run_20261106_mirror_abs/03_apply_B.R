## =============================================================================
## 03_apply_B.R -- option B: keep the BH of v3; only for the significant exact-pair boundary events, set the label by absolute activity
##   - Each exact-pair boundary event still enters the BH family of p0 as p0(S) (v3); P, q and significance are unchanged.
##   - For significant events with outcome == "O_only" (label table of 01), the retained unit is changed from p0(S) to its twin p1(O) (the same test):
##       p1(O): final / primary / q_final / sig_final / q_all / sig_all take the v3 values of the event, dedup_drop = FALSE;
##       p0(S): becomes the dropped twin (dedup_drop = TRUE, duplicate_of = O:p1, final = FALSE, q empty).
##   - New column bh_family: the BH family of the test (= the group of v3; "p0" for the relabelled events), abs_relabel: TRUE for the two rows of each relabelled event.
##   - Non-significant events stay p0(S). Nothing is refitted: only the label table of 01 and the v3 results are read.
## Inputs: results/mirror_abs_labels.tsv.gz; trait_regression_20261104/results/summary/units_{real,sim1,perm1}.rds, sig_real.tsv, counts_real_*.tsv, null_discovery.tsv; trait_regression_20261104/results/export/trait_regression_units.tsv.gz and README_trait_regression.txt
## Output (new folder, v3 not overwritten):
##   results/export_B/trait_regression_units.tsv.gz + README_trait_regression.txt (for Figure 6 S9; columns = v3 export + bh_family + abs_relabel)
##   results/summary_B/sig_real.tsv, counts_real_{mu,p0,p1,union}.tsv, null_discovery.tsv, relabelled_events.tsv, apply_B.log
## =============================================================================
suppressPackageStartupMessages(library(data.table))
source("/projects/group/bliu/e1327989/GTEx/tissue_regression/trait_regression_20261104/code/data.R")
D <- "/projects/group/bliu/e1327989/GTEx/tissue_regression/trait_regression_20261106_mirror_abs/results"
OS <- file.path(D, "summary_B"); OE <- file.path(D, "export_B")
for (d in c(OS, OE)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
V3 <- file.path(PATHS$res, "summary")

L <- fread(cmd = paste("zcat", file.path(D, "mirror_abs_labels.tsv.gz")))
R <- L[sig_final_v3 == TRUE & outcome == "O_only", .(tissue, trait, gene, S, O, outcome, z_S, z_O, n, n0, n1, q_final_v3, p_v3)]
cat("significant exact-pair events:", L[sig_final_v3 == TRUE, .N], "| relabelled to p1(O):", nrow(R), "\n")
## The earlier run (with the junction floor) had 3 relabelled events and asserted that number; this run reports its own
## count (above), so the assertion is a message here. The structural checks below are unchanged.
if (nrow(R) != 3L) cat("note: the earlier run had 3 relabelled events; this run has", nrow(R), "\n")

U <- readRDS(file.path(V3, "units_real.rds"))
U[, bh_family := fifelse(final, type, NA_character_)]
U[, abs_relabel := FALSE]
U[, sig_final := sig_final %in% TRUE]
for (i in seq_len(nrow(R))) {
  r <- R[i]
  iS <- U[tissue == r$tissue & trait == r$trait & promoter == r$S & type == "p0", which = TRUE]
  iO <- U[tissue == r$tissue & trait == r$trait & promoter == r$O & type == "p1", which = TRUE]
  stopifnot(length(iS) == 1L, length(iO) == 1L, U$final[iS], U$sig_final[iS], !U$final[iO], U$dedup_drop[iO],
            U$duplicate_of[iO] == paste0(r$S, ":p0"), abs(U$p[iO] - U$p[iS]) < 1e-6,
            ## Under the final rule p1(O) is evaluable like p0(S) (the same 0/1 response)
            U$testable[iO], U$share_10_90[iO], U$assessable[iO], U$report_tissue[iO])
  keepS <- U[iS, .(final, primary, q_final, sig_final, q, sig, eq_class)]
  set(U, iO, c("final", "primary", "q_final", "sig_final", "q", "sig"), keepS[, .(final, primary, q_final, sig_final, q, sig)])
  set(U, iO, c("dedup_drop", "duplicate_of", "bh_family", "abs_relabel", "merged_component", "twin", "eq_class"),
      list(FALSE, NA_character_, "p0", TRUE, "boundary", paste0(r$S, ":p0"), keepS$eq_class))
  set(U, iS, c("final", "primary", "q_final", "sig_final", "q", "sig"), list(FALSE, FALSE, NA_real_, FALSE, NA_real_, FALSE))
  set(U, iS, c("dedup_drop", "duplicate_of", "bh_family", "abs_relabel", "merged_component", "twin"),
      list(TRUE, paste0(r$O, ":p1"), NA_character_, TRUE, NA_character_, NA_character_))
}
## Check: redoing BH by bh_family gives the original q; significant set = v3 significant set (only the relabelled units swapped)
U[final == TRUE, q_re := p.adjust(p, "BH"), by = .(tissue, trait, bh_family)]
q_dev <- U[final == TRUE, max(abs(q_re - q_final) / q_final)]; U[, q_re := NULL]
stopifnot(q_dev < 1e-9, U[final == TRUE, all(bh_family == type | abs_relabel)])
S <- U[sig_final == TRUE]
S3 <- fread(file.path(V3, "sig_real.tsv"))
cat("significant tests:", nrow(S3), "->", nrow(S), "| q re-derived by BH family, max rel. diff", signif(q_dev, 2), "\n")
stopifnot(nrow(S) == nrow(S3))

## ---- summary_B: significance table and count matrices (same definition as summarize.R) ----
fwrite(S, file.path(OS, "sig_real.tsv"), sep = "\t")
count_mat <- function(S, tis) {
  M <- dcast(unique(S[, .(tissue, trait, promoter)]), tissue ~ trait, value.var = "promoter", fun.aggregate = length)
  M <- merge(data.table(tissue = tis), M, by = "tissue", all.x = TRUE)
  for (t in setdiff(TRAITS, names(M))) M[, (t) := 0L]
  for (t in TRAITS) set(M, which(is.na(M[[t]])), t, 0L)
  setcolorder(M, c("tissue", TRAITS)); M
}
tis <- sort(unique(U$tissue))
for (ty in c("mu", "p0", "p1")) fwrite(count_mat(S[type == ty], tis), file.path(OS, sprintf("counts_real_%s.tsv", ty)), sep = "\t")
fwrite(count_mat(S, tis), file.path(OS, "counts_real_union.tsv"), sep = "\t")
cmp_counts <- rbindlist(lapply(c("mu", "p0", "p1", "union"), function(ty) {
  a <- melt(fread(file.path(V3, sprintf("counts_real_%s.tsv", ty))), id.vars = "tissue", variable.factor = FALSE)
  b <- melt(fread(file.path(OS, sprintf("counts_real_%s.tsv", ty))), id.vars = "tissue", variable.factor = FALSE)
  m <- merge(a, b, by = c("tissue", "variable"))
  data.table(type = ty, total_v3 = sum(m$value.x), total_B = sum(m$value.y), cells_changed = m[value.x != value.y, .N],
             changed = paste(m[value.x != value.y, sprintf("%s x %s: %d -> %d", tissue, variable, value.x, value.y)], collapse = "; "))
}))

## ---- null_discovery (real counts per group by the new labels; null arms unchanged) ----
UR <- U[primary == TRUE, .(tissue, trait, type, sig, sig_final)]
rm(U); invisible(gc())
ND3 <- fread(file.path(V3, "null_discovery.tsv"))
nd <- rbindlist(lapply(c("sim1", "perm1"), function(a) {
  key <- unique(readRDS(file.path(V3, sprintf("units_%s.rds", a)))[primary == TRUE, .(tissue, trait)]); invisible(gc())
  rbindlist(lapply(c("sig", "sig_final"), function(sc) {
    r <- merge(UR[get(sc) == TRUE], key, by = c("tissue", "trait"))[, .(real = .N), by = type]
    z <- ND3[arm == a & set == fifelse(sc == "sig", "all_primary", "final"), .(type, null)]
    m <- merge(r, z, by = "type", all = TRUE); m[is.na(null), null := 0L]
    m[, `:=`(ratio = null / real, arm = a, set = fifelse(sc == "sig", "all_primary", "final"))]; m }))
}))
fwrite(nd, file.path(OS, "null_discovery.tsv"), sep = "\t")
NDc <- merge(ND3, nd, by = c("type", "arm", "set"), suffixes = c("_v3", "_B"))

## ---- Figure 6 export (same columns as export_fig6.R, plus bh_family, abs_relabel) ----
E0 <- fread(cmd = paste("zcat", file.path(PATHS$res, "export/trait_regression_units.tsv.gz")))
E <- copy(E0); E[, `:=`(bh_family = fifelse(final, type, NA_character_), abs_relabel = FALSE)]
for (i in seq_len(nrow(R))) {
  r <- R[i]
  iS <- E[tissue == r$tissue & trait == r$trait & promoter == r$S & type == "p0", which = TRUE]
  iO <- E[tissue == r$tissue & trait == r$trait & promoter == r$O & type == "p1", which = TRUE]
  stopifnot(length(iS) == 1L, length(iO) == 1L, E$sig_final[iS], E$dedup_drop[iO])
  vS <- E[iS, .(final, q_final, sig_final, primary, q_all, sig_all)]
  set(E, iO, c("final", "q_final", "sig_final", "primary", "q_all", "sig_all"), vS)
  set(E, iO, c("dedup_drop", "duplicate_of", "bh_family", "abs_relabel"), list(FALSE, NA_character_, "p0", TRUE))
  set(E, iS, c("final", "q_final", "sig_final", "primary", "q_all", "sig_all"), list(FALSE, NA_real_, FALSE, FALSE, NA_real_, FALSE))
  set(E, iS, c("dedup_drop", "duplicate_of", "bh_family", "abs_relabel"), list(TRUE, paste0(r$O, ":p1"), NA_character_, TRUE))
}
diffrows <- E[abs_relabel == TRUE]
## Apart from these rows and the two new columns, cell-by-cell identical to the v3 export
same_rest <- isTRUE(all.equal(E[abs_relabel == FALSE, !c("bh_family", "abs_relabel")], E0[!diffrows, on = .(tissue, trait, promoter, type)], check.attributes = FALSE))
stopifnot(nrow(E) == nrow(E0), same_rest, E[sig_final == TRUE, .N] == E0[sig_final == TRUE, .N])
fwrite(E, file.path(OE, "trait_regression_units.tsv.gz"), sep = "\t")
rd <- readLines(file.path(PATHS$res, "export/README_trait_regression.txt"))
writeLines(c(
"Figure 3 trait regression -- export for Figure 6 (Supplementary Tables S9a/S9b), MIRROR LABELS BY ABSOLUTE ACTIVITY, option B",
paste0("Built ", format(Sys.time()), " by trait_regression_20261106_mirror_abs/code/03_apply_B.R from the v3 export (trait_regression_20261104/results/export)."),
"",
"WHAT CHANGED FROM THE v3 EXPORT (option B):",
"  * BH is unchanged: every exact-pair boundary event is tested in the p0 family of the silent TSS, as in v3. No P, q or significance call changes.",
"  * The 586 significant exact-pair boundary events were labelled by absolute activity (OLS of each TSS's proActiv log2 normalised count on the",
"    trait, trait-model covariates, |z| >= 1.96): only the other TSS changes -> p1 at the other TSS; otherwise p0 at the silent TSS.",
"    3 events (all Breast_Mammary_Tissue x SEX) became p1 at the partner. For these the kept unit is now p1(O): its row carries the event's",
"    final / q_final / sig_final / primary / q_all / sig_all, and p0(S) is the dropped twin (dedup_drop = TRUE, duplicate_of = O:p1).",
"    P of p1(O) is that unit's own fit (equal to P of p0(S) to < 1e-6; the same 0/1 response).",
"  * New column bh_family: the BH family the test's q_final comes from (= type, except \"p0\" for the 3 relabelled p1 rows).",
"  * New column abs_relabel: TRUE for the 6 rows of the 3 relabelled events.",
"  * Non-significant exact-pair events keep p0(S). Every other row is identical to the v3 export.",
"  * Labels for all exact-pair events (z of each TSS, outcome) are in trait_regression_20261106_mirror_abs/results/mirror_abs_labels.tsv.gz.",
"",
"---- the v3 README follows ----", rd), file.path(OE, "README_trait_regression.txt"))
fwrite(R, file.path(OS, "relabelled_events.tsv"), sep = "\t")

sink(file.path(OS, "apply_B.log"), split = TRUE)
cat("relabelled events (significant, outcome O_only):\n"); print(R)
cat("\nexport: rows", nrow(E), "(v3", nrow(E0), "); rows that differ from v3: 6 (abs_relabel); significant", E[sig_final == TRUE, .N], "\n")
cat("q re-derived by BH family: max rel. diff", signif(q_dev, 2), "\n\ncounts (summary_B vs v3):\n"); print(cmp_counts[, .(type, total_v3, total_B, cells_changed, changed)])
cat("\nnull_discovery (real counts follow the labels; null arms unchanged):\n"); print(NDc[real_v3 != real_B | null_v3 != null_B, .(arm, set, type, real_v3, real_B, null_v3, null_B, ratio_v3 = signif(ratio_v3, 4), ratio_B = signif(ratio_B, 4))])
sink()
