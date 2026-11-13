## units_v3b_from_labels.R (2026-11-06): unit step on the SAVED absolute-activity mirror labels, OPTION B (author 2026-11-06, one
## rule with Figure 3): every boundary event of an exact pair is TESTED as p0 on the silent TSS S (results/final_nogate_v3m, its
## +-500 kb window) and BH is unchanged; the rule label (code/mirror_abs_labels.R -> results/mirror_abs_labels/labels_real.tsv.gz,
## label_applied == TRUE: significant units only) only NAMES the event. No fitting, no BH here.
## For a significant event whose saved label is p1(O) (only O's absolute activity changes):
##   the p1(O) row becomes the primary unit and carries the event's TESTED statistics, copied from p0(S): p, q, leader, n_snp,
##   leader_* and *_at_leader (p0(S) per-SNP test == p1(O) per-SNP test: y_O = 1 - y_S); its own-window values are kept as
##   p_own_window / leader_own_window. The p0(S) row is dropped as its duplicate (duplicate_of = O:p1).
## Every other unit is the v3m unit unchanged (same p and q). Permuted arms: copied unchanged.
## Output: results/final_nogate_v3b/freeze/<tissue>.<arm>.units.rds ; results/final_nogate_v3b/checks.txt. Env LABELS / UNITS_OUT choose the
## label table and output (2026-11-06 refit: LABELS=results/mirror_abs_labels_v3c/labels_real.tsv.gz UNITS_OUT=results/final_nogate_v3c).
suppressPackageStartupMessages(library(data.table))
V2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"
IN <- file.path(V2, "results/final_nogate_v3m/freeze"); OUT <- file.path(V2, Sys.getenv("UNITS_OUT", "results/final_nogate_v3b"), "freeze"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
L <- fread(cmd = paste("zcat", file.path(V2, Sys.getenv("LABELS", "results/mirror_abs_labels/labels_real.tsv.gz"))))
stopifnot(all(L[changed == TRUE, label_applied & sig]))
CH <- L[changed == TRUE]; cat("saved labels:", nrow(L), "; applied (significant):", L[label_applied == TRUE, .N], "; changed:", nrow(CH), "\n")
res <- rbindlist(lapply(list.files(IN, pattern = "units[.]rds$", full.names = TRUE), function(f) {
  U <- as.data.table(readRDS(f)); U[, `:=`(q_v3m = q, dedup_drop_v3m = dedup_drop, p_own_window = p, leader_own_window = leader, tested_as = NA_character_)]
  U[, `:=`(abs_outcome = NA_character_, abs_z_S = NA_real_, abs_z_O = NA_real_, mirror_label = NA_character_, label_applied = NA)]
  if (U$arm[1] == "real") {
    l <- L[tissue == U$tissue[1]]
    for (i in seq_len(nrow(l))) { r <- l[i]
      U[promoter %in% c(r$S, r$O) & type %in% c("p0", "p1") & boundary_exact %in% TRUE,
        `:=`(abs_outcome = r$outcome, abs_z_S = r$z_S, abs_z_O = r$z_O, mirror_label = r$label_unit, label_applied = r$label_applied)] }
    ch <- CH[tissue == U$tissue[1]]
    cols <- intersect(c("p", "q", "leader", "n_snp", grep("^leader_|_at_leader$", names(U), value = TRUE)), names(U))
    cols <- setdiff(cols, "leader_own_window")
    for (i in seq_len(nrow(ch))) { r <- ch[i]
      iS <- U[promoter == r$S & type == "p0", which = TRUE]; iO <- U[promoter == r$O & type == "p1", which = TRUE]
      stopifnot(length(iS) == 1, length(iO) == 1, !U$dedup_drop[iS], U$dedup_drop[iO], U$q[iS] <= 0.05)
      for (cc in cols) set(U, iO, cc, U[[cc]][iS])
      set(U, iO, "dedup_drop", FALSE); set(U, iO, "duplicate_of", NA_character_); set(U, iO, "tested_as", paste0(r$S, ":p0"))
      set(U, iS, "dedup_drop", TRUE); set(U, iS, "duplicate_of", paste0(r$O, ":p1")); set(U, iS, "q", NA_real_) } }
  U[, mirror_rule := Sys.getenv("RULE_TAG", "v3b (abs-activity labels on significant events; tested as p0(S))")]
  tmp <- file.path(OUT, paste0(basename(f), ".tmp")); saveRDS(U, tmp); file.rename(tmp, file.path(OUT, basename(f)))
  U[, .(tissue, arm, promoter, type, p, q, q_v3m, dedup_drop, dedup_drop_v3m, estimable, tested_as)] }))
sink(file.path(dirname(OUT), "checks.txt"), split = TRUE)
cat("units", basename(dirname(OUT)), "(abs-activity mirror labels, option B), labels", Sys.getenv("LABELS", "results/mirror_abs_labels/labels_real.tsv.gz"), ",", format(Sys.time()), "\n")
R <- res[arm == "real"]
cat("\nunits whose primary flag changed: real", R[dedup_drop != dedup_drop_v3m, .N], "; permuted", res[arm != "real" & dedup_drop != dedup_drop_v3m, .N], "\n")
cat("q values: the set of primary q values per tissue is identical to v3m:",
    all(R[, .(a = list(sort(q[!dedup_drop & estimable])), b = list(sort(q_v3m[!dedup_drop_v3m & estimable]))), by = tissue][, mapply(identical, a, b)]), "\n")
cat("\ntested (primary, estimable) units by type, real arm:\n")
print(merge(R[dedup_drop_v3m == FALSE & estimable == TRUE, .(v3m = .N), by = type], R[dedup_drop == FALSE & estimable == TRUE, .(v3b = .N), by = type], by = "type"))
cat("\nsignificant units by type (q <= 0.05), real arm:\n")
S <- merge(R[q_v3m <= 0.05, .(v3m = .N), by = type], R[q <= 0.05, .(v3b = .N), by = type], by = "type", all = TRUE)
print(rbind(S, S[, .(type = "all", v3m = sum(v3m), v3b = sum(v3b))]))
cat("\nrelabelled events:\n"); print(R[!is.na(tested_as)])
sink(); cat("DONE\n")
