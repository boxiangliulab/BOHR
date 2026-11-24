## Mirror rule v3.
## Reason: in exact two-promoter switches, the absolute activity that changes is the minor TSS's; the major TSS
## stays. The zeros identify the switching promoter, so an exact boundary event is labelled p0 on the promoter that goes
## silent (whichever it is), instead of v2's "label on the canonical promoter" (which made most of them p1(canonical)).
## Rule (only the boundary choice changes; pair definitions = mirror table v2, results/mirror/README.txt):
##   boundary_exact pairs (both members, boundary-only pairs included): keep p0(A), p0(B); drop p1(A), p1(B)
##     (p1(A) == p0(B)); duplicate_of for dropped p1(X) = mirror_partner(X):p0
##   mu: unchanged (drop mu(partner) when interior_exact)
##   non-exact two-promoter genes, >= 3-promoter genes: unchanged. Each promoter keeps its own +-500 kb window.
## Estimability rules (units v3) unchanged: mu n_mid >= 80; p0/p1 minor class >= 10% (symmetric under the swap).
## BH within tissue on dedup_drop == FALSE & estimable (q), as before. Old primary q kept as q_mirror_v2.
## Input : results/final_nogate_v3/freeze/*.units.rds (untouched)
## Output: results/final_nogate_v3m/freeze/<tissue>.<arm>.units.rds ; results/mirror_v3/promoter_flags.tsv.gz ;
##         results/mirror_v3/checks.txt
suppressPackageStartupMessages(library(data.table))
V2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"
IN <- file.path(V2, "results/final_nogate_v3/freeze"); OUT <- file.path(V2, "results/final_nogate_v3m/freeze")
MV3 <- file.path(V2, "results/mirror_v3"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE); dir.create(MV3, showWarnings = FALSE)
## Units of the promoters added by the final promoter rules (the v4 delta) enter BH only for the types admitted by a
## pre-specified permutation calibration (diag_v4_calibration.R): mu and p0 admitted, p1 not. File from v4_strata_sizes.R.
V4ADD <- fread(file.path(V2, "work/v4_admission/delta_promoter_features.tsv"))[, paste(tissue, promoter)]; V4_NOT_ADMITTED <- "p1"   # delta (tissue, promoter) pairs; the type not admitted

## ---- mirror table v3: promoter flags with the dropped boundary member named ----
FL <- fread(cmd = paste("zcat", file.path(V2, "results/mirror/promoter_flags.tsv.gz")))
setnames(FL, c("drop_boundary", "dedup_drop"), c("drop_boundary_v2", "dedup_drop_v2"))
FL[, drop_p0 := FALSE]
FL[, drop_p1 := boundary_exact == TRUE]                      # v3: p1 dropped on BOTH members of a boundary_exact pair
FL[, dropped_boundary_member := fifelse(boundary_exact == TRUE, paste0(promoter, ":p1"), NA_character_)]
FL[, kept_as := fifelse(boundary_exact == TRUE, paste0(mirror_partner, ":p0"), NA_character_)]   # p1(X) == p0(partner)
fwrite(FL, file.path(MV3, "promoter_flags.tsv.gz"), sep = "\t")

lab <- c(mu = "mu", p0 = "p1", p1 = "p0")
files <- list.files(IN, pattern = "units[.]rds$", full.names = TRUE)
res <- rbindlist(lapply(files, function(f) {
  U <- as.data.table(readRDS(f))
  U[, `:=`(dedup_drop_v2 = dedup_drop, duplicate_of_v2 = duplicate_of, merged_component_v2 = merged_component, twin_v2 = twin,
           eq_class_v2 = eq_class, q_mirror_v2 = q, q_v2_mirror_v2 = q_v2)]
  be <- U$boundary_exact %in% TRUE
  U[, dedup_drop := fifelse(type == "mu", drop_mu, be & type == "p1")]
  U[, duplicate_of := fifelse(dedup_drop, paste0(mirror_partner, ":", lab[type]), NA_character_)]
  U[, merged_component := fifelse(mirror_role %in% "canonical" & type == "mu" & interior_exact %in% TRUE, "interior",
                         fifelse(type == "p0" & be, "boundary", NA_character_))]
  U[, twin := fifelse(!is.na(merged_component), paste0(mirror_partner, ":", lab[type]), NA_character_)]
  U[, eq_class := fifelse(!is.na(merged_component) | dedup_drop, paste0(mirror_pair_id, "|", fifelse(type == "mu", "interior",
                   fifelse((mirror_role == "canonical" & type == "p0") | (mirror_role == "partner" & type == "p1"), "canonical_off", "canonical_full"))), NA_character_)]
  U[, q_v2 := NA_real_]; U[dedup_drop == FALSE & estimable_mu_rule == TRUE, q_v2 := p.adjust(p, "BH")]
  U[, v4_added := paste(tissue, promoter) %in% V4ADD]; U[, v4_not_admitted := v4_added & type %in% V4_NOT_ADMITTED]   # delta units of a type that was not admitted
  U[v4_not_admitted == TRUE, estimable := FALSE]   # are not tested (they do not enter BH)
  U[, q := NA_real_]; U[dedup_drop == FALSE & estimable == TRUE, q := p.adjust(p, "BH")]
  U[, mirror_rule := "v3"]
  tmp <- file.path(OUT, paste0(basename(f), ".tmp")); saveRDS(U, tmp); file.rename(tmp, file.path(OUT, basename(f)))
  U[, .(tissue, arm, promoter, type, p, q, q_mirror_v2, dedup_drop, dedup_drop_v2, estimable, boundary_exact, interior_exact,
        mirror_role, mirror_partner, leader, n_snp)]
}))

## ---- checks (Figure 5's four + BH edge flips) ----
sink(file.path(MV3, "checks.txt"), split = TRUE)
cat("Mirror rule v3 checks,", format(Sys.time()), "\n\n")
R <- res[arm == "real"]
cat("1. surviving p1 units on boundary_exact pairs (all arms):", res[type == "p1" & boundary_exact %in% TRUE & dedup_drop == FALSE, .N], "\n")
M2 <- res[type == "mu"]
cat("2. mu units: dedup flag changed", M2[dedup_drop != dedup_drop_v2, .N], "(all arms); real-arm significance changed",
    M2[arm == "real" & xor(q <= 0.05 & !is.na(q), q_mirror_v2 <= 0.05 & !is.na(q_mirror_v2)), .N],
    "; max |log10 q change| (BH runs over all types in a tissue, so q moves slightly)", signif(M2[!is.na(q), max(abs(log10(q) - log10(q_mirror_v2)))], 3), "\n")
## tests per tissue
tc <- res[, .(tests_v2 = sum(!dedup_drop_v2 & estimable), tests_v3 = sum(!dedup_drop & estimable)), by = .(tissue, arm)]
cat("3. tissue x arm with a different number of tested (primary, estimable) units:", tc[tests_v2 != tests_v3, .N], "of", nrow(tc), "\n")
print(tc[, .(tests_v2 = sum(tests_v2), tests_v3 = sum(tests_v3)), by = arm])
## significant counts
sg <- function(d, qq) d[get(qq) <= 0.05, .N, by = type]
s2 <- sg(R, "q_mirror_v2"); s3 <- sg(R, "q"); S <- merge(s2, s3, by = "type", suffixes = c("_v2", "_v3"))
S <- rbind(S, S[, .(type = "all", N_v2 = sum(N_v2), N_v3 = sum(N_v3))]); cat("\n4. significant units, real arm (q <= 0.05):\n"); print(S)
moved <- R[type == "p1" & boundary_exact %in% TRUE & q_mirror_v2 <= 0.05]       # v2 p1(canonical) on exact pairs
cat("\n   v2 significant p1 on boundary_exact pairs:", nrow(moved), "\n")
p0v2 <- R[type == "p0" & q_mirror_v2 <= 0.05, .N]
newp0 <- R[type == "p0" & q <= 0.05 & dedup_drop_v2 == TRUE]                       # p0(partner): newly primary and significant
cat("   expected p0 v3 = v2 p0", p0v2, "+ moved", nrow(moved), "=", p0v2 + nrow(moved), "; observed", R[type == "p0" & q <= 0.05, .N], "\n")
## map moved: v2 p1(X) -> v3 p0(partner(X))
mp <- merge(moved[, .(tissue, promoter_v2 = promoter, unit_v2 = paste0(promoter, ":p1"), partner = mirror_partner, p_v2 = p, q_v2 = q_mirror_v2, leader_v2 = leader)],
            R[type == "p0", .(tissue, partner = promoter, p_v3 = p, q_v3 = q, leader_v3 = leader, n_snp_v3 = n_snp)], by = c("tissue", "partner"), all.x = TRUE)
mp[, unit_v3 := paste0(partner, ":p0")]
cat("   moved units with the twin p0 present:", mp[!is.na(p_v3), .N], "; twin significant in v3:", mp[q_v3 <= 0.05, .N],
    "; same leader:", mp[leader_v2 == leader_v3, .N], "; |log10 p diff| > 0.01:", mp[abs(log10(p_v2) - log10(p_v3)) > 0.01, .N], "\n")
np <- R[type == "p0" & dedup_drop_v2 == TRUE & q <= 0.05]
cat("   newly primary p0(partner) significant:", nrow(np), "= twin of a v2-significant p1:", np[paste(tissue, mirror_partner) %in% moved[, paste(tissue, promoter)], .N],
    "+ twin p1 NOT significant in v2 (window shift / BH):", np[!paste(tissue, mirror_partner) %in% moved[, paste(tissue, promoter)], .N], "\n")
cat("\nBH edge flips (real arm), units whose significance changed for reasons other than the label swap:\n")
oth <- R[!(type == "p1" & boundary_exact %in% TRUE) & !(type == "p0" & dedup_drop_v2 == TRUE)]
print(oth[, .(gained = sum(q <= 0.05 & !(q_mirror_v2 <= 0.05) , na.rm = TRUE), lost = sum(q_mirror_v2 <= 0.05 & !(q <= 0.05), na.rm = TRUE)), by = type])
thr <- R[, .(thr_v2 = suppressWarnings(max(p[q_mirror_v2 <= 0.05], na.rm = TRUE)), thr_v3 = suppressWarnings(max(p[q <= 0.05], na.rm = TRUE))), by = tissue]
cat("tissue thresholds: looser", thr[thr_v3 > thr_v2, .N], "stricter", thr[thr_v3 < thr_v2, .N], "same", thr[thr_v3 == thr_v2, .N], "\n")
## permuted arms: realized-FDR numerator with the real-arm thresholds
thr_r <- R[, .(thr = suppressWarnings(max(p[q <= 0.05], na.rm = TRUE)), thr2 = suppressWarnings(max(p[q_mirror_v2 <= 0.05], na.rm = TRUE))), by = tissue]
P <- merge(res[arm != "real"], thr_r, by = "tissue")
pf <- P[, .(perm_v3 = sum(!dedup_drop & estimable & p <= thr), perm_v2 = sum(!dedup_drop_v2 & estimable & p <= thr2)), by = arm]
cat("\nrealized FDR (permuted units passing the real-arm per-tissue threshold, mean of perm0-2 / real significant):\n")
print(pf); cat("v2:", round(mean(pf$perm_v2), 1), "/", S[type == "all", N_v2], "=", round(mean(pf$perm_v2) / S[type == "all", N_v2], 4),
             "; v3:", round(mean(pf$perm_v3), 1), "/", S[type == "all", N_v3], "=", round(mean(pf$perm_v3) / S[type == "all", N_v3], 4), "\n")
sink()
fwrite(mp[, .(tissue, unit_v2, unit_v3, promoter_v2, promoter_v3 = partner, p_v2, p_v3, q_v2, q_v3, leader_v2, leader_v3,
              same_leader = leader_v2 == leader_v3)], file.path(MV3, "unit_map_v2_to_v3.tsv.gz"), sep = "\t")
fwrite(res[arm == "real" & (q <= 0.05 | q_mirror_v2 <= 0.05)], file.path(MV3, "real_sig_v2_or_v3.tsv.gz"), sep = "\t")
cat("DONE\n")
