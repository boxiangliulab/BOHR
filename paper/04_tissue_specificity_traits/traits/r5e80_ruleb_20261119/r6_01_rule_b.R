#!/usr/bin/env Rscript
## Rule (b) for the trait tests of the final version (r5e80): rare binary tests are dropped, then BH is run again.
## It is rule b_drop_rare_binary of the comparison of FDR rules (results/20261119c_fdr_options/):
##   every test whose trait_class is "rare_binary" leaves the final set (all components). trait_class is the run's own label
##   (trait_regression_20261007/code/summarize.R:97-98): rare_binary when n_minor / n < 0.10, where n_minor is the smaller class among
##   the promoter's tested donors and n the number of tested donors; it is set per test (promoter), so a tissue x characteristic can
##   lose some tests and keep others. Continuous characteristics (AGE, BMI) are never rare_binary.
##   Then q = BH over the remaining final tests within the same families (real: tissue x trait x bh_family; nulls: tissue x trait x
##   type, as r5_01_units.R). Nothing is refitted.
## Input:  data_r5e80_20261119/ (r5_01_units.R): the r5e80 real units and the null arms' r5e80 final units.
## Output: data_r5e80_ruleb_20261119/ (promoter-level statistics and aggregates; no donor data):
##   trait_regression_units_r5e80_ruleb.tsv.gz  the r5e80 units with final / q_final / sig_final = rule (b); the r5e80 values kept as
##                                              final_r5e80 / q_final_r5e80 / sig_final_r5e80; drop_rare_binary flag
##   units_{sim1,perm1}_ruleb.rds               the null arms' rule-(b) final units (same columns as the r5e80 files)
##   counts_real_ruleb_{union,mu,p0,p1}.tsv, null_discovery_ruleb.tsv, summary_ruleb.tsv
##   methods_facts_ruleb.tsv, cells_ruleb.tsv, characteristics_ruleb.tsv   (the facts for a Methods sentence)
suppressMessages(library(data.table)); setDTthreads(2L)
F3  <- "/projects/group/bliu/e1327989/GTEx/tissue_regression/figure3_20261007"
R5  <- file.path(F3, "data_r5e80_20261119")
OUT <- file.path(F3, "data_r5e80_ruleb_20261119")
C19 <- file.path(F3, "results/20261119c_fdr_options")
stopifnot(!dir.exists(OUT)); dir.create(OUT)
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
wr  <- function(x, f) fwrite(x, file.path(OUT, f), sep = "\t")

## ---------------- real arm ----------------
say("reading r5e80 units")
U <- fread(cmd = paste("zcat", file.path(R5, "trait_regression_units_r5e80.tsv.gz")), na.strings = c("", "NA"))
say(sprintf("r5e80: %d units, %d final, %d significant", nrow(U), sum(U$final), sum(U$sig_final)))
stopifnot(U[trait %in% c("AGE", "BMI"), all(trait_class == "continuous")], U[final == TRUE, !anyNA(trait_class)])
setnames(U, c("final", "q_final", "sig_final"), c("final_r5e80", "q_final_r5e80", "sig_final_r5e80"))
U[, drop_rare_binary := trait_class %in% "rare_binary"]                                   # RULE (b)
U[, final := final_r5e80 & !drop_rare_binary]                                              # RULE (b)
U[, q_final := NA_real_]; U[final == TRUE, q_final := p.adjust(p, "BH"), by = .(tissue, trait, bh_family)]   # BH re-run
U[, sig_final := final & !is.na(q_final) & q_final <= 0.05]
say(sprintf("rule (b): %d final, %d significant", sum(U$final), sum(U$sig_final)))
fwrite(U, file.path(OUT, "trait_regression_units_r5e80_ruleb.tsv.gz"), sep = "\t", na = "NA")
say("wrote trait_regression_units_r5e80_ruleb.tsv.gz")

## ---- facts for Methods ----
CL <- c("mu", "p0", "p1")
F <- U[final_r5e80 == TRUE]
fx <- rbindlist(lapply(c(CL, "all"), function(ty) { x <- if (ty == "all") F else F[type == ty]
  data.table(component = ty, tests_r5e80 = nrow(x), tests_removed = sum(x$drop_rare_binary), tests_ruleb = sum(x$final),
             significant_r5e80 = sum(x$sig_final_r5e80), significant_removed = sum(x$sig_final_r5e80 & x$drop_rare_binary),
             significant_ruleb = sum(x$sig_final), significant_lost_by_BH_rerun = sum(x$sig_final_r5e80 & !x$drop_rare_binary & !x$sig_final),
             significant_gained_by_BH_rerun = sum(!x$sig_final_r5e80 & x$sig_final)) }))
wr(fx, "methods_facts_ruleb.tsv"); print(fx)
CE <- F[, .(tests = .N, rare = sum(drop_rare_binary), sig_r5e80 = sum(sig_final_r5e80), sig_ruleb = sum(sig_final)), by = .(tissue, trait)]
CE[, status := fifelse(rare == 0L, "kept", fifelse(rare == tests, "removed entirely", "partly removed"))]
wr(CE, "cells_ruleb.tsv")
cat("tissue x characteristic cells with r5e80 final tests:", nrow(CE), "\n"); print(CE[, .N, by = status])
cat("by component (cells with a final test of that component):\n")
print(F[, .(rare = sum(drop_rare_binary), tests = .N), by = .(tissue, trait, type)][
  , .(cells = .N, removed_entirely = sum(rare == tests), partly = sum(rare > 0 & rare < tests)), by = type])
CH <- F[, .(tissues_tested_r5e80 = uniqueN(tissue), tissues_tested_ruleb = uniqueN(tissue[final]), tests_r5e80 = .N, tests_ruleb = sum(final),
            significant_r5e80 = sum(sig_final_r5e80), significant_ruleb = sum(sig_final)), by = trait][order(trait)]
wr(CH, "characteristics_ruleb.tsv"); print(CH)
cat("characteristics with no test left in any tissue:", paste(CH[tests_ruleb == 0, trait], collapse = ", "), "\n")
cat("characteristics with no significant test left:", paste(CH[significant_ruleb == 0, trait], collapse = ", "), "\n")

## ---- counts matrices (as r5_01_units.R) ----
TR <- sort(unique(U$trait)); TS <- sort(unique(U$tissue))
cnt <- function(x, f) {
  M <- dcast(x[, .(n = uniqueN(promoter)), by = .(tissue, trait)], tissue ~ trait, value.var = "n", fill = 0L)
  M <- merge(data.table(tissue = TS), M, by = "tissue", all.x = TRUE)
  for (v in TR) { if (!v %in% names(M)) M[, (v) := 0L]; set(M, which(is.na(M[[v]])), v, 0L) }
  setcolorder(M, c("tissue", TR)); wr(M, f)
}
S <- U[sig_final == TRUE]
cnt(S, "counts_real_ruleb_union.tsv"); for (ty in CL) cnt(S[type == ty], sprintf("counts_real_ruleb_%s.tsv", ty))
REALSIG <- S[, .N, by = type]
SUM <- list(data.table(arm = "real", final_r5e80 = sum(U$final_r5e80), final = sum(U$final), sig_r5e80 = sum(U$sig_final_r5e80), sig = sum(U$sig_final)))
BH_REAL <- S[, .(tissue, trait, promoter, type)]
rm(U, F, S); invisible(gc())

## ---------------- null arms ----------------
ND <- list()
for (arm in c("sim1", "perm1")) {
  N <- readRDS(file.path(R5, sprintf("units_%s_r5e80.rds", arm)))
  n0 <- nrow(N); s0 <- sum(N$sig_final)
  N <- N[!trait_class %in% "rare_binary"]                                                   # RULE (b)
  N[, sig_final := p.adjust(p, "BH") <= 0.05, by = .(tissue, trait, type)]                 # BH re-run
  saveRDS(N, file.path(OUT, sprintf("units_%s_ruleb.rds", arm)))
  say(sprintf("%s: %d r5e80 final -> %d rule (b); significant %d -> %d", arm, n0, nrow(N), s0, sum(N$sig_final)))
  SUM[[length(SUM) + 1L]] <- data.table(arm = arm, final_r5e80 = n0, final = nrow(N), sig_r5e80 = s0, sig = sum(N$sig_final))
  ns <- N[sig_final == TRUE, .(null = .N), by = type]
  ND[[arm]] <- merge(REALSIG[, .(type, real = N)], ns, by = "type", all.x = TRUE)[, `:=`(arm = arm, set = "final_r5e80_ruleb")]
}
ND <- rbindlist(ND); ND[is.na(null), null := 0L]; ND[, ratio := null / real]
wr(ND[, .(type, real, null, ratio, arm, set)], "null_discovery_ruleb.tsv"); print(ND)
wr(rbindlist(SUM), "summary_ruleb.tsv"); print(rbindlist(SUM))

## ---- must equal rule b of the comparison of FDR rules (results/20261119c_fdr_options) ----
S4 <- fread(file.path(C19, "4a_summary_by_rule.tsv"))[rule == "b_drop_rare_binary"]
F3T <- fread(file.path(C19, "3_fdr_by_rule.tsv"))[rule == "b_drop_rare_binary" & cls == "all" & type %in% CL]
ok <- S4$mu == REALSIG[type == "mu", N] && S4$p0 == REALSIG[type == "p0", N] && S4$p1 == REALSIG[type == "p1", N] &&
  all(merge(F3T, dcast(ND, type ~ arm, value.var = "null"), by = "type")[, sim1.x == sim1.y & perm1.x == perm1.y])
say("equal to rule (b) of the FDR-rule comparison (real by component; null discoveries by component in both arms):", ok); stopifnot(ok)
for (f in list.files(OUT)) cat(sprintf("%s\t%s\n", f, tools::md5sum(file.path(OUT, f))))
say("DONE")
