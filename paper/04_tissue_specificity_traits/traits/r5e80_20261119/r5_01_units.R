#!/usr/bin/env Rscript
## The trait tests of the final version, r5e80 (no junction floor; gene median reads >= 5 for the association tests;
## endpoint skip at 80%), for Fig. 3C-F, Extended Data Fig. 3B and Supplementary Tables 3b-f.
## The rule is that of the earlier draft run (endpoint80_20261118/e80_01_traits.R, version nofloor_reads5_E80), unchanged:
##   E80:    a (tissue, trait, promoter) is skipped when p0_raw >= 0.80 or p1_raw >= 0.80 (shares of 0 / 1 before endpoint deletion,
##           as stored by run_trait.R); all its units go.
##   reads5: every (tissue, gene) with median first-intron junction reads < 5 (figures_v4/fig1c/gene_median_all57.tsv) is dropped.
##   BH:     q = BH within tissue x trait x bh_family over the kept final units; significant = q <= 0.05. Nothing is refitted.
## Here the rule is applied to every unit (the draft wrote only the example families and the significant units) and to the two
## null arms of the same run (sim1, perm1: their units carry their own p0_raw / p1_raw; BH within tissue x trait x type, as the
## run's summary, after checking that this reproduces their stored q).
## Outputs (figure3_20261007/data_r5e80_20261119/; promoter-level statistics, no donor data):
##   trait_regression_units_r5e80.tsv.gz  the no-floor export_B units with final / q_final / sig_final replaced by the r5e80 values
##                                        (stored values kept as final_stored / q_final_stored / sig_final_stored; drop flags added)
##   units_sim1_r5e80.rds, units_perm1_r5e80.rds  the null arms' r5e80 final units (tissue, trait, type, trait_class, p, sig_final)
##   counts_real_r5e80_{union,mu,p0,p1}.tsv  tissue x trait counts of significant promoters (as summary/counts_real_*.tsv)
##   null_discovery_r5e80.tsv, summary_r5e80.tsv
suppressMessages(library(data.table)); setDTthreads(2L)
NF    <- "/projects/group/bliu/e1327989/GTEx/nofloor_full/GTEx/tissue_regression"
EXPB  <- file.path(NF, "trait_regression_20261106_mirror_abs/results/export_B/trait_regression_units.tsv.gz")
REAL  <- file.path(NF, "trait_regression_20261007/results/real")
SUMM  <- file.path(NF, "trait_regression_20261104/results/summary")
GMED  <- "/projects/group/bliu/e1327989/GTEx/nofloor_full/figures_v4/fig1c/gene_median_all57.tsv"
DRAFT <- "/projects/group/bliu/e1327989/GTEx/nofloor_full/figures_v4/endpoint80/nofloor_reads5_E80/Figure3/data/traits"
OUT   <- "/projects/group/bliu/e1327989/GTEx/tissue_regression/figure3_20261007/data_r5e80_20261119"
E80 <- 0.80
stopifnot(!dir.exists(OUT)); dir.create(OUT)
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

## gene reads table, keyed as in e80_01_traits.R
G <- fread(GMED)
G[, tkey := tolower(gsub("[^A-Za-z0-9]", "", sub("_n[0-9]+[.]rds$", "", sub("^proActiv_", "", file))))]
G <- G[, .(tkey, gene, median_D)]
stopifnot(!anyDuplicated(G[, .(tkey, gene)]))
reads_drop <- function(U) {
  U[, tkey := tolower(gsub("[^A-Za-z0-9]", "", tissue))]
  stopifnot(length(setdiff(unique(U$tkey), G$tkey)) == 0L)
  U[G, median_D := i.median_D, on = .(tkey, gene)]
  U[, drop_reads := is.na(median_D) | median_D < 5]
  U[, c("tkey", "median_D") := NULL]
  invisible(U)
}

## ---------------- real arm ----------------
say("reading", EXPB)
U <- fread(cmd = paste("zcat", EXPB), na.strings = c("", "NA"))
say(sprintf("export_B: %d units, %d final, %d significant", nrow(U), sum(U$final), sum(U$sig_final, na.rm = TRUE)))
U[, q_chk := NA_real_]; U[final == TRUE, q_chk := p.adjust(p, "BH"), by = .(tissue, trait, bh_family)]
d <- U[final == TRUE, max(abs(q_chk - q_final))]; say(sprintf("stored BH reproduced: max |dq| = %.3g", d)); stopifnot(d < 1e-12)
U[, q_chk := NULL]
R <- rbindlist(lapply(unique(U$tissue), function(t) rbindlist(lapply(unique(U[tissue == t, trait]), function(tr) {
  x <- readRDS(file.path(REAL, t, paste0(tr, ".rds")))$res
  x[is.na(skip), .(tissue = t, trait = tr, promoter, p0_raw, p1_raw)] }))))
stopifnot(!anyDuplicated(R[, .(tissue, trait, promoter)]))
U[R, `:=`(p0_raw = i.p0_raw, p1_raw = i.p1_raw), on = .(tissue, trait, promoter)]
stopifnot(!anyNA(U$p0_raw), !anyNA(U$p1_raw)); rm(R)
U[, drop_endpoint := p0_raw >= E80 | p1_raw >= E80]
reads_drop(U)
U[, keep := !drop_endpoint & !drop_reads]
setnames(U, c("final", "q_final", "sig_final"), c("final_stored", "q_final_stored", "sig_final_stored"))
U[, final := final_stored & keep]
U[, q_final := NA_real_]; U[final == TRUE, q_final := p.adjust(p, "BH"), by = .(tissue, trait, bh_family)]
U[, sig_final := final & !is.na(q_final) & q_final <= 0.05]
say(sprintf("r5e80: %d final (stored final dropped: endpoint %d, reads %d, either %d), %d significant",
            sum(U$final), U[final_stored & drop_endpoint, .N], U[final_stored & drop_reads, .N], U[final_stored & !keep, .N], sum(U$sig_final)))
## must be the draft run's significant units (same rule, same inputs)
DS <- fread(cmd = paste("zcat", file.path(DRAFT, "units_sig_e80_all.tsv.gz")), select = c("tissue", "trait", "promoter", "type", "p", "q_e80"))
cmp <- merge(DS, U[sig_final == TRUE, .(tissue, trait, promoter, type, q_final)], by = c("tissue", "trait", "promoter", "type"), all = TRUE)
say(sprintf("draft significant units %d; here %d; in both %d; max |dq| %.3g", nrow(DS), sum(U$sig_final), cmp[!is.na(q_e80) & !is.na(q_final), .N],
            cmp[, max(abs(q_e80 - q_final), na.rm = TRUE)]))
stopifnot(nrow(cmp) == nrow(DS), nrow(DS) == sum(U$sig_final), cmp[, max(abs(q_e80 - q_final))] < 1e-12)
setcolorder(U, c(setdiff(names(U), c("final_stored", "q_final_stored", "sig_final_stored", "p0_raw", "p1_raw", "drop_endpoint", "drop_reads", "keep")),
                 "final_stored", "q_final_stored", "sig_final_stored", "p0_raw", "p1_raw", "drop_endpoint", "drop_reads", "keep"))
fwrite(U, file.path(OUT, "trait_regression_units_r5e80.tsv.gz"), sep = "\t", na = "NA")
say("wrote trait_regression_units_r5e80.tsv.gz")

## counts matrices (as the run's summary/counts_real_*.tsv: significant promoters per tissue x trait; all tissues of the export)
TR <- sort(unique(U$trait)); TS <- sort(unique(U$tissue))
cnt <- function(x, f) {
  M <- dcast(x[, .(n = uniqueN(promoter)), by = .(tissue, trait)], tissue ~ trait, value.var = "n", fill = 0L)
  M <- merge(data.table(tissue = TS), M, by = "tissue", all.x = TRUE)
  for (v in TR) { if (!v %in% names(M)) M[, (v) := 0L]; set(M, which(is.na(M[[v]])), v, 0L) }
  setcolorder(M, c("tissue", TR)); fwrite(M, file.path(OUT, f), sep = "\t")
}
S <- U[sig_final == TRUE]
cnt(S, "counts_real_r5e80_union.tsv")
for (ty in c("mu", "p0", "p1")) cnt(S[type == ty], sprintf("counts_real_r5e80_%s.tsv", ty))
REALSIG <- S[, .N, by = type]
SUM <- list(data.table(arm = "real", units = nrow(U), final_stored = sum(U$final_stored), final = sum(U$final),
                       sig_stored = sum(U$sig_final_stored, na.rm = TRUE), sig = sum(U$sig_final)))
rm(U, S); invisible(gc())

## ---------------- null arms ----------------
ND <- list()
for (arm in c("sim1", "perm1")) {
  say("reading null arm", arm)
  N <- readRDS(file.path(SUMM, sprintf("units_%s.rds", arm)))
  N <- N[, .(tissue, trait, promoter, gene, type, trait_class, p, final, q_final, sig_final, p0_raw, p1_raw)]
  invisible(gc())
  N[, q_chk := NA_real_]; N[final == TRUE, q_chk := p.adjust(p, "BH"), by = .(tissue, trait, type)]
  d <- N[final == TRUE, max(abs(q_chk - q_final))]; say(sprintf("%s stored BH reproduced: max |dq| = %.3g", arm, d)); stopifnot(d < 1e-12)
  stopifnot(N[final == TRUE, !anyNA(p0_raw) && !anyNA(p1_raw)])
  N[, drop_endpoint := p0_raw >= E80 | p1_raw >= E80]
  reads_drop(N)
  fs <- sum(N$final); ss <- sum(N$sig_final, na.rm = TRUE)
  N <- N[final == TRUE & !drop_endpoint & !drop_reads]
  N[, q_final := p.adjust(p, "BH"), by = .(tissue, trait, type)]
  N[, sig_final := q_final <= 0.05]
  N <- N[, .(tissue, trait, type, trait_class, p, sig_final)]
  saveRDS(N, file.path(OUT, sprintf("units_%s_r5e80.rds", arm)))
  say(sprintf("%s: %d final stored -> %d r5e80 final; significant %d -> %d", arm, fs, nrow(N), ss, sum(N$sig_final)))
  SUM[[length(SUM) + 1L]] <- data.table(arm = arm, units = NA_integer_, final_stored = fs, final = nrow(N), sig_stored = ss, sig = sum(N$sig_final))
  ns <- N[sig_final == TRUE, .(null = .N), by = type]
  ND[[arm]] <- merge(REALSIG[, .(type, real = N)], ns, by = "type", all.x = TRUE)[, `:=`(arm = arm, set = "final_r5e80")]
  rm(N); invisible(gc())
}
ND <- rbindlist(ND); ND[is.na(null), null := 0L]; ND[, ratio := null / real]
fwrite(ND[, .(type, real, null, ratio, arm, set)], file.path(OUT, "null_discovery_r5e80.tsv"), sep = "\t"); print(ND)
fwrite(rbindlist(SUM), file.path(OUT, "summary_r5e80.tsv"), sep = "\t"); print(rbindlist(SUM))
for (f in list.files(OUT)) cat(sprintf("%s\t%s\n", f, tools::md5sum(file.path(OUT, f))))
say("DONE")
