## compare_units.R: the earlier five-read-floor run against the no-floor run (v4); also writes the gene read depth per promoter.
## Earlier run (five-read floor) vs no floor (v4): every significant unit of each version is put in a category by what the other
## version did with the same event, and gets expression, depth and plausibility features. Aggregates only (one row per
## unit or per promoter; no sample-level values are written).
## Event keys (robust to exact mirror pairs, whose p0/p1 label and carrier promoter are conventions):
##   unit  = tissue | promoter | type, with exact-mirror units mapped to the pair's first promoter (p0 <-> p1 swapped
##           when the unit sits on the second promoter; mu stays mu)
##   gc    = tissue | gene | class (class = mu or boundary)
## Significant: q <= 0.05 & estimable & !dedup_drop. Tested: estimable & !dedup_drop & q not NA.
## Usage: Rscript compare_units.R <ncores>    (lab compute node; reads the freeze inputs and the private depth counts)
suppressPackageStartupMessages({ library(data.table); library(parallel) })
G <- "/projects/group/bliu/e1327989/GTEx"; W <- file.path(G, "nofloor_full")
OUT <- file.path(W, "consult/version_choice_20261118/main/out"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
PU <- file.path(G, "QTL_figure/rerun_v2/results/final_nogate_v3c/freeze")
VU <- file.path(W, "GTEx/QTL_figure/rerun_v4/results/final_nogate_v3c/freeze")
VPERM <- file.path(W, "core/v4/units_v3m")
PPERM <- PU
PF <- file.path(G, "QTL_figure/rerun_v2/work/final/freeze"); VF <- file.path(W, "GTEx/QTL_figure/rerun_v4/work/final/freeze")
PRIV <- file.path(W, "core/depth/private")
PINFO <- file.path(G, "v11/Promoter_info_filter_5.txt"); VINFO <- file.path(W, "GTEx/v11/Promoter_info_filter_5.txt")
PCOL <- file.path(G, "QTL_figure/rerun_v2/results/coloc_eqtl_v3c/final_clump_v3c/coloc_eqtl_signals.tsv.gz")
VCOL <- file.path(W, "GTEx/QTL_figure/rerun_v4/results/coloc_eqtl_v3m/final_clump_v3m/coloc_eqtl_signals.tsv.gz")
DEPTH <- file.path(W, "core/depth/v4/unit_tests.tsv.gz")
TPMF <- "/home/e1327989/project/GTEx/TPM_v11/TPM1_genes.txt"
NC <- as.integer(commandArgs(TRUE)[1]); if (is.na(NC)) NC <- 8
fr <- function(f, ...) fread(cmd = paste("zcat", shQuote(f)), ...)
tis <- sub("[.]real[.]units[.]rds$", "", list.files(VU, "[.]real[.]units[.]rds$")); stopifnot(length(tis) == 50)
TT <- Sys.getenv("TEST_TIS"); if (nzchar(TT)) { tis <- strsplit(TT, ",")[[1]]; OUT <- paste0(OUT, "_test"); dir.create(OUT, showWarnings = FALSE) }

## ---- units ----
cols <- c("tissue", "promoter", "gene", "type", "p", "q", "estimable", "dedup_drop", "leader", "leader_logp", "leader_est",
          "n", "n0", "n1", "n_mid", "n_prom_gene", "boundary_exact", "interior_exact", "mirror_partner")
rdu <- function(d, t, arm = "real") { x <- as.data.table(readRDS(file.path(d, sprintf("%s.%s.units.rds", t, arm)))); x[, tissue := t]
  if (!"v4_added" %in% names(x)) x[, v4_added := NA]; if (!"leader_se" %in% names(x)) x[, leader_se := NA_real_]
  x[, c(cols, "v4_added", "leader_se"), with = FALSE] }
P <- rbindlist(mclapply(tis, rdu, d = PU, mc.cores = NC)); V <- rbindlist(mclapply(tis, rdu, d = VU, mc.cores = NC))
for (x in list(P, V)) { x[, sig := q <= 0.05 & estimable %in% TRUE & dedup_drop %in% FALSE & !is.na(q)]
  x[, tested := estimable %in% TRUE & dedup_drop %in% FALSE & !is.na(q)]; x[, cls := fifelse(type == "mu", "mu", "boundary")] }
cat("significant units: paper", P[sig == TRUE, .N], " v4", V[sig == TRUE, .N], "\n")
if (!nzchar(TT)) stopifnot(P[sig == TRUE, .N] == 13811, V[sig == TRUE, .N] == 30366)
swap <- function(ty) fifelse(ty == "p0", "p1", fifelse(ty == "p1", "p0", ty))
evkey <- function(x) { pr <- x$promoter; ty <- x$type; mp <- x$mirror_partner
  ## mu units are one event only when the pair's interior values are exact mirrors (interior_exact);
  ## boundary units when the boundary values are (boundary_exact)
  ex <- fifelse(ty == "mu", x$interior_exact %in% TRUE, x$boundary_exact %in% TRUE)
  m <- ex & !is.na(mp) & mp != ""; sec <- m & mp < pr
  paste(x$tissue, fifelse(sec, mp, pr), fifelse(sec, swap(ty), ty), sep = "|") }
for (x in list(P, V)) { x[, ukey := evkey(x)]; x[, gckey := paste(tissue, gene, cls, sep = "|")]; x[, tgkey := paste(tissue, gene, sep = "|")] }

## ---- promoter sets and promoter-level activity (both versions), gene and promoter reads (v4 counts) ----
prom_stats <- function(t, fold) { td <- readRDS(file.path(fold, paste0(t, ".rds"))); a <- td$act
  data.table(tissue = t, promoter = colnames(a), mean_act = colMeans(a, na.rm = TRUE), f0 = colMeans(a == 0, na.rm = TRUE),
             f1 = colMeans(a == 1, na.rm = TRUE), fna = colMeans(is.na(a))) }
PS <- rbindlist(mclapply(tis, prom_stats, fold = PF, mc.cores = NC)); VS <- rbindlist(mclapply(tis, prom_stats, fold = VF, mc.cores = NC))
depth_stats <- function(t) { td <- readRDS(file.path(VF, paste0(t, ".rds"))); dp <- readRDS(file.path(PRIV, paste0(t, ".rds")))
  stopifnot(identical(as.character(dp$row_id), as.character(td$row_id)))
  pr <- intersect(colnames(td$act), colnames(dp$R)); y <- td$act[, pr, drop = FALSE]; R <- dp$R[, pr, drop = FALSE]
  Dp <- dp$D[, unlist(dp$gene_of[pr]), drop = FALSE]; ok <- !is.na(y); m <- ok & y > 0 & y < 1
  Dp[!ok] <- NA; Rpos <- R; Rpos[!ok | R == 0] <- NA
  data.table(tissue = t, promoter = pr, med_D = apply(Dp, 2, median, na.rm = TRUE), med_R_pos = apply(Rpos, 2, median, na.rm = TRUE),
             mid_share_r1_4 = colSums(m & R <= 4) / pmax(colSums(m), 1), n_mid_v4 = colSums(m)) }
DS <- rbindlist(mclapply(tis, depth_stats, mc.cores = NC, mc.preschedule = FALSE))
fwrite(PS, file.path(OUT, "promoter_stats_paper.tsv.gz"), sep = "\t"); fwrite(VS, file.path(OUT, "promoter_stats_v4.tsv.gz"), sep = "\t")
fwrite(DS, file.path(OUT, "promoter_depth_v4.tsv.gz"), sep = "\t")
## major / minor within gene (v4 activity, promoters of the gene in the tissue's v4 set)
info <- fread(VINFO)[, .(promoter = V1, gene_info = geneId, chr = seqnames, tss = start, strand)]
pinfo <- fread(PINFO)[, .(promoter = V1, gene_info = geneId, chr = seqnames, tss = start, strand)]
cat("promoter annotation: v4", nrow(info), " paper", nrow(pinfo), " paper not in v4", pinfo[!promoter %in% info$promoter, .N], "\n")
in_paper_annot <- pinfo$promoter
info <- rbind(info, pinfo[!promoter %in% info$promoter])
VS <- merge(VS, info[, .(promoter, gene_info)], by = "promoter", all.x = TRUE)
VS[, rank_in_gene := frank(-mean_act, ties.method = "min"), by = .(tissue, gene_info)][, n_in_set := .N, by = .(tissue, gene_info)]

## ---- TPM ----
key <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
TPM <- fread(TPMF)[, .(gene = Name, tk = key(SMTSD), tpm = TPM_mean)]
addtpm <- function(x) { x[, tk := key(tissue)]; x <- merge(x, TPM, by = c("gene", "tk"), all.x = TRUE); x[, tk := NULL]; x }

## ---- categories ----
PSIG <- P[sig == TRUE]; VSIG <- V[sig == TRUE]
inset <- function(x, S) paste(x$tissue, x$promoter) %in% paste(S$tissue, S$promoter)
cat_of <- function(x, other_all, other_sig, other_set, side) {
  ## x: significant units of one version; other_*: the other version's units (all tested rows), its significant units
  ## and its promoter-set table (tissue, promoter)
  ot <- other_all[tested == TRUE]
  c1 <- x$ukey %in% other_sig$ukey
  c2 <- !c1 & x$gckey %in% other_sig$gckey
  c3 <- !c1 & !c2 & x$tgkey %in% other_sig$tgkey
  c4 <- !c1 & !c2 & !c3 & x$ukey %in% ot$ukey
  c5 <- !c1 & !c2 & !c3 & !c4 & inset(x, other_set)
  out <- fifelse(c1, "1 same unit", fifelse(c2, "2 same gene x class", fifelse(c3, "3 class changed",
         fifelse(c4, "4 other version tested it, n.s.", fifelse(c5, "5 promoter in other set, unit not tested", "6 promoter not in other set")))))
  out }
VSIG[, category := cat_of(VSIG, P, PSIG, PS, "v4")]; PSIG[, category := cat_of(PSIG, V, VSIG, VS, "paper")]
## the other version's P of the same event (category 4) for pi1
VSIG <- merge(VSIG, P[tested == TRUE, .(ukey, p_other = p, q_other = q)][!duplicated(ukey)], by = "ukey", all.x = TRUE)
PSIG <- merge(PSIG, V[tested == TRUE, .(ukey, p_other = p, q_other = q)][!duplicated(ukey)], by = "ukey", all.x = TRUE)
## the class of the other version's significant unit for category 3
cc <- function(x, S) { s <- S[, .(other_cls = paste(sort(unique(type)), collapse = "+")), by = tgkey]; merge(x, s, by = "tgkey", all.x = TRUE) }
VSIG <- cc(VSIG, PSIG[, .(tgkey, type)]); PSIG <- cc(PSIG, VSIG[, .(tgkey, type)])

## ---- features ----
feat <- function(x, S, ver) {
  x <- merge(x, S[, .(tissue, promoter, mean_act, f0_set = f0, f1_set = f1, fna)], by = c("tissue", "promoter"), all.x = TRUE)
  x <- merge(x, DS, by = c("tissue", "promoter"), all.x = TRUE)
  x <- merge(x, VS[, .(tissue, promoter, mean_act_v4 = mean_act, rank_in_gene, n_in_set)], by = c("tissue", "promoter"), all.x = TRUE)
  x <- merge(x, info[, .(promoter, chr, tss, strand)], by = "promoter", all.x = TRUE)
  x[, pos := as.numeric(sub("^chr[^_]+_([0-9]+)_.*$", "\\1", leader))]
  x[, tss_dist := fifelse(strand == "-", tss - pos, pos - tss)]
  x[, nlog10p_leader := -leader_logp / log(10)]
  x <- addtpm(x); x[, version := ver]; x }
VSIG <- feat(VSIG, VS, "v4"); PSIG <- feat(PSIG, PS, "paper")

## depth test (v4 units before the p1 admission rule of units_v3m_mirror.R: 30,363 rows; its 3 changes are left NA)
thr <- V[sig == TRUE, .(thr = max(p)), by = tissue]
DT <- fr(DEPTH)[status == "ok"]; DT <- merge(DT, thr, by = "tissue")
DT[, depth_eval := p0 <= thr][, depth_stay := p1 <= thr & sign(b1) == sign(b0)][, strong_eqtl := eqtl_p < 1e-5]
DT <- merge(DT, V[, .(tissue, promoter, type, ukey)], by = c("tissue", "promoter", "type"))
cat("depth tests matched to v4 units:", nrow(DT), "\n")
VSIG <- merge(VSIG, unique(DT[, .(ukey, depth_eval, depth_stay, strong_eqtl, log10D_med)], by = "ukey"), by = "ukey", all.x = TRUE)

## eQTL colocalization of the unit's first signal (k = 1; classes of build_panel_e_signal_v3e.R)
DIST <- c("distinct (H3)", "puQTL-specific (H1)", "different locus (LD windows disjoint)")
## the lead signal is the first signal kept after the backward step (delivery v7 "lead"), not
## always k = 1 (43 v4 and 23 paper leads have k != 1); k = 1 only where the delivery has no lead row
leadk <- function(f) { d <- fr(f, select = c("tissue", "promoter", "type", "k", "lead"))[lead %in% TRUE]; unique(d[, .(tissue, promoter, type, k_lead = k)], by = c("tissue", "promoter", "type")) }
LK <- list(v4 = leadk(file.path(W, "GTEx/QTL_figure/rerun_v4/results/fig5_delivery_v7/signals.tsv.gz")),
           paper = leadk(file.path(G, "QTL_figure/rerun_v2/results/fig5_delivery_v7/signals.tsv.gz")))
colc <- function(f, U, lk) { x <- fr(f, select = c("tissue", "promoter", "type", "k", "egene", "call"))
  x <- merge(x, lk, by = c("tissue", "promoter", "type"), all.x = TRUE); x <- x[k == fifelse(is.na(k_lead), 1L, as.integer(k_lead))]
  eg <- toupper(as.character(x$egene)) == "TRUE"
  x[, eqtl_cls := fifelse(!eg, "no eQTL", fifelse(call == "shared (H4)", "shared", fifelse(call %in% DIST, "distinct", "unresolved")))]
  x <- merge(x, U[, .(tissue, promoter, type, ukey)], by = c("tissue", "promoter", "type"))
  unique(x[, .(ukey, eqtl_cls)], by = "ukey") }
VSIG <- merge(VSIG, colc(VCOL, V, LK$v4), by = "ukey", all.x = TRUE)
PSIG <- merge(PSIG, colc(PCOL, P, LK$paper), by = "ukey", all.x = TRUE)

## sharing: tissues in which the same event (promoter or its exact-mirror pair x type, tissue dropped) is significant
share <- function(x) { x[, ev := sub("^[^|]+[|]", "", ukey)]; x[, n_tis := uniqueN(tissue), by = ev]; x }
VSIG <- share(VSIG); PSIG <- share(PSIG)
VSIG[, why_not_in_paper := fifelse(category != "6 promoter not in other set", NA_character_,
  fifelse(v4_added %in% TRUE, "v4 rule (SD skip dropped)", fifelse(promoter %in% in_paper_annot, "in paper annotation, not in the tissue's set", "not in paper annotation")))]
## kept units: same leader variant and effect direction (non-mirror units only, so no label swap)
KP <- merge(VSIG[category == "1 same unit" & !(boundary_exact %in% TRUE), .(ukey, leader, leader_est, nlog10p_leader)],
            PSIG[!(boundary_exact %in% TRUE), .(ukey, leader_p = leader, est_p = leader_est, nlog10p_p = nlog10p_leader)], by = "ukey")
KC <- KP[, .(units = .N, same_leader = round(mean(leader == leader_p), 3), same_sign = round(mean(sign(leader_est) == sign(est_p)), 3),
             same_sign_same_leader = round(mean(sign(leader_est[leader == leader_p]) == sign(est_p[leader == leader_p])), 3),
             spearman_nlog10p = round(cor(nlog10p_leader, nlog10p_p, method = "spearman"), 3),
             median_nlog10p_v4 = round(median(nlog10p_leader), 1), median_nlog10p_paper = round(median(nlog10p_p), 1))]
fwrite(KC, file.path(OUT, "kept_unit_concordance.tsv"), sep = "\t")
GN <- data.table(what = c("tissue x gene with a significant unit", "genes with a significant unit (any tissue)"),
  paper = c(uniqueN(PSIG$tgkey), uniqueN(PSIG$gene)), v4 = c(uniqueN(VSIG$tgkey), uniqueN(VSIG$gene)),
  both = c(length(intersect(PSIG$tgkey, VSIG$tgkey)), length(intersect(PSIG$gene, VSIG$gene))))
fwrite(GN, file.path(OUT, "gene_counts.tsv"), sep = "\t")
cat("duplicated event keys among significant units: v4", sum(duplicated(VSIG$ukey)), " paper", sum(duplicated(PSIG$ukey)), "\n")
fwrite(VSIG, file.path(OUT, "v4_sig_units_features.tsv.gz"), sep = "\t"); fwrite(PSIG, file.path(OUT, "paper_sig_units_features.tsv.gz"), sep = "\t")

## ---- realized FDR by paper status (v4 permutation arms, v3m units, real-arm threshold of each tissue) ----
pst <- function(x) { st <- fifelse(x$ukey %in% PSIG$ukey | x$gckey %in% PSIG$gckey, "A paper significant (same unit or gene x class)",
  fifelse(x$ukey %in% P[tested == TRUE]$ukey, "B paper tested, n.s.", fifelse(inset(x, PS), "C promoter in paper set, unit not tested",
  "D promoter not in paper set"))); st }
V[, pstat := pst(V)]
perm <- rbindlist(mclapply(tis, function(t) rbindlist(lapply(c("perm0", "perm1", "perm2"), function(a) {
  x <- rdu(VPERM, t, a); x[, cls := fifelse(type == "mu", "mu", "boundary")]; x[, ukey := evkey(x)]; x[, gckey := paste(tissue, gene, cls, sep = "|")]
  x <- merge(x, thr, by = "tissue"); x <- x[estimable %in% TRUE & dedup_drop %in% FALSE & !is.na(p) & p <= thr]
  if (!nrow(x)) return(NULL); x[, arm := a]; x[, pstat := pst(x)]; x[, .(tissue, arm, type, pstat)] })), mc.cores = NC))
RF <- merge(V[sig == TRUE, .(real = .N), by = .(pstat, cls)], perm[, .(perm_pass = .N / 3), by = .(pstat, cls = fifelse(type == "mu", "mu", "boundary"))],
            by = c("pstat", "cls"), all = TRUE)
RF[is.na(perm_pass), perm_pass := 0]; RF[, realized_fdr := round(perm_pass / real, 4)]
RF <- merge(RF, V[tested == TRUE, .(tested = .N), by = .(pstat, cls)], by = c("pstat", "cls"))
fwrite(RF, file.path(OUT, "realized_fdr_by_paper_status.tsv"), sep = "\t")

## same, paper arms, by v4 status (are the paper units the floor makes calibrated?)
PS_v <- function(x) fifelse(x$ukey %in% VSIG$ukey | x$gckey %in% VSIG$gckey, "A v4 significant (same unit or gene x class)",
  fifelse(x$ukey %in% V[tested == TRUE]$ukey, "B v4 tested, n.s.", fifelse(inset(x, VS), "C promoter in v4 set, unit not tested", "D promoter not in v4 set")))
thrP <- P[sig == TRUE, .(thr = max(p)), by = tissue]; P[, vstat := PS_v(P)]
permP <- rbindlist(mclapply(tis, function(t) rbindlist(lapply(c("perm0", "perm1", "perm2"), function(a) {
  f <- file.path(PPERM, sprintf("%s.%s.units.rds", t, a)); if (!file.exists(f)) return(NULL)
  x <- rdu(PPERM, t, a); x[, cls := fifelse(type == "mu", "mu", "boundary")]; x[, ukey := evkey(x)]; x[, gckey := paste(tissue, gene, cls, sep = "|")]
  x <- merge(x, thrP, by = "tissue"); x <- x[estimable %in% TRUE & dedup_drop %in% FALSE & !is.na(p) & p <= thr]
  if (!nrow(x)) return(NULL); x[, arm := a]; x[, vstat := PS_v(x)]; x[, .(tissue, arm, type, vstat)] })), mc.cores = NC))
RFP <- merge(P[sig == TRUE, .(real = .N), by = .(vstat, cls)], permP[, .(perm_pass = .N / uniqueN(permP$arm)), by = .(vstat, cls = fifelse(type == "mu", "mu", "boundary"))],
             by = c("vstat", "cls"), all = TRUE)
RFP[is.na(perm_pass), perm_pass := 0]; RFP[, realized_fdr := round(perm_pass / real, 4)]
fwrite(RFP, file.path(OUT, "paper_realized_fdr_by_v4_status.tsv"), sep = "\t")

## ---- pi1 (Storey, lambda = 0.5) of the other version's P for category-4 units, with a background ----
pi1 <- function(p) { p <- p[is.finite(p)]; if (length(p) < 50) return(NA_real_); round(1 - min(1, mean(p > 0.5) / 0.5), 3) }
bgV <- merge(V[tested == TRUE & sig == FALSE, .(ukey)], P[tested == TRUE & sig == FALSE, .(ukey, p)], by = "ukey")
bgP <- merge(P[tested == TRUE & sig == FALSE, .(ukey)], V[tested == TRUE & sig == FALSE, .(ukey, p)], by = "ukey")
PI <- rbind(
  VSIG[category == "4 other version tested it, n.s.", .(set = "v4 significant, paper tested n.s.: paper P", units = .N, pi1 = pi1(p_other),
      p_lt_0.05 = round(mean(p_other < 0.05), 3), p_lt_0.01 = round(mean(p_other < 0.01), 3)), by = cls],
  bgV[, .(set = "background: n.s. in both versions: paper P", cls = "all", units = .N, pi1 = pi1(p), p_lt_0.05 = round(mean(p < 0.05), 3), p_lt_0.01 = round(mean(p < 0.01), 3))],
  PSIG[category == "4 other version tested it, n.s.", .(set = "paper significant, v4 tested n.s.: v4 P", units = .N, pi1 = pi1(p_other),
      p_lt_0.05 = round(mean(p_other < 0.05), 3), p_lt_0.01 = round(mean(p_other < 0.01), 3)), by = cls],
  bgP[, .(set = "background: n.s. in both versions: v4 P", cls = "all", units = .N, pi1 = pi1(p), p_lt_0.05 = round(mean(p < 0.05), 3), p_lt_0.01 = round(mean(p < 0.01), 3))],
  fill = TRUE)
fwrite(PI, file.path(OUT, "pi1.tsv"), sep = "\t")

## ---- summaries by category ----
summ <- function(x) x[, .(units = .N, mu_share = round(mean(cls == "mu"), 3),
  tpm_median = round(median(tpm, na.rm = TRUE), 1), tpm_lt5 = round(mean(tpm < 5, na.rm = TRUE), 3), tpm_lt10 = round(mean(tpm < 10, na.rm = TRUE), 3),
  gene_reads_median = round(median(med_D, na.rm = TRUE), 0), gene_reads_lt20 = round(mean(med_D < 20, na.rm = TRUE), 3),
  prom_mean_act_v4 = round(median(mean_act_v4, na.rm = TRUE), 3), minor_promoter = round(mean(rank_in_gene > 1, na.rm = TRUE), 3),
  prom_reads_pos_median = round(median(med_R_pos, na.rm = TRUE), 1), mid_r1_4_median = round(median(mid_share_r1_4, na.rm = TRUE), 3),
  zero_share_median = round(median(n0 / n, na.rm = TRUE), 3), nlog10p_median = round(median(nlog10p_leader, na.rm = TRUE), 1),
  q_median = signif(median(q, na.rm = TRUE), 2), tss_10kb = round(mean(abs(tss_dist) <= 1e4, na.rm = TRUE), 3),
  tss_dist_kb_median = round(median(abs(tss_dist), na.rm = TRUE) / 1e3, 1),
  n_tis_median = as.numeric(median(n_tis)), one_tissue = round(mean(n_tis == 1), 3),
  eqtl_shared = round(mean(eqtl_cls == "shared", na.rm = TRUE), 3), eqtl_distinct_or_none = round(mean((eqtl_cls %in% c("distinct", "no eQTL"))[!is.na(eqtl_cls)]), 3),
  coloc_missing = sum(is.na(eqtl_cls))), keyby = .(version, category)]
SV <- summ(VSIG); SP <- summ(PSIG)
SD <- VSIG[depth_eval %in% TRUE, .(depth_units = .N, depth_stay = round(mean(depth_stay, na.rm = TRUE), 3), strong_gene_eqtl = round(mean(strong_eqtl, na.rm = TRUE), 3)), keyby = .(category, cls)]
SC <- rbind(SV, SP, fill = TRUE)
fwrite(SC, file.path(OUT, "summary_by_category.tsv"), sep = "\t"); fwrite(SD, file.path(OUT, "depth_by_category.tsv"), sep = "\t")
CNT <- rbind(VSIG[, .N, keyby = .(version, category, type)], PSIG[, .N, keyby = .(version, category, type)])
fwrite(CNT, file.path(OUT, "counts_by_category.tsv"), sep = "\t")
SW <- rbind(VSIG[category == "3 class changed", .N, keyby = .(version, from_other = other_cls, to = type)],
            PSIG[category == "3 class changed", .N, keyby = .(version, from_other = other_cls, to = type)])
fwrite(SW, file.path(OUT, "class_changes.tsv"), sep = "\t")
## category 3 is gene-level; how many paper boundary units have a significant v4 mu unit at the same promoter
## or at its mirror partner (any mirror flag)?
vmu <- unique(VSIG[type == "mu", paste(tissue, promoter)])
CC3 <- PSIG[category == "3 class changed" & type != "mu", .(units = .N,
  same_promoter = sum(paste(tissue, promoter) %in% vmu),
  promoter_or_partner = sum(paste(tissue, promoter) %in% vmu | (!is.na(mirror_partner) & paste(tissue, mirror_partner) %in% vmu))), by = type]
fwrite(CC3, file.path(OUT, "class_change_same_promoter.tsv"), sep = "\t")
## low expression: v4 units by gene TPM bin and category; depth stay by TPM bin
tb <- function(x) cut(x, c(1, 2, 5, 10, 50, Inf), labels = c("1-2", "2-5", "5-10", "10-50", ">=50"), right = FALSE)
TB <- rbind(VSIG[, .N, keyby = .(version, category, tpm_bin = tb(tpm))], PSIG[, .N, keyby = .(version, category, tpm_bin = tb(tpm))])
TB[, share := round(N / sum(N), 3), by = .(version, category)]
fwrite(TB, file.path(OUT, "tpm_bins_by_category.tsv"), sep = "\t")
DB <- VSIG[depth_eval %in% TRUE, .(units = .N, depth_stay = round(mean(depth_stay, na.rm = TRUE), 3)), keyby = .(tpm_bin = tb(tpm), cls)]
fwrite(DB, file.path(OUT, "depth_stay_by_tpm.tsv"), sep = "\t")
sink(file.path(OUT, "summary.txt"), split = TRUE)
options(width = 250)
cat("significant units: paper", nrow(PSIG), " v4", nrow(VSIG), "\n\nCOUNTS\n"); print(dcast(CNT, version + category ~ type, value.var = "N", fill = 0))
cat("\nCLASS CHANGES (category 3)\n"); print(SW)
cat("\nPAPER BOUNDARY UNITS IN CATEGORY 3: v4 mu unit at the same promoter or its partner\n"); print(CC3)
cat("\nSUMMARY BY CATEGORY\n"); print(SC)
cat("\nDEPTH (v4 units at the tissue threshold; stay = depth-adjusted P still at threshold, same sign)\n"); print(SD)
cat("\nDEPTH STAY BY GENE TPM\n"); print(DB)
cat("\nTPM BINS\n"); print(dcast(TB, version + category ~ tpm_bin, value.var = "share", fill = 0))
cat("\nREALIZED FDR, v4 arms by paper status\n"); print(RF)
cat("\nREALIZED FDR, paper arms by v4 status\n"); print(RFP)
cat("\nPI1\n"); print(PI)
cat("\nKEPT UNITS: leader and direction\n"); print(KC)
cat("\nGENES\n"); print(GN)
cat("\nWHY v4 category-6 promoters were not in the paper set\n"); print(VSIG[category == "6 promoter not in other set", .N, keyby = .(why_not_in_paper, type)])
sink()
