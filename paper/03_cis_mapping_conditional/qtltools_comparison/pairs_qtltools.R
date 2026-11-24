## pairs_qtltools.R: QTLtools on the promoter-variant pairs BOHR tested, with the samples of the original QTLtools run.
## Its rows of promoters without a missing value are the no-fill rows of every arm (fill_qtltools.R).
## QTLtools 1.2 cis --permute on the same promoter-variant pairs as BOHR, with the samples handled as in the original QTLtools
## run: all donors of the tissue, a missing value replaced by the promoter's mean (make_freeze_beds.R), nothing deleted.
##   variants   = exactly the variants BOHR tested for the promoter: each genotype class has at least 5 of BOHR's analysis
##                samples (missing values dropped, rare boundary values deleted: prepare_y), modal imputation of genotypes;
##   phenotypes = as the original inputs: raw (NA -> mean over the donors with a value); rint = rank-normal of that vector over
##                all donors, ties averaged; normal = QTLtools --normal on raw;
##   covariates, window (500 kb), permutations (1000), seed: as the original run.
## So the only difference from the original run is the variant list. Promoters are grouped by BOHR's analysis sample set (it
## fixes the variant list) within a gene; promoters BOHR analysed on all donors go by chromosome with the tissue-level site list.
## Reuse of identical inputs: rows of promoters analysed on all donors come from the exact run when it is finished (same samples,
## values and variants); rows of promoters without missing values where BOHR deleted samples come from the keep run when it is
## finished (its samples were all donors). Everything else, or everything if those runs are not there, is run here.
## Output: pairs/<version>/<mode>/<Tissue>.<arm>.promoters.tsv.gz: promoter, n_s, n_nonmiss, full, n_variants, top_variant, dof,
## p_nominal, p_empirical, p, n_used (samples QTLtools used = all donors), source (run | exact | keep). Group files live in a
## node-local folder and are removed. No genotype or donor-level value is written to the project.
## Usage: Rscript pairs_qtltools.R <v4|paper> <Tissue> <arm> <ncores>
##   PAIRS_MAXGENES=n  test on the first n genes;  PAIRS_NOREUSE=1  run every promoter here;  PAIRS_OUT=dir  write there instead
suppressPackageStartupMessages({ library(data.table); library(VariantAnnotation); library(parallel) })
a <- commandArgs(TRUE); ver <- a[1]; t <- a[2]; arm <- a[3]; NC <- as.integer(a[4]); MAXG <- as.integer(Sys.getenv("PAIRS_MAXGENES", "0")); NOREUSE <- Sys.getenv("PAIRS_NOREUSE", "0") == "1"
G <- "/projects/group/bliu/e1327989/GTEx"; W <- file.path(G, "nofloor_full")
V2 <- c(v4 = file.path(W, "GTEx/QTL_figure/rerun_v4"), paper = file.path(G, "QTL_figure/rerun_v2"))[[ver]]
INFO <- fread(c(v4 = file.path(W, "GTEx/v11/Promoter_info_filter_5.txt"), paper = file.path(G, "v11/Promoter_info_filter_5.txt"))[[ver]]); setnames(INFO, 1, "pid")
VCF <- file.path(G, "Rawdata/v9_vcf_MAF_0.05"); MV <- file.path(W, "consult/fig2e_20261119/matched_variants")
OUT <- Sys.getenv("PAIRS_OUT", file.path(MV, "pairs", ver)); QB <- "/home/e1327989/project/software/QTLtools_1.2_Ubuntu16.04_x86_64/QTLtools_1.2_Ubuntu16.04_x86_64"
HTS <- "/home/e1327989/.local/share/mamba/envs/htsenv/bin"; MIN_GENO <- 5L; FLANK <- 500000L
MODES <- c(raw = "raw", rint = "rint", normal = "raw")                     # mode -> BED kind
TMP <- file.path(Sys.getenv("EXACT_TMP", tempdir()), paste("pairs", ver, t, arm, Sys.getpid(), sep = ".")); dir.create(file.path(TMP, "out"), recursive = TRUE)
on.exit(unlink(TMP, recursive = TRUE), add = TRUE)
fname <- paste(t, arm, "promoters.tsv.gz", sep = "."); final <- file.path(OUT, names(MODES), fname)
if (MAXG == 0 && all(file.exists(final))) { cat("exists:", t, arm, "\n"); quit(save = "no") }
done <- function(dir, logs) { f <- file.path(MV, logs, paste(ver, t, arm, "exit", sep = ".")); all(file.exists(file.path(MV, dir, ver, names(MODES), fname))) && file.exists(f) && identical(trimws(readLines(f, n = 1)), "0") }
have_exact <- !NOREUSE && done("exact", "logs/run"); have_keep <- !NOREUSE && done("exact_keep", "logs/run_keep")

## ---- BOHR's own rules (prepare_y of bohr2_data.R; v4 and paper differ in the SD rule), as in exact_qtltools.R ----
prepare_y <- function(y, p_cov) { y <- y[!is.na(y)]; if (length(y) < 20 || length(y) < p_cov) return(list(skip = "n"))
  p0 <- mean(y == 0); p1 <- mean(y == 1)
  if (p0 >= 0.95 || p1 >= 0.95 || (if (ver == "v4") sd(y) < 0.005 else sd(y) <= 0.05)) return(list(skip = "prop_or_sd"))
  keep <- rep(TRUE, length(y)); if (p0 <= 0.05) keep <- keep & y > 0; if (p1 <= 0.05) keep <- keep & y < 1
  list(skip = NA_character_, y = y[keep]) }
gt_dosage <- function(gt) { D <- matrix(NA_integer_, nrow(gt), ncol(gt), dimnames = dimnames(gt))
  D[gt == "0|0" | gt == "0/0"] <- 0L; D[gt == "0|1" | gt == "1|0" | gt == "0/1" | gt == "1/0"] <- 1L; D[gt == "1|1" | gt == "1/1"] <- 2L; D }
impute_modal <- function(D) { for (i in which(rowSums(is.na(D)) > 0)) { tb <- tabulate(D[i, ] + 1L, nbins = 3L); D[i, is.na(D[i, ])] <- which.max(tb) - 1L }; D }

## ---- the arm's phenotypes and covariates (as run_final2.R / make_freeze_beds.R) ----
td <- readRDS(file.path(V2, "work/final/freeze", paste0(t, ".rds"))); U <- td$U; K <- td$K; nU <- length(U)
if (arm == "real") { A <- td$act[U, , drop = FALSE]; Pe <- td$peer[U, , drop = FALSE] } else {
  m <- fread(file.path(V2, "perm_maps/final/freeze", paste0(t, ".", arm, ".perm_map.tsv")))
  stopifnot(setequal(m$slot, U), setequal(m$donor, U), !anyDuplicated(m$donor))
  d <- setNames(m$donor, m$slot)[U]; A <- td$act[d, , drop = FALSE]; Pe <- td$peer[d, , drop = FALSE]; rownames(A) <- rownames(Pe) <- U }
C <- cbind(td$pcs[U, , drop = FALSE], Pe[, seq_len(K), drop = FALSE]); colnames(C) <- c(paste0("PC", 1:ncol(td$pcs)), if (K > 0) paste0("peer_", seq_len(K)))
pc <- ncol(C)

## ---- promoters: those of the units compared in Fig. 2E, with BOHR's analysis samples (they fix the variant list) ----
fr <- function(f, ...) fread(cmd = paste("zcat", shQuote(f)), ...)
b <- fr(file.path(V2, "results/fig2e_freeze_v3m/bohr", paste0(t, ".real.tsv.gz")), select = c("promoter", "promoters_in_unit"))
X <- fr(file.path(W, "consult/fig2e_20261119/main/out/units_paper_r5e80.tsv.gz"), select = c("set", "tissue", "unit"))
units <- X[tissue == t & set == (if (ver == "v4") "r5e80" else "paper"), unit]
pr <- unique(unlist(strsplit(b[promoter %in% units, promoters_in_unit], ";", fixed = TRUE))); pr <- intersect(pr, colnames(A))
pr <- pr[td$info[pr, "seqnames"] %in% paste0("chr", 1:22)]
P <- rbindlist(lapply(pr, function(p) { y <- A[, p]; names(y) <- U; z <- prepare_y(y, pc)
  if (!is.na(z$skip) || !any(z$y > 0 & z$y < 1)) return(NULL)
  data.table(promoter = p, n_s = length(z$y), n_nonmiss = sum(!is.na(y)), skey = paste(match(names(z$y), U), collapse = ",")) }))
P[, `:=`(chr = td$info[promoter, "seqnames"], tss = as.integer(td$info[promoter, "start"]), gene = td$info[promoter, "geneId"], full = n_s == nU)]
P[, source := fifelse(full & have_exact, "exact", fifelse(!full & n_nonmiss == nU & have_keep, "keep", "run"))]
R <- P[source == "run"]
if (MAXG > 0) R <- R[gene %in% head(unique(R[full == FALSE, gene]), MAXG) | (full == TRUE & chr == "chr22")]
R[, grp := fifelse(full, paste0("full.", chr), paste0("g", as.integer(factor(paste(skey, chr, gene)))))]
Gt <- R[, .(chr = chr[1], from = min(tss), to = max(tss), n = .N, full = full[1], skey = skey[1], n_s = n_s[1]), by = grp]
cat(t, arm, ver, ": donors", nU, "| covariates", pc, "| promoters", nrow(P), "| from the exact run", sum(P$source == "exact"), "| from the keep run", sum(P$source == "keep"),
    "| run here", nrow(R), "(with missing values", R[n_nonmiss < nU, .N], ") | groups", nrow(Gt), "\n")

## ---- BED (raw, rint) as the original inputs (make_freeze_beds.R): all donors, NA -> the promoter's mean; rint over all donors ----
Rw <- matrix(0, nrow(R), nU, dimnames = list(R$promoter, U)); Rn <- Rw
for (i in seq_len(nrow(R))) { y <- A[, R$promoter[i]]; y[is.na(y)] <- mean(y, na.rm = TRUE); stopifnot(all(is.finite(y)), length(unique(y)) > 1)
  Rw[i, ] <- y; Rn[i, ] <- qnorm((rank(y, ties.method = "average") - 0.5) / nU) }
write_bed <- function(M, f) { bd <- data.table(`#Chr` = R$chr, start = R$tss - 1L, end = R$tss, pid = R$promoter, gid = R$promoter, strand = INFO$strand[match(R$promoter, INFO$pid)])
  stopifnot(!anyNA(bd$strand)); bd <- cbind(bd, as.data.table(round(M, 6))); setorder(bd, `#Chr`, start); fwrite(bd, f, sep = "\t", quote = FALSE)
  stopifnot(system2(file.path(HTS, "bgzip"), c("-f", f)) == 0, system2(file.path(HTS, "tabix"), c("-f", "-p", "bed", paste0(f, ".gz"))) == 0) }
write_bed(Rw, file.path(TMP, "raw.bed")); write_bed(Rn, file.path(TMP, "rint.bed")); rm(Rw, Rn)
cv <- data.table(id = colnames(C), t(C)); setnames(cv, c("id", U)); fwrite(cv, file.path(TMP, "cov.txt"), sep = "\t", quote = FALSE)

## ---- per group: phenotype list and the sites BOHR does not test (counted on BOHR's analysis samples); samples = all donors ----
smp <- function(k) U[as.integer(strsplit(k, ",", fixed = TRUE)[[1]])]
one_group <- function(i) { g <- Gt[i]; writeLines(R[grp == g$grp, promoter], file.path(TMP, paste0(g$grp, ".phen")))
  if (g$full) { f <- file.path(MV, "sites", paste0(t, ".", g$chr, ".exclude.txt")); stopifnot(file.exists(f)); return(data.table(grp = g$grp, excl = f, n_region = NA_integer_, n_excluded = NA_integer_)) }
  s <- smp(g$skey)
  v <- readVcf(TabixFile(file.path(VCF, paste0(g$chr, ".vcf.gz"))), "hg38",
               ScanVcfParam(which = GRanges(g$chr, IRanges(max(1L, g$from - FLANK), g$to + FLANK)), samples = U, fixed = NA, info = NA, geno = "GT"))
  f <- file.path(TMP, paste0(g$grp, ".excl"))
  if (nrow(v) == 0) { writeLines(character(0), f); return(data.table(grp = g$grp, excl = f, n_region = 0L, n_excluded = 0L)) }
  D <- impute_modal(gt_dosage(geno(v)$GT[, U, drop = FALSE]))[, s, drop = FALSE]
  bad <- pmin(rowSums(D == 0L), rowSums(D == 1L), rowSums(D == 2L)) < MIN_GENO
  writeLines(unique(rownames(D)[bad]), f); data.table(grp = g$grp, excl = f, n_region = nrow(D), n_excluded = sum(bad)) }
if (nrow(Gt)) { E <- rbindlist(mclapply(seq_len(nrow(Gt)), one_group, mc.cores = NC, mc.preschedule = TRUE)); stopifnot(nrow(E) == nrow(Gt)); Gt <- merge(Gt, E, by = "grp")
  ## ---- QTLtools: as the original run, plus the site list ----
  cmd <- unlist(lapply(names(MODES), function(md) sprintf(
    "%s cis --vcf %s/%s.vcf.gz --bed %s/%s.bed.gz --cov %s/cov.txt --include-phenotypes %s/%s.phen --exclude-sites %s --permute 1000 --window 500000 --seed 123 %s --region %s:%d-%d --out %s/out/%s.%s.txt --silent > %s/out/%s.%s.log 2>&1",
    QB, VCF, Gt$chr, TMP, MODES[[md]], TMP, TMP, Gt$grp, Gt$excl, if (md == "normal") "--normal" else "", Gt$chr, pmax(1L, Gt$from - 1L), Gt$to, TMP, md, Gt$grp, TMP, md, Gt$grp)))
  writeLines(cmd, file.path(TMP, "cmds.txt")); t0 <- Sys.time()
  rc <- system(sprintf("xargs -P %d -I{} bash -c '{}' < %s/cmds.txt", NC, TMP)); cat("QTLtools runs:", length(cmd), "| xargs status", rc, "| minutes", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "\n") }

## ---- collect: rows run here, plus the rows whose inputs another finished run already had ----
COLS <- c("promoter", "n_s", "n_nonmiss", "full", "n_variants", "top_variant", "dof", "p_nominal", "p_empirical", "p")
for (md in names(MODES)) {
  fs <- file.path(TMP, "out", paste0(md, ".", Gt$grp, ".txt")); ok <- file.exists(fs) & file.size(fs) > 0
  d <- rbindlist(lapply(fs[ok], function(f) fread(f, header = FALSE, select = c(1, 6, 8, 12, 16, 18, 19))))
  if (!nrow(d)) d <- data.table(promoter = character(), n_variants = integer(), top_variant = character(), dof = integer(), p_nominal = numeric(), p_empirical = numeric(), p = numeric()) else
    setnames(d, c("promoter", "n_variants", "top_variant", "dof", "p_nominal", "p_empirical", "p"))
  d <- merge(R[, .(promoter, n_s, n_nonmiss, full)], d, by = "promoter", all.x = TRUE)[, source := "run"]
  L <- list(d)
  if (have_exact) L$exact <- fr(file.path(MV, "exact", ver, md, fname), select = COLS)[promoter %in% P[source == "exact", promoter]][, source := "exact"]
  if (have_keep)  L$keep  <- fr(file.path(MV, "exact_keep", ver, md, fname), select = COLS)[promoter %in% P[source == "keep", promoter]][, source := "keep"]
  o <- rbindlist(L, use.names = TRUE)[, n_used := nU][order(promoter)]
  cat(md, ": run here", nrow(d), "| with a result", sum(!is.na(d$p)), "| groups without output", sum(!ok), "| all rows: dof = donors - 2 in", sum(o$dof == nU - 2, na.rm = TRUE), "of", sum(!is.na(o$p)),
      "| by source:", paste(names(table(o$source)), table(o$source), collapse = ", "), "\n")
  if (MAXG == 0) { stopifnot(nrow(o) == nrow(P), !anyDuplicated(o$promoter), all(o$dof == nU - 2, na.rm = TRUE))
    dir.create(file.path(OUT, md), recursive = TRUE, showWarnings = FALSE); fwrite(o, file.path(OUT, md, fname), sep = "\t") } else
    fwrite(o, file.path(MV, "logs", paste("test_pairs", ver, md, t, arm, "tsv", sep = ".")), sep = "\t") }
if (nrow(Gt)) fwrite(Gt[, .(grp, chr, n, full, n_s, n_region, n_excluded)], file.path(MV, "logs", paste("groups_pairs", ver, t, arm, "tsv", sep = ".")), sep = "\t")
cat("DONE\n")
