## fill_qtltools.R: QTLtools on the promoter-variant pairs BOHR tested, with the donors without a value left out (no fill).
## The no-fill tables (pairs_nofill/) are the QTLtools results of Fig. 2E; the zero-fill tables are kept for reference.
## QTLtools 1.2 cis --permute on the same promoter-variant pairs as BOHR (as pairs_qtltools.R), with the two treatments of a
## missing relative activity (the gene has no reads in that donor), in place of the promoter's mean:
##   nofill = the donors with a missing value take no part in that promoter's test (--include-samples = donors with a value);
##            rint = rank-normal over those donors, ties averaged. QTLtools itself always replaces a missing phenotype by the
##            mean (imputePhenotypes, cis mode), so "no fill" has to go through the sample list.
##   zero   = all donors, a missing value set to 0; rint = rank-normal of that vector over all donors, ties averaged.
## Everything else as the pairs run: promoters = those of the units BOHR tested; variants = exactly the variants BOHR tested for
## the promoter (each genotype class has at least 5 of BOHR's analysis samples); covariates, window (500 kb), permutations
## (1000), seed, modes (raw, rint, normal = --normal on raw). BOHR is untouched (it drops missing values).
## Only promoters with a missing value are run here, grouped by (donors with a value, BOHR's analysis samples) within a gene;
## one read of the genotypes per group serves both treatments. Rows of promoters without a missing value are copied from the
## pairs run (same samples, values and variants under either treatment). When the keep run of the tissue and arm is finished
## (exact_keep/: non-missing donors, same variants), its rows are the nofill rows and only zero is run.
## Output: pairs_nofill/<version>/<mode>/<Tissue>.<arm>.promoters.tsv.gz and pairs_zero/<version>/<mode>/...: promoter, n_s,
## n_nonmiss, full, n_variants, top_variant, dof, p_nominal, p_empirical, p, source (run | pairs | keep), n_used (samples
## QTLtools used). Group files live in a node-local folder and are removed. No genotype or donor-level value is written to the
## project.
## Usage: Rscript fill_qtltools.R <v4|paper> <Tissue> <arm> <ncores>
##   FILL_MAXGENES=n  test on the first n genes;  FILL_NOREUSE=1  run nofill here even when the keep run exists;
##   FILL_OUT=dir     write under dir/pairs_nofill and dir/pairs_zero instead
suppressPackageStartupMessages({ library(data.table); library(VariantAnnotation); library(parallel) })
a <- commandArgs(TRUE); ver <- a[1]; t <- a[2]; arm <- a[3]; NC <- as.integer(a[4]); MAXG <- as.integer(Sys.getenv("FILL_MAXGENES", "0")); NOREUSE <- Sys.getenv("FILL_NOREUSE", "0") == "1"
G <- "/projects/group/bliu/e1327989/GTEx"; W <- file.path(G, "nofloor_full")
V2 <- c(v4 = file.path(W, "GTEx/QTL_figure/rerun_v4"), paper = file.path(G, "QTL_figure/rerun_v2"))[[ver]]
INFO <- fread(c(v4 = file.path(W, "GTEx/v11/Promoter_info_filter_5.txt"), paper = file.path(G, "v11/Promoter_info_filter_5.txt"))[[ver]]); setnames(INFO, 1, "pid")
VCF <- file.path(G, "Rawdata/v9_vcf_MAF_0.05"); MV <- file.path(W, "consult/fig2e_20261119/matched_variants")
OUTB <- Sys.getenv("FILL_OUT", MV); OUT <- c(nofill = file.path(OUTB, "pairs_nofill", ver), zero = file.path(OUTB, "pairs_zero", ver))
QB <- "/home/e1327989/project/software/QTLtools_1.2_Ubuntu16.04_x86_64/QTLtools_1.2_Ubuntu16.04_x86_64"
HTS <- "/home/e1327989/.local/share/mamba/envs/htsenv/bin"; MIN_GENO <- 5L; FLANK <- 500000L
MODES <- c(raw = "raw", rint = "rint", normal = "raw")                     # mode -> BED kind
TMP <- file.path(Sys.getenv("EXACT_TMP", tempdir()), paste("fill", ver, t, arm, Sys.getpid(), sep = ".")); dir.create(file.path(TMP, "out"), recursive = TRUE)
fname <- paste(t, arm, "promoters.tsv.gz", sep = "."); final <- unlist(lapply(OUT, function(d) file.path(d, names(MODES), fname)))
if (MAXG == 0 && all(file.exists(final))) { cat("exists:", t, arm, "\n"); quit(save = "no") }
done <- function(dir, logs) { f <- file.path(MV, logs, paste(ver, t, arm, "exit", sep = ".")); all(file.exists(file.path(MV, dir, ver, names(MODES), fname))) && file.exists(f) && identical(trimws(readLines(f, n = 1)), "0") }
if (!done("pairs", "logs/run_pairs")) stop("the pairs run of ", ver, " ", t, " ", arm, " is not finished")
have_keep <- !NOREUSE && done("exact_keep", "logs/run_keep")

## ---- BOHR's own rules (prepare_y of bohr2_data.R; v4 and paper differ in the SD rule), as in pairs_qtltools.R ----
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

## ---- promoters: those of the units compared in Fig. 2E (the pairs run's promoters), with BOHR's analysis samples ----
fr <- function(f, ...) fread(cmd = paste("zcat", shQuote(f)), ...)
b <- fr(file.path(V2, "results/fig2e_freeze_v3m/bohr", paste0(t, ".real.tsv.gz")), select = c("promoter", "promoters_in_unit"))
X <- fr(file.path(W, "consult/fig2e_20261119/main/out/units_paper_r5e80.tsv.gz"), select = c("set", "tissue", "unit"))
units <- X[tissue == t & set == (if (ver == "v4") "r5e80" else "paper"), unit]
pr <- unique(unlist(strsplit(b[promoter %in% units, promoters_in_unit], ";", fixed = TRUE))); pr <- intersect(pr, colnames(A))
pr <- pr[td$info[pr, "seqnames"] %in% paste0("chr", 1:22)]
P <- rbindlist(lapply(pr, function(p) { y <- A[, p]; names(y) <- U; z <- prepare_y(y, pc)
  if (!is.na(z$skip) || !any(z$y > 0 & z$y < 1)) return(NULL)
  data.table(promoter = p, n_s = length(z$y), n_nonmiss = sum(!is.na(y)), skey = paste(match(names(z$y), U), collapse = ","), nmkey = paste(which(!is.na(y)), collapse = ",")) }))
P[, `:=`(chr = td$info[promoter, "seqnames"], tss = as.integer(td$info[promoter, "start"]), gene = td$info[promoter, "geneId"], full = n_s == nU)]
COLS <- c("promoter", "n_s", "n_nonmiss", "full", "n_variants", "top_variant", "dof", "p_nominal", "p_empirical", "p")
PR <- lapply(setNames(names(MODES), names(MODES)), function(md) fr(file.path(MV, "pairs", ver, md, fname), select = COLS))       # the pairs run's rows
stopifnot(setequal(PR$raw$promoter, P$promoter))                              # same promoters and sample counts as the pairs run
chk <- merge(P[, .(promoter, n_s, n_nonmiss)], PR$raw[, .(promoter, a = n_s, b = n_nonmiss)], by = "promoter"); stopifnot(chk$n_s == chk$a, chk$n_nonmiss == chk$b)
R <- P[n_nonmiss < nU]                                                        # a missing value: the two treatments differ from the mean fill here
if (MAXG > 0) R <- R[gene %in% head(unique(gene), MAXG)]
R[, grp := paste0("f", as.integer(factor(paste(nmkey, skey, chr, gene))))]
Gt <- R[, .(chr = chr[1], from = min(tss), to = max(tss), n = .N, skey = skey[1], nmkey = nmkey[1], n_s = n_s[1], n_used = n_nonmiss[1]), by = grp]
KINDS <- c(if (!have_keep) "nofill", "zero")                                  # treatments run here
cat(t, arm, ver, ": donors", nU, "| covariates", pc, "| promoters", nrow(P), "| with missing values (run here)", nrow(R), "| groups", nrow(Gt), "| copied from the pairs run", P[n_nonmiss == nU, .N],
    "| treatments run here:", KINDS, "| nofill rows from the keep run:", have_keep, "\n")

## ---- BED files. nofill: values outside a promoter's donors with a value are placeholders that QTLtools never reads ----
smp <- function(k) U[as.integer(strsplit(k, ",", fixed = TRUE)[[1]])]
Bn <- matrix(0, nrow(R), nU, dimnames = list(R$promoter, U)); Bnr <- Bn; Bz <- Bn; Bzr <- Bn
for (i in seq_len(nrow(R))) { y <- A[, R$promoter[i]]; ok <- !is.na(y); stopifnot(sum(ok) == R$n_nonmiss[i], identical(unname(which(ok)), as.integer(strsplit(R$nmkey[i], ",", fixed = TRUE)[[1]])))
  Bn[i, ] <- mean(y[ok]); Bn[i, ok] <- y[ok]; Bnr[i, ok] <- qnorm((rank(y[ok], ties.method = "average") - 0.5) / sum(ok))
  z <- y; z[!ok] <- 0; stopifnot(all(is.finite(z)), length(unique(z)) > 1); Bz[i, ] <- z; Bzr[i, ] <- qnorm((rank(z, ties.method = "average") - 0.5) / nU) }
write_bed <- function(M, f) { bd <- data.table(`#Chr` = R$chr, start = R$tss - 1L, end = R$tss, pid = R$promoter, gid = R$promoter, strand = INFO$strand[match(R$promoter, INFO$pid)])
  stopifnot(!anyNA(bd$strand)); bd <- cbind(bd, as.data.table(round(M, 6))); setorder(bd, `#Chr`, start); fwrite(bd, f, sep = "\t", quote = FALSE)
  stopifnot(system2(file.path(HTS, "bgzip"), c("-f", f)) == 0, system2(file.path(HTS, "tabix"), c("-f", "-p", "bed", paste0(f, ".gz"))) == 0) }
if (nrow(R)) { if ("nofill" %in% KINDS) { write_bed(Bn, file.path(TMP, "nofill.raw.bed")); write_bed(Bnr, file.path(TMP, "nofill.rint.bed")) }
  write_bed(Bz, file.path(TMP, "zero.raw.bed")); write_bed(Bzr, file.path(TMP, "zero.rint.bed")) }
rm(Bn, Bnr, Bz, Bzr)
cv <- data.table(id = colnames(C), t(C)); setnames(cv, c("id", U)); fwrite(cv, file.path(TMP, "cov.txt"), sep = "\t", quote = FALSE)

## ---- per group: donors with a value (nofill sample list); excluded sites = those BOHR does not test, counted on its analysis samples ----
one_group <- function(i) { g <- Gt[i]; s <- smp(g$skey); nm <- smp(g$nmkey); stopifnot(all(s %in% nm), length(nm) < nU)
  writeLines(nm, file.path(TMP, paste0(g$grp, ".samples"))); writeLines(R[grp == g$grp, promoter], file.path(TMP, paste0(g$grp, ".phen")))
  v <- readVcf(TabixFile(file.path(VCF, paste0(g$chr, ".vcf.gz"))), "hg38",
               ScanVcfParam(which = GRanges(g$chr, IRanges(max(1L, g$from - FLANK), g$to + FLANK)), samples = U, fixed = NA, info = NA, geno = "GT"))
  f <- file.path(TMP, paste0(g$grp, ".excl"))
  if (nrow(v) == 0) { writeLines(character(0), f); return(data.table(grp = g$grp, excl = f, n_region = 0L, n_excluded = 0L)) }
  D <- impute_modal(gt_dosage(geno(v)$GT[, U, drop = FALSE]))[, s, drop = FALSE]
  bad <- pmin(rowSums(D == 0L), rowSums(D == 1L), rowSums(D == 2L)) < MIN_GENO
  writeLines(unique(rownames(D)[bad]), f); data.table(grp = g$grp, excl = f, n_region = nrow(D), n_excluded = sum(bad)) }
if (nrow(Gt)) { E <- rbindlist(mclapply(seq_len(nrow(Gt)), one_group, mc.cores = NC, mc.preschedule = TRUE)); stopifnot(nrow(E) == nrow(Gt)); Gt <- merge(Gt, E, by = "grp")
  ## ---- QTLtools: as the pairs run; nofill adds the sample list ----
  cmd <- unlist(lapply(KINDS, function(kd) unlist(lapply(names(MODES), function(md) sprintf(
    "%s cis --vcf %s/%s.vcf.gz --bed %s/%s.%s.bed.gz --cov %s/cov.txt %s --include-phenotypes %s/%s.phen --exclude-sites %s --permute 1000 --window 500000 --seed 123 %s --region %s:%d-%d --out %s/out/%s.%s.%s.txt --silent > %s/out/%s.%s.%s.log 2>&1",
    QB, VCF, Gt$chr, TMP, kd, MODES[[md]], TMP, if (kd == "nofill") sprintf("--include-samples %s/%s.samples", TMP, Gt$grp) else "", TMP, Gt$grp, Gt$excl, if (md == "normal") "--normal" else "",
    Gt$chr, pmax(1L, Gt$from - 1L), Gt$to, TMP, kd, md, Gt$grp, TMP, kd, md, Gt$grp)))))
  writeLines(cmd, file.path(TMP, "cmds.txt")); t0 <- Sys.time()
  rc <- system(sprintf("xargs -P %d -I{} bash -c '{}' < %s/cmds.txt", NC, TMP)); cat("QTLtools runs:", length(cmd), "| xargs status", rc, "| minutes", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "\n") }

## ---- collect: rows run here (or the keep run's), plus the pairs run's rows of the promoters without a missing value ----
empty <- data.table(promoter = character(), n_variants = integer(), top_variant = character(), dof = integer(), p_nominal = numeric(), p_empirical = numeric(), p = numeric())
for (kd in names(OUT)) for (md in names(MODES)) {
  base <- PR[[md]][promoter %in% P[n_nonmiss == nU, promoter]][, `:=`(source = "pairs", n_used = nU)]; nok <- 0L
  if (kd == "nofill" && have_keep) {
    d <- merge(fr(file.path(MV, "exact_keep", ver, md, fname), select = COLS), R[, .(promoter)], by = "promoter")[, `:=`(source = "keep", n_used = n_nonmiss)]
  } else {
    fs <- file.path(TMP, "out", paste(kd, md, Gt$grp, "txt", sep = ".")); ok <- file.exists(fs) & file.size(fs) > 0; nok <- sum(!ok)
    d <- rbindlist(lapply(fs[ok], function(f) fread(f, header = FALSE, select = c(1, 6, 8, 12, 16, 18, 19))))
    if (!nrow(d)) d <- empty else setnames(d, c("promoter", "n_variants", "top_variant", "dof", "p_nominal", "p_empirical", "p"))
    d <- merge(R[, .(promoter, n_s, n_nonmiss, full)], d, by = "promoter", all.x = TRUE)[, source := "run"]; d[, n_used := if (kd == "zero") nU else n_nonmiss] }
  o <- rbindlist(list(base, d), use.names = TRUE)[order(promoter)]
  x <- merge(o[source != "pairs"], PR[[md]][, .(promoter, nv_pairs = n_variants, p_pairs = p)], by = "promoter")
  cat(kd, md, ": rows", nrow(o), "| by source:", paste(names(table(o$source)), table(o$source), collapse = ", "), "| with a result", sum(!is.na(o$p)), "| groups without output", nok,
      "| dof = samples used - 2 in", sum(o$dof == o$n_used - 2, na.rm = TRUE), "| variants tested as in the pairs run:", sum(x$n_variants == x$nv_pairs, na.rm = TRUE), "of", nrow(x),
      "| P differs from the pairs run in", sum(x$p != x$p_pairs, na.rm = TRUE), "\n")
  if (MAXG == 0) { stopifnot(nrow(o) == nrow(P), !anyDuplicated(o$promoter), all(o$dof == o$n_used - 2, na.rm = TRUE), nok == 0)
    dir.create(file.path(OUT[[kd]], md), recursive = TRUE, showWarnings = FALSE); fwrite(o, file.path(OUT[[kd]], md, fname), sep = "\t") } else
    fwrite(x, file.path(MV, "logs", paste("test_fill", kd, ver, md, t, arm, "tsv", sep = ".")), sep = "\t") }
if (nrow(Gt)) fwrite(Gt[, .(grp, chr, n, n_s, n_used, n_region, n_excluded)], file.path(MV, "logs", paste("groups_fill", ver, t, arm, "tsv", sep = ".")), sep = "\t")
unlink(TMP, recursive = TRUE); cat("DONE\n")
