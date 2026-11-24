## keep_qtltools.R: QTLtools on BOHR's variants and every donor with a value (the second of the matched QTLtools runs).
## Its rows are the no-fill rows of the real arm (fill_qtltools.R) and are reused by pairs_qtltools.R.
## Second input variant for QTLtools 1.2 cis --permute. As exact_qtltools.R except for the samples:
##   samples  = every donor with a non-missing value of the promoter (missing values dropped, as BOHR does). The boundary values
##              that BOHR's prepare_y deletes (a value of 0, or of 1, held by <= 5% of the samples) stay in.
##   variants = exactly the variants BOHR tested for the promoter: each genotype class has at least 5 of BOHR's own analysis
##              samples (after its deletion), missing genotypes imputed to the modal dosage over the tissue's donors.
##   covariates, window, permutations, seed, modes (raw; rint = rank-normal over the samples used, ties averaged; normal): as there.
## The two variants differ only for promoters where BOHR deleted samples. Those promoters are run here, grouped by identical
## (non-missing samples, BOHR analysis samples) within a gene. The rows of every other promoter are copied from the exact run:
## its inputs there (samples, variants, covariates, seed) are the ones this variant asks for.
## Output: exact_keep/<version>/<mode>/<Tissue>.<arm>.promoters.tsv.gz = the columns of the exact run, plus n_used (samples
## QTLtools used) and rerun (TRUE = run here, FALSE = copied). Group files live in a node-local folder and are removed. No
## genotype or donor-level value is written to the project.
## Usage: Rscript keep_qtltools.R <v4|paper> <Tissue> <arm> <ncores>      (KEEP_MAXGENES=n: test on the first n genes)
suppressPackageStartupMessages({ library(data.table); library(VariantAnnotation); library(parallel) })
a <- commandArgs(TRUE); ver <- a[1]; t <- a[2]; arm <- a[3]; NC <- as.integer(a[4]); MAXG <- as.integer(Sys.getenv("KEEP_MAXGENES", "0"))
G <- "/projects/group/bliu/e1327989/GTEx"; W <- file.path(G, "nofloor_full")
V2 <- c(v4 = file.path(W, "GTEx/QTL_figure/rerun_v4"), paper = file.path(G, "QTL_figure/rerun_v2"))[[ver]]
INFO <- fread(c(v4 = file.path(W, "GTEx/v11/Promoter_info_filter_5.txt"), paper = file.path(G, "v11/Promoter_info_filter_5.txt"))[[ver]]); setnames(INFO, 1, "pid")
VCF <- file.path(G, "Rawdata/v9_vcf_MAF_0.05"); MV <- file.path(W, "consult/fig2e_20261119/matched_variants")
EXA <- file.path(MV, "exact", ver); OUT <- file.path(MV, "exact_keep", ver); QB <- "/home/e1327989/project/software/QTLtools_1.2_Ubuntu16.04_x86_64/QTLtools_1.2_Ubuntu16.04_x86_64"
HTS <- "/home/e1327989/.local/share/mamba/envs/htsenv/bin"; MIN_GENO <- 5L; FLANK <- 500000L
MODES <- c(raw = "raw", rint = "rint", normal = "raw")                     # mode -> BED kind
TMP <- file.path(Sys.getenv("EXACT_TMP", tempdir()), paste("keep", ver, t, arm, Sys.getpid(), sep = ".")); dir.create(file.path(TMP, "out"), recursive = TRUE)
on.exit(unlink(TMP, recursive = TRUE), add = TRUE)
fname <- paste(t, arm, "promoters.tsv.gz", sep = "."); final <- file.path(OUT, names(MODES), fname); exf <- setNames(file.path(EXA, names(MODES), fname), names(MODES))
if (MAXG == 0 && all(file.exists(final))) { cat("exists:", t, arm, "\n"); quit(save = "no") }
if (!all(file.exists(exf))) stop("the exact run of ", ver, " ", t, " ", arm, " is not there yet")

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
td <- readRDS(file.path(V2, "work/final/freeze", paste0(t, ".rds"))); U <- td$U; K <- td$K
if (arm == "real") { A <- td$act[U, , drop = FALSE]; Pe <- td$peer[U, , drop = FALSE] } else {
  m <- fread(file.path(V2, "perm_maps/final/freeze", paste0(t, ".", arm, ".perm_map.tsv")))
  stopifnot(setequal(m$slot, U), setequal(m$donor, U), !anyDuplicated(m$donor))
  d <- setNames(m$donor, m$slot)[U]; A <- td$act[d, , drop = FALSE]; Pe <- td$peer[d, , drop = FALSE]; rownames(A) <- rownames(Pe) <- U }
C <- cbind(td$pcs[U, , drop = FALSE], Pe[, seq_len(K), drop = FALSE]); colnames(C) <- c(paste0("PC", 1:ncol(td$pcs)), if (K > 0) paste0("peer_", seq_len(K)))
pc <- ncol(C)

## ---- promoters: those of the units compared in Fig. 2E (the exact run's promoters) ----
fr <- function(f, ...) fread(cmd = paste("zcat", shQuote(f)), ...)
b <- fr(file.path(V2, "results/fig2e_freeze_v3m/bohr", paste0(t, ".real.tsv.gz")), select = c("promoter", "promoters_in_unit"))
X <- fr(file.path(W, "consult/fig2e_20261119/main/out/units_paper_r5e80.tsv.gz"), select = c("set", "tissue", "unit"))
units <- X[tissue == t & set == (if (ver == "v4") "r5e80" else "paper"), unit]
pr <- unique(unlist(strsplit(b[promoter %in% units, promoters_in_unit], ";", fixed = TRUE))); pr <- intersect(pr, colnames(A))
pr <- pr[td$info[pr, "seqnames"] %in% paste0("chr", 1:22)]
P <- rbindlist(lapply(pr, function(p) { y <- A[, p]; names(y) <- U; z <- prepare_y(y, pc)
  if (!is.na(z$skip) || !any(z$y > 0 & z$y < 1)) return(NULL)
  data.table(promoter = p, n_s = length(z$y), n_nonmiss = sum(!is.na(y)), skey = paste(match(names(z$y), U), collapse = ","), nmkey = paste(which(!is.na(y)), collapse = ",")) }))
P[, `:=`(chr = td$info[promoter, "seqnames"], tss = as.integer(td$info[promoter, "start"]), gene = td$info[promoter, "geneId"], full = n_s == length(U))]
E0 <- fr(exf[["raw"]]); stopifnot(setequal(E0$promoter, P$promoter))             # same promoters, same sample counts as the exact run
chk <- merge(P[, .(promoter, n_s, n_nonmiss)], E0[, .(promoter, n_s0 = n_s, n_nm0 = n_nonmiss)], by = "promoter"); stopifnot(chk$n_s == chk$n_s0, chk$n_nonmiss == chk$n_nm0)
R <- P[n_s < n_nonmiss]                                                           # BOHR deleted boundary samples here: run with them kept
if (MAXG > 0) R <- R[gene %in% head(unique(gene), MAXG)]
R[, grp := paste0("k", as.integer(factor(paste(nmkey, skey, chr, gene))))]
Gt <- R[, .(chr = chr[1], from = min(tss), to = max(tss), n = .N, skey = skey[1], nmkey = nmkey[1], n_s = n_s[1], n_used = n_nonmiss[1]), by = grp]
cat(t, arm, ver, ": donors", length(U), "| covariates", pc, "| promoters", nrow(P), "| with deleted samples (run here)", nrow(R), "| groups", nrow(Gt), "| copied", nrow(P) - P[n_s < n_nonmiss, .N], "\n")

## ---- BED (raw, rint) and covariate files; values outside a promoter's non-missing samples are placeholders that QTLtools never reads ----
smp <- function(k) U[as.integer(strsplit(k, ",", fixed = TRUE)[[1]])]
Rw <- matrix(0, nrow(R), length(U), dimnames = list(R$promoter, U)); Rn <- Rw
for (i in seq_len(nrow(R))) { s <- smp(R$nmkey[i]); y <- A[s, R$promoter[i]]; stopifnot(!anyNA(y), length(y) == R$n_nonmiss[i])
  Rw[i, ] <- mean(y); Rw[i, s] <- y; Rn[i, s] <- qnorm((rank(y, ties.method = "average") - 0.5) / length(y)) }
write_bed <- function(M, f) { bd <- data.table(`#Chr` = R$chr, start = R$tss - 1L, end = R$tss, pid = R$promoter, gid = R$promoter, strand = INFO$strand[match(R$promoter, INFO$pid)])
  stopifnot(!anyNA(bd$strand)); bd <- cbind(bd, as.data.table(round(M, 6))); setorder(bd, `#Chr`, start); fwrite(bd, f, sep = "\t", quote = FALSE)
  stopifnot(system2(file.path(HTS, "bgzip"), c("-f", f)) == 0, system2(file.path(HTS, "tabix"), c("-f", "-p", "bed", paste0(f, ".gz"))) == 0) }
write_bed(Rw, file.path(TMP, "raw.bed")); write_bed(Rn, file.path(TMP, "rint.bed")); rm(Rw, Rn)
cv <- data.table(id = colnames(C), t(C)); setnames(cv, c("id", U)); fwrite(cv, file.path(TMP, "cov.txt"), sep = "\t", quote = FALSE)

## ---- per group: samples = non-missing; excluded sites = those BOHR does not test, counted on BOHR's analysis samples ----
one_group <- function(i) { g <- Gt[i]; s <- smp(g$skey); nm <- smp(g$nmkey); stopifnot(all(s %in% nm), length(s) < length(nm))
  writeLines(nm, file.path(TMP, paste0(g$grp, ".samples"))); writeLines(R[grp == g$grp, promoter], file.path(TMP, paste0(g$grp, ".phen")))
  v <- readVcf(TabixFile(file.path(VCF, paste0(g$chr, ".vcf.gz"))), "hg38",
               ScanVcfParam(which = GRanges(g$chr, IRanges(max(1L, g$from - FLANK), g$to + FLANK)), samples = U, fixed = NA, info = NA, geno = "GT"))
  f <- file.path(TMP, paste0(g$grp, ".excl"))
  if (nrow(v) == 0) { writeLines(character(0), f); return(data.table(grp = g$grp, excl = f, n_region = 0L, n_excluded = 0L)) }
  D <- impute_modal(gt_dosage(geno(v)$GT[, U, drop = FALSE]))[, s, drop = FALSE]
  bad <- pmin(rowSums(D == 0L), rowSums(D == 1L), rowSums(D == 2L)) < MIN_GENO
  writeLines(unique(rownames(D)[bad]), f); data.table(grp = g$grp, excl = f, n_region = nrow(D), n_excluded = sum(bad)) }
E <- rbindlist(mclapply(seq_len(nrow(Gt)), one_group, mc.cores = NC, mc.preschedule = TRUE)); stopifnot(nrow(E) == nrow(Gt)); Gt <- merge(Gt, E, by = "grp")

## ---- QTLtools ----
cmd <- unlist(lapply(names(MODES), function(md) sprintf(
  "%s cis --vcf %s/%s.vcf.gz --bed %s/%s.bed.gz --cov %s/cov.txt --include-samples %s/%s.samples --include-phenotypes %s/%s.phen --exclude-sites %s --permute 1000 --window 500000 --seed 123 %s --region %s:%d-%d --out %s/out/%s.%s.txt --silent > %s/out/%s.%s.log 2>&1",
  QB, VCF, Gt$chr, TMP, MODES[[md]], TMP, TMP, Gt$grp, TMP, Gt$grp, Gt$excl, if (md == "normal") "--normal" else "", Gt$chr, pmax(1L, Gt$from - 1L), Gt$to, TMP, md, Gt$grp, TMP, md, Gt$grp)))
writeLines(cmd, file.path(TMP, "cmds.txt")); t0 <- Sys.time()
rc <- system(sprintf("xargs -P %d -I{} bash -c '{}' < %s/cmds.txt", NC, TMP)); cat("QTLtools runs:", length(cmd), "| xargs status", rc, "| minutes", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "\n")

## ---- collect: rerun rows from here, every other row from the exact run ----
for (md in names(MODES)) {
  fs <- file.path(TMP, "out", paste0(md, ".", Gt$grp, ".txt")); ok <- file.exists(fs) & file.size(fs) > 0
  d <- rbindlist(lapply(fs[ok], function(f) fread(f, header = FALSE, select = c(1, 6, 8, 12, 16, 18, 19))))
  setnames(d, c("promoter", "n_variants", "top_variant", "dof", "p_nominal", "p_empirical", "p"))
  d <- merge(R[, .(promoter, n_s, n_nonmiss, full)], d, by = "promoter", all.x = TRUE)[, `:=`(n_used = n_nonmiss, rerun = TRUE)]
  e <- fr(exf[[md]]); x <- merge(d[, .(promoter, n_variants, dof, n_nonmiss, p)], e[, .(promoter, nv_exact = n_variants, p_exact = p)], by = "promoter")
  cat(md, ": run here", nrow(d), "| with a result", sum(!is.na(d$p)), "| groups without output", sum(!ok), "| variants tested equal to the exact run", sum(x$n_variants == x$nv_exact, na.rm = TRUE),
      "| dof = non-missing samples - 2 in", sum(x$dof == x$n_nonmiss - 2, na.rm = TRUE), "| P differs from the exact run in", sum(x$p != x$p_exact, na.rm = TRUE), "\n")
  if (MAXG == 0) { o <- rbind(e[!promoter %in% R$promoter][, `:=`(n_used = n_s, rerun = FALSE)], d, use.names = TRUE)[order(promoter)]; stopifnot(nrow(o) == nrow(P), !anyDuplicated(o$promoter))
    dir.create(file.path(OUT, md), recursive = TRUE, showWarnings = FALSE); fwrite(o, file.path(OUT, md, fname), sep = "\t") } else
    fwrite(x, file.path(MV, "logs", paste("test_keep", ver, md, t, arm, "tsv", sep = ".")), sep = "\t") }
fwrite(Gt[, .(grp, chr, n, n_s, n_used, n_region, n_excluded)], file.path(MV, "logs", paste("groups_keep", ver, t, arm, "tsv", sep = ".")), sep = "\t")
cat("DONE\n")
