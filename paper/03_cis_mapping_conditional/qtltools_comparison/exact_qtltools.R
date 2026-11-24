## exact_qtltools.R: QTLtools on exactly BOHR's inputs (the first of the matched QTLtools runs of Fig. 2E).
## Its rows are reused by keep_qtltools.R and pairs_qtltools.R wherever their inputs are the same.
## Runs QTLtools 1.2 cis --permute on exactly the inputs BOHR analysed, promoter by promoter:
##   samples  = the promoter's analysis samples after BOHR's prepare_y (missing values dropped; a boundary value held by <= 5% of
##              the samples deleted), in the arm's own permutation;
##   variants = the cis variants BOHR keeps for those samples: each genotype class has at least 5 of them (MIN_GENO = 5), with
##              missing genotypes imputed to the modal dosage over the tissue's donors, as in bohr2_data.R;
##   covariates, window (500 kb), permutations (1000) and seed as in the original comparison (make_freeze_beds.R; command in run.sh).
## QTLtools takes one sample list and one site list per run, so promoters are grouped by identical analysis samples (within a
## gene, unless the sample set is the tissue's full donor set, which is run by chromosome with the tissue-level site list of
## make_site_lists.R). Modes: raw; rint (rank-normal over the promoter's own analysis samples, ties averaged, then QTLtools
## regresses the covariates); normal (QTLtools --normal on the raw values).
## Units: v4 = the promoters of the r5e80 units (gene median reads >= 5, neither share of 0 nor of 1 >= 80%); paper = all units.
## Output: exact/<version>/<mode>/<Tissue>.<arm>.promoters.tsv.gz (one row per promoter: samples used, variants tested, top
## variant, nominal, empirical and beta-approximated P). Group files live in a node-local folder and are removed. No genotype
## or donor-level value is written to the project.
## Usage: Rscript exact_qtltools.R <v4|paper> <Tissue> <arm> <ncores>      (EXACT_MAXGENES=n: test on the first n genes)
suppressPackageStartupMessages({ library(data.table); library(VariantAnnotation); library(parallel) })
a <- commandArgs(TRUE); ver <- a[1]; t <- a[2]; arm <- a[3]; NC <- as.integer(a[4]); MAXG <- as.integer(Sys.getenv("EXACT_MAXGENES", "0"))
G <- "/projects/group/bliu/e1327989/GTEx"; W <- file.path(G, "nofloor_full")
V2 <- c(v4 = file.path(W, "GTEx/QTL_figure/rerun_v4"), paper = file.path(G, "QTL_figure/rerun_v2"))[[ver]]
INFO <- fread(c(v4 = file.path(W, "GTEx/v11/Promoter_info_filter_5.txt"), paper = file.path(G, "v11/Promoter_info_filter_5.txt"))[[ver]]); setnames(INFO, 1, "pid")
VCF <- file.path(G, "Rawdata/v9_vcf_MAF_0.05"); MV <- file.path(W, "consult/fig2e_20261119/matched_variants")
OUT <- file.path(MV, "exact", ver); QB <- "/home/e1327989/project/software/QTLtools_1.2_Ubuntu16.04_x86_64/QTLtools_1.2_Ubuntu16.04_x86_64"
HTS <- "/home/e1327989/.local/share/mamba/envs/htsenv/bin"; MIN_GENO <- 5L; FLANK <- 500000L
MODES <- c(raw = "raw", rint = "rint", normal = "raw")                     # mode -> BED kind
TMP <- file.path(Sys.getenv("EXACT_TMP", tempdir()), paste(ver, t, arm, Sys.getpid(), sep = ".")); dir.create(file.path(TMP, "out"), recursive = TRUE)
on.exit(unlink(TMP, recursive = TRUE), add = TRUE)
final <- file.path(OUT, names(MODES), paste(t, arm, "promoters.tsv.gz", sep = "."))
if (MAXG == 0 && all(file.exists(final))) { cat("exists:", t, arm, "\n"); quit(save = "no") }

## ---- BOHR's own rules (prepare_y of bohr2_data.R; v4 and paper differ in the SD rule) ----
SD_MIN <- if (ver == "v4") 0.005 else NA                                    # paper: sd <= 0.05 skip (rerun_v2); set below from its own code
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

## ---- promoters: those of the units compared in Fig. 2E ----
fr <- function(f, ...) fread(cmd = paste("zcat", shQuote(f)), ...)
b <- fr(file.path(V2, "results/fig2e_freeze_v3m/bohr", paste0(t, ".real.tsv.gz")), select = c("promoter", "promoters_in_unit"))
X <- fr(file.path(W, "consult/fig2e_20261119/main/out/units_paper_r5e80.tsv.gz"), select = c("set", "tissue", "unit"))
units <- X[tissue == t & set == (if (ver == "v4") "r5e80" else "paper"), unit]
pr <- unique(unlist(strsplit(b[promoter %in% units, promoters_in_unit], ";", fixed = TRUE))); pr <- intersect(pr, colnames(A))
pr <- pr[td$info[pr, "seqnames"] %in% paste0("chr", 1:22)]
P <- rbindlist(lapply(pr, function(p) { y <- A[, p]; names(y) <- U; z <- prepare_y(y, pc)
  if (!is.na(z$skip) || !any(z$y > 0 & z$y < 1)) return(NULL)
  data.table(promoter = p, n_s = length(z$y), n_nonmiss = sum(!is.na(y)), skey = paste(match(names(z$y), U), collapse = ",")) }))
P[, `:=`(chr = td$info[promoter, "seqnames"], tss = as.integer(td$info[promoter, "start"]), gene = td$info[promoter, "geneId"], full = n_s == length(U))]
if (MAXG > 0) { gs <- head(unique(P[full == FALSE, gene]), MAXG); P <- P[gene %in% gs | (full == TRUE & chr == "chr22")] }
P[, grp := fifelse(full, paste0("full.", chr), paste0("g", as.integer(factor(paste(skey, chr, gene)))))]
Gt <- P[, .(chr = chr[1], from = min(tss), to = max(tss), n = .N, full = full[1], skey = skey[1], n_s = n_s[1]), by = grp]
cat(t, arm, ver, ": donors", length(U), "| covariates", pc, "| promoters", nrow(P), "| full sample set", sum(P$full), "| groups", nrow(Gt), "\n")

## ---- BED (raw, rint) and covariate files; values outside a promoter's analysis samples are placeholders that QTLtools never reads ----
smp <- function(k) U[as.integer(strsplit(k, ",", fixed = TRUE)[[1]])]
Rw <- matrix(0, nrow(P), length(U), dimnames = list(P$promoter, U)); Rn <- Rw
for (i in seq_len(nrow(P))) { s <- smp(P$skey[i]); y <- A[s, P$promoter[i]]
  Rw[i, ] <- mean(y); Rw[i, s] <- y; Rn[i, s] <- qnorm((rank(y, ties.method = "average") - 0.5) / length(y)) }
write_bed <- function(M, f) { bd <- data.table(`#Chr` = P$chr, start = P$tss - 1L, end = P$tss, pid = P$promoter, gid = P$promoter, strand = INFO$strand[match(P$promoter, INFO$pid)])
  stopifnot(!anyNA(bd$strand)); bd <- cbind(bd, as.data.table(round(M, 6))); setorder(bd, `#Chr`, start); fwrite(bd, f, sep = "\t", quote = FALSE)
  stopifnot(system2(file.path(HTS, "bgzip"), c("-f", f)) == 0, system2(file.path(HTS, "tabix"), c("-f", "-p", "bed", paste0(f, ".gz"))) == 0) }
write_bed(Rw, file.path(TMP, "raw.bed")); write_bed(Rn, file.path(TMP, "rint.bed")); rm(Rw, Rn)
cv <- data.table(id = colnames(C), t(C)); setnames(cv, c("id", U)); fwrite(cv, file.path(TMP, "cov.txt"), sep = "\t", quote = FALSE)

## ---- per group: sample list, phenotype list, excluded sites ----
one_group <- function(i) { g <- Gt[i]; s <- smp(g$skey)
  writeLines(s, file.path(TMP, paste0(g$grp, ".samples"))); writeLines(P[grp == g$grp, promoter], file.path(TMP, paste0(g$grp, ".phen")))
  if (g$full) { f <- file.path(MV, "sites", paste0(t, ".", g$chr, ".exclude.txt")); stopifnot(file.exists(f)); return(data.table(grp = g$grp, excl = f, n_region = NA_integer_, n_excluded = NA_integer_)) }
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

## ---- collect ----
for (md in names(MODES)) {
  fs <- file.path(TMP, "out", paste0(md, ".", Gt$grp, ".txt")); ok <- file.exists(fs) & file.size(fs) > 0
  d <- rbindlist(lapply(fs[ok], function(f) fread(f, header = FALSE, select = c(1, 6, 8, 12, 16, 18, 19))))
  setnames(d, c("promoter", "n_variants", "top_variant", "dof", "p_nominal", "p_empirical", "p"))
  d <- merge(P[, .(promoter, grp, n_s, n_nonmiss, full)], d, by = "promoter", all.x = TRUE)
  cat(md, ": promoters with a result", sum(!is.na(d$p)), "of", nrow(P), "| groups without output", sum(!ok), "\n")
  if (MAXG == 0) { dir.create(file.path(OUT, md), recursive = TRUE, showWarnings = FALSE); fwrite(d[, !"grp"], file.path(OUT, md, paste(t, arm, "promoters.tsv.gz", sep = ".")), sep = "\t") } else
    fwrite(d, file.path(MV, "logs", paste("test", ver, md, t, arm, "tsv", sep = ".")), sep = "\t") }
fwrite(Gt[, .(grp, chr, n, full, n_s, n_region, n_excluded)], file.path(MV, "logs", paste("groups", ver, t, arm, "tsv", sep = ".")), sep = "\t")
cat("DONE\n")
