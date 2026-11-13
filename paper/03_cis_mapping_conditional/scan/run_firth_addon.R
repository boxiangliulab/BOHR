## Firth add-on (author, 2026-10-01): the type-proportion thresholds were a false-positive filter of the old pipeline,
## so every promoter with a zero and/or one part gets the Firth test and the thresholds become a summary-time choice.
## For each promoter of the main run's chromosome file that has a zero/one part (family BEINF/BEZI/BEOI, not skipped)
## but no Firth yet (neither elig_p0 nor elig_p1), the main run's inputs are rebuilt with the same code as run_final2.R
## (prepare_y -> orth_design -> promoter_null; read_geno / cis_index / MIN_GENO) and firth_snp_tests (bohr2_addon3.R)
## is run on exactly the SNPs the main run tested. The main run's files are never modified; the plain LRT / Wald /
## estimates of these SNPs stay in the main files, so before/after-Firth comparisons are possible for every SNP.
## FIRTH_ADDON_CHECK=<k>: recompute the first k promoters that already HAVE Firth and report the largest difference to
## the stored values (reproducibility check); output then goes to work/firth_addon_check/.
## Usage: Rscript run_firth_addon.R <def> <tissue> <arm> <chr> <ncores>
## Output (atomic): results/final_firth_addon/<def>/<tissue>/<arm>/<chr>.rds = list(snp, prom, meta)
suppressPackageStartupMessages(library(parallel))
CODE <- Sys.getenv("BOHR2_CODE", "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/code")
source(file.path(CODE, "bohr2_data.R")); source(file.path(CODE, "bohr2_engine.R")); source(file.path(CODE, "bohr2_addon3.R"))
a <- commandArgs(TRUE); def <- a[1]; tis <- a[2]; arm <- a[3]; chr <- a[4]; NC <- as.integer(a[5])
CHECK <- as.integer(Sys.getenv("FIRTH_ADDON_CHECK", "0"))
main_f <- file.path(V2, "results", "final", def, tis, arm, paste0(chr, ".rds"))
out_root <- if (CHECK > 0) file.path(V2, "work", "firth_addon_check") else file.path(V2, "results", "final_firth_addon")
out_dir <- file.path(out_root, def, tis, arm); out <- file.path(out_dir, paste0(chr, ".rds"))
if (file.exists(out)) { cat("exists:", out, "\n"); quit(save = "no") }
if (!file.exists(main_f)) stop("main run file missing: ", main_f)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
t_start <- Sys.time()

M <- readRDS(main_f)
stopifnot(identical(M$meta$def, def), identical(M$meta$tissue, tis), identical(M$meta$arm, arm), identical(M$meta$chr, chr))
PR <- M$prom
has01 <- is.na(PR$skip) & PR$family %in% c("BEINF", "BEZI", "BEOI")
done <- (PR$elig_p0 %in% TRUE) | (PR$elig_p1 %in% TRUE)
todo <- if (CHECK > 0) head(PR$promoter[has01 & done], CHECK) else PR$promoter[has01 & !done]
snp_of <- if (nrow(M$snp)) split(M$snp$snp, M$snp$promoter) else list()

td <- readRDS(file.path(V2, "work", "final", def, paste0(tis, ".rds"))); K <- td$K
stopifnot(K == M$meta$K)
if (arm == "real") { A <- td$act; P <- td$peer } else {
  m <- fread(file.path(V2, "perm_maps", "final", def, paste0(tis, ".", arm, ".perm_map.tsv")))
  stopifnot(setequal(m$slot, td$U), setequal(m$donor, td$U), !anyDuplicated(m$donor))
  d <- setNames(m$donor, m$slot)[td$U]; A <- td$act[d, , drop = FALSE]; P <- td$peer[d, , drop = FALSE]; rownames(A) <- rownames(P) <- td$U }
C_all <- cov_matrix(td, P, K)

FIRTH_COLS <- c("nu_firth_stat", "nu_firth_logp", "nu_firth_est", "nu_firth_se", "nu_firth_conv",
                "tau_firth_stat", "tau_firth_logp", "tau_firth_est", "tau_firth_se", "tau_firth_conv", "firth_rounds")
CHUNK <- 100L
prep_promoter <- function(p) {                       # same steps as run_final2.R::prep_promoter
  t0 <- proc.time()[["elapsed"]]
  z <- prepare_y(A[, p], ncol(C_all) - 1L)
  if (!is.na(z$skip)) stop("skip differs from the main run: ", z$skip)
  s <- names(z$y); y <- unname(z$y); Cs <- orth_design(C_all[s, , drop = FALSE])
  nul <- promoter_null(y, Cs, z$family)
  pos <- td$info[p, "start"]; geno <- read_geno(chr, pos - FLANK, pos + FLANK - 1L, td$U)
  idx <- cis_index(geno, pos); G <- geno$G[s, idx, drop = FALSE]
  cnt <- apply(G, 2, function(x) tabulate(x + 1L, nbins = 3L)); keep <- which(apply(cnt, 2, min) >= MIN_GENO)
  want <- snp_of[[p]]; if (is.null(want)) want <- character(0)
  if (!setequal(colnames(G)[keep], want)) stop("tested SNP set differs from the main run")
  list(promoter = p, family = z$family, y = y, Cs = Cs, nul = nul, G = G[, keep, drop = FALSE], sec = proc.time()[["elapsed"]] - t0)
}
prep_one <- function(p) tryCatch(prep_promoter(p), error = function(e) list(promoter = p, err = conditionMessage(e)))
prep <- mclapply(todo, prep_one, mc.cores = NC, mc.preschedule = FALSE)
dead <- which(!vapply(prep, is.list, TRUE)); for (i in dead) prep[[i]] <- prep_one(todo[i])
incomplete <- sum(!vapply(prep, is.list, TRUE))
prep <- prep[vapply(prep, is.list, TRUE)]
tasks <- list()
for (i in seq_along(prep)) if (!is.null(prep[[i]]$G) && ncol(prep[[i]]$G) > 0) {
  m <- ncol(prep[[i]]$G); w <- c(BEINF = 3, BEZI = 1, BEOI = 1)[[prep[[i]]$family]] * length(prep[[i]]$y) * ncol(prep[[i]]$Cs)^1.5
  for (ch in split(seq_len(m), ceiling(seq_len(m) / CHUNK))) tasks[[length(tasks) + 1]] <- list(i = i, j = ch, w = w * length(ch))
}
tasks <- tasks[order(-vapply(tasks, `[[`, 0, "w"))]
run_chunk <- function(t) {
  P <- prep[[t$i]]; t0 <- proc.time()[["elapsed"]]; fe <- 0L
  rows <- lapply(t$j, function(j) {
    ff <- tryCatch(firth_snp_tests(P$y, P$Cs, P$G[, j], P$nul), error = function(e) NULL)
    f <- setNames(rep(list(NA), length(FIRTH_COLS)), FIRTH_COLS); err <- is.null(ff)
    if (err) { fe <<- fe + 1L; f$nu_firth_conv <- FALSE; f$tau_firth_conv <- FALSE } else f <- ff[FIRTH_COLS]
    c(list(snp = colnames(P$G)[j], firth_error = err), f) })
  list(i = t$i, n = length(t$j), ferr = fe, sec = proc.time()[["elapsed"]] - t0, snp = rbindlist(rows, fill = TRUE))
}
chunk_one <- function(t) tryCatch(run_chunk(t), error = function(e) list(i = t$i, n = 0L, ferr = 0L, sec = NA_real_, snp = NULL, err = conditionMessage(e)))
res <- mclapply(tasks, chunk_one, mc.cores = NC, mc.preschedule = FALSE)
dead <- which(!vapply(res, function(r) is.list(r) && !is.null(r$snp), TRUE))
if (length(dead)) res[dead] <- mclapply(tasks[dead], chunk_one, mc.cores = min(NC, length(dead)), mc.preschedule = FALSE)
incomplete <- incomplete + sum(!vapply(res, function(r) is.list(r) && !is.null(r$snp), TRUE))
res <- res[vapply(res, function(r) is.list(r) && !is.null(r$snp), TRUE)]
snp <- rbindlist(lapply(res, function(r) cbind(promoter = prep[[r$i]]$promoter, r$snp)), fill = TRUE)
prom <- rbindlist(lapply(seq_along(prep), function(i) { P <- prep[[i]]; mine <- res[vapply(res, function(r) r$i == i, TRUE)]
  data.table(promoter = P$promoter, family = if (is.null(P$family)) NA_character_ else P$family,
             n_snp = if (is.null(P$G)) NA_integer_ else ncol(P$G), n_firth = sum(vapply(mine, `[[`, 0L, "n")),
             n_firth_error = sum(vapply(mine, `[[`, 0L, "ferr")), err = if (is.null(P$err)) NA_character_ else P$err,
             sec = (if (is.null(P$sec)) 0 else P$sec) + sum(vapply(mine, function(r) if (is.finite(r$sec)) r$sec else 0, 0))) }), fill = TRUE)
## bookkeeping: every promoter prepared without error, and every tested SNP has a Firth row
if (nrow(prom)) incomplete <- incomplete + prom[!is.na(err) | n_firth != n_snp, .N]
meta <- list(def = def, tissue = tis, arm = arm, chr = chr, K = K, n_todo = length(todo), check = CHECK, cores = NC,
             wall_sec = as.numeric(difftime(Sys.time(), t_start, units = "secs")), host = Sys.info()[["nodename"]], date = format(Sys.time()),
             engine = "bohr2_addon3.R firth_snp_tests (alternating full/restricted search) on promoters not covered by the main run")
meta$incomplete <- incomplete
if (CHECK > 0 && nrow(snp)) {                        # reproducibility against the stored main-run values
  ref <- M$snp[promoter %in% todo, c("promoter", "snp", FIRTH_COLS), with = FALSE]
  cmp <- merge(snp, ref, by = c("promoter", "snp"), suffixes = c("", ".main"))
  for (k in c("nu_firth_stat", "tau_firth_stat", "nu_firth_logp", "tau_firth_logp", "nu_firth_est", "tau_firth_est")) {
    x <- cmp[[k]]; y <- cmp[[paste0(k, ".main")]]; ok <- is.finite(x) & is.finite(y)
    cat(sprintf("CHECK %-15s matched %6d  max|diff| %.3g  NA mismatch %d\n", k, sum(ok), if (any(ok)) max(abs(x[ok] - y[ok])) else NA, sum(is.finite(x) != is.finite(y)))) } }
if (incomplete > 0) { saveRDS(list(snp = snp, prom = prom, meta = meta), sub("\\.rds$", ".incomplete.rds", out))
  cat("INCOMPLETE:", incomplete, "promoters/chunks lost after retry; not published\n"); quit(save = "no", status = 0) }
tmp <- paste0(out, ".tmp"); saveRDS(list(snp = snp, prom = prom, meta = meta), tmp); file.rename(tmp, out)
cat(sprintf("%s %s K%d %s %s: firth add-on %d promoters, %d SNP rows, %.1f min\n", def, tis, K, arm, chr, length(todo), nrow(snp), meta$wall_sec / 60))
