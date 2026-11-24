## merge_v4.R -- scan files of the v4 root (NSCC): rerun_v2's rows for every promoter outside the delta, plus the delta
## scan (rerun_v4_delta), for the main run (results/final) and the Firth add-on (results/final_firth_addon).
## The delta promoters of a chromosome are those of the delta's main file (it lists every promoter it looked at, skipped
## ones included); their rows are dropped from both rerun_v2 files and the delta's rows are added. A chromosome without
## delta promoters is a symbolic link to the rerun_v2 file. rerun_v2 and rerun_v4_delta are only read.
## Checks: same K, def, tissue, arm, chr in both; no lost work (incomplete = 0) in either; every delta promoter appears
## once in the merged main file. Output files are written atomically and never overwritten.
## Usage: Rscript merge_v4.R <arm> <ncores> [tissue ...]   (default: every tissue of the delta root)
suppressPackageStartupMessages({ library(data.table); library(parallel) })
N <- "/scratch/users/nus/e1327989/nofloor_full"; R2 <- file.path(N, "rerun_v2"); RD <- file.path(N, "rerun_v4_delta"); R4 <- file.path(N, "rerun_v4")
a <- commandArgs(TRUE); ARM <- a[1]; NC <- as.integer(a[2])
TIS <- if (length(a) > 2) a[-(1:2)] else sub("[.]rds$", "", list.files(file.path(RD, "work/final/freeze"), "[.]rds$"))
path <- function(root, r, t, c) file.path(root, "results", r, "freeze", t, ARM, paste0("chr", c, ".rds"))
one <- function(t, c) {
  out <- data.table(tissue = t, chr = c, status = NA_character_, delta = NA_integer_, v2_rows_dropped = NA_integer_)
  fm2 <- path(R2, "final", t, c); fmd <- path(RD, "final", t, c); fa2 <- path(R2, "final_firth_addon", t, c); fad <- path(RD, "final_firth_addon", t, c)
  fmo <- path(R4, "final", t, c); fao <- path(R4, "final_firth_addon", t, c)
  if (file.exists(fmo) && file.exists(fao)) { out$status <- "exists"; return(out) }
  if (!all(file.exists(c(fm2, fmd, fa2, fad)))) { out$status <- "input missing"; return(out) }
  M2 <- readRDS(fm2); MD <- readRDS(fmd); A2 <- readRDS(fa2); AD <- readRDS(fad)
  for (x in list(MD, A2, AD)) stopifnot(identical(x$meta$def, M2$meta$def), identical(x$meta$tissue, t), identical(x$meta$arm, ARM), identical(x$meta$chr, paste0("chr", c)))
  stopifnot(M2$meta$K == MD$meta$K, M2$meta$incomplete == 0, MD$meta$incomplete == 0)
  dp <- if (NROW(MD$prom)) as.character(MD$prom$promoter) else character(0)
  dir.create(dirname(fmo), recursive = TRUE, showWarnings = FALSE); dir.create(dirname(fao), recursive = TRUE, showWarnings = FALSE)
  out$delta <- length(dp)
  if (!length(dp)) { if (!file.exists(fmo)) file.symlink(fm2, fmo); if (!file.exists(fao)) file.symlink(fa2, fao); out$status <- "link"; out$v2_rows_dropped <- 0L; return(out) }
  mg <- function(X, D) { P <- as.data.table(X$prom); S <- as.data.table(X$snp)
    list(prom = rbind(P[!promoter %in% dp], as.data.table(D$prom), fill = TRUE),
         snp = rbind(if (nrow(S)) S[!promoter %in% dp] else S, as.data.table(D$snp), fill = TRUE), dropped = P[promoter %in% dp, .N]) }
  m <- mg(M2, MD); ad <- mg(A2, AD)
  stopifnot(m$prom[promoter %in% dp, .N] == length(dp), !anyDuplicated(m$prom$promoter))
  meta_m <- M2$meta; meta_m$n_promoters <- nrow(m$prom); meta_m$v4_merge <- list(delta_promoters = length(dp), v2_promoters_replaced = m$dropped, delta_meta = MD$meta, date = format(Sys.time()))
  meta_a <- A2$meta; meta_a$v4_merge <- list(delta_promoters = length(dp), delta_meta = AD$meta, date = format(Sys.time()))
  w <- function(obj, f) { tmp <- paste0(f, ".tmp"); saveRDS(obj, tmp); file.rename(tmp, f) }
  if (!file.exists(fmo)) w(list(snp = m$snp, prom = m$prom, meta = meta_m), fmo)
  if (!file.exists(fao)) w(list(snp = ad$snp, prom = ad$prom, meta = meta_a), fao)
  out$status <- "merged"; out$v2_rows_dropped <- m$dropped; out }
jobs <- CJ(t = TIS, c = 1:22)
res <- mclapply(seq_len(nrow(jobs)), function(i) tryCatch(one(jobs$t[i], jobs$c[i]), error = function(e) data.table(tissue = jobs$t[i], chr = jobs$c[i], status = paste("ERROR:", conditionMessage(e)))),
                mc.cores = NC, mc.preschedule = FALSE)
R <- rbindlist(res, fill = TRUE)
fwrite(R, file.path(N, "logs", paste0("merge_v4_", ARM, ".tsv")), sep = "\t")
print(R[, .N, by = status]); cat("delta promoters merged:", sum(R$delta, na.rm = TRUE), "; v2 rows replaced:", sum(R$v2_rows_dropped, na.rm = TRUE), "\n")
if (any(!R$status %in% c("merged", "link", "exists"))) { print(R[!status %in% c("merged", "link", "exists")]); quit(status = 1) }
cat("MERGE PASS", ARM, "\n")
