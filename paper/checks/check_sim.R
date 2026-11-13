## run_paper.sh 02.2 and 02.3 (authored for the repository): re-run simulation and fits against the stored ones.
## Usage: Rscript check_sim.R gen <new sim folder> <stored sim folder> <status.tsv>
##        Rscript check_sim.R fit <new rds> <stored rds> <status.tsv>
a <- commandArgs(TRUE)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
if (a[1] == "gen") {
  n <- new.env(); o <- new.env()
  load(file.path(a[2], "sim_data.RData"), envir = n); load(file.path(a[3], "sim_data.RData"), envir = o)
  for (v in ls(o)) add(paste("simulated data:", v), identical(n[[v]], o[[v]]),
                       sprintf("identical to the stored object: %s (%s)", identical(n[[v]], o[[v]]), paste(dim(as.matrix(o[[v]])), collapse = " x ")))
  r <- compare_tables(fread(file.path(a[2], "ground_truth.txt")), fread(file.path(a[3], "ground_truth.txt")), "gene_id", tol = 0)
  add("simulated data: ground truth (true effects of the 2,000 genes)", r$ok, r$detail)
} else {
  N <- readRDS(a[2]); O <- readRDS(a[3])
  hn <- as.data.table(N$header); ho <- as.data.table(O$header)[id %in% hn$id]
  r <- compare_tables(hn, ho, "id", tol = 0); add("fits: promoter rules (skip, family, counts, eligibility)", r$ok, paste(nrow(hn), "promoters of the chunk;", r$detail))
  rn <- as.data.table(N$res); ro <- as.data.table(O$res)[id %in% hn$id]
  r <- compare_tables(rn, ro, "id", tol = 1e-6); add("fits: BOHR and linear-test statistics", r$ok, r$detail)
}
finish(a[4])
