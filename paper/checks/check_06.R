## run_paper.sh 06.1-06.3 (authored for the repository): re-run major/minor-promoter analysis against the stored outputs.
## Usage: Rscript check_06.R encode <new data folder> <stored stage in/fig5_v3data> <status.tsv>
##        Rscript check_06.R stage  <new out folder>  <stored out folder>           <status.tsv>
a <- commandArgs(TRUE)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
if (a[1] == "encode") {
  fs <- c("alltested_tss_elements.tsv", paste0("zscores_", c("dnase", "h3k27ac", "h3k4me3"), "_alltested.tsv.gz"),
          paste0("zscores_", c("dnase", "h3k27ac", "h3k4me3"), "_alltested_columns.tsv"))
  for (f in fs) {
    n <- file.path(a[2], f); o <- file.path(a[3], f)
    if (!file.exists(n)) { add(paste("ENCODE input", f), FALSE, "not written (see the log)"); next }
    same <- unname(tools::md5sum(n)) == unname(tools::md5sum(o))
    if (same) { add(paste("ENCODE input", f), TRUE, "byte-identical to the stage input"); next }
    x <- fr(n); y <- fr(o); x[, .row := .I]; y[, .row := .I]
    r <- compare_tables(x, y, ".row", tol = 1e-6); add(paste("ENCODE input", f), r$ok, paste("not byte-identical;", r$detail))
  }
} else compare_folders(a[2], a[3], "table")
finish(a[4])
