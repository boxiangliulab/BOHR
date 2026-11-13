## run_paper.sh 00.1 (authored for the repository): the BOHR package installed into OUT/lib, and the R packages the
## re-runs need. Usage: Rscript check_setup.R <lib> <status.tsv>
a <- commandArgs(TRUE)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
ok <- suppressWarnings(requireNamespace("BOHR", lib.loc = a[1], quietly = TRUE))
add("BOHR package installed", ok, if (ok) paste("version", as.character(packageVersion("BOHR", lib.loc = a[1])), "in OUT/lib") else "not installed (see the install output in this log)")
need <- c("data.table", "Matrix", "VariantAnnotation", "coloc", "gamlss", "gamlss.dist", "truncnorm", "parallel", "SummarizedExperiment")
have <- vapply(need, requireNamespace, TRUE, quietly = TRUE)
add("R packages of the re-runs", all(have), if (all(have)) paste(sprintf("%s %s", need, vapply(need, function(p) as.character(packageVersion(p)), "")), collapse = ", ")
    else paste("missing:", paste(need[!have], collapse = ", ")))
add("R version", TRUE, R.version.string)
finish(a[2])
