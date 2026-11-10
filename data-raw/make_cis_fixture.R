## Builds tests/testthat/fixtures/cis_ae071ae.rds: the bohr_cis() tables of commit ae071ae (the last version before the
## stage-2b tissue loop and file inputs) on the synthetic configurations of cis_v0_configs() (tests/testthat/helper-cis.R),
## for the single-tissue parity test in test-cis-parity.R. Run it with ae071ae installed, from the current package root:
##   git archive ae071ae | tar -x -C <src>; R CMD INSTALL -l <lib> <src>
##   Rscript data-raw/make_cis_fixture.R <lib>
a <- commandArgs(TRUE)
library(BOHR, lib.loc = a[1])
stopifnot(!"tissue" %in% names(formals(bohr_cis)))                       # the old version
for (f in c("tests/testthat/helper-engine.R", "tests/testthat/helper-cis.R")) source(f)
cfg <- cis_v0_configs()
out <- lapply(cfg, function(args) suppressWarnings(do.call(bohr_cis, args))[c("snp", "units", "promoters", "mirror", "threshold")])
out$meta <- list(commit = "ae071ae", created = format(Sys.Date()), R = R.version.string)
saveRDS(out, "tests/testthat/fixtures/cis_ae071ae.rds", compress = "xz")
cat("fixture:", format(file.size("tests/testthat/fixtures/cis_ae071ae.rds"), big.mark = ","), "bytes\n")
