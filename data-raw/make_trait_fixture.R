## Builds tests/testthat/fixtures/trait_fixture.rds: the trait regression of the ORIGINAL code (ref_trait(): the transcribed
## run_trait.R, data.R, summarize.R, mirror rule v3 and option-B steps, on the original engine and data.R functions) for the
## synthetic donors of sim_trait() (tests/testthat/helper-trait.R), traits SEX, AGE, BMI and MHT2D, with absolute activity.
## The regression test in test-trait.R compares bohr_trait() with it everywhere (no original code needed).
## Run from the package root with the original code available (see helper-engine.R and helper-trait.R):
##   Rscript data-raw/make_trait_fixture.R [lib]
a <- commandArgs(TRUE)
if (length(a)) library(BOHR, lib.loc = a[1]) else library(BOHR)
for (f in c("tests/testthat/helper-engine.R", "tests/testthat/helper-trait.R")) source(f)
stopifnot(!is.null(engine_dir()), !is.null(trait_code_dir()))
x <- sim_trait()
R <- ref_trait(load_trait_code(), x, c("SEX", "AGE", "BMI", "MHT2D"), A = x$absolute)
R$units <- R$units[, c("trait", "promoter", "type", "p", "est", "se", "pe", "final", "q_final", "sig_final", "abs_relabel")]
R$meta <- list(engine = engine_dir(), trait_code = trait_code_dir(), created = format(Sys.Date()), R = R.version.string)
saveRDS(R, "tests/testthat/fixtures/trait_fixture.rds", compress = "xz")
cat("fixture:", format(file.size("tests/testthat/fixtures/trait_fixture.rds"), big.mark = ","), "bytes;", nrow(R$units), "units\n")
