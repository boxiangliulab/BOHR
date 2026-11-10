## Builds tests/testthat/fixtures/conditional_fixture.rds: the conditional analysis of the ORIGINAL pipeline (the transcribed
## run_conditional_v3c.R / merge_conditional_v3c.R steps of ref_conditional(), on the original engine's functions) for the
## significant units of the synthetic tissue sim_cond() (tests/testthat/helper-conditional.R), whose cis units come from
## bohr_cis(). The regression test in test-conditional.R compares bohr_conditional() with it everywhere (no engine needed).
## Run from the package root with the package installed and the original engine available (see helper-engine.R):
##   Rscript data-raw/make_conditional_fixture.R [lib]
a <- commandArgs(TRUE)
if (length(a)) library(BOHR, lib.loc = a[1]) else library(BOHR)
for (f in c("tests/testthat/helper-engine.R", "tests/testthat/helper-cis.R", "tests/testthat/helper-conditional.R")) source(f)
stopifnot(!is.null(engine_dir()))
eng <- load_engine()
x <- sim_cond()
cis <- bohr_cis(x$phenotype, x$genotype, x$variants, x$promoters, covariates = x$covariates, min_n_mid = 60)
U <- cis$units[cis$units$significant, , drop = FALSE]
units <- data.frame(promoter = U$promoter, type = U$component, leader = U$leader, leader_logp = U$leader_log_p, p = U$p, stringsAsFactors = FALSE)
R <- ref_conditional(eng, x, units, cis$threshold)
R$meta <- list(engine = engine_dir(), threshold = cis$threshold, created = format(Sys.Date()), R = R.version.string)
saveRDS(R, "tests/testthat/fixtures/conditional_fixture.rds", compress = "xz")
cat("fixture:", format(file.size("tests/testthat/fixtures/conditional_fixture.rds"), big.mark = ","), "bytes;", nrow(R$signals), "signals,",
    nrow(R$loo), "leave-one-out rows\n")
