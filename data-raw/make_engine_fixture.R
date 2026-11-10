## Builds tests/testthat/fixtures/engine_fixture.rds: small synthetic inputs and the outputs of the ORIGINAL
## analysis engine of the BOHR paper (not of this package), for the regression test test-regression.R.
## Needs the engine (BOHR_ENGINE_DIR or the default storage path) and the package's test helpers.
## Run from the package root:  Rscript data-raw/make_engine_fixture.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
for (f in c("tests/testthat/helper-engine.R", "tests/testthat/helper-cis.R")) source(f)
stopifnot(!is.null(engine_dir()))
eng <- load_engine()
## a small tissue: 80 samples, 5 promoters on chr1 (BEINF, an exact mirror pair, BE, BEZI), variants every 125 kb
x <- sim_tissue(seed = 21, n = 80, spacing = 125000, span = 1.3e6)
keep <- c("PromoterID_101", "PromoterID_102", "PromoterID_103", "PromoterID_106", "PromoterID_107")
x$phenotype <- x$phenotype[keep, ]; x$promoters <- x$promoters[x$promoters$id %in% keep, ]
on <- x$variants$chr == "chr1"; x$variants <- x$variants[on, ]; x$genotype <- x$genotype[on, ]
ref <- ref_tissue(eng, x)
## single promoters: the four families and a separated case
cases <- parity_cases()[c("BEINF_cov", "BEZI_cov", "BEOI_nocov", "BE_cov", "sep_beinf")]
single <- lapply(cases, function(cs) { r <- ref_promoter(eng, cs$y, cs$Z, cs$g)
  list(input = cs, family = r$family, engine = c(r$r[BOHR:::SNP_COLS], if (is.null(r[["f"]])) NULL else r[["f"]][BOHR:::FIRTH_COLS], list(mu_F_logp = r$mu_F_logp))) })
md5 <- tools::md5sum(file.path(engine_dir(), c("bohr2_engine.R", "bohr2_addon3.R", "bohr2_data.R", "run_muF.R", "discovery_nogate.R")))
fx <- list(tissue = x, tissue_engine = ref, single = single,
           meta = list(created = format(Sys.Date()), engine_md5 = setNames(unname(md5), basename(names(md5))), R = R.version.string,
                       note = "Expected values computed with the original engine (tests/testthat/helper-*.R reference pipelines)."))
saveRDS(fx, "tests/testthat/fixtures/engine_fixture.rds", compress = "xz")
cat("fixture:", nrow(ref$snp), "SNP rows,", nrow(ref$units), "units,", length(single), "single-promoter cases;",
    format(file.size("tests/testthat/fixtures/engine_fixture.rds"), big.mark = ","), "bytes\n")
