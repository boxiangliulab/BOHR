## Inputs as objects or files, the tissue loop, out_dir and the VCF reader of bohr_cis().

small_sim <- function(seed = 31, out_dir = NULL, n_tissues = 2) with_seed(seed, simulate_bohr(n_individuals = 90, n_tissues = n_tissues,
                                                                                               n_genes = 4, variant_spacing = 1e5, out_dir = out_dir))
cis_tabs <- function(r) r[c("snp", "units", "promoters", "mirror", "threshold")]
strip_tissue <- function(d) { if (is.null(d)) return(d); d$tissue <- NULL; rownames(d) <- NULL; d }
run_cis <- function(...) suppressWarnings(bohr_cis(..., window = 2e5))

test_that("TSV files are read the same with data.table and with read.delim", {
  d <- withr::local_tempdir()
  sim <- small_sim(out_dir = d)
  f <- file.path(d, c("relative_activity.tsv", "genotypes.tsv", "covariates/T1.tsv"))
  a <- lapply(f, BOHR:::read_matrix, what = "x")
  testthat::local_mocked_bindings(has_pkg = function(pkg) pkg != "data.table")
  b <- lapply(f, BOHR:::read_matrix, what = "x")
  expect_identical(a, b)
  expect_identical(a[[1]], sim$relative); expect_identical(a[[2]], sim$genotype); expect_identical(a[[3]], sim$covariates$T1)
  expect_identical(BOHR:::read_table(file.path(d, "tss.tsv"), "p")$id, sim$tss$id)
})

test_that("bohr_cis() gives the same results from files and from R objects", {
  d <- withr::local_tempdir()
  sim <- small_sim(out_dir = d)
  ro <- run_cis(sim$relative, sim$genotype, sim$variants, sim$tss, covariates = sim$covariates, tissue = sim$samples)
  rf <- run_cis(file.path(d, "relative_activity.tsv"), genotype = file.path(d, "genotypes.tsv"), variants = file.path(d, "variants.tsv"),
                promoters = file.path(d, "tss.tsv"), covariates = file.path(d, "covariates"), tissue = file.path(d, "samples.tsv"))
  expect_identical(cis_tabs(rf), cis_tabs(ro))
  ## one tissue from files: covariates as one file, phenotype restricted to that tissue's samples
  s1 <- sim$samples$sample[sim$samples$tissue == "T1"]
  ph1 <- sim$relative[, s1]; G1 <- sim$genotype[, sim$samples$individual[match(s1, sim$samples$sample)]]; colnames(G1) <- s1
  BOHR:::write_tsv(data.frame(id = rownames(ph1), ph1, check.names = FALSE), file.path(d, "ph1.tsv"))
  BOHR:::write_tsv(data.frame(id = rownames(G1), G1, check.names = FALSE), file.path(d, "g1.tsv"))
  r1o <- run_cis(ph1, G1, sim$variants, sim$tss, covariates = sim$covariates$T1)
  r1f <- run_cis(file.path(d, "ph1.tsv"), file.path(d, "g1.tsv"), file.path(d, "variants.tsv"), file.path(d, "tss.tsv"),
                 covariates = file.path(d, "covariates", "T1.tsv"))
  expect_identical(cis_tabs(r1f), cis_tabs(r1o))
  expect_error(run_cis(file.path(d, "nope.tsv"), G1, sim$variants, sim$tss), "file not found")
  expect_error(run_cis(ph1, G1, NULL, sim$tss), "`variants` is needed")
  expect_error(run_cis(ph1, G1, sim$variants, sim$tss, covariates = file.path(d, "covariates")), "is a folder")
})

test_that("the tissue loop maps each tissue on its own samples, with the pooled mean for mirror pairs", {
  sim <- small_sim(seed = 32)
  r <- run_cis(sim$relative, sim$genotype, sim$variants, sim$tss, covariates = sim$covariates, tissue = sim$samples)
  expect_s3_class(r, "bohr_cis")
  expect_identical(r$settings$tissues, c("T1", "T2"))
  for (k in c("snp", "units", "promoters", "mirror")) expect_identical(names(r[[k]])[1], "tissue")
  expect_named(r$threshold, c("T1", "T2"))
  pooled <- rowMeans(sim$relative, na.rm = TRUE)
  for (tt in c("T1", "T2")) {
    s <- sim$samples[sim$samples$tissue == tt, ]
    G <- sim$genotype[, s$individual]; colnames(G) <- s$sample
    one <- run_cis(sim$relative[, s$sample], G, sim$variants, sim$tss, covariates = sim$covariates[[tt]], major_score = pooled)
    for (k in c("snp", "units", "promoters", "mirror")) expect_identical(strip_tissue(r[[k]][r[[k]]$tissue == tt, , drop = FALSE]), one[[k]], info = paste(tt, k))
    expect_identical(r$threshold[[tt]], one$threshold)
    expect_identical(r$settings$n_samples[[tt]], nrow(s))
  }
  ## q values stay within each tissue
  U <- r$units[r$units$tissue == "T1", ]; prim <- U$estimable & !U$dedup_drop
  expect_equal(U$q[prim], stats::p.adjust(U$p[prim], "BH"))
  expect_output(print(r), "BOHR cis mapping: 2 tissues")
  expect_output(print(r), "q values within each tissue")
})

test_that("one tissue label for every sample equals the single-tissue call", {
  sim <- small_sim(seed = 33, n_tissues = 1)
  a <- run_cis(sim$relative, sim$genotype, sim$variants, sim$tss, covariates = sim$covariates$T1,
               individual = stats::setNames(sim$samples$individual, sim$samples$sample))
  b <- run_cis(sim$relative, sim$genotype, sim$variants, sim$tss, covariates = sim$covariates, tissue = sim$samples)
  for (k in c("snp", "units", "promoters", "mirror")) expect_identical(strip_tissue(b[[k]]), a[[k]], info = k)
  expect_identical(unname(b$threshold), a$threshold)
})

test_that("tissue labels, individuals and covariates can be given in several forms", {
  sim <- small_sim(seed = 34)
  ref <- run_cis(sim$relative, sim$genotype, sim$variants, sim$tss, covariates = sim$covariates, tissue = sim$samples)
  s <- sim$samples
  ## vectors in column order or named by sample, and one covariate table with sample row names
  allcov <- do.call(rbind, unname(sim$covariates))
  v1 <- run_cis(sim$relative, sim$genotype, sim$variants, sim$tss, covariates = allcov, tissue = s$tissue, individual = s$individual)
  v2 <- run_cis(sim$relative, sim$genotype, sim$variants, sim$tss, covariates = allcov[rev(rownames(allcov)), ],
                tissue = stats::setNames(rev(s$tissue), rev(s$sample)), individual = stats::setNames(s$individual, s$sample))
  expect_identical(cis_tabs(v1), cis_tabs(ref)); expect_identical(cis_tabs(v2), cis_tabs(ref))
  ## a factor gives the order of the tissues
  v3 <- run_cis(sim$relative, sim$genotype, sim$variants, sim$tss, covariates = sim$covariates, tissue = factor(s$tissue, c("T2", "T1")),
                individual = s$individual)
  expect_identical(v3$settings$tissues, c("T2", "T1")); expect_identical(unique(v3$units$tissue), c("T2", "T1"))
  ## errors
  expect_error(run_cis(sim$relative, sim$genotype, sim$variants, sim$tss, covariates = sim$covariates["T1"], tissue = s),
               "no element for tissue T2")
  expect_error(run_cis(sim$relative, sim$genotype, sim$variants, sim$tss, covariates = unname(allcov), tissue = s), "row names")
  expect_error(run_cis(sim$relative, sim$genotype, sim$variants, sim$tss, tissue = s$tissue[-1]), "one value per column")
  expect_error(run_cis(sim$relative, sim$genotype, sim$variants, sim$tss, tissue = s[-1, ]), "no row for 1 sample")
  expect_error(run_cis(sim$relative, sim$genotype, sim$variants, sim$tss, tissue = s[, c("sample", "individual")]), "columns sample and tissue")
  expect_error(run_cis(sim$relative, sim$genotype, sim$variants, sim$tss, tissue = s$tissue), "no column for .* sample\\(s\\) of `phenotype`")
  expect_error(run_cis(sim$relative, sim$genotype[, -1], sim$variants, sim$tss, tissue = s), "no column for 1 individual\\(s\\)")
  d <- withr::local_tempdir(); dir.create(file.path(d, "cv"))
  expect_error(run_cis(sim$relative, sim$genotype, sim$variants, sim$tss, covariates = file.path(d, "cv"), tissue = s), "No covariate file for tissue T1")
})

test_that("out_dir writes the four tables, which read back as the returned ones", {
  d <- withr::local_tempdir()
  sim <- small_sim(seed = 35)
  r <- run_cis(sim$relative, sim$genotype, sim$variants, sim$tss, covariates = sim$covariates, tissue = sim$samples,
               out_dir = file.path(d, "res"))
  expect_setequal(list.files(file.path(d, "res")), c("snp.tsv", "units.tsv", "promoters.tsv", "mirror.tsv"))
  u <- utils::read.delim(file.path(d, "res", "units.tsv"), stringsAsFactors = FALSE)
  expect_identical(names(u), names(r$units)); expect_identical(nrow(u), nrow(r$units))
  expect_identical(u$promoter, r$units$promoter); expect_equal(u$q, r$units$q, tolerance = 1e-12)
  sn <- utils::read.delim(file.path(d, "res", "snp.tsv"), stringsAsFactors = FALSE)
  expect_equal(sn$log_p_p0, r$snp$log_p_p0, tolerance = 1e-12); expect_identical(sn$variant, r$snp$variant)
  expect_error(run_cis(sim$relative, sim$genotype, sim$variants, sim$tss, out_dir = c("a", "b")), "out_dir")
})

test_that("bohr_combine() recomputes q values within each tissue", {
  sim <- small_sim(seed = 36)
  args <- list(sim$relative, sim$genotype, sim$variants, sim$tss, covariates = sim$covariates, tissue = sim$samples)
  all <- do.call(run_cis, args)
  a <- do.call(run_cis, c(args, chromosomes = "chr1")); b <- do.call(run_cis, c(args, chromosomes = "chr2"))
  cb <- bohr_combine(a, b)
  o <- function(u) u[order(u$tissue, u$promoter, u$component), , drop = FALSE]
  expect_equal(o(cb$units), o(all$units), ignore_attr = TRUE)
  expect_identical(cb$threshold, all$threshold)
  expect_error(bohr_combine(a, a), "overlap")
})

test_that("read_vcf_dosage() asks for VariantAnnotation when it is missing", {
  d <- withr::local_tempdir()
  sim <- small_sim(seed = 37, out_dir = d)
  testthat::local_mocked_bindings(has_pkg = function(pkg) FALSE)
  expect_error(read_vcf_dosage(file.path(d, "genotypes.vcf")), "needs the Bioconductor packages VariantAnnotation")
  expect_error(read_vcf_dosage(file.path(d, "genotypes.vcf")), "or pass a dosage matrix")
  expect_error(run_cis(sim$relative, file.path(d, "genotypes.vcf"), NULL, sim$tss, tissue = sim$samples), "VariantAnnotation")
  expect_error(read_vcf_dosage(file.path(d, "none.vcf")), "not found")
})

test_that("the VCF written by simulate_bohr() reads back as the dosages, and bohr_cis() takes it", {
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("SummarizedExperiment")
  d <- withr::local_tempdir()
  sim <- small_sim(seed = 38, out_dir = d)
  v <- read_vcf_dosage(file.path(d, "genotypes.vcf"))
  expect_identical(v$genotype, sim$genotype)
  expect_identical(v$variants$id, sim$variants$id)
  expect_identical(v$variants$end, as.numeric(sim$variants$end)); expect_identical(v$variants$chr, sim$variants$chr)
  some <- sim$samples$individual[1:5]
  expect_identical(read_vcf_dosage(file.path(d, "genotypes.vcf"), samples = some)$genotype, sim$genotype[, some])
  expect_error(read_vcf_dosage(file.path(d, "genotypes.vcf"), samples = "nobody"), "no column for 1 sample")
  rv <- run_cis(sim$relative, file.path(d, "genotypes.vcf"), NULL, sim$tss, covariates = sim$covariates, tissue = sim$samples)
  rt <- run_cis(sim$relative, sim$genotype, sim$variants, sim$tss, covariates = sim$covariates, tissue = sim$samples)
  expect_identical(cis_tabs(rv), cis_tabs(rt))
  expect_error(run_cis(sim$relative, file.path(d, "genotypes.vcf"), sim$variants, sim$tss, tissue = sim$samples), "must be NULL")
})
