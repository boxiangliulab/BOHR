## The synthetic input set of simulate_bohr().

test_that("simulate_bohr() returns a consistent input set", {
  sim <- with_seed(51, simulate_bohr(n_individuals = 80, n_tissues = 2, n_genes = 8))
  expect_s3_class(sim, "bohr_sim")
  expect_named(sim, c("relative", "absolute", "samples", "tss", "genotype", "variants", "covariates", "traits", "truth", "files"))
  S <- sim$samples; P <- sim$tss
  expect_identical(colnames(sim$relative), S$sample); expect_identical(colnames(sim$absolute), S$sample)
  expect_false(anyDuplicated(S$sample) > 0)
  expect_true(all(S$individual %in% colnames(sim$genotype)))
  expect_identical(rownames(sim$relative), P$id[P$tested]); expect_identical(rownames(sim$absolute), P$id)
  expect_true(all(sim$relative >= 0 & sim$relative <= 1, na.rm = TRUE))
  expect_true(any(is.na(sim$relative)))
  ## relative activity = absolute activity over the gene's external promoters (internal ones excluded)
  ext <- !P$internal
  den <- rowsum(sim$absolute[ext, ], P$gene[ext], reorder = FALSE)
  rel <- sim$absolute[P$tested, ] / den[P$gene[P$tested], ]
  rel[!is.finite(rel)] <- NA
  expect_equal(sim$relative, rel, tolerance = 1e-11, ignore_attr = TRUE)
  ## gene layouts: the four types, internal promoters, untested promoters in the denominator
  ty <- tapply(P$tested, P$gene, sum)
  expect_setequal(unique(ty), c(2L, 3L))
  expect_true(any(P$internal)); expect_true(any(!P$internal & !P$tested))
  mp <- mirror_pairs(sim$relative, group = P$gene[P$tested])
  gt <- sim$truth$gene_type[match(mp$group, sim$truth$gene)]
  kind <- ifelse(mp$interior_exact, "exact", ifelse(mp$boundary_exact, "boundary_only", "nonexact"))
  types <- rep(c("exact", "three", "boundary_only", "nonexact"), length.out = 8)
  expect_identical(kind, types[match(mp$group, sprintf("G%03d", 1:8))])
  ## genotypes, variants, covariates and traits
  expect_identical(rownames(sim$genotype), sim$variants$id)
  expect_true(all(sim$genotype %in% c(0L, 1L, 2L, NA)))
  expect_true(any(sim$variants$end > sim$variants$pos))
  expect_named(sim$covariates, c("T1", "T2"))
  expect_identical(rownames(sim$covariates$T2), S$sample[S$tissue == "T2"])
  expect_identical(sort(unique(sim$traits$binary_trait)), 0:1)
  expect_null(sim$files)
  expect_output(print(sim), "BOHR synthetic input set")
})

test_that("the truth table follows the gene layouts", {
  sim <- with_seed(52, simulate_bohr(n_individuals = 60, n_tissues = 2, n_genes = 8, prop_causal = 1))
  tr <- sim$truth
  expect_true(all(tr$promoter %in% sim$tss$id[sim$tss$tested]))
  expect_setequal(unique(tr$predictor_type), c("dosage", "binary", "continuous"))
  expect_true(all(tr$predictor[tr$predictor_type == "dosage"] %in% sim$variants$id))
  ## an exact mirror pair: p0 of the partner is p1 of the primary promoter and the reverse; mu changes sign
  ex <- tr[tr$gene_type == "exact" & tr$predictor_type == "dosage", ]
  a <- ex[ex$role == "primary", ]; b <- ex[ex$role == "mirror", ]
  b <- b[match(paste(a$tissue, a$gene), paste(b$tissue, b$gene)), ]
  expect_equal(b$effect_p0, a$effect_p1); expect_equal(b$effect_p1, a$effect_p0); expect_equal(b$effect_mu, -a$effect_mu)
  ## a variant effect on p0 of the primary promoter is twice the per-allele effect (dosage 2 versus 0)
  p0 <- tr[tr$role == "primary" & tr$predictor_type == "dosage" & tr$effect_p0 > 0 & tr$effect_mu == 0, ]
  expect_true(all(abs(p0$effect_p0 - 0.4) < 1e-12))
  ## effects in the first tissue only for alternate causal genes
  expect_true(any(tapply(tr$tissue, tr$gene, function(v) length(unique(v))) == 1L))
  expect_identical(with_seed(52, simulate_bohr(n_individuals = 60, n_tissues = 2, n_genes = 8, prop_causal = 1))$truth, tr)
  z <- with_seed(52, simulate_bohr(n_individuals = 60, n_genes = 4, prop_causal = 0))
  expect_false(any(z$truth$predictor_type == "dosage"))
})

test_that("out_dir writes every input in the formats bohr_cis() reads", {
  d <- withr::local_tempdir()
  sim <- with_seed(53, simulate_bohr(n_individuals = 50, n_tissues = 2, n_genes = 4, out_dir = d))
  expect_true(all(file.exists(sim$files)))
  expect_setequal(list.files(d), c("relative_activity.tsv", "absolute_activity.tsv", "samples.tsv", "tss.tsv", "genotypes.tsv",
                                   "genotypes.vcf", "variants.tsv", "traits.tsv", "truth.tsv", "covariates"))
  expect_setequal(list.files(file.path(d, "covariates")), c("T1.tsv", "T2.tsv"))
  expect_identical(BOHR:::read_matrix(sim$files[["relative"]], "x"), sim$relative)
  expect_identical(BOHR:::read_matrix(sim$files[["absolute"]], "x"), sim$absolute)
  expect_identical(BOHR:::read_matrix(sim$files[["genotype"]], "x"), sim$genotype)
  expect_identical(BOHR:::read_matrix(file.path(d, "covariates", "T2.tsv"), "x"), sim$covariates$T2)
  smp <- utils::read.delim(sim$files[["samples"]], stringsAsFactors = FALSE)
  expect_identical(smp, sim$samples)
  tss <- utils::read.delim(sim$files[["tss"]], stringsAsFactors = FALSE)
  expect_identical(tss$pos, sim$tss$pos); expect_identical(tss$internal, sim$tss$internal)
  ## the VCF, read with base R: one record per variant, GT = number of ALT alleles, ./. for missing
  v <- readLines(sim$files[["vcf"]])
  expect_identical(v[1], "##fileformat=VCFv4.2")
  hdr <- strsplit(v[startsWith(v, "#CHROM")], "\t")[[1]]
  expect_identical(hdr[-(1:9)], colnames(sim$genotype))
  rec <- do.call(rbind, strsplit(v[!startsWith(v, "#")], "\t"))
  expect_identical(rec[, 3], sim$variants$id); expect_identical(as.integer(rec[, 2]), sim$variants$pos)
  gt <- rec[, -(1:9)]
  dos <- matrix(match(gt, c("0/0", "0/1", "1/1")) - 1L, nrow(gt))
  expect_identical(unname(dos), unname(sim$genotype))
  expect_identical(sum(gt == "./."), sum(is.na(sim$genotype)))
})

test_that("simulate_bohr() rejects bad settings", {
  expect_error(simulate_bohr(n_individuals = 5), "n_individuals")
  expect_error(simulate_bohr(n_tissues = 0), "n_tissues")
  expect_error(simulate_bohr(variant_spacing = 10), "variant_spacing")
  expect_error(simulate_bohr(prop_causal = 2), "prop_causal")
  expect_error(simulate_bohr(effects = c(mu = 0.5, p0 = 0.1, p1 = 0.1)), "effects")
  expect_error(simulate_bohr(effects = c(0.1, 0.1, 0.1)), "effects")
})

test_that("bohr_test() recovers the true effects in a large sample", {
  skip_on_cran()
  sim <- with_seed(54, simulate_bohr(n_individuals = 1500, n_tissues = 1, n_genes = 4, prop_causal = 1))
  s <- sim$samples; tr <- sim$truth
  for (r in seq_len(nrow(tr))) {
    y <- sim$relative[tr$promoter[r], s$sample]; ty <- tr$predictor_type[r]
    x <- switch(ty, dosage = sim$genotype[tr$predictor[r], s$individual],
                binary = sim$traits$binary_trait[match(s$individual, sim$traits$individual)],
                continuous = sim$traits$continuous_trait[match(s$individual, sim$traits$individual)])
    b <- bohr_test(y, x, covariates = sim$covariates$T1[s$sample, ], type = ty, contrast = if (ty == "continuous") c(50, 60),
                   control = bohr_control(firth = FALSE))
    for (k in c("mu", "p0", "p1")) {
      e <- b[[paste0("effect_", k)]]; t <- tr[[paste0("effect_", k)]][r]
      if (is.na(e)) expect_true(abs(t) < 1e-12, info = paste(r, k)) else expect_lt(abs(e - t), 0.07)
    }
  }
})
