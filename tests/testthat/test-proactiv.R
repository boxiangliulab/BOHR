## relative_activity() and from_proactiv(): a hand-built set where every promoter's fate is known

## 20 samples, tissues T1 (1-10) and T2 (11-20); absolute activities by promoter
ra_case <- function() {
  s <- paste0("S", 1:20); t1 <- 1:10; t2 <- 11:20
  A <- matrix(0, 17, 20, dimnames = list(sprintf("P%02d", 1:17), s))
  A["P01", ] <- 4; A["P02", ] <- 2; A["P03", ] <- 6                    # G1: P03 internal
  A["P04", ] <- 3                                                         # G2: one annotated promoter
  A["P05", ] <- 3; A["P06", ] <- 0.1                                      # G3: P06 inactive, P05 left alone
  A["P07", ] <- 3; A["P08", ] <- 3; A["P09", 1] <- 10                     # G4: P09 above 0 in 1/10 samples
  A["P10", ] <- 3; A["P11", ] <- 3; A["P11", c(1:5, 11:15)] <- NA         # G5: half missing in both tissues
  A["P12", ] <- 5; A["P13", ] <- 5; A["P14", ] <- 0.3                     # G6: P14 mean relative ~0.03
  A["P15", ] <- 5; A["P16", 1:4] <- 0.1; A["P16", 11:12] <- 50            # G7: P16 passes nonzero in T1, mean in T2
  A["P17", ] <- 1
  gene <- c("G1", "G1", "G1", "G2", "G3", "G3", "G4", "G4", "G4", "G5", "G5", "G6", "G6", "G6", "G7", "G7", "G8")
  internal <- rownames(A) %in% c("P03", "P17")                            # G8: an internal promoter only
  list(A = A, gene = gene, internal = internal, tissue = rep(c("T1", "T2"), each = 10), s = s)
}

test_that("relative_activity() drops each promoter at the expected step", {
  k <- ra_case()
  ra <- relative_activity(k$A, gene = k$gene, internal = k$internal, group = k$tissue)
  expect_s3_class(ra, "bohr_relative")
  why <- stats::setNames(ra$promoters$reason, ra$promoters$id)
  expect_equal(unname(why[c("P03", "P17", "P06", "P04", "P09", "P10", "P11", "P14", "P05")]),
               c("internal", "internal", "inactive", "single_promoter_gene", "nonzero", "missing", "missing", "mean", "few_promoters"))
  expect_equal(rownames(ra$relative), c("P01", "P02", "P07", "P08", "P12", "P13", "P15", "P16"))
  expect_true(all(is.na(why[rownames(ra$relative)])))
  expect_equal(ra$steps$step, c("annotated", "external", "active", "multi_promoter_gene", "nonzero", "missing", "mean", "min_promoters"))
  expect_equal(ra$steps$promoters, c(17, 15, 14, 13, 12, 10, 9, 8))
  expect_equal(ra$steps$genes, c(8, 7, 7, 6, 6, 5, 5, 4))
  expect_equal(ra$steps$genes_multi[1:2], c(6, 6))
  expect_equal(ra$promoters$max_nonzero[ra$promoters$id == "P09"], 0.1)
  expect_equal(ra$promoters$min_missing[ra$promoters$id == "P11"], 0.5)
})

test_that("the denominator holds every external promoter, filtered or not, and no internal one", {
  k <- ra_case()
  ra <- relative_activity(k$A, gene = k$gene, internal = k$internal, group = k$tissue)
  expect_equal(unname(ra$relative["P01", ]), rep(4 / 6, 20))                        # P03 (internal) left out
  expect_equal(unname(ra$relative["P12", ]), rep(5 / 10.3, 20))                     # P14 (dropped) counted
  expect_equal(unname(ra$relative["P07", 1]), 3 / 16)
  all_ext <- relative_activity(k$A, gene = k$gene, internal = k$internal, filter = FALSE)
  expect_equal(nrow(all_ext$relative), 15L)
  expect_equal(all_ext$relative[rownames(ra$relative), ], ra$relative)
  expect_true(all(is.na(all_ext$relative["P10", c(1:5, 11:15)])))                   # a missing value propagates
  B <- k$A; B[c("P05", "P06"), 3] <- 0
  expect_true(is.na(relative_activity(B, gene = k$gene, internal = k$internal, filter = FALSE)$relative["P05", 3]))  # 0 / 0
  expect_equal(all_ext$steps$step, c("annotated", "external"))
})

test_that("criteria may be met in different tissues, and the tissues matter", {
  k <- ra_case()
  ra <- relative_activity(k$A, gene = k$gene, internal = k$internal, group = k$tissue)
  expect_true("P16" %in% rownames(ra$relative))
  one <- relative_activity(k$A, gene = k$gene, internal = k$internal)                 # one group of 20
  expect_identical(one$promoters$reason[one$promoters$id == "P16"], "mean")          # 6/20 above 0 passes; mean 0.095 does not
  expect_equal(one$promoters$max_nonzero[one$promoters$id == "P16"], 0.3)
  ## a threshold is strict: 0.3 is not > 0.3
  strict <- relative_activity(k$A, gene = k$gene, internal = k$internal, min_nonzero = 0.3)
  expect_identical(strict$promoters$reason[strict$promoters$id == "P16"], "nonzero")
  loose <- relative_activity(k$A, gene = k$gene, internal = k$internal, group = k$tissue, min_mean = 0.01, max_missing = 0.6)
  expect_true(all(c("P10", "P11", "P14") %in% rownames(loose$relative)))
  noact <- relative_activity(k$A, gene = k$gene, internal = k$internal, group = k$tissue, require_active = FALSE)
  expect_false("active" %in% noact$steps$step)
  expect_identical(noact$promoters$reason[noact$promoters$id == "P06"], "mean")
  expect_equal(relative_activity(k$A, gene = k$gene, internal = k$internal, group = k$tissue, min_abs = 0.05)$promoters$reason[6], "mean")
})

test_that("relative_activity() takes a TSS table, files and out_dir", {
  k <- ra_case()
  tss <- data.frame(id = rev(rownames(k$A)), gene = rev(k$gene), internal = rev(k$internal))
  a <- relative_activity(k$A, gene = k$gene, internal = k$internal, group = k$tissue)
  b <- relative_activity(k$A, gene = tss, group = stats::setNames(k$tissue, k$s))
  expect_equal(a$relative, b$relative); expect_equal(a$steps, b$steps)
  d <- withr::local_tempdir()
  f <- file.path(d, "abs.tsv")
  write.table(data.frame(id = rownames(k$A), k$A), f, sep = "\t", quote = FALSE, row.names = FALSE)
  c2 <- relative_activity(f, gene = tss, group = data.frame(sample = k$s, tissue = k$tissue), out_dir = file.path(d, "out"))
  expect_equal(c2$relative, a$relative)
  expect_true(all(file.exists(file.path(d, "out", c("relative_activity.tsv", "promoters.tsv", "steps.tsv")))))
  back <- read.delim(file.path(d, "out", "relative_activity.tsv"), check.names = FALSE)
  expect_equal(back$id, rownames(a$relative))
  expect_output(print(a), "8 promoters of 4 genes")
})

test_that("relative_activity() reproduces simulate_bohr()'s relative activity", {
  set.seed(3)
  sim <- simulate_bohr(n_individuals = 40, n_tissues = 2, n_genes = 4)
  r <- relative_activity(sim$absolute, gene = sim$tss, filter = FALSE)$relative
  expect_equal(r[rownames(sim$relative), colnames(sim$relative)], sim$relative, tolerance = 1e-10)
  expect_false(any(sim$tss$id[sim$tss$internal] %in% rownames(r)))
})

test_that("relative_activity() checks its input", {
  k <- ra_case()
  expect_error(relative_activity(k$A, gene = k$gene[-1]), "gene")
  expect_error(relative_activity(k$A, gene = data.frame(id = "P01", gene = "G1")), "no row for 16")
  expect_error(relative_activity(k$A, gene = k$gene, internal = rep(NA, 17)), "internal")
  expect_error(relative_activity(k$A, gene = k$gene, filter = NA), "filter")
  expect_error(relative_activity(k$A, gene = k$gene, min_promoters = 0), "min_promoters")
  expect_error(relative_activity(k$A, gene = k$gene, min_mean = "a"), "min_mean")
  B <- k$A; rownames(B) <- NULL
  expect_error(relative_activity(B, gene = k$gene), "row names")
  expect_error(relative_activity(k$A, gene = k$gene, group = k$tissue[-1]), "tissue")
})

## a proActiv-like result for one tissue: assays and row data named as proActiv() names them
fake_proactiv <- function(k, cols, cond, rel = NULL) {
  A <- k$A[, cols, drop = FALSE]
  if (is.null(rel)) { tot <- rowsum(A, k$gene); rel <- A / tot[k$gene, , drop = FALSE] }   # proActiv: internal in the denominator
  rownames(A) <- rownames(rel) <- NULL
  rd <- S4Vectors::DataFrame(promoterId = 1:17, geneId = k$gene, seqnames = "chr1", start = seq(100, 1700, 100), strand = "+",
                             internalPromoter = k$internal, promoterPosition = 1L)
  mn <- rowMeans(A)
  rd[[paste0(cond, ".class")]] <- ifelse(k$internal, NA, ifelse(mn < 0.25, "Inactive", "Minor"))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(promoterCounts = A, normalizedPromoterCounts = A, absolutePromoterActivity = A, relativePromoterActivity = rel,
                  geneExpression = A),
    rowData = rd, colData = S4Vectors::DataFrame(sampleName = k$s[cols], condition = cond, row.names = k$s[cols]))
}

test_that("from_proactiv() collects proActiv results; relative_activity() uses its classes and relative activity", {
  skip_if_not_installed("SummarizedExperiment")
  k <- ra_case()
  se1 <- fake_proactiv(k, 1:10, "T1"); se2 <- fake_proactiv(k, 11:20, "T2")
  pa <- from_proactiv(list(se1, se2))
  expect_s3_class(pa, "bohr_proactiv")
  expect_equal(dim(pa$absolute), c(15L, 20L))
  expect_equal(rownames(pa$absolute), paste0("PromoterID_", which(!k$internal)))
  expect_equal(colnames(pa$class), c("T1", "T2"))
  expect_equal(unname(pa$condition), k$tissue)
  expect_equal(pa$promoters$internal, k$internal)
  expect_output(print(pa), "17 promoters")
  ## the same calls as from the matrix (proActiv's classes follow the same 0.25 rule here)
  ra <- relative_activity(pa, group = pa$condition)
  m <- relative_activity(k$A, gene = k$gene, internal = k$internal, group = k$tissue)
  expect_equal(sub("PromoterID_", "", rownames(ra$relative)), as.character(as.integer(sub("P", "", rownames(m$relative)))))
  expect_equal(unname(ra$relative), unname(m$relative))
  expect_equal(ra$settings$active_from, "proActiv classes")
  ## classes decide activity: a promoter classified Inactive everywhere goes, one with a missing class too
  se1b <- se1; SummarizedExperiment::rowData(se1b)$T1.class[1] <- "Inactive"
  se2b <- se2; SummarizedExperiment::rowData(se2b)$T2.class[1] <- "Inactive"
  expect_identical(relative_activity(from_proactiv(list(se1b, se2)), group = k$tissue)$promoters$reason[1], NA_character_)
  expect_identical(relative_activity(from_proactiv(list(se1b, se2b)), group = k$tissue)$promoters$reason[1], "inactive")
  se2c <- se2; SummarizedExperiment::rowData(se2c)$T2.class[1] <- NA
  expect_identical(relative_activity(from_proactiv(list(se1, se2c)), group = k$tissue)$promoters$reason[1], "inactive")
  ## criteria on proActiv's relative activity: a promoter whose proActiv value passes is kept although its own does not
  rel <- k$A[, 1:10] / rowsum(k$A[, 1:10], k$gene)[k$gene, ]; rel[14, ] <- 0.5
  pa2 <- from_proactiv(list(fake_proactiv(k, 1:10, "T1", rel = rel), se2))
  r2 <- relative_activity(pa2, group = k$tissue)
  expect_true("PromoterID_14" %in% rownames(r2$relative))
  expect_equal(unname(r2$relative["PromoterID_14", 1]), 0.3 / 10.3)              # the values returned are the external ones
  ## files
  d <- withr::local_tempdir(); f <- file.path(d, c("a.rds", "b.rds")); saveRDS(se1, f[1]); saveRDS(se2, f[2])
  expect_equal(from_proactiv(f)$absolute, pa$absolute)
  expect_error(from_proactiv(list(se1, se2[2:17, ])), "same promoters")
  expect_error(from_proactiv(file.path(d, "none.rds")), "not found")
  expect_error(from_proactiv(list(1)), "SummarizedExperiment")
})
