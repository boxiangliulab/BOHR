## Real-data parity of bohr_trait(), step A (gtexfig R: data.table + Matrix, for the sparse absolute-activity matrix): the inputs of the stored trait regression for one or more
## tissue x trait pairs, as plain objects, plus absolute activity (proActiv, v11/tabs_att_filter_In.rds $abs) for the promoters of
## the exact mirror pairs. With "cis:<tissue>" it also adds the absolute activity of the significant exact-pair p0 units of the
## cis conditional analysis to an input file of 04_extract_conditional_inputs.R. Output stays in the scratch folder.
## Usage: Rscript 07_extract_trait_inputs.R <out_dir> <tissue:trait> [<tissue:trait> ...] [cis:<tissue>:<cond_inputs.rds>]
suppressPackageStartupMessages({ library(data.table); library(Matrix) })
a <- commandArgs(TRUE); OUT <- a[1]; JOBS <- a[-1]
TRG <- "/projects/group/bliu/e1327989/GTEx/tissue_regression/trait_regression_20261104"
V2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"
ABS <- "/home/e1327989/project/GTEx/v11/tabs_att_filter_In.rds"
meta <- readRDS(file.path(TRG, "work/meta.rds"))
MIR <- fread(file.path(TRG, "work/mirror_v3.tsv"))
t0 <- Sys.time(); ab <- readRDS(ABS)$abs; cat("absolute activity loaded:", nrow(ab), "x", ncol(ab), "in", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")
for (j in JOBS) {
  f <- strsplit(j, ":", fixed = TRUE)[[1]]
  if (f[1] == "cis") {                                                   # cis: absolute activity for the conditional label check
    TIS <- f[2]; IN <- f[3]; x <- readRDS(IN)
    td <- readRDS(file.path(V2, "work/final/freeze", paste0(TIS, ".rds")))
    U <- x$units; U <- U[U$component == "p0" & U$boundary_exact %in% TRUE, ]
    need <- unique(c(U$promoter, U$mirror_partner)); stopifnot(all(need %in% colnames(ab)), all(td$row_id %in% rownames(ab)))
    A <- t(as.matrix(ab[td$row_id, need, drop = FALSE])); colnames(A) <- td$U
    x$absolute <- A; saveRDS(x, sub("\\.rds$", "_abs.rds", IN))
    cat(sprintf("%s cis: %d significant exact-pair p0 units, absolute activity of %d promoters x %d samples\n", TIS, nrow(U), nrow(A), ncol(A)))
    next
  }
  TIS <- f[1]; TR <- f[2]
  td <- readRDS(file.path(TRG, "work", paste0(TIS, ".rds")))
  prs <- colnames(td$act)
  m <- meta[td$U, , drop = FALSE]
  cov <- data.frame(SEX = m$SEX, AGE = m$AGE, PC1 = m$PC1, PC2 = m$PC2, PC3 = m$PC3, PC4 = m$PC4, PC5 = m$PC5,
                    DTHHRDY = as.character(m$DTHHRDY), MHABNWBC = m$MHABNWBC, row.names = td$U, stringsAsFactors = FALSE)
  traits <- data.frame(m[, TR, drop = FALSE], row.names = td$U)
  mt <- MIR[tissue == TIS]
  ## the canonical member of each pair as decided in the paper's mirror table (mirror.R: pooled mean, shared table first)
  major <- setNames(rep(0.5, length(prs)), prs); major[mt[mirror_role == "canonical", promoter]] <- 1; major[mt[mirror_role == "partner", promoter]] <- 0
  ex <- unique(mt[boundary_exact == TRUE, promoter]); stopifnot(all(ex %in% colnames(ab)), all(td$row_id %in% rownames(ab)))
  A <- t(as.matrix(ab[td$row_id, ex, drop = FALSE])); colnames(A) <- td$U
  out <- list(phenotype = t(td$act), traits = traits, covariates = cov, promoters = data.frame(id = prs, gene = td$info[prs, "geneId"], stringsAsFactors = FALSE),
              major_score = major, absolute = A, tissue = TIS, trait = TR)
  saveRDS(out, file.path(OUT, sprintf("trait_%s_%s_inputs.rds", TIS, TR)))
  cat(sprintf("%s x %s: %d promoters, %d donors, %d exact-pair promoters with absolute activity; mirror rows %d (boundary_exact %d)\n",
              TIS, TR, length(prs), length(td$U), length(ex), nrow(mt), mt[boundary_exact == TRUE, .N]))
}
