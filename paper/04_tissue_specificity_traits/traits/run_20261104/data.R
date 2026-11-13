## =============================================================================
## data.R -- shared definitions of the trait regression (new BOHR engine) (paths, trait coding, covariates, promoter preprocessing)
## Data: GTEx v11 (the same raw data in every run; folders are named by run date, not by v number, to avoid confusion with GTEx versions)
##
## The only change of this run trait_regression_20261007 relative to trait_regression_20261005 (author 2026-10-07):
##   * Drinking and smoking are now regressed as "ever drinks" MHDRNKSTS and "ever smokes" MHSMKSTS (0/1; blank and 99 are missing), replacing MHDRNKNMB,
##     MHSMKNMB. MHDRNKNMB is the number of drinks within a self-reported period
##     (Day/Week/Month/Infrequent), not comparable across donors, and only drinkers have a value. MHDRNKSTS takes the old position,
##     so the permutation seeds of the other 19 traits are unchanged; their results are linked directly from trait_regression_20261005 (BH is per trait, so unaffected).
##
## The only change of trait_regression_20261005 relative to the previous run trait_regression_20261001 (author 2026-10-05):
##   * The promoter set now uses GTEx v11 TPM: promoters of genes whose mean TPM in the tissue is >= 1 (TPM_v11/TPM1_genes.txt).
##     The previous run reused the v10 TPM filter table from v13; its gene Ensembl versions are 1-3 versions older than the promoter annotation, so matching by
##     versioned ID missed 1,191 genes (about 18% of tissue x promoter pairs, mostly lncRNA / antisense transcripts, also B2M, GSN, ACTA2, etc.).
##     The v11 TPM versions agree with the promoter annotation: 0 version mismatches in all tissues.
##   * The random seed of sim1 is now set by promoter ID (not by matrix position).
##   * Everything else follows the previous run (rules in DECISIONS.md):
##       samples = GTEx v11 analysis freeze (SMAFRZE == "RNASEQ"), one sample per donor (highest SMRIN);
##       engine = rerun_v2 split likelihood (bohr2_engine.R / bohr2_addon3.R copied unchanged, md5 in ENGINE_SOURCE.md5):
##         mu = beta-part LRT + small-sample F correction; p0/p1 = Firth-penalised LRT (computed as usual when cat_complete_sep);
##       endpoint removal = drop when the proportion of 0 or 1 is <= 5% (prepare_y);
##       covariates = SEX, AGE, PC1-5, DTHHRDY (factor; levels with < 10 donors in a tissue are merged into the reference), MHABNWBC; the tested trait is not a covariate;
##       donor phenotypes = dbGaP GRU table (metadata/v10_cov_pca_full_GRU.txt, restoring MH values of 138 donors); PC = PCA computed by the author;
##       binary traits need >= 5 donors per class; mu additionally needs >= 5 per class among the mid-value samples (mu_testable);
##       final rules (applied in the summary script): mu n_mid >= 80, p0/p1 boundary proportion 0.10-0.90, p0/p1 not completely separated, tissue donors >= 50,
##         mirror twins de-duplicated by component (shared table v2), then BH within tissue x trait x component.
## =============================================================================
suppressPackageStartupMessages(library(data.table))

GTEX  <- "/home/e1327989/project/GTEx"
## Run directory: local default; on NSCC the job sets TRG_RUN (the same code runs in both places)
RUN   <- Sys.getenv("TRG_RUN", "/projects/group/bliu/e1327989/GTEx/tissue_regression/trait_regression_20261104")
PATHS <- list(
  act   = file.path(GTEX, "v11/filter_5_recal_att_filter.txt"),
  info  = file.path(GTEX, "v11/Promoter_info_filter_5.txt"),
  meta  = file.path(GTEX, "metadata/v10_cov_pca_full_GRU.txt"),   # Since 2026-10-03 the dbGaP GRU table is used (restores MH values of 138 donors; see metadata/*.changes.tsv)
  attr  = "/home/e1327989/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt",
  tpm11 = file.path(GTEX, "TPM_v11/TPM1_genes.txt"),   # GTEx v11: genes with mean TPM >= 1 in each tissue (long table Name / Description / SMTSD / TPM_mean)
  tissues = file.path(RUN, "work/tissue_list.txt"),            # = copy of v13/tissue_list.txt (the order sets the permutation/simulation seeds)
  work  = file.path(RUN, "work"),
  res   = Sys.getenv("TRG_RES", file.path(RUN, "results")))   # TRG_RES can redirect this elsewhere for a smoke test

## ---- 20 targets (same order as v13) ---------------------------------------
TRAITS <- c("SEX", "AGE", "BMI", "MHDRNKSTS", "MHSMKSTS", "MHASTHMA", "MHCOPD", "MHCVD", "MHHRTDIS", "MHHTN",
            "MHT1D", "MHT2D", "MHDPRSSN", "MHHRTATT", "MHRNLFLR", "MHARTHTS", "MHPNMNIA", "MHCANCERNM",
            "MHLVRDIS", "MHCOCAINE5")
CONT_TRAITS <- c("AGE", "BMI")             # continuous: effect per +1 SD
## all others (including SEX) are binary: effect as unexposed -> exposed (SEX: male -> female)
MH_BINARY <- c("MHDRNKSTS", "MHSMKSTS", "MHASTHMA", "MHCOPD", "MHCVD", "MHHRTDIS", "MHHTN", "MHT1D", "MHT2D", "MHABNWBC", "MHDPRSSN",
               "MHHRTATT", "MHRNLFLR", "MHARTHTS", "MHPNMNIA", "MHCANCERNM", "MHLVRDIS", "MHCOCAINE5")
COVS <- c("SEX", "AGE", paste0("PC", 1:5), "DTHHRDY", "MHABNWBC")
MIN_CLASS <- 5L        # binary trait: tested only with at least 5 donors per class (same as MIN_GENO of the QTL pipeline);
                       # mu also requires at least 5 per class among the intermediate-value samples (n_minor_mid)
MIN_HARDY <- 10L       # DTHHRDY levels with fewer than 10 donors in a tissue are merged into the reference level

## ---- donor metadata: same missing-value coding as v13 ------------------------------------
read_meta <- function() {
  m <- fread(PATHS$meta, select = unique(c("SUBJID", TRAITS, COVS)), data.table = FALSE)
  ## dbGaP dictionary (NSCC .../gtex_v10/Sample_info/phs000424.v10.pht002742.v9.GTEx_Subject_Phenotypes.data_dict.xml):
  ## MH binary 0 = No, 1 = Yes, 98 = Not Reported, 99 = Unknown (there is no 97). v13 handled only 97/99; the data contain only 0/1/99,
  ## so adding 98 changes no result (checked 2026-10-03)
  for (v in MH_BINARY) { m[[v]] <- as.numeric(m[[v]]); m[[v]][m[[v]] %in% c(97, 98, 99)] <- NA }
  ## Note (2026-10-03): MHDRNKNMB is "the number of units within the drinking period MHDRNKPRD"; the period can be Day/Week/Month/Infrequent, so it is not directly comparable across donors;
  ## both have a value 98 (1 / 4 donors) that the dictionary does not define. Author decision 2026-10-07: drinking and smoking use MHDRNKSTS / MHSMKSTS (see file header)
  ## (MHDRNKNMB / MHSMKNMB, the two count variables, are no longer used since 2026-10-07)
  for (v in c("AGE", "BMI", paste0("PC", 1:5), "DTHHRDY")) m[[v]] <- as.numeric(m[[v]])
  m$SEX <- as.numeric(m$SEX == 2)                    # 0 = male, 1 = female
  rownames(m) <- m$SUBJID
  m
}

## ---- covariate design matrix (excluding the tested trait) -----------------------------------
## donors: donors already filtered for non-missing trait / covariates; returns intercept + standardised columns
cov_design <- function(meta, donors, trait) {
  cv <- setdiff(COVS, trait)
  d <- meta[donors, cv, drop = FALSE]
  num <- intersect(cv, c("SEX", "AGE", paste0("PC", 1:5), "MHABNWBC"))
  Z <- as.matrix(d[, num, drop = FALSE])
  if ("DTHHRDY" %in% cv) {                            # factors: the most common level is the reference
    ## Levels with fewer than MIN_HARDY donors in the tissue are merged into the reference level (fixed in advance):
    ## dummies of rare levels are nearly constant with high leverage, a common cause of (quasi-)separation in the categorical part
    tb <- table(d$DTHHRDY); ref <- names(which.max(tb))
    lv <- as.character(d$DTHHRDY); lv[lv %in% names(tb)[tb < MIN_HARDY]] <- ref
    f <- factor(lv); f <- relevel(f, ref = ref)
    if (nlevels(f) >= 2) {                            # if only one level remains after merging (e.g. Bladder) the covariate is constant and is left out
      D <- model.matrix(~ f)[, -1, drop = FALSE]; colnames(D) <- paste0("DTHHRDY_", levels(f)[-1])
      Z <- cbind(Z, D)
    }
  }
  Z <- Z[, apply(Z, 2, function(x) sd(x) > 0), drop = FALSE]   # drop columns that are constant in this donor set
  Z <- scale(Z)
  C <- cbind(`(Intercept)` = 1, Z); rownames(C) <- donors
  C
}

## ---- trait coding: returns g' so that the "dose 0 -> 2" in snp_tests is exactly the target contrast ------
##   continuous: g' = 2 (g - mean) / sd  -> 0 -> 2 = mean -> mean + 1 SD
##   binary:     g' = 2 g             -> 0 -> 2 = unexposed -> exposed
## Test statistics are invariant to linear transformation of g; estimates and SE are multiplied by 2 on output to return to the "per contrast" scale.
trait_contrast <- function(g, trait) {
  if (trait %in% CONT_TRAITS) 2 * (g - mean(g)) / sd(g) else 2 * g
}

## ---- promoter preprocessing and orthogonalisation: copied verbatim from rerun_v2/code/bohr2_data.R -----------------
## (sourcing the original file loads VariantAnnotation; genotypes are not needed here, so only these two functions are taken)
prepare_y <- function(y, p_cov) {
  y <- y[!is.na(y)]
  if (length(y) < 20 || length(y) < p_cov) return(list(skip = "n"))
  p0 <- mean(y == 0); p1 <- mean(y == 1)
  if (p0 >= 0.95 || p1 >= 0.95 || sd(y) <= 0.05) return(list(skip = "prop_or_sd"))
  keep <- rep(TRUE, length(y))
  if (p0 <= 0.05) keep <- keep & y > 0
  if (p1 <= 0.05) keep <- keep & y < 1
  y <- y[keep]
  z <- any(y == 0); o <- any(y == 1)
  list(skip = NA_character_, y = y,
       family = if (z && o) "BEINF" else if (z) "BEZI" else if (o) "BEOI" else "BE")
}
orth_design <- function(C, tol = 1e-7) {
  q <- qr(C, tol = tol); r <- q$rank
  Q <- qr.Q(q)[, seq_len(r), drop = FALSE]
  s <- sign(diag(qr.R(q))[seq_len(r)]); s[s == 0] <- 1
  Q <- sweep(Q, 2, s, "*") * sqrt(nrow(C))
  Q[, 1] <- 1
  dimnames(Q) <- list(rownames(C), c("(Intercept)", if (r > 1) paste0("z", seq_len(r - 1))))
  attr(Q, "dropped") <- colnames(C)[q$pivot[-seq_len(r)]]
  Q
}
