# Real-data parity checks (stages 2, 3a and 3b-1)

Compare the package with the stored results of the paper's original pipeline. GTEx individual-level data never enter
this repository: the scripts read them from the lab's storage and write only to a scratch directory given on the
command line. The comparisons print aggregate numbers only (counts, largest differences, agreement rates).
Paths to the lab storage are set at the top of each script. Results: `BOHR_v4/docs/code_release/STAGE1_2_REPORT.md`
(cis) and `STAGE3A_REPORT.md` (conditional analysis and traits).

## cis mapping, `bohr_cis()` (stage 2)

1. `01_extract_inputs.R <tissue> <chr> <inputs.rds>` (R with data.table): rebuilds the inputs of the stored run
   as plain matrices: activity of the tissue's promoters on the chromosome, GT dosages of every VCF record for the
   analysis samples (missing kept; `bohr_cis()` fills them), PC1-5 and PEER 1..K, and the pooled-mean score
   that picks the canonical member of mirror pairs.
2. `02_run_bohr_cis.R <inputs.rds> <out.rds> <ncores> [max_promoters] [lib]` (R with BOHR installed).
3. `03_compare_stored.R <out.rds> <tissue> <chr> [summary.rds]`: per-SNP statistics, promoter rules, units
   (ACAT P, leader, estimability, mirror flags) and the tissue-wide q values, recomputed with this chromosome's
   units replaced by the package's.

## Independent signals, `bohr_conditional()` (stage 3a)

4. `04_extract_conditional_inputs.R <tissue> <inputs.rds> [units.tsv]` (R with data.table and VariantAnnotation):
   the inputs of the stored conditional runs for the significant units of one tissue (all chromosomes, or the units
   listed in `units.tsv`), with the genotypes of each unit promoter's window read as the original `read_geno()` did,
   and the stored units, threshold and tested pairs (the `cis` input).
5. `05_run_conditional.R stored <inputs.rds> <out.rds> <ncores> [lib]`: `bohr_conditional()` with the stored units as
   `cis`; or `05_run_conditional.R cis <inputs.rds> <cis.rds> <out.rds> <ncores> [lib]`: chained after a `bohr_cis()` run
   of one chromosome (scripts 1-2), with its q values and threshold recomputed over the tissue.
6. `06_compare_conditional.R <run.rds> [summary.rds]`: signal counts per unit, lead variants, conditional and
   leave-one-out statistics, rounds and backward calls against the stored v3c tables (the published 14,800 signals)
   and v3m.

## Donor traits, `bohr_trait()` (stage 3a)

7. `07_extract_trait_inputs.R <out_dir> <tissue:trait> ... [cis:<tissue>:<inputs.rds>]` (R with data.table and Matrix):
   activity, traits, covariates (the Hardy scale as a categorical column), promoters with genes, the paper's
   canonical mirror members and the absolute activity of the exact-pair promoters; `cis:` adds the absolute activity
   for the conditional label check to an input file of script 4.
8. `08_run_trait.R <inputs.rds> <out.rds> <ncores> [lib]`.
9. `09_compare_trait.R <run.rds> [summary.rds]`: P values, estimates, effects, q values, eligibility, significance,
   mirror flags and labels against `export_B/trait_regression_units.tsv.gz`, per-promoter skip reasons against the
   stored rows, and the published examples (PDE5A in breast, CSNK1D in tibial artery).

## The null start of the restricted Firth fits (stage 3b-1)

`bohr_control(firth_null_start = TRUE)` (the default since stage 3b-1) also starts each restricted Firth fit from the null fit.
Scripts 02, 05 and 08 use the package default; `BOHR_FIRTH_NULL_START=FALSE` runs the paper's engine exactly. With
`BOHR_SAVE_SCANS=1`, script 05 also saves every conditional scan (each variant's conditional estimate, SE and log P in every
forward round and leave-one-out scan).

10. `10_compare_null_start.R <cis|cond|trait> <old.rds> <new.rds> [summary.rds]` (base R): compares a run with the paper's
    engine (or of an earlier commit) with a run with the null start: counts and largest differences of the Firth statistics,
    and every changed decision (units, signals, rounds, conditional scans, trait calls). Results: `STAGE3B1_REPORT.md`.

## Tissue means, `tissue_means()` (stage 3b-1)

11. `11_tissue_means_parity.R <out_dir> <ncores> [lib] [tissue ...]` (R with BOHR, or `BOHR_LOAD_ALL=<package dir>`): runs
    `tissue_means()` on the v11 relative activity over the samples of the 944 genotyped donors (the sample set of the stored
    gamlss fits) and compares with `fit/v11_raw_fit_fit_results.txt`: pairs and sample sizes, the sample-mean rules, the means,
    the families, and the 0.10 ROKU effect filter with ROKU's stored calls held fixed. Run it as a PBS job (about 1 GB of input;
    8 cores, 24 GB).
