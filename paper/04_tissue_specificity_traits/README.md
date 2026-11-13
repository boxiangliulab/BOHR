# Tissue specificity and donor traits

**Methods §6; Supplementary Methods §6.** Tissue-specific promoters: one representative mean per promoter and tissue
(intercept-only zero-and-one-inflated beta fit, `raw_fit()`), ROKU calls, and an effect filter of 0.10. Donor traits:
BOHR tests of μ, p0 and p1 against 20 donor traits in each tissue, with simulated and permuted arms for calibration.

**Run order.** Tissue specificity: `tissue_specificity/gamlss_functions_optimized.R` (`raw_fit()`) →
`tissue_specificity/ROKU.R` → `tissue_specificity/roku_effect_filter.R`. Traits: `traits/run_20261005/prep.R`
(per-tissue donor × promoter matrices) → `run_trait.R` of run 20261005 (18 traits) and run 20261007 (drinker and smoker
status), each per tissue and arm → `traits/run_20261104/mirror_v3.R` and `summarize.R` (mirror rule, testability rules,
BH per tissue × trait × component) → `export_fig6.R` → `traits/run_20261106_mirror_abs/01_abs_labels.R` and
`03_apply_B.R` (exact-pair events labelled by absolute activity; the paper's results).

**BOHR package.** `BOHR::tissue_means()` reproduces the tissue means (all 56 tissues, 430,464 pairs; means within about
0.001 where both fits converged; the paper fitted with gamlss) and `BOHR::bohr_trait()` the trait tests (repository
`workflow/parity/`, scripts 07–09 and 11).

`run.sh` lists the commands in the order they were run (paths of the original runs).


## tissue_specificity/

| file | what it does | inputs → outputs |
|---|---|---|
| `tissue_specificity/gamlss_functions_optimized.R` | Provides the published GAMLSS functions, including raw_fit() for the tissue-mean fits used by ROKU (Fig. 3). | sourced library; raw_fit() → fit/v11_raw_fit_fit_results.txt |
| `tissue_specificity/ROKU.R` | Calls tissue-specific promoters with TCC::ROKU on the promoter x tissue matrix of fitted mean activity. | fit/v11_raw_fit_fit_results.txt (tissue-level overall_mean per promoter) → v11_tissue_specific_ROKU_results.rds ($outlier) |
| `tissue_specificity/roku_effect_filter.R` | Keeps a raw ROKU call only when the tissue differs from the median of non-outlier tissues by >= 0.10. | v11_tissue_specific_ROKU_results.rds + fit table → roku_20261007/v11_tissue_specific_ROKU_results_delta010.rds, summary_delta_filter.tsv |

## traits/

| file | what it does | inputs → outputs |
|---|---|---|
| `traits/run_20261005/data.R` | Defines paths, trait coding and read_meta() for run 20261005; sourced by prep.R. | sourced by prep.R (PATHS, read_meta, TRAITS with the drink and smoke count variables) |
| `traits/run_20261005/prep.R` | Builds per-tissue donor x promoter activity matrices (one analysis-freeze sample per donor) and donor metadata. | v11 promoter activity, sample attributes (SMAFRZE, SMRIN), v11 TPM>=1 gene list, donor covariates → work/<tissue>.rds, meta.rds, tissues.tsv |
| `traits/run_20261005/mirror.R` | Marks exact two-promoter mirror pairs per tissue from Figure 4's shared mirror table (the copy that wrote the table used). | Figure 4 mirror/promoter_flags.tsv.gz and mirror_pairs.tsv.gz, work/tissues.tsv, work/<tissue>.rds → work/mirror.tsv, mirror_check.tsv |
| `traits/run_20261005/run_trait.R` | Fits the BOHR trait tests (mu, p0, p1) for one tissue and arm; produced the results of 18 of the 20 donor traits. | work/<tissue>.rds, work/meta.rds, work/tissue_list.txt → results/{real,sim1,perm1}/<tissue>/<trait>.rds (18 traits, linked into runs 20261007 and 20261104) |
| `traits/run_20261007/data.R` | Defines the trait regression's paths, 20 donor traits, covariates, covariate design and outcome screening rules. | sourced by run_trait.R (PATHS, TRAITS, COVS, read_meta, cov_design, trait_contrast, prepare_y) |
| `traits/run_20261007/mirror.R` | Marks exact two-promoter mirror pairs per tissue from Figure 4's shared mirror table (mirror rule v2 pairs). | Figure 4 mirror/promoter_flags.tsv.gz and mirror_pairs.tsv.gz, work/tissues.tsv, work/<tissue>.rds → work/mirror.tsv, mirror_check.tsv |
| `traits/run_20261007/run_trait.R` | Fits the BOHR trait tests (mu, p0, p1) for every promoter of one tissue in one arm (real, simulated or permuted). | work/<tissue>.rds, work/meta.rds, work/tissue_list.txt → results/{real,sim1,perm1}/<tissue>/<trait>.rds |
| `traits/run_20261104/data.R` | Defines paths and trait coding for run 20261104; sourced by summarize.R, mirror_v3.R and the breast expression check. | sourced by summarize.R, mirror_v3.R, expr_check_breast_sex_20261104.R (also calib_strata.R, export_fig6.R, 01_abs_labels.R, 03_apply_B.R) |
| `traits/run_20261104/mirror_v3.R` | Applies mirror rule v3 to exact two-promoter pairs: keeps p0 on both members and drops p1. | work/mirror_v2.tsv (mirror.R output), Figure 4 mirror_v3/promoter_flags.tsv.gz → work/mirror_v3.tsv, mirror_v3_check.tsv |
| `traits/run_20261104/summarize.R` | Applies mirror rule v3 and the final testability rules to the trait-test units; BH per tissue x trait x component. | results/{real,sim1,perm1}/*/*.rds, work/mirror_v3.tsv → results/summary/units_*.rds, sig_*.tsv, counts_*.tsv, null_discovery.tsv |
| `traits/run_20261104/export_fig6.R` | Exports all 6,246,053 testable real-arm trait tests; input of the option-B relabelling and of Fig. 3C. | results/summary/units_real.rds → results/export/trait_regression_units.tsv.gz, README_trait_regression.txt |
| `traits/run_20261106_mirror_abs/01_abs_labels.R` | Labels each exact-pair boundary event by which promoter's absolute activity changes with the trait (option B; run once). | units_real.rds, work/<tissue>.rds, work/meta.rds, v11/tabs_att_filter_In.rds (absolute activity) → results/mirror_abs_labels.tsv.gz |
| `traits/run_20261106_mirror_abs/03_apply_B.R` | Applies the absolute-activity labels to the 586 significant exact-pair events (option B); rewrites export and summary. | mirror_abs_labels.tsv.gz, results/summary/ (units, sig_real, counts, null_discovery), results/export/ → results/export_B/, results/summary_B/ |
| `config/tissue_list.txt` | Lists the 57 tissue names of the trait regression; the order sets the permutation and simulation seeds. | read by prep.R and run_trait.R through data.R PATHS$tissues (RUN/work/tissue_list.txt) |
