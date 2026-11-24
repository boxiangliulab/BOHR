# The BOHR model and simulation benchmarking

**Methods §2–3; Supplementary Methods §2–3.** `engine/` is the BOHR likelihood used by every analysis of the paper: the
zero-and-one-inflated beta model fitted as a beta regression on the values between 0 and 1 and a multinomial logit model
for 0, between and 1; likelihood-ratio tests of μ (with the small-sample F correction) and Firth-penalised
likelihood-ratio tests of p0 and p1. `gamlss_check/` compares the engine with gamlss refits. The simulations measure
calibration and power against QTLtools-style linear tests (Fig. 2A–D, Extended Data Fig. 2, Supplementary Tables 2a–c
and 2f). The comparison with QTLtools on GTEx data (Fig. 2E,F, Supplementary Tables 2d,e) is in
`../03_cis_mapping_conditional/qtltools_comparison/`.

**Run order.** The simulated data were generated once and reused. Simulated activity:
`simulated_activity/newsim_direct_nothreshold_indep.R`. Simulated RNA-seq: `simulated_rnaseq/v12_simulation_nothr.R`
(generator) run by `v12_sim_run_nothr_big.R` for two gene batches (b1, b2), then STAR. For the final version, proActiv was
run again on the unfiltered STAR junction files, with no junction-read floor (`simulated_rnaseq/sim_proactiv_nofloor.R`,
on NSCC, 2 batches × 12 sample sizes); it also writes each gene's junction reads. Fits: `benchmark/bohr2_sim.R` (BOHR) and
`benchmark/baselines_v3.R` (linear tests), for power and for calibration (1,000 genotype permutations of the calibration
genes in `config/`, drawn as in `INLINE_STEPS.md`), with the final promoter rules (`SKIP_FRAC=0.80`, `SD_SKIP=0`). The
calibration runs were repeated with four more permutation seed bases (Supplementary Table 2b). Summaries, all with
`DATASET=E80_nofloor_reads5` (simulated RNA-seq genes with a median of fewer than five gene reads are dropped):
`benchmark/combined_p_and_p0p1_v3.R` (null FPR and QQ points of the p0 and p1 tests); then, with `COMBINE=fisher_sidak`
(the combined BOHR P value of Methods §2), `benchmark/prepare_data_v3.R` (power, thresholds, FPR, QQ, effect recovery),
`benchmark/prepare_power_by_effect.R` (per-promoter calls at N = 300), `benchmark/prepare_power_by_OR.R` (power by
per-allele odds ratio, run with the main bins and with `OR_BINS=lin5`) and `benchmark/multiseed_final_v3.R` (genomic
inflation over five permutation sets).

**BOHR package.** `BOHR::bohr_fit()`, `BOHR::bohr_test()` and `BOHR::simulate_response()` implement this engine; the
package's tests compare it with stored engine fits (`tests/testthat/test-engine-parity.R`). The defaults of
`BOHR::bohr_control()` (`R/control.R`: `min_sd = 0.05`, `max_boundary = 0.95`) are the rules of the earlier version. The
final analyses skipped a response with SD < 0.005 or with at least 80% of its values at one boundary (the package skips
SD <= `min_sd`, and shares >= `max_boundary`).

**Notes.** `engine/bohr2_data.R` holds the input paths and the promoter rules used by the cis-mapping scripts of
`../03_cis_mapping_conditional/`, which source the three `bohr2_*.R` files from one folder (`BOHR2_CODE`). Its
`prepare_y()` holds the final SD rule: skip only SD < 0.005 (the earlier version skipped SD <= 0.05). It still skips a
share of 0 or of 1 of at least 95% at the scan; the 80% boundary-share rule and the gene-reads rule (gene median of at least
five first-intron junction reads) of the final version are applied to the units, in
`../03_cis_mapping_conditional/units/build_r5e80_units.R`. In the simulation, `bohr2_sim.R` applies the 80% rule inside
`prepare_y()` (`SKIP_FRAC`) and `prepare_data_v3.R` the gene-reads rule (`gene_reads.tsv`).

The switches of the summary and fit scripts (`SKIP_FRAC`, `SD_SKIP`, `DATASET`, `COMBINE`, `OR_BINS`) default to the
earlier version, which they reproduce. `combined_p_and_p0p1_v3.R`, `prepare_power_by_effect.R` and
`multiseed_final_v3.R` evaluate the first part of `prepare_data_v3.R` when they run. `combined_p_and_p0p1_v3.R` ran before
the `COMBINE` switch was added to `prepare_data_v3.R`; at its default (`simes`) the switch changes nothing, and the p0 and
p1 tests do not depend on it. Its output folder (`data_v3_E80_nofloor_reads5/diag_tail`) was linked into the folder of the
final summaries (`data_v3_E80_nofloor_reads5_fishersidak/diag_tail`), which Extended Data Fig. 2C and Supplementary Table
2c read. Scripts that only draw the figures or format the Supplementary Tables are not included.

`run.sh` lists the commands in the order they were run (paths of the original runs).


## engine/

| file | what it does | inputs → outputs |
|---|---|---|
| `engine/bohr2_engine.R` | Implements the BOHR likelihood: split zero/one-inflated beta fits, null fits and per-SNP LRT and Wald tests. | sourced library (base R only) |
| `engine/bohr2_addon3.R` | Adds Firth-penalised p0/p1 likelihood-ratio tests (multi-start, alternating search, BFGS polish). | sourced library (needs bohr2_engine.R) |
| `engine/bohr2_data.R` | Defines BOHR data handling: input paths, promoter skips (SD < 0.005, share at one boundary >= 0.95), endpoint deletion and family choice, cis window, genotypes, covariates. | sourced library; reads v11 activity matrix, promoter info, genotype PCs, PEER factors, v9 VCF |
| `engine/pe_delta.R` | Computes probability-scale effects (dosage 0 to 2) with delta-method SEs from the fits behind the tests. | sourced library (after bohr2_engine.R, bohr2_addon3.R) |

## gamlss_check/

| file | what it does | inputs → outputs |
|---|---|---|
| `gamlss_check/prep_pilot.R` | Caches the pilot tissue inputs that the gamlss equivalence check reads (Supplementary Methods 2, check against gamlss). | v11 relative activity, genotype PCs, PEER factors, promoter lists → work/tissue/<T>.rds, results/stage0_prep_checks.tsv |
| `gamlss_check/equiv_test.R` | Checks the BOHR engine against gamlss refits (null, full, reduced models) on real and simulated data. | work/tissue/<T>.rds (pilot cache), v9 VCF → results/<OUTSUB>/{equiv_summary.txt, equiv_models.tsv, equiv_lrt.tsv} |
| `gamlss_check/bohr2_engine.stage1b.R` | Holds the BOHR engine version that ran the gamlss equivalence check of Supplementary Methods 2 (before amendment 4). | sourced as code/bohr2_engine.R by equiv_test.R in run 2 → results/stage1_equiv_run2/equiv_summary.txt |

## simulated_activity/ and simulated_rnaseq/

| file | what it does | inputs → outputs |
|---|---|---|
| `simulated_activity/newsim_direct_nothreshold_indep.R` | Generates simulated-activity data: 2,000 genes, BEINF draws, one genotype vector per gene, true effects. | args N 2000 2; env TRUNC_L=0.005 GEN_ONLY=1 → m1_nothr_g2000/sim_N<N>_G2000_realfloor/{sim_data.RData, ground_truth.txt} |
| `simulated_rnaseq/v12_simulation_nothr.R` | Defines the simulated RNA-seq generator: two promoters per gene, BEINF usage, genotypes and polyester reads. | Promoter_info_filter_5.RData, GRCh38 cDNA FASTA → test_output/{sample_*_{1,2}.fasta, sim_tx_info_extended.txt, genotypes.vcf} |
| `simulated_rnaseq/v12_sim_run_nothr_big.R` | Runs the simulated RNA-seq generator for one N and batch: 1,200 genes, four effect types at 0.25 each. | args N sim_dir gene_seed data_seed prob → method2_polyester/sim_<N>_nothr_big{,_b2}/test_output/ |
| `simulated_rnaseq/sim_proactiv_nofloor.R` | Runs proActiv on the unfiltered simulated junction files (no read floor), recalculates relative activity over external promoters and writes each gene's reads (raw counts of its non-internal promoters summed per sample). | <sim_dir>/star_output/*_SJ.out.tab, GENCODE v47 GTF → <sim_dir>/nofloor/{proactiv_result,recal_att_f,promoter_info}.RData, gene_reads.tsv, count_check.txt |

## benchmark/

| file | what it does | inputs → outputs |
|---|---|---|
| `benchmark/bohr2_sim.R` | Fits BOHR to every simulated promoter (power) or to 1,000 genotype permutations of the calibration set (calib); env SKIP_FRAC and SD_SKIP set the promoter skips. | args mode method dataset N out ncores [K]; env PERM_IDS, PERM_SEED_BASE, SKIP_FRAC, SD_SKIP; simulated data (datasets nothr2000, nf_b1, nf_b2) → results/e80_20261119/nofloor/<m1,b1,b2>/{power,calib*}/N<N>.rds |
| `benchmark/baselines_v3.R` | Computes QTLtools' nominal tests in R (raw, --normal, RINT-first) on BOHR's samples and permutations; takes prepare_y() and load_sim() from bohr2_sim.R, so the same rules apply. | args mode method dataset N out; env PERM_IDS, PERM_SEED_BASE → results/e80_20261119/nofloor/baselines_<M1,M2b1,M2b2>{,_seed<K>}/{power,calib}/N<N>.rds |
| `benchmark/combined_p_and_p0p1_v3.R` | Computes the per-promoter null FPR and QQ points of the BOHR p0 and p1 tests (Extended Data Fig. 2C, Supplementary Table 2c), and null tails of promoter-level combined P values. | permuted nulls (first part of prepare_data_v3.R; seed 1 for the E80 datasets) → data_v3_E80_nofloor_reads5/diag_tail/{p0p1_*,combined_p_*}.csv |
| `benchmark/prepare_data_v3.R` | Builds the Fig. 2A–D Source Data: n_mid >= 5 filter, mirror de-duplication, gene-reads filter, F correction, the combined BOHR P (COMBINE), calibrated thresholds, power, FPR, QQ, effect recovery. | BOHR and QTLtools result rds, design tables, recal_att_f.RData and gene_reads.tsv (DATASET=E80_nofloor_reads5, COMBINE=fisher_sidak) → data_v3_E80_nofloor_reads5_fishersidak/{power,thresholds,fpr_perpromoter,qq,tail_*,effect}.csv |
| `benchmark/prepare_power_by_effect.R` | Per-promoter calls and power by simulated effect size at N = 300, from the same runs, tests and calibrated thresholds as Fig. 2A. | first part of prepare_data_v3.R, thresholds.csv → data_v3_E80_nofloor_reads5_fishersidak/power_by_effect_N300/{per_promoter.csv, power_by_effect_bin.csv, CHECKS.txt} |
| `benchmark/prepare_power_by_OR.R` | Power by the per-allele odds ratio of the simulated effect (Fig. 2B, Supplementary Table 2f): replays the generators' parameters at each dosage and re-bins the per-promoter calls. | power_by_effect_N300/per_promoter.csv, ground truth and design tables (env OR_BINS) → power_by_OR_N300{,_lin5}/{per_promoter_OR.csv, power_by_OR_bin.csv, CHECKS.txt} |
| `benchmark/multiseed_final_v3.R` | Computes lambda and tail ratios for each of five permutation sets and pooled (Supplementary Table 2b). | calibration results of sets 1-5 (via the first part of prepare_data_v3.R; env MS_DATASET, COMBINE) → data_v3_E80_nofloor_reads5_fishersidak/seed_replicates/final_stats/lambda_tail_by_seed.csv |

## config/

| file | what it does | inputs → outputs |
|---|---|---|
| `INLINE_STEPS.md` | Records the two one-line R commands that drew the calibration genes (simulated activity; RNA-seq batch 1). | design tables → ids_calib_random150.txt, b1/ids_calib_random150genes.txt (config rows) |
| `config/m1_nothr_g2000/ids_calib_random150.txt` | Lists the 150 simulated-activity calibration genes (seed 20261001). | INLINE_STEPS.md section 1 → PERM_IDS of the simulated-activity calibration runs |
| `config/m2_nothr/b1/ids_calib_random150genes.txt` | Lists both promoters of the 150 RNA-seq batch-1 calibration genes (seed 20261001). | INLINE_STEPS.md section 2 → PERM_IDS of the batch-1 calibration runs |
| `config/m2_nothr/b2/ids_calib_random150genes.txt` | Lists both promoters of the 150 RNA-seq batch-2 calibration genes (seed 20261002). | runner_m2nothr_local.sh lines 10-11 → PERM_IDS of the batch-2 calibration runs |
