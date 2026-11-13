# The BOHR model and simulation benchmarking

**Methods §2–3; Supplementary Methods §2–3.** `engine/` is the BOHR likelihood used by every analysis of the paper: the
zero-and-one-inflated beta model fitted as a beta regression on the values between 0 and 1 and a multinomial logit model
for 0, between and 1; likelihood-ratio tests of μ (with the small-sample F correction) and Firth-penalised
likelihood-ratio tests of p0 and p1. `gamlss_check/` compares the engine with gamlss refits. The simulations measure
calibration and power against QTLtools-style linear tests (Fig. 2, Extended Data Fig. 2, Supplementary Table 2).

**Run order.** Simulated activity: `simulated_activity/newsim_direct_nothreshold_indep.R`. Simulated RNA-seq:
`simulated_rnaseq/v12_simulation_nothr.R` (generator) run by `v12_sim_run_nothr_big.R`, then STAR, the five-read junction
filter and `sim_proactiv_realjunc.R` (proActiv). Fits: `benchmark/bohr2_sim.R` (BOHR) and `benchmark/baselines_v3.R`
(linear tests); calibration genes in `config/`, drawn as in `INLINE_STEPS.md`. Summaries: `benchmark/prepare_data_v3.R`
(power, thresholds, FPR), `benchmark/multiseed_final_v3.R` (genomic inflation), `benchmark/combined_p_and_p0p1_v3.R`
(null FPR of the p0 and p1 tests) and `benchmark/prepare_power_by_effect.R` (power by effect size, Fig. 2B).

**BOHR package.** `BOHR::bohr_fit()`, `BOHR::bohr_test()` and `BOHR::simulate_response()` implement this engine; the
package's tests compare it with stored engine fits (`tests/testthat/test-engine-parity.R`).

**Note.** `engine/bohr2_data.R` holds the input paths and the promoter rules used by the cis-mapping scripts of
`../03_cis_mapping_conditional/`, which source the three `bohr2_*.R` files from one folder (`BOHR2_CODE`).

`run.sh` lists the commands in the order they were run (paths of the original runs).


## engine/

| file | what it does | inputs → outputs |
|---|---|---|
| `engine/bohr2_engine.R` | Implements the BOHR likelihood: split zero/one-inflated beta fits, null fits and per-SNP LRT and Wald tests. | sourced library (base R only) |
| `engine/bohr2_addon3.R` | Adds Firth-penalised p0/p1 likelihood-ratio tests (multi-start, alternating search, BFGS polish). | sourced library (needs bohr2_engine.R) |
| `engine/bohr2_data.R` | Defines BOHR data handling: input paths, endpoint rules and family choice, cis window, genotypes, covariates. | sourced library; reads v11 activity matrix, promoter info, genotype PCs, PEER factors, v9 VCF |
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
| `simulated_rnaseq/sim_proactiv_realjunc.R` | Runs proActiv on the re-filtered simulated junctions and recalculates relative activity over external promoters. | star_output/*SJ.out_realfilter.tab, GENCODE v47 GTF → realjunc/{proactiv_result,recal_att_f,promoter_info}.RData, count_check.txt |

## benchmark/

| file | what it does | inputs → outputs |
|---|---|---|
| `benchmark/bohr2_sim.R` | Fits BOHR to every simulated promoter (power) or to 1,000 genotype permutations of the calibration set (calib). | args mode method dataset N out ncores [K]; env PERM_IDS, PERM_SEED_BASE; simulated data → results/<set>/{power,calib*}/N<N>.rds |
| `benchmark/baselines_v3.R` | Computes QTLtools' nominal tests in R (raw, --normal, RINT-first) on BOHR's samples and permutations. | args mode method dataset N out; env PERM_IDS, PERM_SEED_BASE → results/baselines_v3*/<tag>/{power,calib}/N<N>.rds |
| `benchmark/prepare_data_v3.R` | Builds the Fig. 2 Source Data: n_mid >= 5 filter, mirror de-duplication, F correction, Simes P, power and FPR. | BOHR and QTLtools result rds, design tables, recal_att_f.RData (DATASET=nothr_realjunc) → data_v3_nothr_realjunc/{power,thresholds,fpr_perpromoter,qq,tail_*,effect}.csv |
| `benchmark/multiseed_final_v3.R` | Computes lambda and tail ratios for each of five permutation sets and pooled (Supplementary Table 2b). | calibration results of sets 1-5 (via prepare_data_v3.R lines 1-82) → seed_replicates/final_stats/lambda_tail_by_seed.csv |
| `benchmark/combined_p_and_p0p1_v3.R` | Computes the per-promoter null FPR and QQ points of the BOHR p0 and p1 tests (Extended Data Fig. 2B, Supplementary Table 2c). | permuted nulls of sets 1-5 (header of prepare_data_v3.R) → diag_tail/p0p1_{fpr_perpromoter,qq}_seed1.csv, combined_p_*.csv |
| `benchmark/prepare_power_by_effect.R` | Power by simulated effect size at N = 300 (Fig. 2B), from the final runs and the permutation-calibrated thresholds of prepare_data_v3.R. |  |

## config/

| file | what it does | inputs → outputs |
|---|---|---|
| `INLINE_STEPS.md` | Records the two one-line R commands that drew the calibration genes (simulated activity; RNA-seq batch 1). | design tables → ids_calib_random150.txt, b1/ids_calib_random150genes.txt (config rows) |
| `config/m1_nothr_g2000/ids_calib_random150.txt` | Lists the 150 simulated-activity calibration genes (seed 20261001). | INLINE_STEPS.md section 1 → PERM_IDS of the simulated-activity calibration runs |
| `config/m2_nothr/b1/ids_calib_random150genes.txt` | Lists both promoters of the 150 RNA-seq batch-1 calibration genes (seed 20261001). | INLINE_STEPS.md section 2 → PERM_IDS of the batch-1 calibration runs |
| `config/m2_nothr/b2/ids_calib_random150genes.txt` | Lists both promoters of the 150 RNA-seq batch-2 calibration genes (seed 20261002). | runner_m2nothr_local.sh lines 10-11 → PERM_IDS of the batch-2 calibration runs |
