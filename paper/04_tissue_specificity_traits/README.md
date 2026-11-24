# Tissue specificity and donor traits

**Methods §6; Supplementary Methods §6.** Tissue-specific promoters: one representative mean per promoter and tissue
over the GTEx v11 analysis-freeze samples of the genotyped donors (`BOHR::tissue_means()`), ROKU calls in 56 tissues, and
an effect filter of 0.10 (Fig. 3A,B, Extended Data Fig. 3A, Supplementary Tables 3g,h). Donor traits: BOHR tests of μ, p0
and p1 against 20 donor traits in each tissue, with simulated and permuted arms for calibration (Fig. 3C–F, Extended Data
Fig. 3B, Supplementary Tables 3b–f).

**Run order, tissue specificity.** One script, `tissue_specificity/roku_freeze51.R`, computes the tissue means, the ROKU
calls and the 0.10 filter. It makes four runs: `freeze56` (the paper's calls: analysis-freeze samples, 56 tissues),
`freeze51`, `all56` and `all51` (the 51-tissue runs leave out the three stomach LCM sub-tissues and the two cell lines;
the `all` runs use every sample of the genotyped donors, as the earlier run did). Before the runs, check C0 runs ROKU on
the stored means of the earlier run and must give its stored calls, and check C1 runs `tissue_means()` on the earlier
run's samples and must give its stored means. In the final run C0 gave identical calls (22,939 promoters × 56 tissues),
C1 gave the 1,028,471 stored means to within 7 × 10⁻¹⁶, and `freeze56` used 19,626 samples and gave 13,634 promoters with
a call (8,815 with an increase, 5,085 with a decrease).

**Run order, donor traits.** All 20 traits were fitted in one run folder, `trait_regression_20261007`.
1. `traits/run_20261005/prep.R` (per-tissue donor × promoter matrices and donor metadata; it sources the `data.R` of run
   20261007 and writes into that run folder) → `traits/run_20261007/mirror.R` (mirror table of the cis-mapping units).
2. `traits/run_20261007/run_trait.R` for every tissue and arm (real, sim1, perm1): NSCC, one job per tissue running its
   three arms in turn, 32 cores.
3. Summary run folder `trait_regression_20261104` (links to the fits and inputs of run 20261007; `work/mirror_v2.tsv` is a
   copy of its `work/mirror.tsv`): `traits/run_20261104/mirror_v3.R` → `summarize.R real sim1 perm1` (mirror rule,
   testability rules, BH per tissue × trait × component) → `export_fig6.R`.
4. `traits/run_20261106_mirror_abs/01_abs_labels.R` and `03_apply_B.R` (exact-pair events labelled by absolute activity).
5. `traits/r5e80_20261119/r5_01_units.R`: the test set of the final version. A test is dropped when its promoter has at
   least 80% of donors at 0 or at 1 before endpoint deletion (`p0_raw`, `p1_raw` of `run_trait.R`), or when its gene's
   median first-intron junction reads in the tissue are below 5; then BH again per tissue × trait × family. Nothing is
   refitted.
6. `traits/r5e80_ruleb_20261119/r6_01_rule_b.R`: rare binary tests (minor class < 10% of the promoter's tested donors, the
   `trait_class` of `summarize.R`) are dropped, then BH again. Its output is the paper's result: 4,577 significant tests
   (μ 2,993, p0 1,084, p1 500).

**The trait run of the final version.** Steps 1–4 ran once on the no-floor data, in the run whose cis-mapping units drop
the SD ≤ 0.05 skip (run root `rerun_v4`; `mirror.R` and `mirror_v3.R` read its mirror tables). The final version (r5e80)
did not rerun them: each promoter is tested on its own, so its rules only remove tests, and BH is redone (steps 5 and 6).
The trait tests keep their own outcome screen (`run_20261007/data.R`, `prepare_y()`: a test is skipped when at least 95%
of the values are 0 or 1, or their SD is 0.05 or less).

**Changes for the final run.** Two scripts differ from the earlier run:
- `traits/run_20261005/prep.R` (line 11) sources the `data.R` of run 20261007 instead of run 20261005, so that one run
  folder holds all 20 traits with the drinking and smoking coding of run 20261007 (MHDRNKSTS, MHSMKSTS).
- `traits/run_20261106_mirror_abs/03_apply_B.R` (lines 24–26): the earlier run asserted its 3 relabelled exact-pair events
  (`stopifnot(nrow(R) == 3L)`); the final run reports its count instead (405 significant exact-pair events, 5 relabelled).
  The structural checks are unchanged. The text that the script writes into `export_B/README_trait_regression.txt` and
  `apply_B.log` still gives the earlier counts (586 events, 3 relabelled, 6 rows).

**One deviation from the code as run.** As run, line 20 of `tissue_specificity/roku_freeze51.R` sourced a copy of
`R/ROKU.R` of TCC 1.52.0 (Bioconductor), because TCC is not installed on the lab cluster. That file is under the GPL and
is not included here; the line now reads `ROKU <- TCC::ROKU`. Check C0 of the same run showed that the copy gives the calls
of `TCC::ROKU` exactly (the earlier run's calls, made on NSCC with `TCC::ROKU`). The C0 message names `TCC::ROKU`.

**Paths.** The trait scripts of steps 1–4 keep the paths of the earlier runs. In the final run they ran as copies in which
only folder names were changed:

| path in the scripts | folder in the final run |
|---|---|
| `/projects/group/bliu/e1327989/GTEx`, `/home/e1327989/project/GTEx` | `/projects/group/bliu/e1327989/GTEx/nofloor_full/GTEx` (lab) |
| the same, in the copies that ran on NSCC (`run_20261007/data.R`, `run_trait.R`) | `/scratch/users/nus/e1327989/nofloor_full/GTEx` and `…/nofloor_full/tissue_regression` |
| `…/GTEx/QTL_figure/rerun_v2` (mirror tables read by `mirror.R` and `mirror_v3.R`) | `…/nofloor_full/GTEx/QTL_figure/rerun_v4` |

On NSCC the jobs set `TRG_RUN=/scratch/users/nus/e1327989/nofloor_full/tissue_regression/trait_regression_20261007_v4`.
The scripts of steps 5–6 and `roku_freeze51.R` were written for the no-floor data and carry its paths
(`…/GTEx/nofloor_full/…`); `roku_freeze51.R` loads the BOHR package from the repository root with `pkgload::load_all()`.

**Inputs without a script here.** `r5_01_units.R` reads the gene medians of first-intron junction reads,
`nofloor_full/figures_v4/fig1c/gene_median_all57.tsv`, which `nofloor_full/scripts/gene_depth_all57.R` computed on NSCC
over all samples of each tissue's proActiv run. The check inputs of `roku_freeze51.R` (C0, C1) are the outputs of the
earlier tissue-specificity run on all samples (`nofloor_full/GTEx/fit/`); the check of `r5_01_units.R` (lines 66–71)
reads the significant tests of an earlier draft run of the same rule, and that of `r6_01_rule_b.R` (lines 94–99) a
comparison of FDR rules (`results/20261119c_fdr_options`).

**Counts written into the scripts.** `r5_01_units.R` and `r6_01_rule_b.R` assert equality with the check inputs above.
`03_apply_B.R` prints and writes the earlier counts noted above. No other script of this step asserts a count.

**BOHR package.** The paper's tissue means are computed with `BOHR::tissue_means()` itself (`roku_freeze51.R`).
`BOHR::bohr_trait()` reproduces the trait tests of `run_trait.R` (repository `workflow/parity/`, scripts 07–09).

`run.sh` lists the commands in the order they were run.


## tissue_specificity/

| file | what it does | inputs → outputs |
|---|---|---|
| `tissue_specificity/roku_freeze51.R` | Computes one mean per promoter and tissue with BOHR::tissue_means(), the ROKU calls (TCC::ROKU, upper.limit = 0.25, sort = TRUE) and the 0.10 effect filter, for four sample and tissue sets; checks C0 and C1 first. | no-floor relative activity (v11/filter_5_recal_att_filter.RData), v9 VCF header (genotyped donors), sample attributes (SMAFRZE), the earlier run's fit/ outputs (checks) → <out>/{freeze56,freeze51,all56,all51}/v11_raw_fit_fit_results.txt, v11_tissue_specific_ROKU_results.rds, roku_20261007/{v11_tissue_specific_ROKU_results_delta010.rds, summary_delta_filter.tsv}; summary_runs.tsv, samples_per_tissue.tsv |

## traits/

| file | what it does | inputs → outputs |
|---|---|---|
| `traits/run_20261005/prep.R` | Builds per-tissue donor × promoter activity matrices (one analysis-freeze sample per donor) and donor metadata, in the run folder of 20261007. | v11 promoter activity, sample attributes (SMAFRZE, SMRIN), v11 TPM ≥ 1 gene list, donor covariates → work/<tissue>.rds, meta.rds, tissues.tsv |
| `traits/run_20261007/data.R` | Defines the trait regression's paths, 20 donor traits, covariates, covariate design and outcome screening rules. | sourced by prep.R, mirror.R and run_trait.R (PATHS, TRAITS, COVS, read_meta, cov_design, trait_contrast, prepare_y) |
| `traits/run_20261007/mirror.R` | Marks exact two-promoter mirror pairs per tissue from the mirror tables of cis mapping. | QTL_figure mirror/promoter_flags.tsv.gz and mirror_pairs.tsv.gz, work/tissues.tsv, work/<tissue>.rds → work/mirror.tsv, mirror_check.tsv |
| `traits/run_20261007/run_trait.R` | Fits the BOHR trait tests (mu, p0, p1) for every promoter of one tissue in one arm (real, simulated or permuted), all 20 traits. | work/<tissue>.rds, work/meta.rds, work/tissue_list.txt → results/{real,sim1,perm1}/<tissue>/<trait>.rds |
| `traits/run_20261104/data.R` | Defines paths and trait coding of the summary run folder; sourced by its scripts. | sourced by mirror_v3.R, summarize.R, export_fig6.R, 01_abs_labels.R, 03_apply_B.R |
| `traits/run_20261104/mirror_v3.R` | Applies mirror rule v3 to exact two-promoter pairs: keeps p0 on both members and drops p1. | work/mirror_v2.tsv (copy of mirror.R output), QTL_figure mirror_v3/promoter_flags.tsv.gz → work/mirror_v3.tsv, mirror_v3_check.tsv |
| `traits/run_20261104/summarize.R` | Applies mirror rule v3 and the testability rules to the trait-test units; BH per tissue × trait × component; labels rare binary tests (trait_class). | results/{real,sim1,perm1}/*/*.rds, work/mirror_v3.tsv → results/summary/units_*.rds, sig_*.tsv, counts_*.tsv, null_discovery.tsv |
| `traits/run_20261104/export_fig6.R` | Exports every testable real-arm trait test. | results/summary/units_real.rds → results/export/trait_regression_units.tsv.gz, README_trait_regression.txt |
| `traits/run_20261106_mirror_abs/01_abs_labels.R` | Labels each exact-pair boundary event by which promoter's absolute activity changes with the trait. | units_real.rds, work/<tissue>.rds, work/meta.rds, v11/tabs_att_filter_In.rds (absolute activity) → results/mirror_abs_labels.tsv.gz |
| `traits/run_20261106_mirror_abs/03_apply_B.R` | Applies the absolute-activity labels to the significant exact-pair events; rewrites export and summary. | mirror_abs_labels.tsv.gz, results/summary/ (units, sig_real, counts, null_discovery), results/export/ → results/export_B/, results/summary_B/ |
| `traits/r5e80_20261119/r5_01_units.R` | Drops the tests of the final version's rules (80% of donors at 0 or 1; gene median reads < 5) and reruns BH, in the real arm and both null arms. | export_B/trait_regression_units.tsv.gz, results/real/*/*.rds (p0_raw, p1_raw), results/summary/units_{sim1,perm1}.rds, gene_median_all57.tsv → data_r5e80_20261119/{trait_regression_units_r5e80.tsv.gz, units_{sim1,perm1}_r5e80.rds, counts_real_r5e80_*.tsv, null_discovery_r5e80.tsv, summary_r5e80.tsv} |
| `traits/r5e80_ruleb_20261119/r6_01_rule_b.R` | Drops rare binary tests and reruns BH (real arm per tissue × trait × family, null arms per tissue × trait × type); the paper's trait results. | data_r5e80_20261119/ → data_r5e80_ruleb_20261119/{trait_regression_units_r5e80_ruleb.tsv.gz, units_{sim1,perm1}_ruleb.rds, counts_real_ruleb_*.tsv, null_discovery_ruleb.tsv, summary_ruleb.tsv, methods_facts_ruleb.tsv, cells_ruleb.tsv, characteristics_ruleb.tsv} |
| `config/tissue_list.txt` | Lists the 57 tissue names of the trait regression; the order sets the permutation and simulation seeds. | read by prep.R and run_trait.R through data.R PATHS$tissues (RUN/work/tissue_list.txt) |
