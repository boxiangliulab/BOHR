# cis mapping, comparison with QTLtools and independent signals

**Methods §4, §5 and §7; Supplementary Methods §4, §5 and §7.** For each of 50 tissues: inputs (one sample per donor,
five genotype PCs, PEER factors capped by sample size), a per-SNP BOHR scan of every promoter within 500 kb (μ, p0, p1),
units (promoter × component) with window P values combined by ACAT, estimability and mirror rules, BH per tissue, then
forward conditional selection with a backward step. The final run gives 28,302 puQTLs with 30,043 independent signals.

**Final rules.** The final run used the activity of `../01_promoter_quantification/` (every junction, no read floor).
A promoter enters the scan when, on the samples of the tissue's genotyped donors, its SD is ≥ 0.005, its share of
missing values < 0.4 and its gene's mean TPM ≥ 1, and when it is not skipped (n ≥ 20 and n ≥ 5 + k; shares at 0 and at 1
each < 0.95; SD ≥ 0.005, `prepare_y()` in `../02_bohr_model_and_simulation/engine/bohr2_data.R`). The earlier rule also
skipped SD ≤ 0.05; the final run does not. Two rules are then applied to the units: a promoter is not tested in a tissue
when at least 80% of its values among the tissue's analysed samples lie at one boundary, or when its gene has a median of
fewer than five first-intron junction reads among those samples (`units/build_r5e80_units.R`, with the exclusions of
`units/make_version_inputs.R`); BH is then redone within each tissue. The estimability rules (μ: at least 80 samples
between 0 and 1; p0 and p1: a boundary share of 10–90%) are unchanged.

**Run order.**
1. `inputs/`: `peer/peer_analysis.R` (PEER factors per tissue on NSCC, from the PEER input written by
   `../01_promoter_quantification/export_txt.R`; the number of factors per tissue is the chromosome 22 scan of the earlier
   run, `config/tissue_list_best_peer_table.txt`, made with `peer/gamlss_QTL_peer_test.R` and `peer/peer_plot.R`) →
   `promoter_sets_rule.R` (the per-tissue promoter sets: first with the SD ≤ 0.05 skip, then with `NF_SD_SKIP=none`, the
   final sets) → `prep_freeze_nofloor.R` (per-tissue input objects; a copy of `prep_all.R` that takes the sets from the
   rule) → `build_v4.R` (checks the final input objects against the first ones and writes the delta: the promoters new in
   the final sets and the promoters that the SD ≤ 0.05 skip removed on the scan samples).
2. `scan/`: `run_final2.R` and `run_firth_addon.R` per tissue, arm (real, perm0–perm2) and chromosome on NSCC: first the
   first sets (engine with the SD ≤ 0.05 skip), then the delta with the final engine; `merge_v4.R` joins the two scans per
   arm; `run_muF.R` adds the F-corrected μ P values.
3. `units/`: `mirror_pairs.R`, `discovery_nogate.R` (units, ACAT, BH; four arms), `units_v3_estimable.R`,
   `units_v3m_mirror.R` (mirror rule; the delta units enter BH only for the admitted types, μ and p0). The admission was
   decided by rules fixed before any permuted delta result existed: `v4_strata_sizes.R` (stratum sizes),
   `diag_v4_fit.R` (fit diagnostics) and `diag_v4_calibration.R` (permutation calibration). Then the final units:
   `depth_extract.R` (raw junction reads per promoter and gene, NSCC) → `compare_units.R` (median gene reads per tissue and
   promoter, `promoter_depth_v4.tsv.gz`; the rest of the script compares the earlier and the no-floor run and feeds no
   result) → `make_version_inputs.R` (the promoters each rule set excludes; the final set is `nofloor_reads5_E80`) →
   `build_r5e80_units.R` (units v3m of the final rules).
4. `conditional/`: a fresh run from the leader of every significant unit: `cond_fresh_list.R` (unit list and the stubs the
   merges need) → `run_conditional_v3m.R` per (tissue, chromosome), through `run_cond_one.sh` →
   `merge_conditional_v3m.R` (backward step) → `merge_conditional_v3c.R` with `LABELS=none` → `units/mirror_abs_labels.R`
   (labels of exact-pair boundary events by absolute activity) → `units/units_v3b_from_labels.R` (units v3c, P and q of
   v3m) → `merge_conditional_v3c.R` with the labels (the signal table). Then the per-signal tables read downstream (delivery
   v6): `fig5_delivery_v4.R`, `fig5_proxies.R`, `fig5_pertype.R`, `fig5_pe_moments_v6.R`, `crosstype_conditional_v6.R`,
   `crosstype_conditional_v6x.R`, `fig5_signals_partner.R`; and `fig5_delivery_v7.R`, which marks each puQTL's lead signal
   (the first signal kept after the backward step) for Fig. 4 and Supplementary Tables 4–5.
5. `ld/`: `ld_tool.R` (the LD standard: PLINK 1.9 r on each tissue's donors) and `ld_v3.R` (donor lists, signal LD).
6. `qtltools_comparison/` (Fig. 2E and 2F; Supplementary Tables 2d and 2e): on the final run's root, `make_freeze_beds.R`,
   QTLtools v1.2 `cis --permute 1000` (raw, rank-normal first, `--normal`) on every tissue and arm, and
   `fig2e_deliver_freeze_v3m_pairs.R`, which forms the comparison units; `fig2e_r5e80_units.R` keeps the units that pass the
   final promoter rules. QTLtools was then run again on the promoter–variant pairs that BOHR tested (`make_site_lists.R`,
   `exact_qtltools.R`, `keep_qtltools.R`, `pairs_qtltools.R`, `fill_qtltools.R`); the no-fill tables, in which the donors
   without a value are left out, are the QTLtools results of Fig. 2E. `class_calib.R` (`UNIT_MODE=merge_all`),
   `dd2_arms_units.R`, `dd11_fisher_rules.R` (`PART=sidak`) and `prepare_panel_e_classes.R` (`PART=sidak`) give BOHR's calls
   by rule BF, the classes of Fig. 2E, and the permuted calls and chance range of Fig. 2F. The perm1 and perm2 maps were
   drawn by `make_perm_maps_freeze_extra.R` in the earlier run and reused.

**Paths and roots.** The scripts keep the paths of the earlier run. The final run used copies in which only the folder
prefixes and the run root differ: `/projects/group/bliu/e1327989/GTEx` and `/home/e1327989/project/GTEx` became
`/projects/group/bliu/e1327989/GTEx/nofloor_full/GTEx` on the lab and `/scratch/users/nus/e1327989/nofloor_full/GTEx`
on NSCC, and the run root `QTL_figure/rerun_v2` became one of three roots: `rerun_v2` (first promoter sets, their scans),
`rerun_v4` (final promoter sets, merged scans and the units of all four arms) and `rerun_r5e80` (units of the final rules,
the conditional analysis and the per-signal tables; it links to the scans, unit tables and inputs of `rerun_v4`). The
scripts written for the final run (`prep_freeze_nofloor.R`, `build_v4.R`, `merge_v4.R`, the scripts of the final units
and of the admission, `cond_fresh_list.R`, and those of `qtltools_comparison/` from `fig2e_r5e80_units.R` on) keep the paths
they ran with. `../02_bohr_model_and_simulation/engine/bohr2_data.R` sets the input paths (`BOHR2_GTEX`, `BOHR2_ROOT`,
`BOHR2_VCF`); scripts source the engine from `BOHR2_CODE`, and write under `BOHR2_ROOT/results/` unless `FINAL_ROOT` or
`COND_ROOT` is set.

**Code changes in the copies that ran** (kept in the files here): `units/discovery_nogate.R` makes the Firth columns
numeric before filling them from the add-on (an all-missing column is logical, and the fill would store TRUE) and carries
the lead variant's standard error (`leader_se`); `conditional/fig5_pertype.R` has the same numeric fix;
`units/units_v3m_mirror.R` reads the delta list and keeps the delta's p1 units out of BH; `units/mirror_abs_labels.R`
runs its comparison with a check table of the earlier run only when that table exists; `ld/ld_v3.R` reads the delivery v6
signals and stops before a comparison with a PLINK table of the earlier run, after writing every output;
`conditional/fig5_delivery_v7.R` reports the earlier run's counts instead of asserting them.

**Carried over from the earlier run.** The permutation maps (perm0 from `inputs/prep_all.R`, perm1 and perm2 from
`qtltools_comparison/make_perm_maps_freeze_extra.R`; the donors are the same and `prep_freeze_nofloor.R` checks it), the
PEER factor count per tissue, and two tables of pair IDs only (`work/crosstype_pair_ids_v1.tsv`,
`results/fig5_delivery_v3/crosstype_conditional_pairs.tsv.gz`), which keep the pair IDs of the earlier deliveries.
`inputs/prep_all.R` is the earlier run's input builder: it took the promoter sets from the result files of the GAMLSS scan
(`QTL/new_v11_LRT`), which applied the rule written out in `inputs/promoter_sets_rule.R` (on the five-read data the rule
gives the earlier sets in all 50 tissues, 131,361 tissue × promoter inputs).

**BOHR package.** `BOHR::bohr_cis()` reproduces steps 2–3 and `BOHR::bohr_conditional()` step 4 in one call each: on the
stored inputs of the earlier run the per-SNP statistics, unit P values, q values, signal counts and lead variants agree
with this pipeline (repository `workflow/parity/`, scripts 01–06). The package's default skips (`R/control.R`:
`min_sd = 0.05`, `max_boundary = 0.95`) are the earlier rules, not the final ones.

`run.sh` lists the commands in the order they were run (paths of the original runs).


## inputs/

| file | what it does | inputs → outputs |
|---|---|---|
| `inputs/peer/0104_peer.R` | Records how the PEER input was built (commented-out lines 48-57: donor means within each tissue, missing set to 0). | filter_5_recal_att.txt (text export of 0103's recalculated activity) → filter_5_recal_att_filter_SMTSD_mean.txt (lines 48-57, commented out); PEER part → v10_* files (not used) |
| `inputs/peer/peer_analysis.R` | Estimates 60 PEER factors for one tissue from the donor-mean relative-activity matrix (final run: one NSCC job per tissue). | tissue, filter_5_recal_att_filter_SMTSD_mean.txt, output dir → v11_filter_5_peer_factors_<tissue>.txt |
| `inputs/peer/gamlss_QTL_peer_test.R` | Runs the earlier gamlss LRT scan for one tissue on chr22 (+-10 kb) with k PEER factors and 5 genotype PCs (sets K0, earlier run). | v11 matrix, Promoter_info_filter_5.txt, genes > 1 TPM, genotype PCs, PEER factors, chr22 VCF → <name>.<tissue>.chr22.factor_peer_<k>_pca_5.results.rds |
| `inputs/peer/peer_plot.R` | Counts, per tissue and k, the promoters significant at q <= 0.05 after ACAT over their chr22 variants. | ../tissue_list_size20_table.txt; <tissue>/sensitivity_test_10k.*.rds → <tissue>_peer_test_count.txt, <tissue>_peer_test_count.pdf |
| `inputs/gamlss_QTL_optimized.R` | Runs the published GAMLSS puQTL scan; its promoter filter (SD, NA share, TPM, lines 115-137) is the filter of the promoter-set rule. | v11 activity, promoter info, PCs, PEER, TPM1_genes.txt, VCF → QTL/new_v11_LRT/<T>/chr<N>/*.results.rds |
| `inputs/gamlss_functions_optimized.R` | Provides the published GAMLSS functions that `gamlss_QTL_optimized.R` sources; its lines 548–557 hold the skip rules that the promoter-set rule reproduces. (Earlier in `04_tissue_specificity_traits/`, where its `raw_fit()` gave the tissue means of the earlier version.) | sourced library |
| `inputs/promoter_sets_rule.R` | Writes the per-tissue promoter sets of cis mapping by the rule above; `NF_SD_SKIP=none` drops the SD <= 0.05 skip (final sets); `NF_TISSUES` runs a subset. | GTEx root (best-PEER table, PCs, TPM genes), relative activity (.txt or .RData), PEER folder → tissue, promoter table |
| `inputs/prep_freeze_nofloor.R` | Builds the per-tissue BOHR input objects of the final run (freeze samples, one per donor, promoter set, PEER cap K) and copies the earlier run's permutation maps after checking the donors. | promoter-set table, v11 activity and promoter info, sample attributes, PCs, PEER, best-PEER table → work/final/freeze/<T>.rds, work/final/tissues.tsv, perm_maps/final/freeze/ |
| `inputs/build_v4.R` | Checks the final input objects against the first ones (same donors, K, covariates, activity) and writes the delta to scan (new promoters, promoters skipped for SD <= 0.05 on the scan samples). | rerun_v2 and rerun_v4 work/final/freeze/<T>.rds → rerun_v4_delta/work/final/freeze/<T>.rds, perm maps, core/v4/delta_promoters.tsv, delta_summary.tsv |
| `inputs/prep_all.R` | Earlier run: builds the per-tissue BOHR inputs and the perm0 permutation maps (reused by the final run) from the published promoter lists. | v11 activity, sample attributes, PCs, PEER, best-PEER table, published promoter lists → work/final/<def>/<T>.rds, work/final/tissues.tsv, perm0 maps |
| `inputs/config/tissue_list_size20_table.txt` | Lists the 57 GTEx v11 tissues with at least 20 samples and their sample counts. | hand-made → read by peer_plot.R line 9 |
| `inputs/config/tissue_list_size20_table_add2.txt` | Lists the four small tissues of the last PEER-scan submission and their sample counts. | hand-made → read by submit_peer_jobs.sh line 4 |
| `inputs/config/tissue_list_best_peer_table.txt` | Gives the best PEER factor count per tissue (chr22 scan); prep_freeze_nofloor.R and promoter_sets_rule.R read it as k. | read from GTEx/QTL/tissue_list_best_peer_table.txt |

## scan/

| file | what it does | inputs → outputs |
|---|---|---|
| `scan/run_final2.R` | Runs the main per-SNP BOHR scan (mu, p0, p1; Firth p0/p1 when f0 or f1 >= 0.35) per tissue, arm and chr. | work/final/<def>/<T>.rds, perm_maps/final/<def>/<T>.<arm>.perm_map.tsv, v9 VCF → results/final/<def>/<T>/<arm>/chr<N>.rds |
| `scan/run_firth_addon.R` | Adds the Firth p0/p1 tests for every other promoter with a 0/1 part, on the SNPs of the main scan. | results/final/<def>/<T>/<arm>/chr*.rds, work/final/<def>/<T>.rds, perm maps → results/final_firth_addon/<def>/<T>/<arm>/chr<N>.rds |
| `scan/merge_v4.R` | Joins, per tissue, arm and chromosome, the scan of the first promoter sets (rows outside the delta) and the scan of the delta, for the main run and the Firth add-on. | rerun_v2 and rerun_v4_delta results/{final,final_firth_addon} → rerun_v4/results/{final,final_firth_addon}, logs/merge_v4_<arm>.tsv |
| `scan/run_muF.R` | Adds the small-sample F-corrected mu P value per SNP from the stored beta-part LRT. | results/final/<def>/<T>/<arm>/chr*.rds, work/final/<def>/<T>.rds → results/final_muF/<def>/<T>/<arm>/chr<N>.rds |

## units/

| file | what it does | inputs → outputs |
|---|---|---|
| `units/mirror_pairs.R` | Builds the exact mirror-pair table (per-component exactness) and promoter flags for the 50 tissues. | work/final/freeze/*.rds → results/mirror/{mirror_pairs.tsv.gz, promoter_flags.tsv.gz} |
| `units/discovery_nogate.R` | Forms units (promoter x mu/p0/p1), combines window P values by ACAT, drops mirror twins, runs BH; carries the lead variant's standard error. | results/{final,final_muF,final_firth_addon}/<def>/<T>/<arm>/, results/mirror/promoter_flags.tsv.gz → results/final_nogate/<def>/<T>.<arm>.units.rds |
| `units/units_v3_estimable.R` | Applies the estimability rules (mu n_mid >= 80; p0/p1 share 10-90%) before BH and recomputes q per tissue. | results/final_nogate/<def>/*.units.rds → results/final_nogate_v3/<def>/*.units.rds, freeze_thresholds_v2_v3.tsv |
| `units/units_v3m_mirror.R` | Applies mirror rule v3 (exact boundary events held as p0), marks the delta's p1 units as not tested, reruns BH per tissue and writes checks (realized FDR on perm0-2). | results/final_nogate_v3/freeze/*.units.rds, results/mirror/promoter_flags.tsv.gz, work/v4_admission/delta_promoter_features.tsv → results/final_nogate_v3m/freeze/, results/mirror_v3/ |
| `units/v4_strata_sizes.R` | Counts, from the phenotypes only, the units the delta adds per type and stratum, and writes the per-promoter features of the delta that units_v3m_mirror.R reads. | rerun_v4_delta freeze objects, core/v4/delta_promoters.tsv → core/v4/strata_sizes.tsv, core/v4/delta_promoter_features.tsv |
| `units/diag_v4_fit.R` | Fit diagnostics of the delta units (real arm): share of converged SNP tests and null models, per type. | rerun_v4_delta scans, rerun_v2 scans (chr20-22, reference) → v4_diag/fit_*.tsv |
| `units/diag_v4_calibration.R` | Applies the three admission rules to the delta units of each type (fit diagnostics, permutation calibration on perm0-2, extreme tail) and writes the decision. | units v3m of rerun_v4 (four arms), delta_promoter_features.tsv, fit_diagnostics.tsv → core/v4/calibration_decision.tsv, calibration_all.tsv, real_added_counts.tsv |
| `units/depth_extract.R` | Extracts the raw proActiv counts of each tissue's promoters and the gene totals over all non-internal promoters, for the scan samples (sample-level, kept private). | per-tissue proActiv rds of the final run, Promoter_info_all.RData, rerun_v4 freeze objects → depth/<T>.rds, depth/extract_checks.tsv |
| `units/compare_units.R` | Writes the median gene reads per tissue and promoter (depth_stats(), the only output used downstream); the rest compares the units of the earlier run and of the no-floor run. | freeze objects of both runs, depth/<T>.rds, units v3c of both runs → out/promoter_depth_v4.tsv.gz, comparison tables |
| `units/make_version_inputs.R` | Lists, for three rule sets, the promoters excluded by the 80% boundary-share rule and the gene-reads rule, and the units kept, new and lost after BH again; nofloor_reads5_E80 is the final set. | freeze objects, promoter_depth_v4.tsv.gz, units v3c, fig2e_freeze_v3m → figures_v4/inputs/<version>/{promoters_excluded, units_keep, units_new_topup, signals_keep, fig2e/} |
| `units/build_r5e80_units.R` | Marks the excluded promoters' units as not tested and recomputes BH within each tissue, giving the units v3m of the final rules. | rerun_v4 results/final_nogate_v3m/freeze, promoters_excluded.tsv.gz → rerun_r5e80/results/final_nogate_v3m/freeze/<T>.<arm>.units.rds, build_r5e80_units.tsv |
| `units/mirror_abs_labels.R` | Labels each significant exact-pair boundary event by the TSS whose absolute activity changes. | final_nogate_v3m, conditional_final_v3c signals (env COND_FINAL), v11 tabs_att_filter_In.rds → results/mirror_abs_labels_v3c/labels_real.tsv.gz |
| `units/units_v3b_from_labels.R` | Adds the saved label columns to the v3m units without refitting or BH, giving units v3c (P and q identical to v3m). | final_nogate_v3m, labels_real.tsv.gz (env LABELS, UNITS_OUT) → results/final_nogate_v3c/freeze/<T>.<arm>.units.rds, checks.txt |

## conditional/

| file | what it does | inputs → outputs |
|---|---|---|
| `conditional/cond_fresh_list.R` | Lists every significant unit of units v3m for a run from its leader, the (tissue, chromosome) jobs, and the stubs that let the merges run without an earlier conditional analysis. | results/final_nogate_v3m/freeze/*.real.units.rds, work/final/freeze → work/cond_v3m/{units_to_run.tsv, jobs.tsv}, results/conditional_final_v3/ (stub), work/cond_v3c/units_to_run.tsv (header only) |
| `conditional/run_cond_one.sh` | Runs run_conditional_v3m.R for one (tissue, chromosome) with its own log. | tissue, chromosome; env RS, CODE, NC, LOGDIR → as run_conditional_v3m.R |
| `conditional/run_conditional_v3m.R` | Runs forward conditional selection and leave-one-out scans for the listed units, per tissue and chromosome. | work/cond_v3m/units_to_run.tsv, final_nogate_v3m, results/final, work/final → results/conditional_v3m/<def>/<T>/real/<chr>.rds |
| `conditional/merge_conditional_v3m.R` | Merges the forward runs, recomputes the backward step and writes the signal tables. | results/conditional_v3m, results/conditional_final_v3 (stub), final_nogate_v3m → results/conditional_final_v3m/freeze_real_*.tsv.gz, loo/ |
| `conditional/merge_conditional_v3c.R` | Applies the saved mirror labels to the signal tables (and replaces re-run units, none in the final run), giving the 30,043 independent signals. | conditional_final_v3m, conditional_v3c, labels (env LABELS) → results/conditional_final_v3c/{freeze_real_units,signals,rounds}.tsv.gz, loo/ |
| `conditional/fig5_delivery_v4.R` | Builds the per-signal table (signals.tsv.gz) that the Figure 4 data builders and Supplementary Tables 4-5 read. | env COND_FINAL, UNITS_DIR, FIG5_OUT; work/final, results/final → <FIG5_OUT>/{units,signals,genotype_pheno,crosstype_ld,boundary_only_partner_pheno}.tsv.gz |
| `conditional/fig5_proxies.R` | Per-genotype phenotype at the LD proxies (r2 >= 0.8, PLINK rule) of every signal, for the per-signal tables read by the major/minor analysis. | <FIG5_OUT>/signals.tsv.gz, work/final, v9 VCF → <FIG5_OUT>/genotype_pheno_proxies.tsv.gz |
| `conditional/fig5_pertype.R` | Per-component marginal statistics (mu, p0, p1) at every signal SNP and LD proxy, looked up from the main-run SNP tables. | <FIG5_OUT> signals and proxies, results/{final,final_muF,final_firth_addon} → <FIG5_OUT>/{signals,proxies}_pertype.tsv.gz, boundary_only_partner_pheno_proxies.tsv.gz |
| `conditional/fig5_pe_moments_v6.R` | Probability-scale effects (dosage 0 to 2, delta-method SEs; engine/pe_delta.R) and interior-only moments by genotype for every signal. | <FIG5_OUT> signals, work/final, results/final → <FIG5_OUT> effect and moment tables |
| `conditional/crosstype_conditional_v6.R` | Conditional cross-type tests for pairs of mu and boundary signals of the same promoter (the pair table that 71_A1_location.R reads). | <FIG5_OUT> signals and crosstype_ld, work/crosstype_pair_ids_v1.tsv (IDs) → <FIG5_OUT>/crosstype_conditional_{snps,pairs}.tsv.gz |
| `conditional/crosstype_conditional_v6x.R` | The same tests for cross-member pairs of exact mirror pairs (mu of one promoter x p0 of its partner), and the joined pair table crosstype_conditional_pairs_all. | <FIG5_OUT> tables, results/fig5_delivery_v3/crosstype_conditional_pairs.tsv.gz (IDs) → <FIG5_OUT>/crosstype_conditional_crossmember_*.tsv.gz, crosstype_conditional_pairs_all.tsv.gz |
| `conditional/fig5_signals_partner.R` | Adds SNP position, TSS, strand and mirror-pair columns per signal (signals_partner.tsv.gz) for Fig. 4C and Supplementary Tables 4-5. | <FIG5_OUT>/signals.tsv.gz, mirror_v3 promoter_flags (env FLAGS), v11 Promoter_info_filter_5.txt → <FIG5_OUT>/signals_partner.tsv.gz |
| `conditional/fig5_delivery_v7.R` | Marks each puQTL's lead signal (the first signal kept after the backward step) and its marginal P; delivery v7 = delivery v6 plus these two columns. | results/fig5_delivery_v6 → results/fig5_delivery_v7/ (signals.tsv.gz with lead, lead_logp_marginal; other files copied), lead_changes_v6_to_v7.tsv |

## ld/

| file | what it does | inputs → outputs |
|---|---|---|
| `ld/ld_tool.R` | Defines the shared LD standard (PLINK 1.9 signed r on tissue donors) as an R library and command-line tool. | library/CLI: (tissue, snp_a, snp_b) table, v9 VCF, results/ld_v3/donors → r, r2, n_pair per pair |
| `ld/ld_v3.R` | Writes each tissue's freeze donor list (eQTL N, donors of Fig. 4B and Supplementary Table 4a, PLINK keep lists) and signal LD tables. | work/final/freeze/<T>.rds, fig5_delivery_v6 signals, v9 VCF → results/ld_v3/{donors/<T>.txt, signal_proxies.tsv.gz, signal_pairs.tsv.gz} |

## qtltools_comparison/

| file | what it does | inputs → outputs |
|---|---|---|
| `qtltools_comparison/make_perm_maps_freeze_extra.R` | Draws the perm1 and perm2 permutation maps (fixed seeds) for the freeze tissues. Drawn in the earlier run; the final run reused them for every permuted arm. | work/final/freeze/<T>.rds → perm_maps/final/freeze/<T>.perm{1,2}.perm_map.tsv, perm_maps/MANIFEST.tsv |
| `qtltools_comparison/make_freeze_beds.R` | Builds the QTLtools BEDs (raw, rank-normal first) and covariate files from BOHR's freeze objects and permutation maps: all donors, a missing value set to the promoter's mean. | work/final/freeze/<T>.rds, perm maps → QTLtools_linear_compare/freeze_all/inputs/<T>.<arm>.{raw,rint}.bed.gz, <T>.<arm>.cov.txt |
| `qtltools_comparison/fig2e_deliver_freeze_v3m_pairs.R` | Forms the comparison units (an exact mirror pair counts once) and compares BOHR with this first QTLtools run. Its bohr/ tables give each unit's promoters to the steps below. | final_nogate_v3m units, QTLtools prm files → results/fig2e_freeze_v3m/{counts.tsv, counts_display.tsv, bohr/, qtltools/} |
| `qtltools_comparison/fig2e_r5e80_units.R` | Keeps the units that pass the final promoter rules (gene median reads >= 5; neither share of 0 nor of 1 >= 80%) and tabulates the stored calls of the earlier and the final version. Later steps read only its unit list. | results/fig2e_freeze_v3m of both versions, freeze objects, promoter_depth_v4.tsv.gz (`../units/compare_units.R`) → out/units_paper_r5e80.tsv.gz, summary tables |
| `qtltools_comparison/make_site_lists.R` | Lists, per tissue and chromosome, the variants that fail BOHR's genotype rule (at least five donors in each genotype class) among the tissue's donors, for QTLtools `--exclude-sites`. | freeze objects, v9 VCF (MAF >= 0.05) → sites/<T>.chr<N>.exclude.txt, sites/summary.tsv |
| `qtltools_comparison/exact_qtltools.R` | Runs QTLtools on BOHR's exact inputs per promoter: its analysis samples and the variants it tested; modes raw, rank-normal first and `--normal`. The next two runs copy its rows where their inputs are the same. | freeze objects, perm maps, unit list, fig2e_freeze_v3m/bohr, sites/, VCF → exact/v4/<mode>/<T>.<arm>.promoters.tsv.gz |
| `qtltools_comparison/keep_qtltools.R` | Runs QTLtools on every donor with a value and on BOHR's variants, for the promoters where BOHR removed rare boundary samples; other rows are copied from the exact run. Real arm only. | as exact_qtltools.R, exact/ → exact_keep/v4/<mode>/<T>.real.promoters.tsv.gz |
| `qtltools_comparison/pairs_qtltools.R` | Runs QTLtools on the promoter–variant pairs that BOHR tested, with the samples of the first run (all donors, mean fill); rows with identical inputs are copied from the exact and keep runs. | as exact_qtltools.R, exact/, exact_keep/ → pairs/v4/<mode>/<T>.<arm>.promoters.tsv.gz |
| `qtltools_comparison/fill_qtltools.R` | Gives the no-fill results: promoters with a missing value are run on the donors with a value (real arm: the keep run's rows), all other rows come from the pairs run. Also writes a zero-fill reference. | pairs/, exact_keep/, freeze objects, perm maps, VCF → pairs_nofill/v4/<mode>/<T>.<arm>.promoters.tsv.gz, pairs_zero/ |
| `qtltools_comparison/class_calib.R` | Builds the Fig. 2E units (`UNIT_MODE=merge_all`: a partial mirror pair is one unit; every computed component test counts) and the calls of BOHR and QTLtools per unit under several schemes. | final_nogate_v3m units (four arms), fig2e_freeze_v3m/bohr, unit list, pairs/, pairs_nofill/ → out/K_{units,map,summary,schemes,thresholds}_r5e80_merge_all |
| `qtltools_comparison/dd2_arms_units.R` | Rebuilds the same units for the real and the three permuted arms and counts the calls of BOHR and QTLtools, with BH inside each arm and at the real arm's thresholds. | final_nogate_v3m units, K_map, K_units, pairs/, pairs_nofill/ → dd2_units.tsv.gz, dd2_arms.tsv, dd2_checks.txt |
| `qtltools_comparison/dd11_fisher_rules.R` | Calls BOHR's units by rule BF: a component test passes BH over all component tests of the tissue, or both parts' smallest P < 0.05 and the Fisher combination of the Šidák-corrected parts passes BH over the units. Other rules are computed for comparison. | final_nogate_v3m units, unit list, dd2_units.tsv.gz → dd11_units_sidak.tsv.gz, dd11_arms_sidak.tsv, dd11_fisher_rules_sidak.txt |
| `qtltools_comparison/prepare_panel_e_classes.R` | Counts Fig. 2E (BOHR's calls by class against QTLtools raw and rank-normal first, no fill) and Fig. 2F (calls in the permuted arms, and the range expected by chance from simulated uniform P values). | dd11_units_sidak, dd11_arms_sidak, dd2_arms, dd2_units, K_summary → panel_e_classes_20261120_sidak/{classes.csv, permuted_by_arm.csv, null_band.csv, CHECKS.txt} |

The first three scripts of `qtltools_comparison/` keep the root of the earlier run; the final run used them on `rerun_v4`
(the mapping above). The scripts from `fig2e_r5e80_units.R` on keep the paths with which they ran; `v4` as their version
argument means the final run, `paper` the earlier one. Fig. 2E and 2F are drawn, and Supplementary Tables 2d and 2e
formatted, from these outputs by scripts that are not included.
