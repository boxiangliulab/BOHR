# cis mapping, comparison with QTLtools and independent signals

**Methods §4, §5 and §7; Supplementary Methods §4, §5 and §7.** For each of 50 tissues: inputs (one sample per donor,
five genotype PCs, PEER factors capped by sample size), a per-SNP BOHR scan of every promoter within 500 kb (μ, p0, p1),
units (promoter × component) with window P values combined by ACAT, estimability and mirror rules, BH per tissue, then
forward conditional selection with a backward step, giving the 14,800 independent signals.

**Run order.**
1. `inputs/`: `peer/peer_analysis.R` (PEER factors; input built as recorded in `peer/0104_peer.R`),
   `peer/gamlss_QTL_peer_test.R` and `peer/peer_plot.R` (the chromosome 22 scan that set the PEER count per tissue,
   `config/tissue_list_best_peer_table.txt`), `gamlss_QTL_optimized.R` (lines 115–137 set the tested promoters), then
   `prep_all.R` (per-tissue input objects and permutation maps).
2. `scan/run_final2.R` per tissue, arm and chromosome; `scan/run_firth_addon.R` and `scan/run_muF.R` add the remaining
   Firth tests and the F-corrected μ P values.
3. `units/`: `mirror_pairs.R`, `discovery_nogate.R` (units, ACAT, BH), `units_v2_estimable.R`, `units_v3_estimable.R`
   (estimability rules), `units_v3m_mirror.R` (mirror rule), `mirror_abs_labels.R` and `units_v3b_from_labels.R`
   (labels by absolute activity; P and q unchanged). The paper's units are the last output (units v3c).
4. `conditional/`: the forward runs were made incrementally as the unit rules were settled, and the final table takes
   each unit from its latest run: `run_conditional_ng.R` (all units v1), `run_conditional_ng2.R` and
   `run_conditional_ng3.R` (units new or continued under v2 and v3), `merge_conditional_ng3.R` and `loo_est_patch.R`;
   then `cond_v3m_select.R` and `run_conditional_v3m.R` (units new or continued under the mirror rule),
   `merge_conditional_v3m.R` (backward step), and `run_conditional_v3c.R` with `merge_conditional_v3c.R` (two units refitted
   on 2026-11-06; 14,800 signals). `fig5_delivery.R` and `fig5_delivery_v4.R` write the per-signal tables read downstream;
   `fig5_proxies.R`, `fig5_pertype.R`, `fig5_pe_moments_v6.R`, `crosstype_conditional_v6.R` and `fig5_signals_partner.R`
   add the per-signal tables that the major/minor analysis reads (`../06_major_minor_promoters/`).
5. `ld/`: `ld_tool.R` (the LD standard: PLINK 1.9 r on each tissue's donors) and `ld_v3.R` (donor lists, signal LD).
6. `qtltools_comparison/`: `make_perm_maps_freeze_extra.R`, `make_freeze_beds.R`, QTLtools v1.2 `cis --permute 1000`
   (raw, `--normal`, RINT-first; command in `run.sh`), `fig2e_deliver_freeze_v3m_pairs.R` (Fig. 2E counts).

**BOHR package.** `BOHR::bohr_cis()` reproduces steps 2–3 and `BOHR::bohr_conditional()` step 4 in one call each: on the
stored inputs the per-SNP statistics, unit P values, q values, signal counts and lead variants agree with this pipeline
(repository `workflow/parity/`, scripts 01–06). The subset test (`../test_subset/`) re-runs `run_final2.R` and
`run_conditional_v3m.R` as they were run.

**Paths.** `../02_bohr_model_and_simulation/engine/bohr2_data.R` sets the input paths (`BOHR2_GTEX`, `BOHR2_ROOT`,
`BOHR2_VCF`); scripts source the engine from `BOHR2_CODE`, and write under `BOHR2_ROOT/results/` unless `FINAL_ROOT` or
`COND_ROOT` is set.

`run.sh` lists the commands in the order they were run (paths of the original runs).


## inputs/

| file | what it does | inputs → outputs |
|---|---|---|
| `inputs/peer/0104_peer.R` | Records how the PEER input was built (commented-out lines 48-57: donor means within each tissue, missing set to 0). | filter_5_recal_att.txt (text export of 0103's recalculated activity) → filter_5_recal_att_filter_SMTSD_mean.txt (lines 48-57, commented out); PEER part → v10_* files (not used) |
| `inputs/peer/peer_analysis.R` | Estimates 60 PEER factors for one tissue from the donor-mean relative-activity matrix. | tissue, filter_5_recal_att_filter_SMTSD_mean.txt, output dir → v11_filter_5_peer_factors_<tissue>.txt |
| `inputs/peer/gamlss_QTL_peer_test.R` | Runs the earlier gamlss LRT scan for one tissue on chr22 (+-10 kb) with k PEER factors and 5 genotype PCs. | v11 matrix, Promoter_info_filter_5.txt, genes > 1 TPM, genotype PCs, PEER factors, chr22 VCF → <name>.<tissue>.chr22.factor_peer_<k>_pca_5.results.rds |
| `inputs/peer/peer_plot.R` | Counts, per tissue and k, the promoters significant at q <= 0.05 after ACAT over their chr22 variants. | ../tissue_list_size20_table.txt; <tissue>/sensitivity_test_10k.*.rds → <tissue>_peer_test_count.txt, <tissue>_peer_test_count.pdf |
| `inputs/gamlss_QTL_optimized.R` | Runs the published GAMLSS puQTL scan; its promoter filter (SD, NA share, TPM) sets the tested promoter set. | v11 activity, promoter info, PCs, PEER, TPM1_genes.txt, VCF → QTL/new_v11_LRT/<T>/chr<N>/*.results.rds |
| `inputs/prep_all.R` | Builds per-tissue BOHR inputs (freeze samples, donors, promoter set, PEER cap K) and perm0 permutation maps. | v11 activity, sample attributes, PCs, PEER, best-PEER table, published promoter lists → work/final/<def>/<T>.rds, work/final/tissues.tsv, perm0 maps |
| `inputs/config/tissue_list_size20_table.txt` | Lists the 57 GTEx v11 tissues with at least 20 samples and their sample counts. | hand-made → read by peer_plot.R line 9 |
| `inputs/config/tissue_list_size20_table_add2.txt` | Lists the four small tissues of the last PEER-scan submission and their sample counts. | hand-made → read by submit_peer_jobs.sh line 4 |
| `inputs/config/tissue_list_best_peer_table.txt` | Gives the best PEER factor count per tissue (chr22 scan); prep_all.R reads it as k_now. | read by prep_all.R (GTEx/QTL/tissue_list_best_peer_table.txt) |

## scan/

| file | what it does | inputs → outputs |
|---|---|---|
| `scan/run_final2.R` | Runs the main per-SNP BOHR scan (mu, p0, p1; Firth p0/p1 when f0 or f1 >= 0.35) per tissue, arm and chr. | work/final/<def>/<T>.rds, perm_maps/final/<def>/<T>.<arm>.perm_map.tsv, v9 VCF → results/final/<def>/<T>/<arm>/chr<N>.rds |
| `scan/run_firth_addon.R` | Adds the Firth p0/p1 tests for every other promoter with a 0/1 part, on the SNPs of the main scan. | results/final/<def>/<T>/<arm>/chr*.rds, work/final/<def>/<T>.rds, perm maps → results/final_firth_addon/<def>/<T>/<arm>/chr<N>.rds |
| `scan/run_muF.R` | Adds the small-sample F-corrected mu P value per SNP from the stored beta-part LRT. | results/final/<def>/<T>/<arm>/chr*.rds, work/final/<def>/<T>.rds → results/final_muF/<def>/<T>/<arm>/chr<N>.rds |

## units/

| file | what it does | inputs → outputs |
|---|---|---|
| `units/mirror_pairs.R` | Builds the exact mirror-pair table (per-component exactness) and promoter flags for the 50 tissues. | work/final/freeze/*.rds → results/mirror/{mirror_pairs.tsv.gz, promoter_flags.tsv.gz} |
| `units/discovery_nogate.R` | Forms units (promoter x mu/p0/p1), combines window P values by ACAT, drops mirror twins, runs BH. | results/{final,final_muF,final_firth_addon}/<def>/<T>/<arm>/, results/mirror/promoter_flags.tsv.gz → results/final_nogate/<def>/<T>.<arm>.units.rds |
| `units/units_v2_estimable.R` | Applies the mu n_mid >= 80 rule to units v1 and recomputes BH, giving units v2 for the second conditional run. | results/final_nogate/<def>/*.units.rds → results/final_nogate_v2/<def>/<T>.<arm>.units.rds |
| `units/units_v3_estimable.R` | Applies the estimability rules (mu n_mid >= 80; p0/p1 share 10-90%) before BH and recomputes q per tissue. | results/final_nogate/<def>/*.units.rds → results/final_nogate_v3/<def>/*.units.rds, freeze_thresholds_v2_v3.tsv |
| `units/units_v3m_mirror.R` | Applies mirror rule v3 (exact boundary events held as p0), reruns BH per tissue and writes checks. | results/final_nogate_v3/freeze/*.units.rds, results/mirror/promoter_flags.tsv.gz → results/final_nogate_v3m/freeze/, results/mirror_v3/ |
| `units/mirror_abs_labels.R` | Labels each significant exact-pair boundary event by the TSS whose absolute activity changes (0 relabels after the refit). | final_nogate_v3m, conditional_final_v3c signals (env COND_FINAL), v11 tabs_att_filter_In.rds → results/mirror_abs_labels_v3c/labels_real.tsv.gz |
| `units/units_v3b_from_labels.R` | Adds the saved label columns to the v3m units without refitting or BH, giving units v3c (P and q identical to v3m). | final_nogate_v3m, labels_real.tsv.gz (env LABELS, UNITS_OUT) → results/final_nogate_v3c/freeze/<T>.<arm>.units.rds, checks.txt |

## conditional/

| file | what it does | inputs → outputs |
|---|---|---|
| `conditional/run_conditional_ng.R` | Runs forward conditional selection and leave-one-out scans for every significant unit of units v1, per tissue, arm and chromosome. | results/final_nogate (units v1), results/final, work/final, perm maps → results/conditional_nogate/<def>/<T>/<arm>/<chr>.rds |
| `conditional/run_conditional_ng2.R` | Runs the conditional supplement for units that are new or continue under units v2. | results/final_nogate_v2, results/conditional_nogate, results/final → results/conditional_nogate_supp/<def>/<T>/real/<chr>.rds |
| `conditional/run_conditional_ng3.R` | Runs the conditional supplement for units that are new or continue under units v3. | results/final_nogate_v3, conditional_nogate, conditional_nogate_supp, results/final → results/conditional_nogate_supp3/<def>/<T>/real/<chr>.rds |
| `conditional/loo_est_patch.R` | Refits the p0/p1 leave-one-out estimates that a naming bug left as NA in the conditional runs. | <COND_FINAL>/freeze_real_signals.tsv.gz, work/final, genotypes → <COND_FINAL>/loo_est_patch.tsv |
| `conditional/merge_conditional_ng3.R` | Merges the three forward runs of units v3 with one backward pass at the units-v3 tissue thresholds. | conditional_nogate, _supp, _supp3, final_nogate_v3, loo_est_patch.tsv → results/conditional_final_v3/ (env COND_FINAL) |
| `conditional/fig5_delivery.R` | Builds the 2026-10-02 per-signal tables of the units-v3 conditional analysis, whose signal list ld_v3.R reads. | conditional_final_v3, final_nogate_v3 (env COND_FINAL, UNITS_DIR, FIG5_OUT) → results/fig5_delivery_v3/{units,signals,...}.tsv.gz |
| `conditional/cond_v3m_select.R` | Selects the units to re-run under mirror rule v3 (new, continued and truncated units) and writes their job list. | conditional_final_v3, final_nogate_v3m, work/final/freeze → work/cond_v3m/{units_to_run.tsv, jobs.tsv} |
| `conditional/run_conditional_v3m.R` | Runs forward conditional selection and leave-one-out scans for the 2,660 units re-run under mirror rule v3. | work/cond_v3m/units_to_run.tsv, final_nogate_v3m, results/final, work/final → results/conditional_v3m/<def>/<T>/real/<chr>.rds |
| `conditional/merge_conditional_v3m.R` | Merges re-run and carried units, recomputes the backward step and writes the final signal tables. | results/conditional_v3m, results/conditional_final_v3, final_nogate_v3m → results/conditional_final_v3m/freeze_real_*.tsv.gz, loo/ |
| `conditional/units_to_run_v3c.tsv` | Lists the two units re-run in the 2026-11-06 conditional refit. | read by run_conditional_v3c.R, merge_conditional_v3c.R, c1_analysis_set_pairs_v3c.R and cond_v3c.pbs |
| `conditional/run_conditional_v3c.R` | Re-runs forward selection and leave-one-out scans for the two units refitted on 2026-11-06 (Putamen 47842 p0, Aorta 38390 p1). | v3c/units_to_run.tsv, final_nogate_v3m, results/final, work/final → results/conditional_v3c/<def>/<T>/real/<chr>.rds |
| `conditional/merge_conditional_v3c.R` | Replaces the two refitted units in the v3m conditional tables and keeps the v3m backward step, giving the 14,800 independent signals. | conditional_final_v3m, conditional_v3c, labels (env LABELS) → results/conditional_final_v3c/{freeze_real_units,signals,rounds}.tsv.gz, loo/ |
| `conditional/fig5_delivery_v4.R` | Builds the per-signal table (signals.tsv.gz) that the Figure 4 data builders and Supplementary Tables 4-5 read. | env COND_FINAL, UNITS_DIR, FIG5_OUT; work/final, results/final → <FIG5_OUT>/{units,signals,genotype_pheno,crosstype_ld,boundary_only_partner_pheno}.tsv.gz |
| `conditional/fig5_proxies.R` | Per-genotype phenotype at the LD proxies (r2 >= 0.8, PLINK rule) of every signal, for the per-signal tables read by the major/minor analysis. |  |
| `conditional/fig5_pertype.R` | Per-component marginal statistics (mu, p0, p1) at every signal SNP and LD proxy, looked up from the main-run SNP tables. |  |
| `conditional/fig5_pe_moments_v6.R` | Probability-scale effects (dosage 0 to 2, delta-method SEs; engine/pe_delta.R) and interior-only moments by genotype for every signal. |  |
| `conditional/crosstype_conditional_v6.R` | Conditional cross-type tests for pairs of mu and boundary signals of the same promoter (the pair table that 71_A1_location.R reads). |  |
| `conditional/fig5_signals_partner.R` | Adds SNP position, TSS, strand and mirror-pair columns per signal (signals_partner.tsv.gz) for Fig. 4C and Supplementary Tables 4-5. | <FIG5_OUT>/signals.tsv.gz, mirror_v3 promoter_flags (env FLAGS), v11 Promoter_info_filter_5.txt → <FIG5_OUT>/signals_partner.tsv.gz |

## ld/

| file | what it does | inputs → outputs |
|---|---|---|
| `ld/ld_tool.R` | Defines the shared LD standard (PLINK 1.9 signed r on tissue donors) as an R library and command-line tool. | library/CLI: (tissue, snp_a, snp_b) table, v9 VCF, results/ld_v3/donors → r, r2, n_pair per pair |
| `ld/ld_v3.R` | Writes each tissue's freeze donor list (eQTL N, donors of Fig. 4B and Supplementary Table 4a, PLINK keep lists) and signal LD tables. | work/final/freeze/<T>.rds, fig5_delivery_v3 signals, v9 VCF → results/ld_v3/{donors/<T>.txt, signal_proxies.tsv.gz, signal_pairs.tsv.gz} |

## qtltools_comparison/

| file | what it does | inputs → outputs |
|---|---|---|
| `qtltools_comparison/make_perm_maps_freeze_extra.R` | Draws the perm1 and perm2 permutation maps (fixed seeds) for the freeze tissues. | work/final/freeze/<T>.rds → perm_maps/final/freeze/<T>.perm{1,2}.perm_map.tsv, perm_maps/MANIFEST.tsv |
| `qtltools_comparison/make_freeze_beds.R` | Builds QTLtools BEDs (raw, RINT-first) and covariate files from BOHR freeze tissue objects and perm maps. | work/final/freeze/<T>.rds, perm maps → QTLtools_linear_compare/freeze_all/inputs/<T>.<arm>.{raw,rint}.bed.gz, <T>.<arm>.cov.txt |
| `qtltools_comparison/fig2e_deliver_freeze_v3m_pairs.R` | Compares BOHR and QTLtools calls (real and 3 permuted arms, 50 tissues) for Fig. 2E and Supp. Tables 2d-e. | final_nogate_v3m units, mirror_v3 flags, QTLtools prm files → results/fig2e_freeze_v3m/{counts.tsv, counts_display.tsv, bohr/, qtltools/} |
