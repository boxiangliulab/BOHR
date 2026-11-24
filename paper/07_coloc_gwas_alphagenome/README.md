# Colocalization with GWAS traits, and variant effect prediction at MRAS

**Methods §10–11; Supplementary Methods §10–11** (§11 and §12 in the code comments, which predate the removal of a
section). Each independent puQTL signal is tested by coloc against 51 GWAS traits in the 50 tissues, over the union of the
puQTL LD window and the LD window of each clumped GWAS lead (P < 10⁻⁷). The GWAS loci are clumped with an
ancestry-matched LD reference: the 799 European donors for the 32 European or mainly European GWAS files, all 943 donors
for the 19 East Asian files. A colocalization is high-confidence at H4 ≥ 0.8 (2,181 of 90,886 tests, in 367 genes). At
MRAS, AlphaGenome predicts the effect of the candidate variants (72 in the final version) on CAGE and PRO-cap signal at
the four promoters (Fig. 6, Extended Data Fig. 6, Supplementary Tables 7–9).

**Run order.**
1. GWAS inputs (earlier run; the final run reads its outputs): `gwas_inputs/prep_magic_eur.py` and `prep_ra_adipo_eur.R`
   (the other files were used as downloaded; `config/` lists the files and sample sizes) → `gwas_coloc/10_gwas_sig_snps.R`
   per GWAS file. These steps do not depend on the puQTL units. `gwas_coloc/ld_v3_eur_mras.R` (earlier run) writes the
   European donor lists (self-reported RACE == 3; `ALL.EUR.txt`, 799 of the 943 donors); its second part checks the LD of
   the MRAS variant rs11713141 (v784) with the hypertension GWAS lead in every tissue.
2. puQTL side (stage S6.1–S6.2 of the r5e80 run): `gwas_coloc/01_build_analysis_set.R` → `12a_puqtl_windows.R` (windows of
   the earlier LD-window table: 5,134 of the 30,168 analysis rows) → `12b_extra_puqtl_windows_nscc.sh` (NSCC: PLINK
   windows of the other 22,137 tissue × lead pairs). The r5e80 run then ran steps 12–20 with the 943-donor clumps for every
   GWAS file; the run of item 5 replaced it.
3. MRAS candidates of the final version: `alphagenome/20_candidates.R` (variants with r² ≥ 0.8, over each tissue's donors,
   to the lead of any significant MRAS signal: 72, of which 71 were among the 170 variants already predicted; it also
   rebuilds the 161 candidates of the earlier version) → `alphagenome_perbase_MRAS_retry.py` for the one new candidate
   (`config/r5e80_new_candidate_to_score.txt`) → `score_r5e80_candidates.py` (all 171 variants on the fixed CAGE and
   PRO-cap peaks; the 170 earlier variants reproduce the earlier tables exactly).
4. GWAS clumping (step 11) on the lab cluster: `gwas_coloc/11_clump_local.R` for every GWAS file, once with all 943 donors
   (`clump_ALL.pbs`) and once with the 799 European donors (`clump_EUR.pbs`). With the 943 donors it reproduces the NSCC
   clumping of the earlier run byte for byte (`leads.tsv` and `lead_partners.tsv` of all 51 GWAS files). `make_eurmix.R` then takes the European clumping for
   the 32 files whose ancestry group in Supplementary Table 7a is European (31) or Multi-ancestry (1), and the 943-donor
   clumping for the 19 East Asian files. 96 leads of the European files are monomorphic in the 799 donors; they are kept
   as leads without partners, so their window is the lead ± 50 kb.
5. Steps 12–20 with the ancestry-matched clumps (`gwas_coloc/coloc_eur.pbs`): `12_pairs.R` → `13_subset_gwas.R` per GWAS
   file → `14_coloc_tissue.R` per tissue → `15_merge.R` → `16_make_canonical.R` (the paper's table, GWAS-lead tier 10⁻⁷)
   → `20_mras_panel_ldwin.R` (MRAS tables of Fig. 6B, C, E). The coloc functions are those of
   `../05_coloc_eqtl/v3c/coloc_final_lib.R`.

The AlphaGenome predictions of the earlier version are inputs of item 3: `alphagenome_perbase_MRAS.py` (95 variants) and
`alphagenome_perbase_MRAS_retry.py` (75 added candidates) → `mras_peak_finalize.py` (CAGE peaks) →
`mras_score_variants.py` (PRO-cap peaks) → `score_v3_candidates.py` (scores of the 170 variants). The AlphaGenome scripts
read the API key from the environment variable `ALPHAGENOME_API_KEY`; no key is stored.

**Changes for the final run.** The scripts of steps 01, 10, 12a, 12b and 12–20 are unchanged. Step 11 ran as
`11_clump_local.R`: the NSCC script of the earlier run with its three paths (PLINK, reference folder, donor list) read
from the environment variables `CLUMP_PLINK`, `CLUMP_REF` and `CLUMP_DONORS`, whose defaults are the NSCC paths; the rest
of the code is identical. The earlier script and its NSCC array job are not included.

**The coloc library of steps 12–20.** Steps 01 and 12a sourced the run copy of `../05_coloc_eqtl/v3c/coloc_final_lib.R`
that the eQTL chain ran (`nofloor_full/run_code/lab_r5e80/coloc/coloc_final_lib.R`). Steps 12–20 of the final run
sourced a separate copy, `nofloor_full/consult/freeze_ld_20261124/gwas/run_code_eur/coloc/coloc_final_lib.R`: the same code
with the roots of the final GWAS run (`NDIR`, the sample-size and trait lists, the GWAS folder) and one more line after
the `UNITSET` check, `if (UNITSET != "v3c") stop(...)`. That line only stops a run with another unit set; it changes no
result. `../05_coloc_eqtl/v3c/coloc_final_lib.R` does not have it.

**Not included.** The job that ran the AlphaGenome prediction of the new candidate read the API key from a local file,
so it is not included; `run.sh` gives its command (1-Mb interval, CAGE, PRO-cap, DNase and ATAC). It ran a copy of
`alphagenome_perbase_MRAS_retry.py` that differs from it only in its docstring and in the comment of the key line. The
scripts that build Supplementary Tables 7–9 and the figures are not included.

**Paths.** The scripts of steps 01, 10, 12a, 12b and 12–20 keep the paths of the earlier runs. In the final run they ran
as copies in which only folder names were changed:

| path in the scripts | folder in the final run |
|---|---|
| `/projects/group/bliu/e1327989/GTEx`, `/home/e1327989/project/GTEx` | `/projects/group/bliu/e1327989/GTEx/nofloor_full/GTEx` |
| `…/GTEx/QTL_figure/rerun_v2` | `…/nofloor_full/GTEx/QTL_figure/rerun_r5e80` |
| `…/GTEx/new_GWAS/coloc_final_nogate_v3/scripts/coloc_final_lib.R` (sourced) | steps 01 and 12a: `…/nofloor_full/run_code/lab_r5e80/coloc/coloc_final_lib.R`; steps 12–20: `…/nofloor_full/consult/freeze_ld_20261124/gwas/run_code_eur/coloc/coloc_final_lib.R` (above) |
| `NDIR` of `coloc_final_lib.R`, `…/GTEx/new_GWAS/coloc_final_nogate_v3`: steps 01, 12a, 12b | `…/nofloor_full/GTEx/new_GWAS/coloc_final_nogate_v3_r5e80` |
| the same, steps 12–20 | `…/nofloor_full/consult/freeze_ld_20261124/gwas/coloc_final_nogate_v3_r5e80_eur` |
| NSCC `/scratch/users/nus/e1327989/fig6_ldwin` (12b) | `/scratch/users/nus/e1327989/nofloor_full/fig6_ldwin_r5e80` |

The root of steps 12–20 links to the r5e80 root for the puQTL side (analysis set, `puqtl_windows.tsv`,
`puqtl_extra_windows.tsv`) and to the earlier run for the GWAS files, `gwas_sig` and the GWAS subsets; its
`coloc_ldunion_20261012/nscc_out/clump` is a link to `clump_EURmix`. The library reads `UNITSET` (`v3c`). The scripts
added for the final run (`11_clump_local.R`, the three lab jobs, `make_eurmix.R`, `20_candidates.R`,
`score_r5e80_candidates.py`) carry its paths; `ld_v3_eur_mras.R` carries those of the earlier run.

**Inputs without a script here.** The LD reference `LD_ref_v9_plink19` (PLINK 1.9 bfiles of the GTEx v9 genotypes, built
on NSCC and copied to the lab cluster, md5-checked) holds the donor lists in the two-column form that PLINK reads:
`donors/ALL.txt` has the 943 donors of `ld_v3.R` (`../03_cis_mapping_conditional/ld/`) and `donors_EUR/ALL.EUR.txt` the 799
of `ld_v3_eur_mras.R` (checked: the same donor sets); the step that wrote them was not kept. `make_eurmix.R` reads the
ancestry group of each GWAS file from Supplementary Table 7a. `20_candidates.R` reads the earlier Supplementary Tables 9a
and 9d and an rsID map (dbSNP 156); in its labels and output columns, "paper" means the earlier version.

**Counts written into the scripts.** `gwas_coloc/15_merge.R` asserts 50 per-tissue files (line 17).
`score_r5e80_candidates.py` asserts 171 predicted variants, 72 candidates and 273 CAGE and 6 PRO-cap sense tracks;
`20_candidates.R` asserts that the rebuilt earlier set equals the earlier Supplementary Table 9a; `make_eurmix.R` that
its 51 files are those of Supplementary Table 7a. `coloc_eur.pbs` stops unless the 943-donor clumps equal the NSCC
outputs of the earlier run for all 51 files and the 799-donor clumping has finished for all 51.

`run.sh` lists the commands in the order they were run.


## gwas_inputs/ and config/

| file | what it does | inputs → outputs |
|---|---|---|
| `gwas_inputs/prep_magic_eur.py` | Converts the MAGIC 2021 European fasting glucose, fasting insulin and 2-h glucose files to the pipeline format. | GWAS Catalog harmonised MAGIC files (raw/) → gwas_eur_20261009/prepared/*_hg38.txt |
| `gwas_inputs/prep_ra_adipo_eur.R` | Prepares the Okada RA European and ADIPOGen adiponectin files: GRCh38 positions, log-OR and SE. | RA_GWASmeta_European_v2.txt.gz, ADIPOGen discovery file, harmonised MAGIC FG file, UCSC chain → prepared/ RA and adiponectin files |
| `config/config_gwas_sample_sizes_v3eur.txt` | Gives the GWAS sample size supplied to coloc for each of the 51 files. | read by coloc_final_lib.R gwas_N() |
| `config/traits_51_basenames.txt` | Lists the 51 GWAS file base names. | the clumping jobs and make_eurmix.R read an identical copy (traits.txt) |
| `config/tissues_50.txt` | Lists the 50 tissues of the colocalization, one per line (the freeze tissues of work/final/tissues.tsv). | read by ../03_cis_mapping_conditional/run.sh and ../05_coloc_eqtl/run.sh |
| `config/trait_maps.R` | Maps GWAS file names to plotted trait labels, labels to eight categories and categories to colours. | sourced by 15_merge.R, 16_make_canonical.R and the figure scripts, by its stored path (`figure6/scripts/trait_maps.R`; this copy is byte-identical) |

## gwas_coloc/

| file | what it does | inputs → outputs |
|---|---|---|
| `gwas_coloc/01_build_analysis_set.R` | Builds the coloc analysis set: significant units x independent signals, plus the non-significant MRAS units. | units v3c and conditional signals v3c → results_v3c/analysis_set.rds, analysis_windows.rds |
| `gwas_coloc/10_gwas_sig_snps.R` | Selects one GWAS file's SNPs with P < 1e-5 that are typed in the GTEx v9 genotypes, as clumping input. | prepared GWAS file, v9 bim files (full MAF >= 0.005 set and panel) → gwas_sig/<trait>.tsv.gz, <trait>.counts.tsv |
| `gwas_coloc/ld_v3_eur_mras.R` | Writes the European donor lists (RACE == 3) of every tissue and of all donors; checks the LD of v784 with the hypertension GWAS lead. | ld_v3 donor lists, metadata/v10_cov_pca_full_GRU.txt, v9 VCF → results/ld_v3/donors/*.EUR.txt, mras_v784_gwas_lead.tsv |
| `gwas_coloc/11_clump_local.R` | Clumps one GWAS file's significant SNPs recursively on PLINK r2 and finds each lead's LD partners, with the donors of CLUMP_DONORS. | gwas_sig/<trait>.tsv.gz, v9 PLINK bfiles (MAF >= 0.005 full set, MAF >= 0.05 panel), donor list → <outdir>/<trait>/{leads,lead_partners}.tsv, DONE |
| `gwas_coloc/clump_ALL.pbs` | Lab job: 11_clump_local.R for the 51 GWAS files with all 943 donors, two at a time. | traits.txt, gwas_sig/ → clump_ALL/ |
| `gwas_coloc/clump_EUR.pbs` | Lab job: the same with the 799 European donors. | traits.txt, gwas_sig/ → clump_EUR/ |
| `gwas_coloc/make_eurmix.R` | Assembles the ancestry-matched clumps: European clumping for European and mainly European files, 943-donor clumping for East Asian files; leads without an r2 row get no partners. | clump_ALL/, clump_EUR/, Supplementary Table 7a (ancestry_group) → clump_EURmix/<trait>/, clump_EURmix_sources.tsv |
| `gwas_coloc/12a_puqtl_windows.R` | Attaches the earlier LD-window table's window to every puQTL analysis signal and lists the leads without one. | QTL_figure coloc_eqtl_ldwin/lead_windows.tsv and coloc_eqtl_v3c analysis set, results_v3c/analysis_set.rds → run_v3c/puqtl_windows.tsv, puqtl_leads_missing.tsv |
| `gwas_coloc/12b_extra_puqtl_windows_nscc.sh` | Computes LD windows (PLINK 1.9 r2 >= 0.1 within 1 Mb, tissue donors) for the puQTL leads without one. | puqtl_leads_missing.tsv, v9 panel bfiles, tissue donor lists → extra_win/puqtl_extra_windows.tsv (read as nscc_out/puqtl_extra_windows.tsv) |
| `gwas_coloc/coloc_eur.pbs` | Lab job of steps 12–20 with the ancestry-matched clumps; checks the clumps first. | clump_ALL/, clump_EUR/, clump_EURmix/ → the outputs of 12–20 in coloc_final_nogate_v3_r5e80_eur |
| `gwas_coloc/12_pairs.R` | Forms candidate puQTL signal x GWAS locus pairs whose LD windows overlap; records intersection and union. | puqtl_windows.tsv, puqtl_extra_windows.tsv, nscc_out/clump/*/{leads,lead_partners}.tsv → run_v3c/pairs/, gwas_loci.tsv.gz, union_intervals/ |
| `gwas_coloc/13_subset_gwas.R` | Cuts one GWAS file to the union regions of its candidate pairs. | prepared GWAS file, run_v3c/union_intervals/ → run_v3c/gwas_subset/<trait>.rds |
| `gwas_coloc/14_coloc_tissue.R` | Runs puQTL x GWAS coloc for one tissue on each candidate pair's union window (marginal or LOO scans). | analysis set, pairs, GWAS subsets, per-SNP statistics and LOO scans of cis mapping → run_v3c/per_tissue/coloc_ldwin_<tissue>.rds |
| `gwas_coloc/15_merge.R` | Merges the per-tissue coloc pairs and keeps the best tested pair per signal x trait at each GWAS-lead tier. | run_v3c/per_tissue/coloc_ldwin_<tissue>.rds → run_v3c/results/coloc_ldwin_pairs.rds, coloc_ldwin_signal_trait_<tier>.rds, summary_totals_by_tier.tsv |
| `gwas_coloc/16_make_canonical.R` | Builds the canonical coloc master and its summary totals at one GWAS-lead tier. | run_v3c/results/coloc_ldwin_signal_trait_<tier>.rds → results_v3c/coloc_ldwin_merged_1e7_v3c_20261106.rds, run_v3c/results/canonical_1e7/summary_*.tsv, gene_rank_* |
| `gwas_coloc/20_mras_panel_ldwin.R` | Builds the MRAS x cardiovascular GWAS coloc tables over all tissues and the rs11713141 direction per MRAS unit. | results_v3c master, analysis set, pairs, GWAS subsets, per-SNP statistics and LOO scans → results_v3c/mras/*.tsv (inputs of Fig. 6B, C, E) |

## alphagenome/

| file | what it does | inputs → outputs |
|---|---|---|
| `alphagenome/alphagenome_perbase_MRAS.py` | Calls the AlphaGenome API per variant (1-Mb interval) and saves per-base REF/ALT tracks over the four MRAS promoters. | variant list → per_snp/<variant>.npz (REF/ALT per modality: CAGE, PRO-cap, DNase, ATAC), meta/<MOD>_tracks.parquet |
| `alphagenome/alphagenome_perbase_MRAS_retry.py` | Copy of the prediction script with retry on quota and server errors and per-variant memory release (added variants). | added-variant lists → alphagenome_results/perbase_MRAS_v3add/ and perbase_MRAS_r5e80add/ per_snp/*.npz |
| `alphagenome/perbase_MRAS.pbs` | Job of the first AlphaGenome run (95 variants; 1-Mb interval; CAGE, PROCAP, DNASE, ATAC). The key comes from the environment. | MRAS_allpromoter_snplist.txt → perbase_MRAS/ |
| `alphagenome/perbase_MRAS_v3add.pbs` | Job of the first add-on attempt (74 variants). | MRAS_v3_new_candidates_to_score.txt → perbase_MRAS_v3add/ |
| `alphagenome/perbase_MRAS_v3add_retry.pbs` | Job of the add-on with the retry copy (75 variants; resumes, skipping variants already predicted); parts of the list ran as copies of this job. | MRAS_v3_new_candidates_to_score_ALL75.txt → perbase_MRAS_v3add/ |
| `alphagenome/mras_peak_finalize.py` | Calls the common predicted CAGE peak set on the adrenal sense-strand track averaged over the 95 REF predictions. | perbase_MRAS/per_snp/*.npz, meta/CAGE_tracks.parquet, exon table → tables/predicted_cage_peaks_COMMON.tsv, fig_peak_common.{pdf,png} |
| `alphagenome/mras_score_variants.py` | Scores the first 95 variants on the fixed CAGE and PRO-cap peaks and calls the common PRO-cap peak set. | perbase_MRAS/per_snp/*.npz, track metadata, predicted_cage_peaks_COMMON.tsv, old proxy tables → tables/predicted_procap_peaks_COMMON.tsv, score_*_peaks.tsv |
| `alphagenome/score_v3_candidates.py` | Scores the 170 variants of the earlier version on the fixed CAGE and PRO-cap peaks (log2FC, P45112 share). | per_snp/*.npz of the 95- and 75-variant runs, predicted_*_peaks_COMMON.tsv, 161-candidate list → alphagenome_addon/tables/score_*_v3.tsv, ranking_v3_candidates.tsv |
| `alphagenome/20_candidates.R` | Builds the MRAS candidate set of the final version (r2 >= 0.8 to the lead of any significant MRAS signal) and compares it with the earlier set and the predicted variants. | r5e80 analysis set, r5e80 and earlier ld_v3/signal_proxies.tsv.gz, earlier Supplementary Tables 9a and 9d, rsID map → checks/candidates_r5e80_vs_paper.tsv, mras_signals_r5e80.tsv |
| `alphagenome/score_r5e80_candidates.py` | Scores the 171 predicted variants on the fixed CAGE and PRO-cap peaks and flags the 72 candidates; checks that the 170 earlier variants reproduce the earlier tables. | per_snp/*.npz of the three runs, predicted_*_peaks_COMMON.tsv, candidates_r5e80_vs_paper.tsv, 161-candidate list, score_*_v3.tsv → alphagenome/tables/score_*_r5e80.tsv, ranking_r5e80_candidates.tsv |
| `alphagenome/config/MRAS_allpromoter_snplist.txt` | Lists the 95 MRAS-locus variants of the first AlphaGenome run. | read by perbase_MRAS.pbs → alphagenome_perbase_MRAS.py |
| `alphagenome/config/MRAS_v3_new_candidates_to_score.txt` | Lists the added MRAS candidate variants (74) of the first add-on job. | read by perbase_MRAS_v3add.pbs (first attempt) |
| `alphagenome/config/MRAS_v3_new_candidates_to_score_ALL75.txt` | Lists the added MRAS candidate variants (75) predicted by the add-on job. | read by perbase_MRAS_v3add_retry.pbs and score_v3_candidates.py |
| `alphagenome/config/MRAS_v3_candidate_set_LDstandard_PLINK.tsv` | Lists the 161 MRAS candidate variants of the earlier version (r2 >= 0.8 proxies of its 105 significant MRAS signals) with their signals and maximum r2. | read by score_v3_candidates.py and score_r5e80_candidates.py (in_v3_candidate_set) |
| `alphagenome/config/r5e80_new_candidate_to_score.txt` | Lists the one candidate of the final version that had no AlphaGenome prediction (chr3_138437434_C_T_b38). | read by the prediction job of the final version (run.sh) → perbase_MRAS_r5e80add/ |
