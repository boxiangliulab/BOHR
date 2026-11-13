# Colocalization with GWAS traits, and variant effect prediction at MRAS

**Methods §10–11; Supplementary Methods §10–11** (§11 and §12 in the code comments, which predate the removal of a
section). Each independent puQTL signal is tested by coloc against 51 GWAS traits in the 50 tissues, over the union of the
puQTL LD window and the LD window of each clumped GWAS lead (P < 10⁻⁷); a colocalization is high-confidence at H4 ≥ 0.8
(1,125 in 193 genes). At MRAS, AlphaGenome predicts the effect of the candidate variants on CAGE and PRO-cap signal at the
four promoters (Fig. 6, Extended Data Fig. 6, Supplementary Tables 7–9).

**Run order.** GWAS inputs: `gwas_inputs/prep_magic_eur.py` and `prep_ra_adipo_eur.R` (the other files were used as
downloaded; `config/` lists the files and sample sizes). `gwas_coloc/`: `01_build_analysis_set.R` → `10_gwas_sig_snps.R` →
`11_clump_nscc.R` (per trait) → `12a_puqtl_windows.R` and `12b_extra_puqtl_windows_nscc.sh` → `12_pairs.R` →
`13_subset_gwas.R` → `14_coloc_tissue.R` (per tissue) → `15_merge.R` → `16_make_canonical.R` (the paper's table, GWAS-lead
tier 10⁻⁷) → `20_mras_panel_ldwin.R` (MRAS tables of Fig. 6B, C, E). The coloc functions are those of
`../05_coloc_eqtl/v3c/coloc_final_lib.R`. `alphagenome/`: `alphagenome_perbase_MRAS.py` (95 variants) and
`alphagenome_perbase_MRAS_retry.py` (the added candidates) → `mras_peak_finalize.py` (CAGE peaks) →
`mras_score_variants.py` (PRO-cap peaks) → `score_v3_candidates.py` (scores of the 170 predicted variants). The
AlphaGenome scripts read the API key from the environment variable `ALPHAGENOME_API_KEY`; no key is stored.

`run.sh` lists the commands in the order they were run (paths of the original runs).


## gwas_inputs/ and config/

| file | what it does | inputs → outputs |
|---|---|---|
| `gwas_inputs/prep_magic_eur.py` | Converts the MAGIC 2021 European fasting glucose, fasting insulin and 2-h glucose files to the pipeline format. | GWAS Catalog harmonised MAGIC files (raw/) → gwas_eur_20261009/prepared/*_hg38.txt |
| `gwas_inputs/prep_ra_adipo_eur.R` | Prepares the Okada RA European and ADIPOGen adiponectin files: GRCh38 positions, log-OR and SE. | RA_GWASmeta_European_v2.txt.gz, ADIPOGen discovery file, harmonised MAGIC FG file, UCSC chain → prepared/ RA and adiponectin files |
| `config/config_gwas_sample_sizes_v3eur.txt` | Gives the GWAS sample size supplied to coloc for each of the 51 files. | read by coloc_final_lib.R gwas_N() |
| `config/traits_51_basenames.txt` | Lists the 51 GWAS file base names indexed by the NSCC clumping array job. | read by 11_clump_nscc.pbs (array 1-51) |
| `config/tissues_50.txt` | Lists the 50 tissues of the GWAS colocalization, one per line (array index of the per-tissue coloc job). | read by v3m_14_coloc_tissue.pbs (and v3c_14_coloc_tissue.pbs) |
| `config/trait_maps.R` | Maps GWAS file names to plotted trait labels, labels to eight categories and categories to colours. | sourced by 15_merge.R, 16_make_canonical.R, 26_text_numbers.R, build_S7_to_S10.R and the figure scripts |

## gwas_coloc/

| file | what it does | inputs → outputs |
|---|---|---|
| `gwas_coloc/01_build_analysis_set.R` | Builds the coloc analysis set: significant units x independent signals, plus the non-significant MRAS units. | Figure 4 v3m units and conditional signals → results_v3m/analysis_set.rds, analysis_windows.rds |
| `gwas_coloc/10_gwas_sig_snps.R` | Selects one GWAS file's SNPs with P < 1e-5 that are typed in the GTEx v9 genotypes, as clumping input. | prepared GWAS file, v9 bim files (full MAF >= 0.005 set and panel) → gwas_sig/<trait>.tsv.gz, <trait>.counts.tsv |
| `gwas_coloc/11_clump_nscc.R` | Clumps one GWAS file's significant SNPs recursively on PLINK r2 and finds each lead's LD partners. | gwas_sig/<trait>.tsv.gz, v9 PLINK bfiles (MAF >= 0.005 full set, MAF >= 0.05 panel), all-donor list → nscc_out/clump/<trait>/{leads,lead_partners}.tsv |
| `gwas_coloc/12a_puqtl_windows.R` | Attaches Figure 4's LD window to every puQTL analysis signal and lists the leads without one. | Figure 4 coloc_eqtl_ldwin/lead_windows.tsv and analysis set, results_v3m/analysis_set.rds → run_v3m/puqtl_windows.tsv, puqtl_leads_missing.tsv |
| `gwas_coloc/12b_extra_puqtl_windows_nscc.sh` | Computes LD windows (PLINK 1.9 r2 >= 0.1 within 1 Mb) for the 122 MRAS puQTL leads missing from Figure 4's table. | puqtl_leads_missing.tsv, v9 panel bfiles, tissue donor lists → extra_win/puqtl_extra_windows.tsv (read as nscc_out/puqtl_extra_windows.tsv) |
| `gwas_coloc/12_pairs.R` | Forms candidate puQTL signal x GWAS locus pairs whose LD windows overlap; records intersection and union. | puqtl_windows.tsv, extra windows (12b; Figure 4 v3m NSCC run), nscc_out/clump/*/{leads,lead_partners}.tsv → run_v3m/pairs/, gwas_loci.tsv.gz, union_intervals/ |
| `gwas_coloc/13_subset_gwas.R` | Cuts one GWAS file to the union regions of its candidate pairs (step 6 of the coloc chain). | prepared GWAS file, run_<set>/union_intervals/ → run_<set>/gwas_subset/<trait>.rds |
| `gwas_coloc/14_coloc_tissue.R` | Runs puQTL x GWAS coloc for one tissue on each candidate pair's union window (marginal or LOO scans). | analysis set, pairs, GWAS subsets, Figure 4 per-SNP statistics and LOO scans → run_v3m/per_tissue/coloc_ldwin_<tissue>.rds |
| `gwas_coloc/15_merge.R` | Merges the per-tissue coloc pairs and keeps the best tested pair per signal x trait at each GWAS-lead tier. | run_v3m/per_tissue/coloc_ldwin_<tissue>.rds → run_v3m/results/coloc_ldwin_pairs.rds, coloc_ldwin_signal_trait_<tier>.rds, summary_totals_by_tier.tsv |
| `gwas_coloc/16_make_canonical.R` | Builds the canonical coloc master and its summary totals at one GWAS-lead tier; live copy that also names the v3c output file. | run_<set>/results/coloc_ldwin_signal_trait_<tier>.rds → results_v3c/coloc_ldwin_merged_1e7_v3c_20261106.rds, run_v3c/results/canonical_1e7/summary_*.tsv, gene_rank_* |
| `gwas_coloc/20_mras_panel_ldwin.R` | Builds the MRAS x cardiovascular GWAS coloc tables over all tissues and the rs11713141 direction per MRAS unit. | results_v3m master, analysis set, pairs, GWAS subsets, Figure 4 per-SNP statistics and LOO scans → results_v3m/mras/*.tsv (inputs of Fig. 6B, C, E) |

## alphagenome/

| file | what it does | inputs → outputs |
|---|---|---|
| `alphagenome/alphagenome_perbase_MRAS.py` | Calls the AlphaGenome API per variant (1-Mb interval) and saves per-base REF/ALT tracks over the four MRAS promoters. | variant list → per_snp/<variant>.npz (REF/ALT per modality: CAGE, PRO-cap, DNase, ATAC), meta/<MOD>_tracks.parquet |
| `alphagenome/alphagenome_perbase_MRAS_retry.py` | Copy of the prediction script with retry on quota and server errors and per-variant memory release (add-on variants). | add-on variant lists → alphagenome_results/perbase_MRAS_v3add/per_snp/*.npz |
| `alphagenome/mras_peak_finalize.py` | Calls the common predicted CAGE peak set on the adrenal sense-strand track averaged over the 95 REF predictions. | perbase_MRAS/per_snp/*.npz, meta/CAGE_tracks.parquet, exon table → tables/predicted_cage_peaks_COMMON.tsv, fig_peak_common.{pdf,png} |
| `alphagenome/mras_score_variants.py` | Scores the first 95 variants on the fixed CAGE and PRO-cap peaks and calls the common PRO-cap peak set. | perbase_MRAS/per_snp/*.npz, track metadata, predicted_cage_peaks_COMMON.tsv, old proxy tables → tables/predicted_procap_peaks_COMMON.tsv, score_*_peaks.tsv |
| `alphagenome/score_v3_candidates.py` | Scores all 170 predicted MRAS variants on the fixed CAGE and PRO-cap peaks (log2FC, P45112 share). | per_snp/*.npz of the 95- and 75-variant runs, predicted_*_peaks_COMMON.tsv, 161-candidate list → alphagenome_addon/tables/score_*_v3.tsv, ranking_v3_candidates.tsv |
| `alphagenome/config/MRAS_allpromoter_snplist.txt` | Lists the 95 MRAS-locus variants of the first AlphaGenome run (Supplementary Table 9, first-run variants). | read by perbase_MRAS.pbs → alphagenome_perbase_MRAS.py |
| `alphagenome/config/MRAS_v3_new_candidates_to_score.txt` | Lists the add-on MRAS candidate variants (74) predicted by one AlphaGenome add-on job (Supplementary Table 9). | read by perbase_MRAS_v3add.pbs (first attempt) |
| `alphagenome/config/MRAS_v3_new_candidates_to_score_ALL75.txt` | Lists the add-on MRAS candidate variants (75) predicted by one AlphaGenome add-on job (Supplementary Table 9). | read by perbase_MRAS_v3add_retry.pbs and score_v3_candidates.py |
| `alphagenome/config/MRAS_v3_candidate_set_LDstandard_PLINK.tsv` | Lists the 161 MRAS candidate variants (r2 >= 0.8 proxies of the 105 significant MRAS signals) with their signals and maximum r2 (Fig. 6D, Supplementary Table 9a). | read by score_v3_candidates.py (in_v3_candidate_set) and build_S7_to_S10.R (check) |
