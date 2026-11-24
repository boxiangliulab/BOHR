# Major and minor promoters

**Methods §9; Supplementary Methods §9.** In genes with two tested promoters, the promoter with the higher mean usage
pooled over tissues is the major one (in genes with three or more, the highest-usage promoter; the rest are minor).
The analysis asks which promoter's absolute activity each signal changes, how that relates to the GTEx eQTL of the gene,
and how ENCODE chromatin marks differ at major and minor TSSs (Fig. 5, Extended Data Fig. 5, Supplementary Table 6).

**Run order.** In the final run:
1. `encode/187a_alltested_tss_elements_v3c.R` (the TSS cCRE of every promoter tested in units v3c) →
   `encode/187b_extract_encode_zscores_alltested_v3c.py` (ENCODE SCREEN DNase, H3K27ac and H3K4me3 Z-scores of these
   elements) → `encode/187c_check_zscores_v3c.R` (the Z-scores equal the earlier extracts on the elements in both). This
   ran next to the eQTL colocalization of `../05_coloc_eqtl/`.
2. After `../05_coloc_eqtl/`: `prep/00_build_analysis_set.R` → `prep/01_extract_eqtl.py` per tissue (the GTEx v11 eQTLs
   within ± 200 kb of every analysis signal, which 72 and 124* read).
3. The stage folders (`F5S`; commands in `run.sh`, no script makes them). `in/` links the per-signal tables of
   `../03_cis_mapping_conditional/conditional/fig5_delivery_v4.R` and the other outputs of the final run, and `code/`
   holds `00_helpers.R`, `theme_nature.R` and the copies of 122a and 134 that other scripts read (table below).
4. `stage_v6/` (two-promoter genes): `71_A1_location.R` (where each signal's variants lie) → `72_A2_eqtl_direction.R`
   (roles and GTEx eQTL direction) → `73b_B1v2_partner_conditional.R` (own, partner and gene-output absolute activity on
   dosage) → `74b_B1v2_patterns.R` (which promoter changes) → `97_I_continuous_effects.R` (μ/boundary signal pairs on the
   same promoter, read by 101) → `101_K_overlap_multi_location.R` (per-signal pattern table) → `104_M_genotype_usage.R`
   (genotype-group usage states). Chromatin: `122a_tissue_matched_chromatin_v4.R` → `132_marks_by_usage_role.R` (Extended Data Fig. 5D),
   `134_major_minor_pooled_chromatin.R` (Fig. 5A) and `184_fig5a_major_minor_pvalues.R` (its P values). Example loci:
   `124_example_candidates_eqtl.R`, `124b_…` and `124c_…` rank candidates by fixed rules (Fig. 5F–H).
5. `stage_v6c/101_K_overlap_multi_location.R` in the second stage folder, with the major promoter of genes with three or
   more promoters set by usage pooled over tissues (the paper's rule), with the stage's `00_helpers.R` and `122a`.

Script 97 was written for the burst analysis, which the paper no longer reports; it is kept because 101 reads its pair
table. The final run also ran two side analyses, 122c (122a without the exact mirror pairs) and 129 (chromatin states by
role and puQTL class); no figure or table of the paper reads their outputs, so they are not included. The analysis uses the BOHR results; it has no counterpart in
the BOHR package.

**Stage folders.** `in/` of the v6 stage holds links; the v6c stage links the same `in/` and the v6 stage's outputs of
71–104 and of the chromatin chain (`out/A1_location`, `A2_eqtl_sqtl`, `B1v2_multipromoter`, `I_continuous_effects`,
`M_genotype_usage`, `PROMOTER_LEVEL_v4`).

| in the stage | from |
|---|---|
| `in/v11/` | `Promoter_info_all.RData`, `Promoter_info_filter_5.txt`, `tabs_att_filter_In.rds` of the promoter quantification (`../01_promoter_quantification/`) |
| `in/freeze/` | the per-tissue cis-mapping inputs `work/final/freeze/<tissue>.rds` (individual-level) and `work/final/tissues.tsv` (= `config/tissues.tsv`) |
| `in/freeze_v3c` | `results/final_nogate_v3c/freeze` (units v3c) |
| `in/delivery` | `results/fig5_delivery_v6` (per-signal tables of `../03_cis_mapping_conditional/conditional/`) |
| `in/mirror/` | `results/mirror/mirror_pairs.tsv.gz`, `results/mirror_v3/promoter_flags.tsv.gz` |
| `in/coloc/` | `eqtl_subset/` and `eqtl_windows.tsv` (step 2), `coloc_eqtl_signals.tsv.gz` (the final eQTL coloc table, `../05_coloc_eqtl/v3m/06_final_v3m.R`) |
| `in/coloc_final_clump/` | the same `coloc_eqtl_signals.tsv.gz` |
| `in/fig5_v3data/` | the outputs of 187a and 187b under the names 122a reads (`alltested_tss_elements.tsv`, `zscores_<mark>_alltested.tsv.gz`, `zscores_<mark>_alltested_columns.tsv`), and `gencode.v47.annotation.gtf.gz` |
| `in/gtf`, `vcf`, `encode`, `external`, `jaspar`, `tata`, `tpm`, `snp_level`, `fig5_old`, `alphagenome` | the earlier stage's inputs, which do not depend on promoter activity |
| `code/` (v6) | `stage_v6/00_helpers.R`, `stage_v6/theme_nature.R`, `stage_v6/122a_tissue_matched_chromatin_v4.R`, `stage_v6/134_major_minor_pooled_chromatin.R` |
| `code/` (v6c) | `stage_v6c/00_helpers.R`, `stage_v6/theme_nature.R`, `stage_v6c/122a_tissue_matched_chromatin_v4.R` |

In the earlier run, 72 read the coloc table of the signal ± 200 kb eQTL run; in the final run
`in/coloc/coloc_eqtl_signals.tsv.gz` is the final eQTL coloc table. Only the descriptive coloc columns of
`A2_signal_eqtl.tsv` depend on it.

**Changes for the final run.**
- 71, 72 and 73b asserted 14,823 per-signal rows (the earlier run's count); the final run asserts a non-empty table.
- 72 reads the eQTL coloc table only when it exists (otherwise an empty table with the same columns) and skips its two
  coloc summaries when the table is empty. The table was present in the final run, so nothing changed there.
- 187a and 187b replace `25a_alltested_tss_elements.R` and `22b2_extract_encode_zscores_alltested.py` of the earlier
  run. 187a takes the tested promoters from units v3c (25a took them from units v3, so promoters tested only in later
  unit sets had no element); 187b reads the accession list of 187a and skips a mark only when the extracted matrix is
  missing (the archives that 22b2 checked for had been moved). Element rule, cCRE file, window, matrix reading and
  output format are those of 25a and 22b2.
- `prep/00_build_analysis_set.R` is `../05_coloc_eqtl/v3/00_build_analysis_set.R` with three changes: signals from
  `conditional_final_v3c`, output folder `results/coloc_eqtl_v3_fig5` (the eQTL chain of `../05_coloc_eqtl/` owns
  `results/coloc_eqtl_v3`), and a stop before the final comparison with Figure 6's analysis set when that file is
  absent. `prep/01_extract_eqtl.py` is `../05_coloc_eqtl/v3/01_extract_eqtl.py` with the same output folder.

**Paths.** The scripts keep the paths of the earlier runs. In the final run they ran as copies in which only folder
names were changed: `/projects/group/bliu/e1327989/GTEx` and `/home/e1327989/project/GTEx` →
`/projects/group/bliu/e1327989/GTEx/nofloor_full/GTEx` (links to the inputs that do not depend on promoter activity, new
folders for those that do); `QTL_figure/rerun_v2` → `nofloor_full/GTEx/QTL_figure/rerun_r5e80`;
`TATA_Enhancer/Figure5_rerun_v4/b1_tss_elements_v3c_20261113` →
`nofloor_full/GTEx/TATA_Enhancer/Figure5_rerun_r5e80/b1_tss_elements_v3c_20261113`; the coloc library sourced by
`prep/00_build_analysis_set.R` → the run copy of `../05_coloc_eqtl/v3c/coloc_final_lib.R` (read with `UNITSET=v3c`).
The stage scripts read `F5S` (the stage folder; `nofloor_full/f5_stage_nofloor_r5e80_v6` and `_v6c`); their defaults
name the earlier stage folders. 187a and 187c compare with the earlier element table and Z-score extracts
(`/home/e1327989/f5_stage_20261104_v4/in/fig5_v3data`, copies of the outputs of 25a and 22b2).

**Counts written into the scripts.** `stage_v6c/101_K_overlap_multi_location.R` (line 37) asserts 50 tissues, and 184
asserts B = 2,000 and the seed of 122a and 134. These hold in the final run.

`run.sh` lists the commands in the order they were run.


## encode/ and prep/

| file | what it does | inputs → outputs |
|---|---|---|
| `encode/187a_alltested_tss_elements_v3c.R` | Assigns each promoter tested in units v3c the overlapping ENCODE cCRE (TSS ± 200 bp) whose midpoint is closest to the TSS, and compares with the earlier table. | final_nogate_v3c freeze *real.units.rds, Promoter_info_filter_5.txt, ENCODE GRCh38-cCREs.bed → data/alltested_tss_elements_v3c.tsv, data/zscore_accessions_alltested_v3c.txt, data/promoters_new_in_v3c.tsv |
| `encode/187b_extract_encode_zscores_alltested_v3c.py` | Extracts the ENCODE SCREEN core-collection DNase, H3K27ac and H3K4me3 Z-scores of the TSS cCREs chosen by 187a. | SCREEN human_cCRE_<mark>_core_collection_zscore.npy, GRCh38-cCREs.Rows.bed, biosample columns, zscore_accessions_alltested_v3c.txt → data/zscores_{dnase,h3k27ac,h3k4me3}_alltested_v3c.tsv.gz, *_columns.tsv |
| `encode/187c_check_zscores_v3c.R` | Checks that the new Z-score extracts equal the earlier ones on every shared accession and that the column files are identical; stops otherwise. | data/zscores_*_alltested_v3c*, earlier extracts → log |
| `prep/00_build_analysis_set.R` | Builds the analysis set (every backward-kept signal of every significant unit) and the ± 200 kb window of each signal per tissue and gene. | final_nogate_v3c (read_units_final, UNITSET v3c), conditional_final_v3c, GTEx v11 eGenes, ld_v3 donors → results/coloc_eqtl_v3_fig5/{analysis_set.rds, eqtl_windows.tsv} |
| `prep/01_extract_eqtl.py` | Cuts the GTEx v11 eQTL all-pairs to the genes and windows of one tissue. | coloc_eqtl_v3_fig5/eqtl_windows.tsv, GTEx v11 all-pairs parquet → coloc_eqtl_v3_fig5/eqtl_subset/<tissue>.tsv.gz |

## stage_v6/

| file | what it does | inputs → outputs |
|---|---|---|
| `stage_v6/00_helpers.R` | Sets the stage paths and the shared helpers (gene-clustered bootstrap ratio tables, P values, file readers) for the stage-v6 scripts. | F5S stage folder → IN, OUTROOT, ratio_table(), boot_weights(), fr() |
| `stage_v6/theme_nature.R` | Provides the ggplot2 theme that 122a sources for its plots; copied into the `code/` folder of both stages. | sourced from F5S/code/ |
| `stage_v6/71_A1_location.R` | Locates each signal's candidate variants relative to all promoters of the gene and lists within-unit μ/boundary pairs (inputs of 74b, 97 and 101). | delivery v6 signals(_partner), genotype_pheno_proxies, crosstype_conditional_pairs_all, signals_pertype, promoter_flags, Promoter_info_filter_5, freeze, GENCODE v26 → A1_location/A1_signal_location.tsv, A1_within_unit_pairs.tsv |
| `stage_v6/72_A2_eqtl_direction.R` | Assigns major/minor roles and the regulatory route of each signal in two-promoter genes, with the GTEx eQTL direction. | delivery v6 signals, signals_partner, genotype_pheno(_proxies), mirror promoter_flags, eQTL subset, eQTL coloc table, GENCODE v26 → A2_eqtl_sqtl/A2_signal_eqtl.tsv, A2_splice_region_by_signal.tsv |
| `stage_v6/73b_B1v2_partner_conditional.R` | Fits per-signal OLS models of own, partner and other-promoter absolute activity and gene output on lead dosage, conditional on the unit's other signals. | delivery v6 signals, mirror promoter_flags, proActiv absolute activity, Promoter_info_all, freeze donors and covariates, VCF → B1v2_multipromoter/B1v2_fits.tsv.gz |
| `stage_v6/74b_B1v2_patterns.R` | Codes which promoter's activity changes with each retained signal (minor only, major only, both, neither), oriented to the allele raising the minor promoter. | B1v2_fits, A2_signal_eqtl, A2_splice_region_by_signal, A1_signal_location, delivery v6 signals, promoter_flags → B1v2_multipromoter/B1v2_signal_patterns.tsv.gz |
| `stage_v6/97_I_continuous_effects.R` | Finds μ/boundary signal pairs on the same switching promoter (read by 101) and fits linear share effects. | B1v2_fits, B1v2_signal_patterns, A1_within_unit_pairs, delivery v6 signals, absolute activity, freeze, VCF → I_continuous_effects/I1_F2_pairs.tsv.gz, I1/I2 tables |
| `stage_v6/101_K_overlap_multi_location.R` | Builds the per-signal pattern table (two-promoter and >= 3-promoter frames, shared-lead flags); >= 3-promoter genes still use the usage < 0.5 rule here. | B1v2_signal_patterns, A1_signal_location, delivery v6 signals and genotype_pheno_proxies, Promoter_info_filter_5, I1_F2_pairs → K_overlap_multi_location/K_signals.tsv.gz, K3_* |
| `stage_v6/104_M_genotype_usage.R` | Computes per-signal genotype-group states of the minor promoter (unused, exclusive, interior usage) that Fig. 5C and Supplementary Table 6e average. | delivery v6 genotype_pheno and signals_pe_moments, B1v2_signal_patterns, K_signals, A1, absolute activity, freeze, VCF → M_genotype_usage/M1_signal_genotype_states.tsv.gz and M1/M2 summaries |
| `stage_v6/122a_tissue_matched_chromatin_v4.R` | Matches ENCODE biosamples to GTEx tissues and calls DNase, H3K27ac and H3K4me3 at each tested TSS (Fig. 5A); 129, 132 and 134 evaluate its lines 1-83. | SCREEN Z-score matrices and TSS elements (fig5_v3data), units v3c (freeze_v3c) → PROMOTER_LEVEL_v4/28v4_tissue_matched_tss_chromatin.tsv, 28v4_tissue_map_used.tsv; objects MAP, MAPb, TP, state_for |
| `stage_v6/132_marks_by_usage_role.R` | Computes the share of matched ENCODE samples with each mark at major and minor TSSs by puQTL class (Extended Data Fig. 5D, Supplementary Table 6c). | 122a objects (lines 1-83 evaluated), mirror_pairs → PROMOTER_LEVEL_v4/132_marks_by_usage_role.tsv |
| `stage_v6/134_major_minor_pooled_chromatin.R` | Computes DNase, H3K27ac and H3K4me3 shares and TSS states at major versus minor promoters with gene-bootstrap CIs (Fig. 5A, Supplementary Table 6b). | 122a objects (lines 1-83 evaluated), mirror_pairs → PROMOTER_LEVEL_v4/134_major_minor_pooled_chromatin.tsv |
| `stage_v6/184_fig5a_major_minor_pvalues.R` | Computes the P values of the major-versus-minor chromatin comparison of Fig. 5A (gene-clustered bootstrap and Welch test). | 134 objects (evaluated up to its bootstrap function) → PROMOTER_LEVEL_v4/184_fig5a_major_minor_pvalues.tsv |
| `stage_v6/124_example_candidates_eqtl.R` | Ranks example-locus candidates for Fig. 5F–H and Supplementary Table 6f (rules h, i and j of the script) with each lead's gene-output P and GTEx v11 eQTL P. | K_signals, A1, B1v2_signal_patterns, proxies_pertype, GTEx eQTL subset, absolute activity, freeze, VCF → FIG5_example_candidates_v4/candidates_v4_eqtl.tsv, contact sheets |
| `stage_v6/124b_example_i_minor_on.R` | Ranks μ major-only example candidates (rule i) with the minor promoter active in at least 80% of donors. | as 124 → FIG5_example_candidates_v4/candidates_v4_eqtl_i_minor_on.tsv |
| `stage_v6/124c_example_j_strong.R` | Ranks μ redistribution example candidates (rule j) with larger opposite effects, with gene-output and GTEx eQTL P. | as 124 → FIG5_example_candidates_v4/candidates_v4_eqtl_j_strong.tsv |

## stage_v6c/

| file | what it does | inputs → outputs |
|---|---|---|
| `stage_v6c/00_helpers.R` | Sets the stage paths and the shared helpers for the stage-v6c scripts (same code as the stage-v6 helpers). | F5S stage folder → IN, OUTROOT, ratio_table(), boot_weights(), fr() |
| `stage_v6c/101_K_overlap_multi_location.R` | Builds the per-signal pattern table with the major promoter of genes with >= 3 tested promoters set by usage pooled over tissues (Fig. 5B, Supplementary Table 6d). | as stage_v6/101, plus freeze activity matrices and mirror_pairs (pooled-usage check) → K_overlap_multi_location/K_signals.tsv.gz, K2_pooled_major_three_plus.tsv.gz |
| `stage_v6c/122a_tissue_matched_chromatin_v4.R` | Provides the ENCODE biosample map and TSS mark matrices (lines 1-83) for Supplementary Table 6a; staged in `code/` of the v6c stage. | fig5_v3data, freeze_v3c units → objects MAP, MAPb, TP, state_for (in memory) |

## config/

| file | what it does | inputs → outputs |
|---|---|---|
| `config/tissues.tsv` | Lists the 50 tissues with sample counts, the number of PEER factors (K_C3) used as covariates and the promoters tested, as in the final run. | -> read as in/freeze/tissues.tsv by 71, 73b, 97 and stage_v6c/101 (rows with def == 'freeze') |
