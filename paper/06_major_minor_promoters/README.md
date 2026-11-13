# Major and minor promoters

**Methods §9; Supplementary Methods §9.** In genes with two tested promoters, the promoter with the higher mean usage
pooled over tissues is the major one (in genes with three or more, the highest-usage promoter; the rest are minor).
The analysis asks which promoter's absolute activity each signal changes, how that relates to the GTEx eQTL of the gene,
and how ENCODE chromatin marks differ at major and minor TSSs (Fig. 5, Extended Data Fig. 5, Supplementary Table 6).

**Run order.** The scripts run in a stage folder (`F5S`, set in `00_helpers.R`) that holds the per-signal tables of
`../03_cis_mapping_conditional/conditional/fig5_delivery_v4.R`. `stage_v6/` (two-promoter genes):
`71_A1_location.R` (where each signal's variants lie) → `72_A2_eqtl_direction.R` (roles and GTEx eQTL direction) →
`73b_B1v2_partner_conditional.R` (own, partner and gene-output absolute activity on dosage) → `74b_B1v2_patterns.R`
(which promoter changes) → `97_I_continuous_effects.R` (μ/boundary signal pairs on the same promoter, read by 101) →
`101_K_overlap_multi_location.R` (per-signal pattern table) → `104_M_genotype_usage.R` (genotype-group usage states). Chromatin: `encode/25a_…` and `encode/22b2_…`
(TSS cCREs, ENCODE SCREEN Z-scores) → `122a_tissue_matched_chromatin_v4.R` → `134_major_minor_pooled_chromatin.R` (Fig. 5A),
`184_fig5a_major_minor_pvalues.R` (its P values) and `132_marks_by_usage_role.R` (Extended Data Fig. 5D). Example loci:
`124_example_candidates_eqtl.R`, `124b_…` and `124c_…` rank candidates by fixed rules (Fig. 5F–H). `stage_v6c/`
re-runs `101_K_overlap_multi_location.R` with the major promoter of genes with three or more promoters set by usage
pooled over tissues (the paper's rule), with the stage's `00_helpers.R` and `122a`. Script 97 was written for the burst
analysis, which the paper no longer reports; it is kept because 101 reads its pair table.

The analysis uses the BOHR results; it has no counterpart in the BOHR package.

`run.sh` lists the commands in the order they were run (paths of the original runs).


## stage_v6/

| file | what it does | inputs → outputs |
|---|---|---|
| `stage_v6/00_helpers.R` | Sets the stage paths and the shared helpers (gene-clustered bootstrap ratio tables, P values, file readers) for the stage-v6 scripts. | F5S stage folder → IN, OUTROOT, ratio_table(), boot_weights(), fr() |
| `stage_v6/71_A1_location.R` | Locates each signal's candidate variants relative to all promoters of the gene and lists within-unit μ/boundary pairs (inputs of 74b, 97 and 101). | delivery v6 signals(_partner), genotype_pheno_proxies, crosstype_conditional_pairs_all, signals_pertype, promoter_flags, Promoter_info_filter_5, freeze, GENCODE v26 → A1_location/A1_signal_location.tsv, A1_within_unit_pairs.tsv |
| `stage_v6/72_A2_eqtl_direction.R` | Assigns major/minor roles and the regulatory route of each signal in two-promoter genes, with the GTEx eQTL direction, on Figure 4 delivery v6. | delivery v6 signals, signals_partner, genotype_pheno(_proxies), mirror promoter_flags, eQTL subset, GENCODE v26 → A2_eqtl_sqtl/A2_signal_eqtl.tsv, A2_splice_region_by_signal.tsv |
| `stage_v6/73b_B1v2_partner_conditional.R` | Fits per-signal OLS models of own, partner and other-promoter absolute activity and gene output on lead dosage, conditional on the unit's other signals (delivery v6). | delivery v6 signals, mirror promoter_flags, proActiv absolute activity, Promoter_info_all, freeze donors and covariates, VCF → B1v2_multipromoter/B1v2_fits.tsv.gz |
| `stage_v6/74b_B1v2_patterns.R` | Codes which promoter's activity changes with each retained signal (minor only, major only, both, neither), oriented to the allele raising the minor promoter. | B1v2_fits, A2_signal_eqtl, A2_splice_region_by_signal, A1_signal_location, delivery v6 signals, promoter_flags → B1v2_multipromoter/B1v2_signal_patterns.tsv.gz |
| `stage_v6/97_I_continuous_effects.R` | Finds μ/boundary signal pairs on the same switching promoter (the 407 pairs behind Fig. 5F and Supplementary Table 6f) and fits linear share effects. | B1v2_fits, B1v2_signal_patterns, A1_within_unit_pairs, delivery v6 signals, absolute activity, freeze, VCF → I_continuous_effects/I1_F2_pairs.tsv.gz, I1/I2 tables |
| `stage_v6/101_K_overlap_multi_location.R` | Builds the per-signal pattern table (two-promoter and >= 3-promoter frames, shared-lead flags) that Fig. 5 r18 reads; >= 3-promoter genes still use the usage < 0.5 rule here. | B1v2_signal_patterns, A1_signal_location, delivery v6 signals and genotype_pheno_proxies, Promoter_info_filter_5, I1_F2_pairs → K_overlap_multi_location/K_signals.tsv.gz, K3_* |
| `stage_v6/104_M_genotype_usage.R` | Computes per-signal genotype-group states of the minor promoter (unused, exclusive, interior usage) that Fig. 5C and Supplementary Table 6e average. | delivery v6 genotype_pheno and signals_pe_moments, B1v2_signal_patterns, K_signals, A1, absolute activity, freeze, VCF → M_genotype_usage/M1_signal_genotype_states.tsv.gz and M1/M2 summaries |
| `stage_v6/122a_tissue_matched_chromatin_v4.R` | Matches ENCODE biosamples to GTEx tissues and calls DNase, H3K27ac and H3K4me3 at each tested TSS (Fig. 5A); 132, 134 and 147 evaluate its lines 1-83. | SCREEN Z-score matrices and TSS elements (fig5_v3data), Figure 4 v3c units (freeze_v3c) → PROMOTER_LEVEL_v4/28v4_tissue_matched_tss_chromatin.tsv, 28v4_tissue_map_used.tsv; objects MAP, MAPb, TP, state_for |
| `stage_v6/134_major_minor_pooled_chromatin.R` | Computes DNase, H3K27ac and H3K4me3 shares and TSS states at major versus minor promoters with gene-bootstrap CIs (Fig. 5A, Supplementary Table 6b). | 122a objects (lines 1-83 evaluated), mirror_pairs → PROMOTER_LEVEL_v4/134_major_minor_pooled_chromatin.tsv |
| `stage_v6/184_fig5a_major_minor_pvalues.R` | P values of the major-versus-minor chromatin comparison of Fig. 5A (gene-clustered bootstrap and Welch test). |  |
| `stage_v6/132_marks_by_usage_role.R` | Computes the share of matched ENCODE samples with each mark at major and minor TSSs by puQTL class (Extended Data Fig. 5B, Supplementary Table 6c). | 122a objects (lines 1-83 evaluated), mirror_pairs → PROMOTER_LEVEL_v4/132_marks_by_usage_role.tsv |
| `stage_v6/124_example_candidates_eqtl.R` | Ranks example-locus candidates for Fig. 5G-I (rules h, i, j) with each lead's gene-output P and GTEx v11 eQTL P. | K_signals, A1, B1v2_signal_patterns, proxies_pertype, GTEx eQTL subset, absolute activity, freeze, VCF → FIG5_example_candidates_v4/candidates_v4_eqtl.tsv, contact sheets |
| `stage_v6/124b_example_i_minor_on.R` | Ranks μ major-only example candidates with the minor promoter active in at least 80% of donors (sheet for Fig. 5H). | as 124 → FIG5_example_candidates_v4/candidates_v4_eqtl_i_minor_on.tsv |
| `stage_v6/124c_example_j_strong.R` | Ranks μ redistribution example candidates with larger opposite effects (sheet for Fig. 5I, TNK2), with gene-output and GTEx eQTL P. | as 124 → FIG5_example_candidates_v4/candidates_v4_eqtl_j_strong.tsv |

## stage_v6c/

| file | what it does | inputs → outputs |
|---|---|---|
| `stage_v6c/00_helpers.R` | Sets the stage paths and the shared helpers for the stage-v6c scripts (same code as the stage-v6 helpers). | F5S stage folder → IN, OUTROOT, ratio_table(), boot_weights(), fr() |
| `stage_v6c/101_K_overlap_multi_location.R` | Builds the per-signal pattern table with the major promoter of genes with >= 3 tested promoters set by usage pooled over tissues (Supplementary Table 6d, Extended Data Fig. 5D). | as code_v6stage/101, plus freeze activity matrices and mirror_pairs (pooled-usage check) → K_overlap_multi_location/K_signals.tsv.gz, K2_pooled_major_three_plus.tsv.gz |
| `stage_v6c/122a_tissue_matched_chromatin_v4.R` | Provides the ENCODE biosample map and TSS mark matrices that the stage-v6c 147 evaluates (lines 1-83) for Supplementary Table 6a. | fig5_v3data, freeze_v3c units → objects MAP, MAPb, TP, state_for (in memory) |

## encode/ and config/

| file | what it does | inputs → outputs |
|---|---|---|
| `encode/25a_alltested_tss_elements.R` | Assigns each tested promoter the overlapping ENCODE cCRE (TSS ±200 bp) whose midpoint is closest to the TSS. | Figure 4 v3 freeze *real.units.rds, Promoter_info_filter_5.txt, ENCODE GRCh38-cCREs.bed → data/alltested_tss_elements.tsv, data/zscore_accessions_alltested.txt |
| `encode/22b2_extract_encode_zscores_alltested.py` | Extracts ENCODE SCREEN core-collection DNase, H3K27ac and H3K4me3 Z-scores for the TSS cCREs chosen by 25a. | SCREEN human_*_core_collection_zscores.tar.gz, zscore_accessions_alltested.txt → data/zscores_{dnase,h3k27ac,h3k4me3}_alltested.tsv.gz, *_columns.tsv |
| `config/tissues.tsv` | Lists the 50 freeze tissues with sample counts and the number of PEER factors (K_C3) used as covariates. | -> read as in/freeze/tissues.tsv by 71, 73b and 97 (rows with def == 'freeze') |
