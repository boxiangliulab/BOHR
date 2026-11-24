# Promoter quantification and filtering

**Methods §1; Supplementary Methods §1.** Junction reads give each promoter's activity in every sample (proActiv), and
the promoters analysed are the external promoters that pass the filter of Supplementary Methods §1: 92,557 external
promoters in 52,037 genes come down to 22,939 promoters in 9,371 genes.

**Run order.** The final run used every junction of the GTEx v11 STAR junction tables, with no minimum read count:
`run_proActiv_arm.R` (proActiv per tissue with at least 20 samples, 57 tissues, one NSCC job per tissue; arm `nofloor`) →
`0103_proActiv_data_process_tissue.R` (merge, filter, relative activity over the gene's external promoters) →
`rebuild_tabs_att.R` (rebuilds the absolute-activity file of 0103 by sample name) → `export_txt.R` (the text exports that
later steps read, and the PEER input). `run_proActiv_arm.R` repeats the proActiv call of `0102_run_proActiv_tissue.R`
with the same GENCODE v47 annotation, samples, tissue labels and sample order; it reads the sample table
`att_sample.info.RData` that 0102 wrote. The proActiv loop of 0102 itself read junction files with at least five uniquely
mapped reads and is not part of the final run. `release_comparison/` builds the GENCODE v19 and v26 annotations and
compares the promoter catalogues of GTEx releases (Fig. 1A,B counts; not repeated for the final run). Its GTEx v8 input
was made by `0101_transform.sh` (v8 junction counts to SJ.out.tab files) and `0101_junction_filter.sh` (junctions with at
least five uniquely mapped reads); the v11 junction tables are GTEx's own STAR files and pass through neither script.

**BOHR package.** `BOHR::relative_activity()` (with `BOHR::from_proactiv()`) reproduces step 0103. It was checked on the
stored proActiv outputs of the earlier five-read run: every filtering count, the final set and all 10,826 × 22,865 values
are identical (repository `workflow/parity/12_relative_activity_parity.R`).

**Notes.** The 57 per-tissue proActiv result files of the final run stay on NSCC scratch; the 0103 outputs were copied to
the lab (`GTEx/v11_nofloor/`, md5 checked). 0103 writes a misaligned absolute-activity file (lines 195–198 bind by
position); `rebuild_tabs_att.R` rebuilt it by sample name, and only the rebuilt file was used. The relative activity was
never affected. In the final run `rebuild_tabs_att.R` ran in the NSCC folder of the 0103 outputs: only the folder in
`setwd()` differed. 0103 writes RData only; later steps read text exports (`filter_5_recal_att_filter.txt`,
`Promoter_info_filter_5.txt`) and PEER reads a donor-mean matrix (`filter_5_recal_att_filter_SMTSD_mean.txt`).
`export_txt.R` was written for this analysis and ran in the final run, on NSCC for the cis-mapping inputs and PEER and on
the lab copy of the 0103 outputs; on the outputs of the earlier five-read run it reproduces the stored files (the two
exports exactly, the PEER input to 1.1 × 10⁻¹⁵).

**Exception.** `run_proActiv_arm.R` reads the annotation object `promoterAnnotation.gencode.v47.rds`, which was built
once by the first lines of a separate junction-count script that is not included here. That script makes the same
`preparePromoterAnnotation()` call on the same GENCODE v47 GTF as `0102_run_proActiv_tissue.R` lines 11–14.

`run.sh` lists the commands in the order they were run (paths of the original runs).


## Files

| file | what it does | inputs → outputs |
|---|---|---|
| `0101_transform.sh` | Converts per-sample GTEx v8 junction-count GCT files to STAR SJ.out.tab format, with the count in column 7 (release comparison). | split_junctions/<sample>.junction.gct.gz → star_junctions/<sample>.SJ.out.tab |
| `0101_junction_filter.sh` | Keeps junctions with at least five uniquely mapped reads (column 7 >= 5) in every per-sample SJ.out.tab file (release comparison; not used by the final run). | star_junctions/*.SJ.out.tab → filtered_junctions/*.SJ.out.tab |
| `0102_run_proActiv_tissue.R` | Builds the GENCODE v47 promoter annotation and the sample table, and runs proActiv per GTEx v11 tissue with at least 20 samples on the five-read junction files (earlier run). The final run reused its sample table. | GENCODE v47 GTF, GTEx sample attributes, filtered SJ.out.tab per sample → promoterAnnotation.gencode.v47.RData, att_sample.info.RData, 57 proActiv_<tissue>_n<N>.rds |
| `run_proActiv_arm.R` | Runs proActiv for one tissue (index 1–57) on every junction of the GTEx v11 STAR tables (arm `nofloor`), with the annotation, samples and tissue labels of 0102. | tissue index, cores, arm; promoterAnnotation.gencode.v47.rds, att_sample.info.RData, Star_out_tab/<sample>.SJ.out.tab.gz → junction_floor/nofloor/proActiv/proActiv_<tissue>_n<N>.rds |
| `0103_proActiv_data_process_tissue.R` | Merges the 57 proActiv runs, filters to 22,939 external promoters and recalculates relative activity. | 57 proActiv rds, att_sample.info.RData → filter_5_rdata, Promoter_info_all, trel/tabs_att_filter_In.rds, Promoter_info_filter_5, filter_5_recal_att_filter.RData |
| `rebuild_tabs_att.R` | Rebuilds the absolute-activity file filter_5_tabs_att.RData by sample name, fixing 0103's misaligned output. | tabs_att_filter_In.rds, filter_5_trel_att.RData → filter_5_tabs_att.RData (old file kept as .BAD.RData) |
| `export_txt.R` | Writes the text exports of the 0103 outputs and the PEER input (mean over a donor's samples within each tissue, missing set to 0; the recipe of `0104_peer.R` lines 48–57). | filter_5_recal_att_filter.RData, Promoter_info_filter_5.RData → filter_5_recal_att_filter.txt, Promoter_info_filter_5.txt, filter_5_recal_att_filter_SMTSD_mean.txt |

## release_comparison/

| file | what it does | inputs → outputs |
|---|---|---|
| `release_comparison/0102_run_proActiv_v6.R` | Builds the GENCODE v19 promoter annotation for GTEx v6p (lines 20-23); its proActiv loop runs Whole Blood only. | GENCODE v19 GTF, GTEx v7 sample attributes, v6 filtered junctions → promoterAnnotation.gencode.v19.RData, att_sample_info_v6.RData, one proActiv rds |
| `release_comparison/run_proActiv_tissue.R` | Builds the GENCODE v26 promoter annotation for GTEx v8 (lines 12-15) and runs proActiv per v8 tissue. | GENCODE v26 GTF, GTEx v8 sample attributes, v8 filtered junctions → att_sample.info.RData, proActiv_<tissue>_n<N>.rds |
| `release_comparison/version_comparison_analysis.R` | Compares proActiv promoter annotations across GTEx releases; its Plot 8 and Plot 9C hold the Fig. 1A,B counts. | GTEx_promoter.csv, Promoter_info_all_v{7,8,10,11}.RData → Fig1-Fig9 comparison plots (PDF, PNG) |
