# Promoter quantification and filtering

**Methods §1; Supplementary Methods §1.** Junction reads give each promoter's activity in every sample (proActiv), and
the promoters analysed are the external promoters that pass the filter of Supplementary Methods §1: 92,557 external
promoters in 52,037 genes come down to 10,826 promoters in 4,744 genes.

**Run order.** `0101_transform.sh` and `0101_junction_filter.sh` (junction files, five-read floor) →
`0102_run_proActiv_tissue.R` (GENCODE v47 annotation; proActiv per tissue with at least 20 samples; 57 tissues) →
`0103_proActiv_data_process_tissue.R` (merge, filter, relative activity over the gene's external promoters) →
`rebuild_tabs_att.R` (rebuilds the absolute-activity file of 0103 by sample name) → `export_txt.R` (the text exports that
later steps read, and the PEER input; authored, see Notes). `release_comparison/` builds the
GENCODE v19 and v26 annotations and compares the promoter catalogues of GTEx releases (Fig. 1A,B counts).

**BOHR package.** `BOHR::relative_activity()` (with `BOHR::from_proactiv()`) reproduces step 0103: on the stored proActiv
outputs every filtering count of Supplementary Methods §1, the final set and all 10,826 × 22,865 values are identical
(repository `workflow/parity/12_relative_activity_parity.R`). The subset test (`../test_subset/`) repeats this on one
chromosome.

**Notes.** The 57 per-tissue proActiv result files that 0103 reads are no longer stored; its outputs are. 0103 writes a
misaligned absolute-activity file (lines 195–198 bind by position); `rebuild_tabs_att.R` rebuilt it by sample name, and
only the rebuilt file was used. The relative activity was never affected. 0103 writes RData only; later steps read text
exports (`filter_5_recal_att_filter.txt`, `Promoter_info_filter_5.txt`) and PEER reads a donor-mean matrix
(`filter_5_recal_att_filter_SMTSD_mean.txt`). These were written outside the scripts kept here; `export_txt.R` was
written for this repository and, on the paper's 0103 outputs, reproduces the stored files (the two exports exactly, the
PEER input to 1.1 × 10⁻¹⁵).

`run.sh` lists the commands in the order they were run (paths of the original runs).


## Files

| file | what it does | inputs → outputs |
|---|---|---|
| `0101_transform.sh` | Converts per-sample GTEx junction-count GCT files to STAR SJ.out.tab format, with the count in column 7. | split_junctions/<sample>.junction.gct.gz → star_junctions/<sample>.SJ.out.tab |
| `0101_junction_filter.sh` | Keeps junctions with at least five uniquely mapped reads (column 7 >= 5) in every per-sample SJ.out.tab file. | star_junctions/*.SJ.out.tab → filtered_junctions/*.SJ.out.tab |
| `0102_run_proActiv_tissue.R` | Builds the GENCODE v47 promoter annotation and runs proActiv per GTEx v11 tissue with at least 20 samples. | GENCODE v47 GTF, GTEx sample attributes, filtered SJ.out.tab per sample → promoterAnnotation.gencode.v47.RData, att_sample.info.RData, 57 proActiv_<tissue>_n<N>.rds |
| `0103_proActiv_data_process_tissue.R` | Merges the 57 proActiv runs, filters to 10,826 external promoters and recalculates relative activity. | 57 proActiv rds, att_sample.info.RData → filter_5_rdata, Promoter_info_all, trel/tabs_att_filter_In.rds, Promoter_info_filter_5, filter_5_recal_att_filter.RData |
| `rebuild_tabs_att.R` | Rebuilds the absolute-activity file filter_5_tabs_att.RData by sample name, fixing 0103's misaligned output. | tabs_att_filter_In.rds, filter_5_trel_att.RData → filter_5_tabs_att.RData (old file kept as .BAD.RData) |
| `export_txt.R` | Authored: writes the text exports of the 0103 outputs and the PEER input (mean over a donor's samples within each tissue, missing set to 0; the recipe of `0104_peer.R` lines 48–57). | filter_5_recal_att_filter.RData, Promoter_info_filter_5.RData → filter_5_recal_att_filter.txt, Promoter_info_filter_5.txt, filter_5_recal_att_filter_SMTSD_mean.txt |

## release_comparison/

| file | what it does | inputs → outputs |
|---|---|---|
| `release_comparison/0102_run_proActiv_v6.R` | Builds the GENCODE v19 promoter annotation for GTEx v6p (lines 20-23); its proActiv loop runs Whole Blood only. | GENCODE v19 GTF, GTEx v7 sample attributes, v6 filtered junctions → promoterAnnotation.gencode.v19.RData, att_sample_info_v6.RData, one proActiv rds |
| `release_comparison/run_proActiv_tissue.R` | Builds the GENCODE v26 promoter annotation for GTEx v8 (lines 12-15) and runs proActiv per v8 tissue. | GENCODE v26 GTF, GTEx v8 sample attributes, v8 filtered junctions → att_sample.info.RData, proActiv_<tissue>_n<N>.rds |
| `release_comparison/version_comparison_analysis.R` | Compares proActiv promoter annotations across GTEx releases; its Plot 8 and Plot 9C hold the Fig. 1A,B counts. | GTEx_promoter.csv, Promoter_info_all_v{7,8,10,11}.RData → Fig1-Fig9 comparison plots (PDF, PNG) |
