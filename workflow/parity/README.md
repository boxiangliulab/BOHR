# Real-data parity check of `bohr_cis()` (stage 2)

Compares `BOHR::bohr_cis()` with the stored results of the paper's original pipeline
(`rerun_v2/results/final`, `final_muF`, `final_firth_addon` and `final_nogate_v3m`) on one GTEx tissue and
chromosome. GTEx individual-level data never enter this repository: the scripts read them from the lab's
storage and write only to a scratch directory given on the command line. The comparison prints aggregate
numbers only (counts, largest differences, agreement rates).

1. `01_extract_inputs.R <tissue> <chr> <inputs.rds>` (R with data.table): rebuilds the inputs of the stored run
   as plain matrices: activity of the tissue's promoters on the chromosome, GT dosages of every VCF record for the
   analysis samples (missing kept; `bohr_cis()` fills them), PC1-5 and PEER 1..K, and the pooled-mean score
   that picks the canonical member of mirror pairs.
2. `02_run_bohr_cis.R <inputs.rds> <out.rds> <ncores> [max_promoters] [lib]` (R with BOHR installed).
3. `03_compare_stored.R <out.rds> <tissue> <chr> [summary.rds]`: per-SNP statistics, promoter rules, units
   (ACAT P, leader, estimability, mirror flags) and the tissue-wide q values, recomputed with this chromosome's
   units replaced by the package's.

Paths to the lab storage are set at the top of scripts 1 and 3. Results: `BOHR_v4/docs/code_release/STAGE1_2_REPORT.md`.
