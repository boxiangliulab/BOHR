#!/bin/bash
# Authored summary (not an original script): the commands of this step in the order they were run, with the
# arguments and environment of the original jobs. Paths are those of the original runs (lab cluster and NSCC); set them
# for your system before running anything. Loops over tissues, arms and chromosomes ran as PBS array or feeder jobs.
set -euo pipefail
cd 01_promoter_quantification
# 1. junction files (GTEx per-sample junction counts -> STAR SJ.out.tab format), then the five-read floor
bash 0101_transform.sh
bash 0101_junction_filter.sh
# 2. GENCODE v47 promoter annotation and proActiv per tissue (57 tissues with >= 20 samples; ran on NSCC)
Rscript 0102_run_proActiv_tissue.R
# 3. merge the 57 runs, filter promoters (SI 1), relative activity over the gene's external promoters
#    (run in the folder that holds ./proActiv/*.rds and att_sample.info.RData)
Rscript 0103_proActiv_data_process_tissue.R
# 4. rebuild the absolute-activity file of step 3 by sample name
Rscript rebuild_tabs_att.R
# release comparison (Fig. 1A,B): annotations of GTEx v6p (GENCODE v19) and v8 (v26), then the comparison
Rscript release_comparison/0102_run_proActiv_v6.R
Rscript release_comparison/run_proActiv_tissue.R
Rscript release_comparison/version_comparison_analysis.R
