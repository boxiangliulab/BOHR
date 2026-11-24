#!/bin/bash
# Authored summary (not an original script): the commands of this step in the order they were run, with the
# arguments and environment of the original jobs. Paths are those of the original runs (lab cluster and NSCC); set them
# for your system before running anything. Loops over tissues, arms and chromosomes ran as PBS array or feeder jobs.
set -euo pipefail
cd 01_promoter_quantification
# 0. earlier run, outputs reused: GENCODE v47 promoter annotation and the sample table att_sample.info.RData
#    (its proActiv loop on the five-read junction files is not part of the final run)
Rscript 0102_run_proActiv_tissue.R
# 1. proActiv per tissue on every junction of the GTEx v11 STAR tables (no read floor); NSCC array 1-57, 16 cores each,
#    run in the folder that holds the annotation object (built with the call of 0102 lines 11-14, see README)
for k in $(seq 1 57); do Rscript run_proActiv_arm.R "$k" 16 nofloor; done
# 2. merge the 57 runs, filter promoters (SI 1), relative activity over the gene's external promoters
#    (NSCC, run in /scratch/users/nus/e1327989/junction_floor/nofloor, which holds ./proActiv/*.rds and att_sample.info.RData)
Rscript 0103_proActiv_data_process_tissue.R
# 3. rebuild the absolute-activity file of step 2 by sample name (setwd() pointed at the folder of step 2)
Rscript rebuild_tabs_att.R
#    the outputs of steps 2-3 were then copied to the lab (GTEx/v11_nofloor/, md5 checked)
# 4. text exports read by later steps (relative activity, promoter table) and the PEER input of
#    ../03_cis_mapping_conditional/inputs/peer/: on NSCC into the run tree, and again on the lab copy
Rscript export_txt.R /scratch/users/nus/e1327989/junction_floor/nofloor \
  /scratch/users/nus/e1327989/nofloor_full/GTEx/v11 /scratch/users/nus/e1327989/nofloor_full/GTEx/peer/v11
Rscript export_txt.R /projects/group/bliu/e1327989/GTEx/nofloor_full/GTEx/v11 /projects/group/bliu/e1327989/GTEx/nofloor_full/GTEx/v11 \
  /projects/group/bliu/e1327989/GTEx/nofloor_full/fork_checks/peer_input_lab
# release comparison (Fig. 1A,B; earlier run, not repeated): GTEx v8 junction counts to SJ.out.tab files, junctions with
# at least five uniquely mapped reads, annotations of GTEx v6p (GENCODE v19) and v8 (v26), then the comparison
bash 0101_transform.sh
bash 0101_junction_filter.sh
Rscript release_comparison/0102_run_proActiv_v6.R
Rscript release_comparison/run_proActiv_tissue.R
Rscript release_comparison/version_comparison_analysis.R
