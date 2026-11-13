#!/bin/bash
# Authored summary (not an original script): the commands of this step in the order they were run, with the
# arguments and environment of the original jobs. Paths are those of the original runs (lab cluster and NSCC); set them
# for your system before running anything. Loops over tissues, arms and chromosomes ran as PBS array or feeder jobs.
set -euo pipefail
cd 06_major_minor_promoters
export F5S=/path/to/stage_v6           # stage folder: in/ holds the per-signal tables of fig5_delivery_v4.R (delivery v6)
# ENCODE TSS elements and Z-scores
Rscript encode/25a_alltested_tss_elements.R
python3 encode/22b2_extract_encode_zscores_alltested.py
# two-promoter genes
for s in 71_A1_location 72_A2_eqtl_direction 73b_B1v2_partner_conditional 74b_B1v2_patterns 97_I_continuous_effects \
         101_K_overlap_multi_location 104_M_genotype_usage 134_major_minor_pooled_chromatin 184_fig5a_major_minor_pvalues \
         132_marks_by_usage_role 124_example_candidates_eqtl 124b_example_i_minor_on 124c_example_j_strong; do
  Rscript stage_v6/$s.R
done
# the pattern table with the paper's rule for genes with three or more promoters (stage v6c)
F5S=/path/to/stage_v6c Rscript stage_v6c/101_K_overlap_multi_location.R
