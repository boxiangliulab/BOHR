#!/bin/bash
# Authored summary (not an original script): the commands of this step in the order they were run, with the
# arguments and environment of the original jobs. Paths are those of the original runs (lab cluster and NSCC); set them
# for your system before running anything. Loops over tissues, arms and chromosomes ran as PBS array or feeder jobs.
set -euo pipefail
cd 04_tissue_specificity_traits
# tissue specificity: raw_fit() of gamlss_functions_optimized.R gives fit/v11_raw_fit_fit_results.txt (one mean per
# promoter and tissue); then ROKU and the 0.10 effect filter
Rscript tissue_specificity/ROKU.R
Rscript tissue_specificity/roku_effect_filter.R
# donor traits: inputs, then one job per tissue and arm (real, sim1, perm1); run 20261007 adds drinker and smoker status
Rscript traits/run_20261005/prep.R
Rscript traits/run_20261005/mirror.R
for T in $(cat config/tissue_list.txt); do for ARM in real sim1 perm1; do
  Rscript traits/run_20261005/run_trait.R "$T" "$ARM" 8
  for TRAIT in drinker smoker; do Rscript traits/run_20261007/run_trait.R "$T" "$ARM" 8 "$TRAIT"; done
done; done
# mirror rule, testability rules and BH; export; exact-pair events labelled by absolute activity (the paper's results)
Rscript traits/run_20261104/mirror_v3.R
Rscript traits/run_20261104/summarize.R
Rscript traits/run_20261104/export_fig6.R
Rscript traits/run_20261106_mirror_abs/01_abs_labels.R
Rscript traits/run_20261106_mirror_abs/03_apply_B.R
