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
# donor traits, run 20261005: inputs (work/<tissue>.rds, meta.rds), mirror table, then one job per tissue and arm
# (real, sim1, perm1) for all traits; 18 of its 20 traits are the paper's
Rscript traits/run_20261005/prep.R
Rscript traits/run_20261005/mirror.R
for T in $(cat config/tissue_list.txt); do for ARM in real sim1 perm1; do
  Rscript traits/run_20261005/run_trait.R "$T" "$ARM" 8
done; done
# run 20261007: drinking (MHDRNKSTS) and smoking (MHSMKSTS) recoded and refitted. Its work/<tissue>.rds and meta.rds are
# links to run 20261005's work/, and its results of the other 18 traits are links to run 20261005's results.
export TRG_RUN=/path/to/trait_regression_20261007     # run root of 20261007 (code/, work/, results/)
Rscript traits/run_20261007/mirror.R                   # work/mirror.tsv (mirror table v2 pairs)
for T in $(cat config/tissue_list.txt); do for ARM in real sim1 perm1; do
  for TRAIT in MHDRNKSTS MHSMKSTS; do Rscript traits/run_20261007/run_trait.R "$T" "$ARM" 8 "$TRAIT"; done
done; done
# run 20261104 refits nothing: its work/ and results/{real,sim1,perm1} are links to run 20261007, and its
# work/mirror_v2.tsv is a link to run 20261007's work/mirror.tsv (read by mirror_v3.R).
# mirror rule, testability rules and BH; export; exact-pair events labelled by absolute activity (the paper's results)
Rscript traits/run_20261104/mirror_v3.R
Rscript traits/run_20261104/summarize.R
Rscript traits/run_20261104/export_fig6.R
Rscript traits/run_20261106_mirror_abs/01_abs_labels.R
Rscript traits/run_20261106_mirror_abs/03_apply_B.R
