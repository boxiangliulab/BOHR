#!/bin/bash
# Authored summary (not an original script): the commands of this step in the order they were run, with the
# arguments and environment of the original jobs. The trait scripts of part A carry the paths of the earlier runs; in the
# final run the folders were those of the no-floor tree (README.md, "Paths"). Set them for your system before running
# anything. Loops over tissues and arms ran as PBS array jobs.
set -euo pipefail
cd 04_tissue_specificity_traits

## A. Donor traits: one run folder for all 20 traits (lab cluster; the fits on NSCC)
T7=/projects/group/bliu/e1327989/GTEx/nofloor_full/GTEx/tissue_regression/trait_regression_20261007
T4=/projects/group/bliu/e1327989/GTEx/nofloor_full/GTEx/tissue_regression/trait_regression_20261104
# A1. inputs (work/<tissue>.rds, meta.rds; prep.R sources the data.R of run 20261007) and the mirror table (8 cores)
export TRG_RUN=$T7
Rscript traits/run_20261005/prep.R
Rscript traits/run_20261007/mirror.R
# A2. the fits: one NSCC job per tissue with a work/<tissue>.rds (57), its three arms in turn, 32 cores. work/ and the code
#     (data.R, run_trait.R, ../02_bohr_model_and_simulation/engine/bohr2_engine.R and bohr2_addon3.R) were copied to the
#     NSCC run folder, and results/{real,sim1,perm1} were copied back to $T7/results.
export TRG_RUN=/scratch/users/nus/e1327989/nofloor_full/tissue_regression/trait_regression_20261007_v4
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
for T in $(grep -v '^#' config/tissue_list.txt); do for ARM in real sim1 perm1; do
  Rscript traits/run_20261007/run_trait.R "$T" "$ARM" 32
done; done
# A3. summary run folder 20261104: links to the inputs and fits of run 20261007; work/mirror_v2.tsv = its mirror table
for f in $T7/work/*; do ln -sfn "$f" "$T4/work/$(basename "$f")"; done
rm -f $T4/work/mirror.tsv; cp $T7/work/mirror.tsv $T4/work/mirror_v2.tsv
for a in real sim1 perm1; do ln -sfn $T7/results/$a $T4/results/$a; done
export TRG_RUN=$T4
Rscript traits/run_20261104/mirror_v3.R
Rscript traits/run_20261104/summarize.R real sim1 perm1
Rscript traits/run_20261104/export_fig6.R
Rscript traits/run_20261106_mirror_abs/01_abs_labels.R
Rscript traits/run_20261106_mirror_abs/03_apply_B.R

## B. Final version (r5e80): its test set with BH again, then the rare binary tests dropped with BH again (2 cores each)
Rscript traits/r5e80_20261119/r5_01_units.R
Rscript traits/r5e80_ruleb_20261119/r6_01_rule_b.R

## C. Tissue specificity on the analysis-freeze samples: tissue means, ROKU and the 0.10 filter (16 cores);
##    the paper uses <out>/freeze56. Needs the BOHR package source (repository root) and TCC.
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
Rscript tissue_specificity/roku_freeze51.R /projects/group/bliu/e1327989/GTEx/nofloor_full/consult/freeze_ld_20261124/roku/out 16
