#!/bin/bash
# Authored summary (not an original script): the commands of this step in the order they were run, with the
# arguments and environment of the original jobs. Paths are those of the original runs (lab cluster and NSCC); set them
# for your system before running anything. Loops over tissues, arms and chromosomes ran as PBS array or feeder jobs.
set -euo pipefail
export BOHR2_CODE=/path/to/code            # a folder holding ../02_bohr_model_and_simulation/engine/bohr2_*.R and these scripts
export BOHR2_ROOT=/path/to/run_root        # work/, perm_maps/, results/
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
TISSUES=$(cut -f1 tissues_50.txt)        # the 50 tissues (work/final/tissues.tsv, def == freeze)
# 1. inputs: PEER factors per tissue (57 jobs), the chromosome 22 PEER scan, then the per-tissue input objects
for T in $TISSUES; do Rscript inputs/peer/peer_analysis.R "$T" filter_5_recal_att_filter_SMTSD_mean.txt peer/v11; done
Rscript inputs/prep_all.R
# 2. per-SNP scan, per tissue and arm (real, perm0 to perm2), chromosome by chromosome; Firth add-on; F-corrected mu
for T in $TISSUES; do for ARM in real perm0; do for c in $(seq 1 22); do
  Rscript scan/run_final2.R freeze "$T" "$ARM" chr$c 64
  Rscript scan/run_firth_addon.R freeze "$T" "$ARM" chr$c 64
done; done; done
Rscript scan/run_muF.R 8
# 3. units, ACAT and BH; estimability and mirror rules; labels by absolute activity (units v3c)
Rscript units/mirror_pairs.R
NC=8 Rscript units/discovery_nogate.R freeze real $TISSUES
Rscript units/units_v2_estimable.R
Rscript units/units_v3_estimable.R
Rscript units/units_v3m_mirror.R
# 4. independent signals (see README: incremental runs, the latest per unit is kept)
for T in $TISSUES; do for c in $(seq 1 22); do Rscript conditional/run_conditional_ng.R freeze "$T" real chr$c 64; done; done
#    run_conditional_ng2.R and run_conditional_ng3.R for the (tissue, chr) lists of new and continued units, same arguments
COND_FINAL=results/conditional_final_v3 Rscript conditional/merge_conditional_ng3.R
COND_FINAL=results/conditional_final_v3 Rscript conditional/loo_est_patch.R
COND_FINAL=results/conditional_final_v3 Rscript conditional/merge_conditional_ng3.R
Rscript conditional/cond_v3m_select.R
#    for each (tissue, chr) of work/cond_v3m/units_to_run.tsv:
#    Rscript conditional/run_conditional_v3m.R freeze <tissue> real <chr> 8
Rscript conditional/merge_conditional_v3m.R
#    refit of two units (conditional/units_to_run_v3c.tsv):
#    Rscript conditional/run_conditional_v3c.R freeze <tissue> real <chr> 8
LABELS=none Rscript conditional/merge_conditional_v3c.R
COND_FINAL=results/conditional_final_v3c LABEL_DIR=results/mirror_abs_labels_v3c NC=10 Rscript units/mirror_abs_labels.R
LABELS=results/mirror_abs_labels_v3c/labels_real.tsv.gz UNITS_OUT=results/final_nogate_v3c Rscript units/units_v3b_from_labels.R
LABELS=results/mirror_abs_labels_v3c/labels_real.tsv.gz Rscript conditional/merge_conditional_v3c.R
#    per-signal tables (delivery v6, read by ../06_major_minor_promoters/); environment COND_FINAL, UNITS_DIR and FIG5_OUT as above
export COND_FINAL=results/conditional_final_v3c UNITS_DIR=results/final_nogate_v3c FIG5_OUT=results/fig5_delivery_v6
NC=10 Rscript conditional/fig5_delivery_v4.R
NC=10 Rscript conditional/fig5_proxies.R
NC=10 Rscript conditional/fig5_pertype.R
NC=10 Rscript conditional/fig5_pe_moments_v6.R
NC=10 Rscript conditional/crosstype_conditional_v6.R
FLAGS=results/mirror_v3/promoter_flags.tsv.gz Rscript conditional/fig5_signals_partner.R
# 5. LD standard: donor lists and signal LD
NC=12 Rscript ld/ld_v3.R
# 6. QTLtools comparison (Fig. 2E): permutation maps, BED and covariate files, QTLtools v1.2, counts
Rscript qtltools_comparison/make_perm_maps_freeze_extra.R
for T in $TISSUES; do Rscript qtltools_comparison/make_freeze_beds.R "$T"; done
#    for each tissue, arm (real, perm0-2), mode (raw, rint; --normal on raw) and chromosome:
#    QTLtools cis --vcf chr$c.vcf.gz --bed $T.$arm.$kind.bed.gz --cov $T.$arm.cov.txt --permute 1000 \
#      --window 500000 --seed 123 [--normal] --region chr$c --out prm.$arm.chr$c.txt
FIG2E_TIS=ALL Rscript qtltools_comparison/fig2e_deliver_freeze_v3m_pairs.R
