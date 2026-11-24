#!/bin/bash
# Authored summary (not an original script): the commands of this step in the order they were run, with the
# arguments and environment of the original jobs. Paths are those of the original runs (lab cluster and NSCC); set them
# for your system before running anything. Loops over tissues, arms and chromosomes ran as PBS array or feeder jobs.
set -euo pipefail
cd 02_bohr_model_and_simulation
# check against gamlss (Supplementary Methods 2): pilot inputs, then refits with gamlss; run 2 used gamlss_check/bohr2_engine.stage1b.R
Rscript gamlss_check/prep_pilot.R
OUTSUB=stage1_equiv_run2 Rscript gamlss_check/equiv_test.R
NS="50 100 150 200 250 300 350 400 450 500 550 600"
# simulated data, generated once and reused by the final version
#   simulated activity (M1): 2,000 genes per sample size
for N in $NS; do TRUNC_L=0.005 GEN_ONLY=1 Rscript simulated_activity/newsim_direct_nothreshold_indep.R $N 2000 2; done
#   simulated RNA-seq (M2): polyester reads (batch 1 seeds 42/123, batch 2 seeds 43/223), then STAR (NSCC)
for N in $NS; do Rscript simulated_rnaseq/v12_sim_run_nothr_big.R $N $SIM_DIR 42 123 prob; done
# 1. proActiv on the unfiltered STAR junction files, no read floor (NSCC, one job per batch and N, 4 cores);
#    SIM_DIR = method2_polyester/sim_<N>_nothr_big (batch 1) or sim_<N>_nothr_big_b2 (batch 2); output in <SIM_DIR>/nofloor/
for N in $NS; do for SUFFIX in nothr_big nothr_big_b2; do
  Rscript simulated_rnaseq/sim_proactiv_nofloor.R method2_polyester/sim_${N}_${SUFFIX}
done; done
#    copied to the lab as newmethod_bohr2/m2_nofloor/<b1|b2>/sim_<N>/ (recal_att_f.RData, gene_reads.tsv), with the design table
#    and genotypes of the same run (test_output/sim_tx_info_extended.txt, genotypes.vcf.gz)
# 2. fits with the final promoter rules (80% boundary-share skip, SD < 0.005 skip): BOHR and the linear tests, power and
#    calibration (1,000 genotype permutations of the calibration genes, permutation seed base 20261001)
export SKIP_FRAC=0.80 SD_SKIP=0
O=results/e80_20261119/nofloor
for N in $NS; do
  Rscript benchmark/bohr2_sim.R power M1 nothr2000 $N $O/m1/power/N$N.rds 8
  Rscript benchmark/baselines_v3.R power M1 nothr2000 $N $O/baselines_M1/power/N$N.rds
  PERM_IDS=config/m1_nothr_g2000/ids_calib_random150.txt Rscript benchmark/bohr2_sim.R calib M1 nothr2000 $N $O/m1/calib/N$N.rds 8 1000
  PERM_IDS=config/m1_nothr_g2000/ids_calib_random150.txt Rscript benchmark/baselines_v3.R calib M1 nothr2000 $N $O/baselines_M1/calib/N$N.rds
  for B in b1 b2; do   # simulated RNA-seq batches, quantified without a read floor
    Rscript benchmark/bohr2_sim.R power M2 nf_$B $N $O/$B/power/N$N.rds 8
    Rscript benchmark/baselines_v3.R power M2 nf_$B $N $O/baselines_M2$B/power/N$N.rds
    PERM_IDS=config/m2_nothr/$B/ids_calib_random150genes.txt Rscript benchmark/bohr2_sim.R calib M2 nf_$B $N $O/$B/calib/N$N.rds 8 1000
    PERM_IDS=config/m2_nothr/$B/ids_calib_random150genes.txt Rscript benchmark/baselines_v3.R calib M2 nf_$B $N $O/baselines_M2$B/calib/N$N.rds
  done
done
# 3. calibration with four more permutation sets (Supplementary Table 2b): seed bases 20262025, 20282025, 20292025, 20302025
for TAG in 2 3 4 5; do
  case $TAG in 2) BASE=20262025 ;; *) BASE=$((20252025 + TAG * 10000)) ;; esac
  for N in $NS; do
    PERM_SEED_BASE=$BASE PERM_IDS=config/m1_nothr_g2000/ids_calib_random150.txt Rscript benchmark/bohr2_sim.R calib M1 nothr2000 $N $O/m1/calib_seed$TAG/N$N.rds 8 1000
    PERM_SEED_BASE=$BASE PERM_IDS=config/m1_nothr_g2000/ids_calib_random150.txt Rscript benchmark/baselines_v3.R calib M1 nothr2000 $N $O/baselines_M1_seed$TAG/calib/N$N.rds
    for B in b1 b2; do
      PERM_SEED_BASE=$BASE PERM_IDS=config/m2_nothr/$B/ids_calib_random150genes.txt Rscript benchmark/bohr2_sim.R calib M2 nf_$B $N $O/$B/calib_seed$TAG/N$N.rds 8 1000
      PERM_SEED_BASE=$BASE PERM_IDS=config/m2_nothr/$B/ids_calib_random150genes.txt Rscript benchmark/baselines_v3.R calib M2 nf_$B $N $O/baselines_M2${B}_seed$TAG/calib/N$N.rds
    done
  done
done
unset SKIP_FRAC SD_SKIP
# 4. summaries. DATASET=E80_nofloor_reads5 reads the runs above and drops simulated RNA-seq genes with a median of fewer than
#    five gene reads. A first pass with the default COMBINE (Simes) wrote data_v3_E80_nofloor_reads5/; of it, only the p0 and p1
#    calibration files of combined_p_and_p0p1_v3.R (diag_tail/) enter the final figures and tables.
DATASET=E80_nofloor_reads5 Rscript benchmark/prepare_data_v3.R
DATASET=E80_nofloor_reads5 Rscript benchmark/combined_p_and_p0p1_v3.R
#    final pass: the combined BOHR P value of Methods, Section 2 (Fisher, Sidak-corrected 0/1 part) -> data_v3_E80_nofloor_reads5_fishersidak/
export DATASET=E80_nofloor_reads5 COMBINE=fisher_sidak
Rscript benchmark/prepare_data_v3.R
EFFECT_N=300 Rscript benchmark/prepare_power_by_effect.R
EFFECT_N=300 Rscript benchmark/prepare_power_by_OR.R
EFFECT_N=300 OR_BINS=lin5 Rscript benchmark/prepare_power_by_OR.R
#    the p0 and p1 calibration files do not depend on COMBINE; the first pass's folder was linked into the final one
ln -s ../data_v3_E80_nofloor_reads5/diag_tail figures_bohr2/data_v3_E80_nofloor_reads5_fishersidak/diag_tail
MS_DATASET=E80_nofloor_reads5 Rscript benchmark/multiseed_final_v3.R
