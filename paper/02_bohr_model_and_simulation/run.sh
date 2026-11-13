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
# simulated activity (M1): 2,000 genes per sample size
for N in $NS; do TRUNC_L=0.005 GEN_ONLY=1 Rscript simulated_activity/newsim_direct_nothreshold_indep.R $N 2000 2; done
# simulated RNA-seq (M2): polyester reads (batch 1 seeds 42/123, batch 2 seeds 43/223), STAR, five-read junction floor, proActiv
for N in $NS; do Rscript simulated_rnaseq/v12_sim_run_nothr_big.R $N $SIM_DIR 42 123 prob; done
Rscript simulated_rnaseq/sim_proactiv_realjunc.R $SIM_DIR
# fits: BOHR and the linear tests (power; calibration = 1,000 genotype permutations of the 150 calibration genes)
for N in $NS; do
  Rscript benchmark/bohr2_sim.R power M1 nothr2000 $N results/m1nothr_g2000/power/N$N.rds 8
  PERM_IDS=config/m1_nothr_g2000/ids_calib_random150.txt Rscript benchmark/bohr2_sim.R calib M1 nothr2000 $N results/m1nothr_g2000/calib/N$N.rds 8 1000
  Rscript benchmark/baselines_v3.R power M1 nothr2000 $N results/baselines_v3/nothr_M1/power/N$N.rds
  PERM_IDS=config/m1_nothr_g2000/ids_calib_random150.txt Rscript benchmark/baselines_v3.R calib M1 nothr2000 $N results/baselines_v3/nothr_M1/calib/N$N.rds
  for B in b1 b2; do   # RNA-seq batches; same calls with M2, dataset nothr_realjunc and config/m2_nothr/$B/ids_calib_random150genes.txt
    Rscript benchmark/bohr2_sim.R power M2 nothr_realjunc $N results/m2nothr_realjunc/$B/power/N$N.rds 8
  done
done
# summaries (Fig. 2, Extended Data Fig. 2, Supplementary Table 2)
DATASET=nothr_realjunc Rscript benchmark/prepare_data_v3.R
MS_DATASET=nothr_realjunc Rscript benchmark/multiseed_final_v3.R
Rscript benchmark/combined_p_and_p0p1_v3.R
EFFECT_N=300 Rscript benchmark/prepare_power_by_effect.R
