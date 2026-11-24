#!/bin/bash
# Authored summary (not an original script): the commands of this step in the order they were run, with the
# arguments and environment of the original jobs. Part A is the earlier run, whose outputs the final run reads; part B
# is the final run. The scripts of steps 1, 2 and 5 carry the paths of the earlier runs; in part B the folders are those
# of the final run (README.md, "Paths"). Set them for your system before running anything. Loops over GWAS files and
# tissues ran as PBS array or xargs jobs.
set -euo pipefail
cd 07_coloc_gwas_alphagenome
export UNITSET=v3c                        # read by coloc_final_lib.R (../05_coloc_eqtl/v3c/coloc_final_lib.R), which every coloc script sources

## A. Earlier run: GWAS inputs, significant GWAS SNPs, European donor lists, AlphaGenome predictions at MRAS
python3 gwas_inputs/prep_magic_eur.py
Rscript gwas_inputs/prep_ra_adipo_eur.R
for i in $(seq 1 51); do Rscript gwas_coloc/10_gwas_sig_snps.R $i; done        # coloc_ldunion_20261012/gwas_sig/
Rscript gwas_coloc/ld_v3_eur_mras.R                                            # results/ld_v3/donors/*.EUR.txt (RACE == 3)
# AlphaGenome at MRAS (needs ALPHAGENOME_API_KEY in the environment; jobs alphagenome/perbase_MRAS*.pbs, AG_SEQLEN=1048576)
python -u alphagenome/alphagenome_perbase_MRAS.py alphagenome/config/MRAS_allpromoter_snplist.txt perbase_MRAS CAGE,PROCAP,DNASE,ATAC          # 95 variants
python -u alphagenome/alphagenome_perbase_MRAS.py alphagenome/config/MRAS_v3_new_candidates_to_score.txt perbase_MRAS_v3add CAGE,PROCAP,DNASE,ATAC   # first add-on attempt
#    the added variants were completed with the retry copy, over the full list and parts of it (one job per part):
python -u alphagenome/alphagenome_perbase_MRAS_retry.py alphagenome/config/MRAS_v3_new_candidates_to_score_ALL75.txt perbase_MRAS_v3add CAGE,PROCAP,DNASE,ATAC
python alphagenome/mras_peak_finalize.py
python alphagenome/mras_score_variants.py
python alphagenome/score_v3_candidates.py

## B. Final run
R5=/projects/group/bliu/e1327989/GTEx/nofloor_full/GTEx/QTL_figure/rerun_r5e80                 # rerun_v2 in the scripts
K=/projects/group/bliu/e1327989/GTEx/nofloor_full/consult/freeze_ld_20261124/gwas                     # clumping and steps 12-20
# B1. puQTL side (r5e80 run, stage S6.1-S6.2; NDIR = nofloor_full/GTEx/new_GWAS/coloc_final_nogate_v3_r5e80)
Rscript gwas_coloc/01_build_analysis_set.R
Rscript gwas_coloc/12a_puqtl_windows.R
bash gwas_coloc/12b_extra_puqtl_windows_nscc.sh        # NSCC job, 8 cores; W = /scratch/users/nus/e1327989/nofloor_full/fig6_ldwin_r5e80
#    extra_win/puqtl_extra_windows.tsv copied to coloc_ldunion_20261012/nscc_out/. The r5e80 run then ran steps 12-20 with
#    the 943-donor clumps of the earlier run for every GWAS file; B4 replaced that run.
# B2. MRAS candidates of the final version; AlphaGenome prediction of the one new candidate (conda env alphagenome,
#     API key in ALPHAGENOME_API_KEY); scores of all 171 predicted variants
Rscript alphagenome/20_candidates.R
export AG_SEQLEN=1048576
python -u alphagenome/alphagenome_perbase_MRAS_retry.py alphagenome/config/r5e80_new_candidate_to_score.txt perbase_MRAS_r5e80add CAGE,PROCAP,DNASE,ATAC
python -u alphagenome/score_r5e80_candidates.py
# B3. GWAS clumping (step 11) on the lab cluster, with all 943 donors and with the 799 European donors
#     (jobs gwas_coloc/clump_ALL.pbs and clump_EUR.pbs, 8 cores each; GWAS files whose output has DONE are skipped)
export CLUMP_PLINK=/home/e1327989/.local/share/mamba/envs/gtexfig/bin/plink CLUMP_REF=/projects/group/bliu/e1327989/GTEx/Rawdata/LD_ref_v9_plink19
SIG=/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/coloc_ldunion_20261012/gwas_sig
for v in ALL EUR; do
  if [ $v = ALL ]; then export CLUMP_DONORS=$CLUMP_REF/donors/ALL.txt; else export CLUMP_DONORS=$CLUMP_REF/donors_EUR/ALL.EUR.txt; fi
  cat config/traits_51_basenames.txt | xargs -P 2 -I{} Rscript gwas_coloc/11_clump_local.R {} $SIG $K/clump_$v 4
done
Rscript gwas_coloc/make_eurmix.R                        # clump_EURmix/: European clumping for 32 files, 943 donors for 19
# B4. steps 12-20 with the ancestry-matched clumps (job gwas_coloc/coloc_eur.pbs, 16 cores; root
#     $K/coloc_final_nogate_v3_r5e80_eur, whose coloc_ldunion_20261012/nscc_out/clump is a link to clump_EURmix)
Rscript gwas_coloc/12_pairs.R
seq 1 51 | xargs -P 8 -I{} Rscript gwas_coloc/13_subset_gwas.R {}
awk -F'\t' 'NR > 1 && $2 == "freeze" {print $1}' $R5/work/final/tissues.tsv | xargs -P 4 -I{} Rscript gwas_coloc/14_coloc_tissue.R {} 4
Rscript gwas_coloc/15_merge.R
TIER=1e-7 Rscript gwas_coloc/16_make_canonical.R
TMAIN=1e-7 Rscript gwas_coloc/20_mras_panel_ldwin.R
