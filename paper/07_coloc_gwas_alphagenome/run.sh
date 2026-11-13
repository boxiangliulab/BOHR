#!/bin/bash
# Authored summary (not an original script): the commands of this step in the order they were run, with the
# arguments and environment of the original jobs. Paths are those of the original runs (lab cluster and NSCC); set them
# for your system before running anything. Loops over tissues, arms and chromosomes ran as PBS array or feeder jobs.
set -euo pipefail
cd 07_coloc_gwas_alphagenome
export UNITSET=v3c                        # read by coloc_final_lib.R, which every script sources (copy: ../05_coloc_eqtl/v3c/coloc_final_lib.R)
# GWAS inputs that needed conversion
python3 gwas_inputs/prep_magic_eur.py
Rscript gwas_inputs/prep_ra_adipo_eur.R
# coloc with 51 GWAS traits in 50 tissues
Rscript gwas_coloc/01_build_analysis_set.R
for i in $(seq 1 51); do Rscript gwas_coloc/10_gwas_sig_snps.R $i; done
for TR in $(cat config/traits_51_basenames.txt); do Rscript gwas_coloc/11_clump_nscc.R "$TR" sig out 8; done   # NSCC: gwas_coloc/11_clump_nscc.pbs
Rscript gwas_coloc/12a_puqtl_windows.R
bash gwas_coloc/12b_extra_puqtl_windows_nscc.sh
#    12_pairs.R also reads run_<UNITSET>/extra_windows*.tsv: config/extra_windows_v3m_fig4nscc.tsv, the windows of the 25
#    leads new to units v3m, taken from the v3m LD windows of ../05_coloc_eqtl (v3m/03_ld_clump_nscc.pbs; no script)
Rscript gwas_coloc/12_pairs.R
for i in $(seq 1 51); do Rscript gwas_coloc/13_subset_gwas.R $i; done
for T in $(cat config/tissues_50.txt); do Rscript gwas_coloc/14_coloc_tissue.R "$T" 8; done   # NSCC: v3m_14_ / v3c_14_coloc_tissue.pbs
Rscript gwas_coloc/15_merge.R
TIER=1e-7 Rscript gwas_coloc/16_make_canonical.R
TMAIN=1e-7 Rscript gwas_coloc/20_mras_panel_ldwin.R
# AlphaGenome at MRAS (needs ALPHAGENOME_API_KEY in the environment; jobs alphagenome/perbase_MRAS*.pbs). The candidate set
# alphagenome/config/MRAS_v3_candidate_set_LDstandard_PLINK.tsv is the union of r2 >= 0.8 PLINK proxies (lead included)
# of the 105 significant v3 MRAS signals under the LD standard of ../03_cis_mapping_conditional/ld/; no script was kept.
python -u alphagenome/alphagenome_perbase_MRAS.py alphagenome/config/MRAS_allpromoter_snplist.txt perbase_MRAS CAGE,PROCAP,DNASE,ATAC          # 95 variants
python -u alphagenome/alphagenome_perbase_MRAS.py alphagenome/config/MRAS_v3_new_candidates_to_score.txt perbase_MRAS_v3add CAGE,PROCAP,DNASE,ATAC   # first add-on attempt
#    the add-on variants were completed with the retry copy, over the full list and parts of it (one job per part):
python -u alphagenome/alphagenome_perbase_MRAS_retry.py alphagenome/config/MRAS_v3_new_candidates_to_score_ALL75.txt perbase_MRAS_v3add CAGE,PROCAP,DNASE,ATAC
python alphagenome/mras_peak_finalize.py
python alphagenome/mras_score_variants.py
python alphagenome/score_v3_candidates.py
