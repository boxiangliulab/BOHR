#!/bin/bash
# Authored summary (not an original script): the commands of this step in the order they were run, with the
# arguments and environment of the original jobs. Paths are those of the original runs (lab cluster and NSCC); set them
# for your system before running anything. Loops over tissues, arms and chromosomes ran as PBS array or feeder jobs.
set -euo pipefail
cd 05_coloc_eqtl
TISSUES=$(cat ../07_coloc_gwas_alphagenome/config/tissues_50.txt)   # the 50 tissues
# 0. units v3 signals (analysis set), and the earlier signal +-200 kb run, superseded by the LD windows but read by
#    ldwin/07_merge_classify.R (its all-shared rule tables)
Rscript v3/00_build_analysis_set.R
for T in $TISSUES; do python3 v3/01_extract_eqtl.py "$T"; Rscript v3/02_coloc_tissue.R "$T" 8; done
Rscript v3/03_merge_classify.R
Rscript v3/04_classify_all_shared.R
# 1. LD windows (PLINK 1.9 on each tissue's donors), pairs with the GTEx top eQTL variant, coloc
Rscript ldwin/00_leads.R
#    PLINK: --r2 --ld-window-r2 0.1 --ld-window-kb 1000 around each lead (ld/ld_tool.R standard): NSCC job
#    ldwin/01_ldwin_nscc.pbs; eQTL leads absent from the MAF >= 0.05 panel: NSCC job ldwin/rarelead.pbs (unfiltered VCF)
#    panel/panel_snps.tsv.gz (read by 03, 08 and v3m/06) = columns 1, 2 and 4 of the panel's PLINK .bim files:
for c in $(seq 1 22); do cut -f1,2,4 LD_ref_v9_plink19/bfile/chr$c.bim; done | gzip > panel/panel_snps.tsv.gz
Rscript ldwin/02_windows.R
Rscript ldwin/03_rarelead_windows.R
Rscript ldwin/04_pairs.R
cut -f1 eqtl_windows.tsv | tail -n +2 | sort -u | xargs -P 6 -I{} python3 ldwin/05_extract_eqtl.py {}
for T in $TISSUES; do Rscript ldwin/06_coloc_tissue.R "$T" 8; done
Rscript ldwin/07_merge_classify.R
# 2. clumped eQTL leads of each eGene, coloc against every lead, best call per signal (2026-10-12 table)
for T in $TISSUES; do python3 ldwin/clumped_leads/s1_sig_pairs.py "$T"; Rscript ldwin/clumped_leads/s2_clump_nscc.R "$T" sig out 16; done   # NSCC: s2_clump_nscc.pbs
Rscript ldwin/clumped_leads/s3_pairs.R
for T in $TISSUES; do python3 ldwin/clumped_leads/s4_extract_eqtl.py "$T"; Rscript ldwin/clumped_leads/s5_coloc_tissue_v2.R "$T" 8; done
Rscript ldwin/clumped_leads/s6_merge_compare.R     # signals_sens.tsv.gz, read by 08_final_clump.R
Rscript ldwin/08_final_clump.R
# 3. signals changed by the mirror rule (v3m)
Rscript v3m/01_analysis_set_pairs.R
Rscript v3m/02_new_inputs.R
python3 v3m/02b_sig_pairs.py
#    NSCC job v3m/03_ld_clump_nscc.pbs: LD windows of the new leads and clumping of the new eGenes
Rscript v3m/04_pairs_new_all.R
for K in top clump; do tail -n +2 eqtl_windows_${K}_new.tsv | cut -f1 | sort -u | xargs -P 4 -I{} python3 v3m/04b_extract_eqtl.py {} $K; done
cat tissues_new.txt | xargs -P 4 -I{} Rscript v3m/05_coloc_tissue_v3m.R {} 4
Rscript v3m/06_final_v3m.R
# 4. the two refitted units (v3c): the paper's table
Rscript v3c/c1_analysis_set_pairs_v3c.R
for T in Artery_Aorta Brain_Putamen_basal_ganglia; do Rscript v3c/c5_coloc_tissue_v3c.R $T 2; done
Rscript v3c/c6_final_v3c.R
Rscript v3c/c7_textclean_v3c.R
