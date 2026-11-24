#!/bin/bash
# Authored summary (not an original script): the commands of this step in the order they were run, with the
# arguments and environment of the original jobs. Part A is the earlier LD-window run, whose LD windows, clumped eQTL
# leads and eQTL subsets the final run reads; part B is the final run. The scripts carry the paths of the earlier runs;
# in part B the folders are those of the final run (README.md, "Paths"). Set them for your system before running
# anything. Loops over tissues ran as PBS array or feeder jobs.
set -euo pipefail
cd 05_coloc_eqtl

## A. Earlier run (units v3, 2026-10): analysis set, LD windows, clumped eQTL leads, eQTL subsets, coloc
TISSUES=$(cat ../07_coloc_gwas_alphagenome/config/tissues_50.txt)   # the 50 tissues
# A0. units v3 signals (analysis set), and the signal +-200 kb run, read by ldwin/07_merge_classify.R (all-shared rule)
Rscript v3/00_build_analysis_set.R
for T in $TISSUES; do python3 v3/01_extract_eqtl.py "$T"; Rscript v3/02_coloc_tissue.R "$T" 8; done
Rscript v3/03_merge_classify.R
Rscript v3/04_classify_all_shared.R
# A1. LD windows (PLINK 1.9 on each tissue's donors), pairs with the GTEx top eQTL variant, coloc
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
# A2. clumped eQTL leads of each eGene, coloc against every lead, best call per signal (the earlier final table)
for T in $TISSUES; do python3 ldwin/clumped_leads/s1_sig_pairs.py "$T"; Rscript ldwin/clumped_leads/s2_clump_nscc.R "$T" sig out 16; done   # NSCC: s2_clump_nscc.pbs
Rscript ldwin/clumped_leads/s3_pairs.R
for T in $TISSUES; do python3 ldwin/clumped_leads/s4_extract_eqtl.py "$T"; Rscript ldwin/clumped_leads/s5_coloc_tissue_v2.R "$T" 8; done
Rscript ldwin/clumped_leads/s6_merge_compare.R     # signals_sens.tsv.gz, read by 08_final_clump.R
Rscript ldwin/08_final_clump.R

## B. Final run: the v3m chain once, on the final units and signals (v3c), every signal new
R5=/projects/group/bliu/e1327989/GTEx/nofloor_full/GTEx/QTL_figure/rerun_r5e80   # the run root (rerun_v2 in the scripts)
O=$R5/results/coloc_eqtl_v3m
NW=/scratch/users/nus/e1327989/nofloor_full/eqtl_ld_r5e80                        # NSCC work folder (W of 03_ld_clump_nscc.pbs)
export UNITSET=v3c                                                               # coloc_final_lib.R: units v3c, conditional_final_v3c
# B0. set-up (no script): the earlier LD-window run as a read-only cache, an empty earlier analysis set (no coloc result
#     is carried), the output folders
ln -s /projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/results/coloc_eqtl_ldwin $R5/results/coloc_eqtl_ldwin
mkdir -p $R5/results/coloc_eqtl_v3 $O/new_inputs/sig
Rscript -e "library(data.table); saveRDS(data.table(tissue = character(), promoter = character(), type = character(), k = integer(), signal_snp = character(), stat = character()), '$R5/results/coloc_eqtl_v3/analysis_set.rds')"
# B1. analysis set and pairs, the inputs still missing (lab job, 4 cores)
Rscript v3m/01_analysis_set_pairs.R && Rscript v3m/02_new_inputs.R && python3 v3m/02b_sig_pairs.py
# B2. NSCC (128 cores): LD windows of the new leads and clumping of the new eGenes. Copied to NW: new_inputs/snps/*.txt
#     (to snps/), new_inputs/sig/*.tsv.gz (to sig/), ldwin/clumped_leads/s2_clump_nscc.R and v3m/03_ld_clump_nscc.pbs
#     (in the copy of the job that ran, W and the job log pointed to NW)
mkdir -p $NW/snps $NW/sig $NW/logs
cp $O/new_inputs/snps/*.txt $NW/snps/; cp $O/new_inputs/sig/*.tsv.gz $NW/sig/; cp ldwin/clumped_leads/s2_clump_nscc.R v3m/03_ld_clump_nscc.pbs $NW/
(cd $NW && qsub 03_ld_clump_nscc.pbs)
#     after the job: its outputs back to the lab
mkdir -p $O/nscc; cp -r $NW/out $NW/clump_out $NW/logs $O/nscc/
# B3. all pairs and the eQTL subsets still needed (lab job, 16 cores)
Rscript v3m/04_pairs_new_all.R
for K in top clump; do tail -n +2 $O/eqtl_windows_${K}_new.tsv | cut -f1 | sort -u | xargs -P 8 -I{} python3 v3m/04b_extract_eqtl.py {} $K; done
# B4. coloc per tissue, 4 tissues at a time with 4 cores each (lab job, 16 cores)
Rscript -e "cat(unique(readRDS('$O/analysis_set.rds')\$tissue), sep = '\n')" | xargs -P 4 -I{} Rscript v3m/05_coloc_tissue_v3m.R {} 4
# B5. final signal table and gene x tissue classes; per-signal classes (Fig. 4E)
Rscript v3m/06_final_v3m.R
Rscript v3m/eqtl_signal_classes.R $O/final_clump_v3m/coloc_eqtl_signals.tsv.gz $O/final_clump_v3m/signal_classes
# B6. the names read downstream (GWAS step 12a, figure code)
ln -sfn coloc_eqtl_v3m $R5/results/coloc_eqtl_v3c
ln -sfn final_clump_v3m $O/final_clump_v3c
