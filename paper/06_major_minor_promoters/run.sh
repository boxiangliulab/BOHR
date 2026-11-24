#!/bin/bash
# Authored summary (not an original script): the commands of this step in the order they were run in the final run,
# with the arguments and environment of the original jobs. The scripts carry the paths of the earlier runs; the folders
# below are those of the final run (README.md, "Paths"). Set them for your system before running anything. Each stage
# script ran as its own lab PBS job (8 cores).
set -euo pipefail
cd 06_major_minor_promoters
HERE=$PWD
G=/projects/group/bliu/e1327989/GTEx/nofloor_full/GTEx                          # /projects/group/bliu/e1327989/GTEx in the scripts
R5=$G/QTL_figure/rerun_r5e80                                                    # QTL_figure/rerun_v2 in the scripts
B=$G/TATA_Enhancer/Figure5_rerun_r5e80/b1_tss_elements_v3c_20261113/data        # outputs of 187a and 187b
E=$R5/results/coloc_eqtl_v3_fig5                                                # outputs of prep/
FC=$R5/results/coloc_eqtl_v3m/final_clump_v3m                                   # final eQTL coloc table (../05_coloc_eqtl)
OLD=/home/e1327989/f5_stage_20261103/in                                         # earlier stage inputs that do not depend on promoter activity
S6=/projects/group/bliu/e1327989/GTEx/nofloor_full/f5_stage_nofloor_r5e80_v6    # stage folders (F5S)
S6C=/projects/group/bliu/e1327989/GTEx/nofloor_full/f5_stage_nofloor_r5e80_v6c
export UNITSET=v3c

# 1. TSS elements and ENCODE Z-scores of every promoter tested in units v3c (lab job, 4 cores; ran next to ../05_coloc_eqtl)
mkdir -p $B
Rscript encode/187a_alltested_tss_elements_v3c.R && python3 encode/187b_extract_encode_zscores_alltested_v3c.py && Rscript encode/187c_check_zscores_v3c.R

# 2. after ../05_coloc_eqtl: GTEx v11 eQTLs within +-200 kb of every analysis signal, read by 72 and 124* (lab job, 8 cores)
mkdir -p $E/eqtl_subset
Rscript prep/00_build_analysis_set.R
tail -n +2 $E/eqtl_windows.tsv | cut -f1 | sort -u | xargs -P 8 -I{} python3 prep/01_extract_eqtl.py {}

# 3. stage folders (no script)
mkdir -p $S6/in/v11 $S6/in/freeze $S6/in/mirror $S6/in/coloc $S6/in/coloc_final_clump $S6/in/fig5_v3data $S6/code $S6/out $S6/logs
for f in Promoter_info_all.RData Promoter_info_filter_5.txt tabs_att_filter_In.rds; do ln -sfn $G/v11/$f $S6/in/v11/$f; done
for f in $R5/work/final/freeze/*.rds; do ln -sfn $f $S6/in/freeze/; done
ln -sfn $R5/work/final/tissues.tsv $S6/in/freeze/tissues.tsv                     # = config/tissues.tsv
ln -sfn $R5/results/final_nogate_v3c/freeze $S6/in/freeze_v3c
ln -sfn $R5/results/fig5_delivery_v6 $S6/in/delivery
ln -sfn $R5/results/mirror/mirror_pairs.tsv.gz $S6/in/mirror/mirror_pairs.tsv.gz
ln -sfn $R5/results/mirror_v3/promoter_flags.tsv.gz $S6/in/mirror/promoter_flags.tsv.gz
ln -sfn $E/eqtl_subset $S6/in/coloc/eqtl_subset; ln -sfn $E/eqtl_windows.tsv $S6/in/coloc/eqtl_windows.tsv
ln -sfn $FC/coloc_eqtl_signals.tsv.gz $S6/in/coloc/coloc_eqtl_signals.tsv.gz
ln -sfn $FC/coloc_eqtl_signals.tsv.gz $S6/in/coloc_final_clump/coloc_eqtl_signals.tsv.gz
ln -sfn $B/alltested_tss_elements_v3c.tsv $S6/in/fig5_v3data/alltested_tss_elements.tsv
for m in h3k4me3 h3k27ac dnase; do
  ln -sfn $B/zscores_${m}_alltested_v3c.tsv.gz $S6/in/fig5_v3data/zscores_${m}_alltested.tsv.gz
  ln -sfn $B/zscores_${m}_alltested_v3c_columns.tsv $S6/in/fig5_v3data/zscores_${m}_alltested_columns.tsv
done
ln -sfn /home/e1327989/f5_stage_20261104_v4/in/fig5_v3data/gencode.v47.annotation.gtf.gz $S6/in/fig5_v3data/gencode.v47.annotation.gtf.gz
for d in gtf vcf encode external jaspar tata tpm snp_level fig5_old alphagenome; do ln -sfn $OLD/$d $S6/in/$d; done
cp stage_v6/theme_nature.R stage_v6/122a_tissue_matched_chromatin_v4.R stage_v6/00_helpers.R stage_v6/134_major_minor_pooled_chromatin.R $S6/code/
mkdir -p $S6C/code $S6C/out $S6C/logs; ln -sfn $S6/in $S6C/in
cp stage_v6/theme_nature.R stage_v6c/122a_tissue_matched_chromatin_v4.R stage_v6c/00_helpers.R $S6C/code/
for d in A1_location A2_eqtl_sqtl B1v2_multipromoter I_continuous_effects M_genotype_usage PROMOTER_LEVEL_v4; do ln -sfn $S6/out/$d $S6C/out/$d; done

# 4. two-promoter genes, chromatin and example loci (stage v6); the run also ran 122c and 129 after 122a, side analyses
#    whose outputs the paper does not use (not included)
export F5S=$S6
for s in 71_A1_location 72_A2_eqtl_direction 73b_B1v2_partner_conditional 74b_B1v2_patterns 97_I_continuous_effects \
         101_K_overlap_multi_location 104_M_genotype_usage 122a_tissue_matched_chromatin_v4 \
         132_marks_by_usage_role 134_major_minor_pooled_chromatin 184_fig5a_major_minor_pvalues \
         124_example_candidates_eqtl 124b_example_i_minor_on 124c_example_j_strong; do
  (cd $S6 && Rscript $HERE/stage_v6/$s.R)
done

# 5. the pattern table with the paper's rule for genes with three or more promoters (stage v6c)
export F5S=$S6C
(cd $S6C && Rscript $HERE/stage_v6c/101_K_overlap_multi_location.R)
