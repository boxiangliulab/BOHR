#!/bin/bash
# Authored summary (not an original script): the commands of this step in the order they were run, with the
# arguments and environment of the original jobs. Paths are those of the original runs (lab cluster and NSCC); set them
# for your system before running anything. Loops over tissues, arms and chromosomes ran as PBS array or feeder jobs.
# Final run: roots rerun_v2 (first promoter set), rerun_v4 (final promoter set, units of four arms) and rerun_r5e80 (final
# units onwards; it links to the scans, S2 tables and inputs of rerun_v4). The scripts' own paths name the earlier run;
# the copies that ran differed only in the folder prefixes and the run root (README, 'Paths and roots').
set -euo pipefail
NF=/projects/group/bliu/e1327989/GTEx/nofloor_full          # lab tree of the final run
N=/scratch/users/nus/e1327989/nofloor_full                   # NSCC tree of the final run
STORED=/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2 # earlier run (permutation maps, ID tables)
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
TISSUES=$(cat ../07_coloc_gwas_alphagenome/config/tissues_50.txt)   # the 50 tissues
ARMS="real perm0 perm1 perm2"

# 1. inputs
#    PEER factors per tissue (NSCC array of 57 jobs, 32 cores; the 7 tissues outside the analysis were stopped) on the
#    PEER input of ../01_promoter_quantification/export_txt.R. The factor count per tissue (K0) is the chromosome 22 scan
#    of the earlier run (inputs/config/tissue_list_best_peer_table.txt; inputs/peer/gamlss_QTL_peer_test.R, peer_plot.R).
for T in $(cut -f2 $N/GTEx/peer/v11/filter_5_recal_att_filter_SMTSD_mean.txt | tail -n +2 | sort -u); do
  Rscript inputs/peer/peer_analysis.R "$T" $N/GTEx/peer/v11/filter_5_recal_att_filter_SMTSD_mean.txt $N/GTEx/peer/v11; done
#    first promoter sets (skip SD <= 0.05) and input objects, lab; two batches: 49 tissues, then Whole_Blood, whose PEER
#    ran longer. The earlier run's permutation maps are copied.
export BOHR2_CODE=$NF/code BOHR2_GTEX=$NF/GTEx BOHR2_ROOT=$NF/GTEx/QTL_figure/rerun_v2
TIS49=$(grep -vx Whole_Blood ../07_coloc_gwas_alphagenome/config/tissues_50.txt | paste -sd,)
NF_TISSUES=$TIS49 NF_PARTIAL=1 Rscript inputs/promoter_sets_rule.R $NF/GTEx $NF/GTEx/v11/filter_5_recal_att_filter.RData $NF/GTEx/peer/v11 $NF/promoter_sets_nofloor_a.tsv
NF_TISSUES=$TIS49 NF_PARTIAL=1 Rscript inputs/prep_freeze_nofloor.R $NF/promoter_sets_nofloor_a.tsv $STORED
NF_TISSUES=Whole_Blood NF_PARTIAL=1 Rscript inputs/promoter_sets_rule.R $NF/GTEx $NF/GTEx/v11/filter_5_recal_att_filter.RData $NF/GTEx/peer/v11 $NF/promoter_sets_nofloor_wb.tsv
NF_TISSUES=Whole_Blood NF_PARTIAL=1 NF_APPEND=1 Rscript inputs/prep_freeze_nofloor.R $NF/promoter_sets_nofloor_wb.tsv $STORED
#    final promoter set (no SD <= 0.05 skip), its input objects, and the delta to scan (lab)
NF_SD_SKIP=none Rscript inputs/promoter_sets_rule.R $NF/GTEx $NF/GTEx/v11/filter_5_recal_att_filter.RData $NF/GTEx/peer/v11 $NF/promoter_sets_nofloor_v4.tsv
BOHR2_ROOT=$NF/GTEx/QTL_figure/rerun_v4 Rscript inputs/prep_freeze_nofloor.R $NF/promoter_sets_nofloor_v4.tsv $STORED
Rscript inputs/build_v4.R                                     # -> rerun_v4_delta (delta inputs), core/v4/delta_promoters.tsv
Rscript units/v4_strata_sizes.R                               # sizes of the admission strata, before any permuted delta scan

# 2. per-SNP scan (NSCC): every tissue, arm and chromosome of the first set (root rerun_v2, engine of the earlier SD rule;
#    feeder jobs of 64-128 cores, chromosome by chromosome), then the delta (root rerun_v4_delta, engine code of rerun_v4:
#    ../02_bohr_model_and_simulation/engine/bohr2_data.R with the SD < 0.005 skip), then one merge per arm
export BOHR2_VCF=/home/users/nus/e1327989/YuhaoDong/GTEx/v9_vcf_MAF_0.05 BOHR2_GTEX=/home/users/nus/e1327989/YuhaoDong/GTEx
#    (each job ran the copies in $BOHR2_CODE, with the run root as working directory)
export BOHR2_ROOT=$N/rerun_v2 BOHR2_CODE=$N/rerun_v2/code
for T in $TISSUES; do for A in $ARMS; do for c in $(seq 1 22); do
  Rscript scan/run_final2.R freeze "$T" "$A" chr$c 64
  Rscript scan/run_firth_addon.R freeze "$T" "$A" chr$c 64
done; done; done
export BOHR2_ROOT=$N/rerun_v4_delta BOHR2_CODE=$N/rerun_v4/code
for A in $ARMS; do for T in $TISSUES; do for c in $(seq 1 22); do
  Rscript scan/run_final2.R freeze "$T" "$A" chr$c 64
  Rscript scan/run_firth_addon.R freeze "$T" "$A" chr$c 64
done; done; done
for A in $ARMS; do Rscript scan/merge_v4.R "$A" 48; done          # -> rerun_v4/results/{final,final_firth_addon}
Rscript units/diag_v4_fit.R 32                                    # fit diagnostics of the delta units, real arm (rule 1)

# 3. units (NSCC, root rerun_v4, one job of 128 cores; ROOT = CODE root, BOHR2_GTEX=$N/GTEx)
export BOHR2_ROOT=$N/rerun_v4 BOHR2_CODE=$N/rerun_v4/code BOHR2_GTEX=$N/GTEx; cd $N/rerun_v4
Rscript scan/run_muF.R 128 $(for T in $TISSUES; do for A in $ARMS; do printf 'freeze/%s/%s ' $T $A; done; done)
Rscript units/mirror_pairs.R
for A in $ARMS; do NC=128 Rscript units/discovery_nogate.R freeze "$A" $TISSUES; done
Rscript units/units_v3_estimable.R
Rscript units/units_v3m_mirror.R            # first pass: every delta unit admitted
#    admission of the delta units (rules fixed before any permuted delta result): on the lab, from the units of the four
#    arms; then the decision table and the delta list go to rerun_v4/work/v4_admission/ and units_v3m_mirror.R runs again
Rscript units/diag_v4_calibration.R $NF/core/v4/units_v3m      # -> core/v4/calibration_decision.tsv (mu, p0 admitted; p1 not)
cp $NF/core/v4/delta_promoter_features.tsv $NF/core/v4/calibration_decision.tsv $N/rerun_v4/work/v4_admission/
Rscript units/units_v3m_mirror.R            # final units v3m of the final promoter set (delta p1 units not tested)
#    final rule set (r5e80): gene read depth (NSCC; sample-level counts stay private) -> median gene reads per tissue and
#    promoter (lab) -> promoters each rule excludes -> units v3m of the final rule set, BH again within tissue (lab)
Rscript units/depth_extract.R 8             # NSCC, run from $N; output $N/depth/<tissue>.rds, fetched to the lab (core/depth/private)
Rscript units/compare_units.R 12            # writes promoter_depth_v4.tsv.gz (med_D); the rest compares the two read rules
Rscript units/make_version_inputs.R         # -> figures_v4/inputs/nofloor_reads5_E80/promoters_excluded.tsv.gz
Rscript units/build_r5e80_units.R           # -> rerun_r5e80/results/final_nogate_v3m (copied to the NSCC root, md5 checked)

# 4. independent signals and per-signal tables (NSCC, root rerun_r5e80)
export BOHR2_ROOT=$N/rerun_r5e80 BOHR2_CODE=$N/rerun_r5e80/code; cd $N/rerun_r5e80
Rscript conditional/cond_fresh_list.R       # every significant unit from its leader; stubs; work/cond_v3m/jobs.tsv
#    every (tissue, chr) of work/cond_v3m/jobs.tsv: an array of 16 jobs (128 cores), 16 runs at a time with 8 cores each
tail -n +2 work/cond_v3m/jobs.tsv | cut -f1,2 | RS=Rscript CODE=$BOHR2_CODE NC=8 LOGDIR=logs xargs -P 16 -L 1 conditional/run_cond_one.sh
Rscript conditional/merge_conditional_v3m.R
LABELS=none Rscript conditional/merge_conditional_v3c.R
COND_FINAL=results/conditional_final_v3c LABEL_DIR=results/mirror_abs_labels_v3c NC=64 Rscript units/mirror_abs_labels.R
LABELS=results/mirror_abs_labels_v3c/labels_real.tsv.gz UNITS_OUT=results/final_nogate_v3c Rscript units/units_v3b_from_labels.R
LABELS=results/mirror_abs_labels_v3c/labels_real.tsv.gz Rscript conditional/merge_conditional_v3c.R
#    per-signal tables (delivery v6, read by ../05 and ../06); ID tables of the earlier run are copied first
cp $STORED/work/crosstype_pair_ids_v1.tsv work/; cp $STORED/results/fig5_delivery_v3/crosstype_conditional_pairs.tsv.gz results/fig5_delivery_v3/
export COND_FINAL=results/conditional_final_v3c UNITS_DIR=results/final_nogate_v3c FIG5_OUT=results/fig5_delivery_v6 NC=128
Rscript conditional/fig5_delivery_v4.R
Rscript conditional/fig5_proxies.R
Rscript conditional/fig5_pertype.R
Rscript conditional/fig5_pe_moments_v6.R
Rscript conditional/crosstype_conditional_v6.R
Rscript conditional/crosstype_conditional_v6x.R
FLAGS=results/mirror_v3/promoter_flags.tsv.gz Rscript conditional/fig5_signals_partner.R
# 5. LD standard: donor lists and signal LD (NSCC, same job)
Rscript ld/ld_v3.R
#    the tables were then copied to the lab root rerun_r5e80; the lead-signal flag (delivery v7, lab)
Rscript conditional/fig5_delivery_v7.R

# 6. QTLtools comparison (Fig. 2E, 2F; Supplementary Tables 2d, 2e), on the lab cluster, root of the final run.
#    The perm1 and perm2 maps were drawn by qtltools_comparison/make_perm_maps_freeze_extra.R in the earlier run and reused.
export BOHR2_ROOT=$NF/GTEx/QTL_figure/rerun_v4 BOHR2_CODE=$NF/GTEx/QTL_figure/rerun_v4/code BOHR2_GTEX=$NF/GTEx BOHR2_VCF=$NF/GTEx/Rawdata/v9_vcf_MAF_0.05
#    6a. QTLtools inputs, QTLtools v1.2 on all variants and all donors, and the comparison units
for T in $TISSUES; do ARMS='real perm0 perm1 perm2' Rscript qtltools_comparison/make_freeze_beds.R "$T"; done
#    one job per tissue and arm group (real and perm0 first, then perm1 and perm2), 16 runs at a time; for each arm, mode
#    (raw and normal read the raw BED, rint the rint BED; normal adds --normal) and chromosome:
#    QTLtools cis --vcf v9_vcf_MAF_0.05/chr$c.vcf.gz --bed inputs/$T.$arm.$kind.bed.gz --cov inputs/$T.$arm.cov.txt \
#      --permute 1000 --window 500000 --seed 123 [--normal] --region chr$c --out results/$mode/$T/prm.$arm.chr$c.txt --silent
FIG2E_TIS=ALL Rscript qtltools_comparison/fig2e_deliver_freeze_v3m_pairs.R
#    6b. the units of Fig. 2E under the final promoter rules (reads promoter_depth_v4.tsv.gz of units/compare_units.R)
Rscript qtltools_comparison/fig2e_r5e80_units.R 14
#    6c. QTLtools on the promoter-variant pairs that BOHR tested. Version argument v4 = the final run. One lab job per tissue
#        and arm (14-16 cores; EXACT_TMP = a node-local folder; OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1),
#        in this order, because each run copies the rows of the earlier ones where the inputs are the same.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript qtltools_comparison/make_site_lists.R 14
for T in $TISSUES; do for ARM in real perm0; do Rscript qtltools_comparison/exact_qtltools.R v4 "$T" $ARM 16; done; done
#    perm1 and perm2 of the exact run were run for 31 tissues; for the other tissues pairs_qtltools.R ran every promoter of
#    those arms itself
EXACT31="Adrenal_Gland Artery_Coronary Bladder Brain_Amygdala Brain_Anterior_cingulate_cortex_BA24 Brain_Caudate_basal_ganglia
  Brain_Cerebellar_Hemisphere Brain_Cerebellum Brain_Cortex Brain_Frontal_Cortex_BA9 Brain_Hippocampus Brain_Hypothalamus
  Brain_Nucleus_accumbens_basal_ganglia Brain_Putamen_basal_ganglia Brain_Spinal_cord_cervical_c-1 Brain_Substantia_nigra
  Cells_EBV-transformed_lymphocytes Esophagus_Gastroesophageal_Junction Kidney_Cortex Liver Minor_Salivary_Gland Ovary Pancreas
  Pituitary Prostate Small_Intestine_Terminal_Ileum Spleen Stomach Testis Uterus Vagina"
for T in $EXACT31; do for ARM in perm1 perm2; do Rscript qtltools_comparison/exact_qtltools.R v4 "$T" $ARM 16; done; done
for T in $TISSUES; do Rscript qtltools_comparison/keep_qtltools.R v4 "$T" real 14; done
for T in $TISSUES; do for ARM in real perm0 perm1 perm2; do Rscript qtltools_comparison/pairs_qtltools.R v4 "$T" $ARM 14; done; done
for T in $TISSUES; do for ARM in real perm0 perm1 perm2; do Rscript qtltools_comparison/fill_qtltools.R v4 "$T" $ARM 14; done; done
#    6d. BOHR's calls, the classes of Fig. 2E and the permuted calls of Fig. 2F
UNIT_MODE=merge_all Rscript qtltools_comparison/class_calib.R v4 6
R_DATATABLE_NUM_THREADS=1 Rscript qtltools_comparison/dd2_arms_units.R 3
PART=sidak R_DATATABLE_NUM_THREADS=1 Rscript qtltools_comparison/dd11_fisher_rules.R 3
PART=sidak R_DATATABLE_NUM_THREADS=1 Rscript qtltools_comparison/prepare_panel_e_classes.R
