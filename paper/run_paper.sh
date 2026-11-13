#!/bin/bash
# Whole-paper test run (authored for the repository; not part of the original runs).
#
#   bash paper/run_paper.sh paper/config_lab.sh
#
# Runs the analysis of the paper, step by step (folders 01-07 of paper/), with the as-run scripts, and compares every
# result with the stored results and the published numbers. Each step is split into numbered sub-steps of two kinds:
#   rerun   an as-run script is run again on a small subset (two tissues, one chromosome, two donor characteristics, one
#           simulated sample size) and its output is compared with the stored output of the original run;
#   verify  where a re-run is impossible (inputs no longer stored, a remote model) or too large (all 50 tissues), the
#           step's numbers are recomputed from the stored outputs and compared with the published ones.
# For every sub-step the run prints what it does and where its log is, writes the log under OUT/logs/, and prints one
# line per check: PASS, FAIL or SKIP with the reason; after a FAIL it names the log to read. OUT/summary.tsv collects all
# checks, OUT/check.txt is the printed report. Step P (paper/with_package/) runs the same subset with the BOHR package.
#
# Settings come from the config file and can be overridden on the command line, e.g.
#   STEPS="03 05" KEEP_GOING=1 OUT=/scratch/$USER/bohr_paper bash paper/run_paper.sh paper/config_lab.sh
# The script submits itself as one PBS job unless SUBMIT=0 or it already runs inside a job; heavy work never runs on
# the login node. The inputs are individual-level GTEx data: do not copy OUT elsewhere, and delete it when done.
set -o pipefail
CONFIG=$(readlink -f "${1:-${CONFIG:?usage: run_paper.sh <config>}}")
source "$CONFIG"
PAPER=${PAPER_DIR:-$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)}   # PBS runs a copy of this script
VARS="OUT STEPS ONLY KEEP_GOING RESUME TISSUES CHR TRAITS SIM_N SIM_CHUNK NCORES RSCRIPT PYTHON SIM_RLIB BOHR_SRC GTEX VCF RUN_ROOT
      SIM_ROOT TRG_ROOT GWAS_ROOT F5S_V6 F5S_V6C F5_ENCODE AG_ROOT SUPP"

if [ -z "$PBS_JOBID" ] && [ "$SUBMIT" = 1 ] && command -v qsub >/dev/null; then
  mkdir -p "$OUT"
  for v in $VARS; do printf '%s=%q\n' "$v" "${!v}"; done > "$OUT/config_resolved.sh"   # the job reads these values
  extra=(); [ -n "$PBS_PROJECT" ] && extra=(-P "$PBS_PROJECT")
  qsub -N bohr_paper -q "$PBS_QUEUE" -l "$PBS_RESOURCES" -l walltime="$PBS_WALLTIME" "${extra[@]}" -j oe \
       -o "$OUT/run_paper.log" -v CONFIG="$OUT/config_resolved.sh",SUBMIT=0,PAPER_DIR="$PAPER" "$PAPER/run_paper.sh"
  echo "submitted; progress: $OUT/check.txt; log: $OUT/run_paper.log; summary: $OUT/summary.tsv"
  exit 0
fi

[ -f "$PAPER/02_bohr_model_and_simulation/engine/bohr2_engine.R" ] || { echo "paper/ not found at $PAPER"; exit 1; }
mkdir -p "$OUT"/{logs,status,done} || exit 1
cd "$OUT" || exit 1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
RS=$RSCRIPT; CK=$PAPER/checks; LIB=$OUT/lib
[ -z "$BOHR_SRC" ] && BOHR_SRC=$(dirname "$PAPER")
read -r -a TIS <<< "$TISSUES"; read -r -a TRT <<< "$TRAITS"
SUMMARY=$OUT/summary.tsv; REPORT=$OUT/check.txt
[ -f "$SUMMARY" ] || printf 'substep\tlevel\tcheck\tresult\tdetail\tminutes\tlog\n' > "$SUMMARY"
say() { echo "$@" | tee -a "$REPORT"; }
say "== run_paper.sh  host $(hostname)  start $(date '+%F %T')  steps: $STEPS${ONLY:+  sub-steps: $ONLY}"
say "   tissues: $TISSUES  chromosome: $CHR  characteristics: $TRAITS  simulated N: $SIM_N  cores: $NCORES"
N_FAIL=0

## status rows written directly by the driver (a SKIP with its reason, or a FAIL when a step stopped with an error)
row() { printf '%s\t%s\t%s\n' "$2" "$3" "$4" >> "$OUT/status/$1.tsv"; }
want() { local id=$1; [[ " $STEPS " == *" ${id%%.*} "* ]] || return 1; [ -z "$ONLY" ] || [[ " $ONLY " == *" $id "* ]]; }

## run one sub-step: substep <id> <level> <description> <function>
substep() {
  local id=$1 level=$2 desc=$3 fn=$4 log=$OUT/logs/$1.log st=$OUT/status/$1.tsv
  want "$id" || return 0
  if [ "$RESUME" = 1 ] && [ -f "$OUT/done/$id" ] && [ -f "$st" ]; then
    say "[$id] $desc — passed in an earlier run (RESUME=1); log: $log"
  else
    say "[$id] $desc — log: $log"
    rm -f "$st" "$OUT/done/$id"; printf 'check\tresult\tdetail\n' > "$st"
    awk -F'\t' -v id="$id" '$1 != id' "$SUMMARY" > "$SUMMARY.tmp" && mv "$SUMMARY.tmp" "$SUMMARY"   # rows of an earlier attempt
    local t0=$SECONDS rc
    ( set -e; echo "[$id] $desc"; echo "start $(date '+%F %T') on $(hostname)"; "$fn" ) > "$log" 2>&1; rc=$?
    echo "end $(date '+%F %T'); exit status $rc" >> "$log"
    local min; min=$(awk -v s=$((SECONDS - t0)) 'BEGIN { printf "%.1f", s / 60 }')
    [ "$(tail -n +2 "$st" | wc -l)" -eq 0 ] && [ $rc -eq 0 ] && row "$id" "$desc" FAIL "the sub-step wrote no check"
    [ $rc -ne 0 ] && row "$id" "$desc: run" FAIL "stopped with an error (exit status $rc)"
    tail -n +2 "$st" | awk -F'\t' -v id="$id" -v lv="$level" -v m="$min" -v lg="logs/$id.log" 'BEGIN{OFS="\t"} {print id, lv, $1, $2, $3, m, lg}' >> "$SUMMARY"
    grep -qP '\tFAIL\t' "$st" || touch "$OUT/done/$id"
    say "      ($min min)"
  fi
  local nf; nf=$(tail -n +2 "$st" | awk -F'\t' '$2 == "FAIL"' | wc -l)
  tail -n +2 "$st" | awk -F'\t' '{ printf "      %-4s  %s  %s\n", $2, $1, $3 }' | tee -a "$REPORT"
  if [ "$nf" -gt 0 ]; then
    N_FAIL=$((N_FAIL + nf)); say "      -> read $log"
    [ "$KEEP_GOING" = 1 ] || { finish_run; exit 1; }
  fi
}
finish_run() {
  say ""
  say "$(awk -F'\t' 'NR > 1 { n[$4]++ } END { printf "totals: %d PASS, %d FAIL, %d SKIP", n["PASS"], n["FAIL"], n["SKIP"] }' "$SUMMARY")  (summary: $SUMMARY)"
  say "end $(date '+%F %T')"
}
## R and Python with the package library first on the path
R_() { R_LIBS=$LIB${R_LIBS:+:$R_LIBS} "$RS" "$@"; }
## a copy of an as-run script with fixed output folders pointed into OUT (checks/redirect.py lists the changed lines)
redirect() { "$PYTHON" "$CK/redirect.py" "$@"; }

#### 00 setup ##########################################################################################################
s00_1() {
  mkdir -p "$LIB"
  "$(dirname "$RS")/R" CMD INSTALL --no-test-load --library="$LIB" "$BOHR_SRC"
  R_LIBS=$LIB:$SIM_RLIB "$RS" "$CK/check_setup.R" "$LIB" "$OUT/status/00.1.tsv"
}
substep 00.1 setup "install the BOHR package of this repository into OUT/lib and check the software" s00_1

#### 01 promoter quantification ########################################################################################
s01_1() {
  row 01.1 "rerun of 0101-0102 (junction files, proActiv per tissue)" SKIP "needs the per-sample GTEx junction files and 57 proActiv runs (NSCC); the per-tissue proActiv results are no longer stored, so 01.2 repeats the filter from the stored outputs of 0103"
  "$RS" "$CK/verify_01.R" filter "$GTEX/v11" "$SUPP" "$OUT/status/01.1.tsv"
}
s01_2() {
  R_ "$PAPER/test_subset/step1_filter.R" "$GTEX/v11" "$CHR" step1
  CHECK_PARTS=1 CHECK_OUT=$OUT/status/01.2.tsv GTEX=$GTEX "$RS" "$PAPER/test_subset/check_subset.R" "$OUT" "$RUN_ROOT" "$CHR" "${TIS[@]}"
}
s01_3() { "$RS" "$CK/verify_01.R" release "$GTEX" "$SUPP" "$OUT/status/01.3.tsv"; }
substep 01.1 verify "promoter annotation and filter of SI §1: stored outputs of 0103 against the published counts" s01_1
substep 01.2 rerun  "promoter filter on $CHR with BOHR::relative_activity() on the stored proActiv outputs (as 0103)" s01_2
substep 01.3 verify "Fig. 1A: single- and multi-promoter genes per GTEx release from the stored annotations (Supplementary Table 1a)" s01_3

#### 02 BOHR model and simulations ######################################################################################
SIM=$OUT/sim
s02_1() {
  row 02.1 "rerun of the gamlss check (equiv_test.R)" SKIP "419 gamlss refits on the cached pilot inputs; the stored outputs of the check are verified instead"
  "$RS" "$CK/verify_02.R" gamlss "$RUN_ROOT" "$SUPP" "$OUT/status/02.1.tsv"
}
s02_2() {
  mkdir -p "$SIM/gen" && cd "$SIM/gen"
  TRUNC_L=0.005 GEN_ONLY=1 R_LIBS=$SIM_RLIB "$RS" "$PAPER/02_bohr_model_and_simulation/simulated_activity/newsim_direct_nothreshold_indep.R" "$SIM_N" 2000 2
  "$RS" "$CK/check_sim.R" gen "$SIM/gen/sim_N${SIM_N}_G2000_realfloor" "$SIM_ROOT/newmethod_bohr2/m1_nothr_g2000/sim_N${SIM_N}_G2000_realfloor" "$OUT/status/02.2.tsv"
}
s02_3() {
  mkdir -p "$SIM/fit"
  local E=$PAPER/02_bohr_model_and_simulation/engine
  ## the fit sources the engine copy of its run folder; the copy here sources the engine of paper/ (identical md5)
  redirect "$PAPER/02_bohr_model_and_simulation/benchmark/bohr2_sim.R" "$SIM/fit/bohr2_sim.R" \
    'file.path(NM, "code_copy", "bohr2_engine.R")' "\"$E/bohr2_engine.R\"" 'file.path(NM, "code_copy", "bohr2_addon3.R")' "\"$E/bohr2_addon3.R\""
  "$RS" "$SIM/fit/bohr2_sim.R" power M1 nothr2000 "$SIM_N" "$SIM/fit/power_N${SIM_N}.rds" "$NCORES" 1000 $SIM_CHUNK
  "$RS" "$CK/check_sim.R" fit "$SIM/fit/power_N${SIM_N}.rds" "$SIM_ROOT/newmethod_bohr2/results/m1nothr_g2000/power/N${SIM_N}.rds" "$OUT/status/02.3.tsv"
}
s02_4() {
  row 02.4 "rerun of the full benchmark (12 sample sizes, RNA-seq simulation, 1,000 permutations)" SKIP "12 sample sizes in two frameworks with 1,000 permutations each, run as PBS arrays; 02.2-02.3 re-run one sample size, and the stored summaries are checked here"
  "$RS" "$CK/verify_02.R" summaries "$SIM_ROOT/figures_bohr2/data_v3_nothr_realjunc" "$SUPP" "$OUT/status/02.4.tsv"
}
substep 02.1 verify "BOHR engine against gamlss (Supplementary Methods §2): stored summary of the check" s02_1
substep 02.2 rerun  "simulated activity, N = $SIM_N, 2,000 genes (newsim_direct_nothreshold_indep.R)" s02_2
substep 02.3 rerun  "BOHR fits of simulated promoters, N = $SIM_N, chunk ${SIM_CHUNK/ / of } (benchmark/bohr2_sim.R)" s02_3
substep 02.4 verify "Fig. 2 and Supplementary Table 2: power and calibrated thresholds from the stored benchmark outputs" s02_4

#### 03 cis mapping and independent signals ############################################################################
## code folder (the scripts source each other from one folder) and a shadow run root: stored inputs linked, new results here
setup_03() {
  mkdir -p code shadow/results/final shadow/results/conditional_v3m
  for f in "$PAPER"/02_bohr_model_and_simulation/engine/bohr2_*.R "$PAPER"/03_cis_mapping_conditional/scan/run_final2.R \
           "$PAPER"/03_cis_mapping_conditional/conditional/run_conditional_v3m.R; do ln -sf "$f" code/; done
  ln -sfn "$RUN_ROOT/work" shadow/work; ln -sfn "$RUN_ROOT/perm_maps" shadow/perm_maps
  ln -sfn "$RUN_ROOT/results/final_nogate_v3m" shadow/results/final_nogate_v3m
}
export_03() { export BOHR2_CODE=$OUT/code BOHR2_ROOT=$OUT/shadow BOHR2_GTEX=$GTEX BOHR2_VCF=$VCF; }
s03_1() {
  setup_03; export_03
  for T in "${TIS[@]}"; do echo "== run_final2.R $T $CHR"; "$RS" code/run_final2.R freeze "$T" real "$CHR" "$NCORES"; done
  CHECK_PARTS=2 CHECK_OUT=$OUT/status/03.1.tsv GTEX=$GTEX "$RS" "$PAPER/test_subset/check_subset.R" "$OUT" "$RUN_ROOT" "$CHR" "${TIS[@]}"
}
s03_2() {
  setup_03; export_03
  for T in "${TIS[@]}"; do echo "== run_conditional_v3m.R $T $CHR"; "$RS" code/run_conditional_v3m.R freeze "$T" real "$CHR" "$NCORES"; done
  CHECK_PARTS=3 CHECK_OUT=$OUT/status/03.2.tsv GTEX=$GTEX "$RS" "$PAPER/test_subset/check_subset.R" "$OUT" "$RUN_ROOT" "$CHR" "${TIS[@]}"
}
s03_3() {
  row 03.3 "rerun of the scan and conditional analysis in all 50 tissues" SKIP "50 tissues x 22 chromosomes x 4 arms (real and permuted), run as PBS arrays on NSCC; 03.1-03.2 re-run the subset, and the stored tables of all tissues are checked here"
  "$RS" "$CK/verify_03.R" "$RUN_ROOT" "$SUPP" "$OUT/status/03.3.tsv"
}
substep 03.1 rerun  "per-SNP BOHR scan of every tested promoter on $CHR in each tissue (scan/run_final2.R)" s03_1
substep 03.2 rerun  "independent signals of the re-run units on $CHR (conditional/run_conditional_v3m.R)" s03_2
substep 03.3 verify "units, significant units and the 14,800 signals of all 50 tissues against Supplementary Tables 4a, 4d and 5d" s03_3

#### 04 tissue specificity and donor characteristics ###################################################################
s04_1() {
  row 04.1 "rerun of ROKU (tissue_specificity/ROKU.R)" SKIP "ran on NSCC with TCC::ROKU (not installed here); its stored output is the input of the effect filter re-run below"
  mkdir -p "$OUT/roku/out"
  redirect "$PAPER/04_tissue_specificity_traits/tissue_specificity/roku_effect_filter.R" "$OUT/roku/roku_effect_filter.R" \
    'OUT <- file.path(FIT, "roku_20261007")' "OUT <- \"$OUT/roku/out\""
  "$RS" "$OUT/roku/roku_effect_filter.R"
  "$RS" "$CK/check_04.R" roku "$OUT/roku/out" "$GTEX/fit/roku_20261007" "$OUT/status/04.1.tsv"
}
s04_2() {
  local R=$OUT/traits/run_20261005; mkdir -p "$R/code" "$OUT/traits/results"
  for f in data.R run_trait.R; do ln -sf "$PAPER/04_tissue_specificity_traits/traits/run_20261005/$f" "$R/code/"; done
  for f in bohr2_engine.R bohr2_addon3.R; do ln -sf "$PAPER/02_bohr_model_and_simulation/engine/$f" "$R/code/"; done
  ln -sfn "$TRG_ROOT/trait_regression_20261005/work" "$R/work"
  for T in "${TIS[@]}"; do
    echo "== run_trait.R $T real ${TRT[*]}"
    TRG_RUN=$R TRG_OUT=$OUT/traits/results "$RS" "$R/code/run_trait.R" "$T" real "$NCORES" "${TRT[@]}"
  done
  "$RS" "$CK/check_04.R" traits "$OUT/traits/results/real" "$TRG_ROOT/trait_regression_20261005/results/real" "$OUT/status/04.2.tsv" "${TIS[@]}"
}
s04_3() {
  row 04.3 "rerun of the 20 characteristics x 3 arms in all tissues" SKIP "20 characteristics x 3 arms (real, simulated, permuted) in every tissue, run as PBS arrays; 04.2 re-runs the subset, and the stored final tables are checked here"
  "$RS" "$CK/verify_04.R" "$TRG_ROOT" "$GTEX/fit/roku_20261007" "$SUPP" "$OUT/status/04.3.tsv"
}
substep 04.1 rerun  "tissue-specific promoters: the 0.10 effect filter on the stored ROKU calls (roku_effect_filter.R)" s04_1
substep 04.2 rerun  "donor characteristics ${TRT[*]} in each tissue, real arm (traits/run_20261005/run_trait.R)" s04_2
substep 04.3 verify "final characteristic tests and ROKU calls against Supplementary Tables 3b, 3g and 3h" s04_3

#### 05 colocalization with eQTLs #######################################################################################
s05_1() {
  local S=$OUT/coloc_eqtl/results; mkdir -p "$S/coloc_eqtl_ldwin/per_tissue"
  ln -sfn "$RUN_ROOT/results/coloc_eqtl_v3" "$S/coloc_eqtl_v3"
  ln -sf "$RUN_ROOT/results/coloc_eqtl_ldwin/pairs.tsv.gz" "$S/coloc_eqtl_ldwin/"
  ln -sfn "$RUN_ROOT/results/coloc_eqtl_ldwin/eqtl_subset" "$S/coloc_eqtl_ldwin/eqtl_subset"
  redirect "$PAPER/05_coloc_eqtl/ldwin/06_coloc_tissue.R" "$OUT/coloc_eqtl/06_coloc_tissue.R" \
    'source("/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/scripts/coloc_final_lib.R")' "source(\"$PAPER/05_coloc_eqtl/coloc_final_lib.R\")" \
    'V2L <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/results"' "V2L <- \"$S\""
  for T in "${TIS[@]}"; do echo "== 06_coloc_tissue.R $T"; "$RS" "$OUT/coloc_eqtl/06_coloc_tissue.R" "$T" "$NCORES"; done
  "$RS" "$CK/check_05_07.R" eqtl "$S/coloc_eqtl_ldwin/per_tissue" "$RUN_ROOT/results/coloc_eqtl_ldwin/per_tissue" "$OUT/status/05.1.tsv" "${TIS[@]}"
}
s05_2() {
  row 05.2 "rerun of the eQTL colocalization in all tissues" SKIP "50 tissues with clumped eQTL leads (NSCC); 05.1 re-runs the subset, and the stored final table is checked here"
  "$RS" "$CK/check_05_07.R" eqtl_final "$RUN_ROOT/results/coloc_eqtl_v3c/final_clump_v3c" "$SUPP" "$OUT/status/05.2.tsv"
}
substep 05.1 rerun  "puQTL x eQTL colocalization over the LD windows in each tissue (ldwin/06_coloc_tissue.R)" s05_1
substep 05.2 verify "Fig. 4E classes: 8,297 gene-tissue pairs, 1,600 distinct, 3,263 distinct or without eQTL (Supplementary Tables 5a, 5c)" s05_2

#### 06 major and minor promoters #######################################################################################
s06_1() {
  local W=$OUT/encode/W E=$OUT/encode/E; mkdir -p "$W/data" "$E"
  redirect "$PAPER/06_major_minor_promoters/encode/25a_alltested_tss_elements.R" "$OUT/encode/25a.R" \
    'W <- "/home/e1327989/project/GTEx/TATA_Enhancer/Figure5_rerun_v3"' "W <- \"$W\""
  "$RS" "$OUT/encode/25a.R"
  ## 22b2 reads the extracted SCREEN matrices and only checks that the downloaded archive exists (moved to NSCC); the
  ## shadow folder links the extracted matrices and holds empty stand-ins for the archives
  for m in dnase h3k27ac h3k4me3; do ln -sfn "$F5_ENCODE/../encode_v4_zscores/$m" "$E/$m"; : > "$E/human_${m}_core_collection_zscores.tar.gz"; done
  redirect "$PAPER/06_major_minor_promoters/encode/22b2_extract_encode_zscores_alltested.py" "$OUT/encode/22b2.py" \
    'E = "/home/e1327989/project/GTEx/TATA_Enhancer/encode_v4_zscores"' "E = \"$E\"" \
    'W = "/home/e1327989/project/GTEx/TATA_Enhancer/Figure5_rerun_v3"' "W = \"$W\""
  "$PYTHON" "$OUT/encode/22b2.py"
  "$RS" "$CK/check_06.R" encode "$W/data" "$F5S_V6/in/fig5_v3data" "$OUT/status/06.1.tsv"
}
F5_CHAIN="71_A1_location 72_A2_eqtl_direction 73b_B1v2_partner_conditional 74b_B1v2_patterns 97_I_continuous_effects
          101_K_overlap_multi_location 104_M_genotype_usage 122a_tissue_matched_chromatin_v4 134_major_minor_pooled_chromatin
          184_fig5a_major_minor_pvalues 132_marks_by_usage_role 124_example_candidates_eqtl 124b_example_i_minor_on 124c_example_j_strong"
s06_2() {
  local S=$OUT/f5_v6; mkdir -p "$S/code" "$S/out"
  ln -sfn "$F5S_V6/in" "$S/in"
  for f in "$PAPER"/06_major_minor_promoters/stage_v6/*.R; do ln -sf "$f" "$S/code/"; done
  ln -sf "$F5S_V6/code/theme_nature.R" "$S/code/"          # figure theme of the stage (plots are not compared)
  export F5S=$S FONTCONFIG_FILE=$F5S_V6/code/fonts_nohint.conf
  cd "$S"
  for s in $F5_CHAIN; do echo "== $s.R"; "$RS" "code/$s.R"; done
  "$RS" "$CK/check_06.R" stage "$S/out" "$F5S_V6/out" "$OUT/status/06.2.tsv"
}
s06_3() {
  local S=$OUT/f5_v6c; mkdir -p "$S/code" "$S/out/K_overlap_multi_location"
  ln -sfn "$F5S_V6C/in" "$S/in"
  for f in "$PAPER"/06_major_minor_promoters/stage_v6c/*.R; do ln -sf "$f" "$S/code/"; done
  for d in "$F5S_V6C"/out/*/; do d=$(basename "$d"); [ "$d" = K_overlap_multi_location ] || ln -sfn "$F5S_V6C/out/$d" "$S/out/$d"; done
  export F5S=$S; cd "$S"
  "$RS" code/101_K_overlap_multi_location.R
  "$RS" "$CK/check_06.R" stage "$S/out/K_overlap_multi_location" "$F5S_V6C/out/K_overlap_multi_location" "$OUT/status/06.3.tsv"
}
substep 06.1 rerun  "TSS cCREs and ENCODE Z-scores of all tested promoters (encode/25a, 22b2)" s06_1
substep 06.2 rerun  "two-promoter genes: routes, patterns, genotype usage, chromatin and example loci (stage_v6, 14 scripts)" s06_2
substep 06.3 rerun  "per-signal pattern table with the paper's rule for genes with three or more promoters (stage_v6c/101)" s06_3

#### 07 colocalization with GWAS traits; AlphaGenome at MRAS ############################################################
s07_1() {
  local N=$OUT/gwas/ndir L; L=$N/coloc_ldunion_20261012/run_v3c; mkdir -p "$L/per_tissue"
  for d in scripts gwas_eur_20261009 results_v3c; do ln -sfn "$GWAS_ROOT/$d" "$N/$d"; done
  for d in pairs gwas_subset; do ln -sfn "$GWAS_ROOT/coloc_ldunion_20261012/run_v3c/$d" "$L/$d"; done
  redirect "$PAPER/05_coloc_eqtl/v3c/coloc_final_lib.R" "$OUT/gwas/coloc_final_lib.R" \
    'NDIR <- "/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3"' "NDIR <- \"$N\""
  redirect "$PAPER/07_coloc_gwas_alphagenome/gwas_coloc/14_coloc_tissue.R" "$OUT/gwas/14_coloc_tissue.R" \
    'source("/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/scripts/coloc_final_lib.R")' "source(\"$OUT/gwas/coloc_final_lib.R\")"
  for T in "${TIS[@]}"; do echo "== 14_coloc_tissue.R $T"; UNITSET=v3c "$RS" "$OUT/gwas/14_coloc_tissue.R" "$T" "$NCORES"; done
  "$RS" "$CK/check_05_07.R" gwas "$L/per_tissue" "$GWAS_ROOT/coloc_ldunion_20261012/run_v3c/per_tissue" "$OUT/status/07.1.tsv" "${TIS[@]}"
}
s07_2() {
  row 07.2 "rerun of the GWAS colocalization in all tissues" SKIP "51 traits x 50 tissues with clumping on NSCC; 07.1 re-runs the subset, and the stored final table is checked here"
  "$RS" "$CK/check_05_07.R" gwas_final "$GWAS_ROOT/results_v3c" "$SUPP" "$OUT/status/07.2.tsv"
}
s07_3() {
  row 07.3 "AlphaGenome predictions (alphagenome_perbase_MRAS*.py)" SKIP "calls the remote AlphaGenome model (needs ALPHAGENOME_API_KEY and network access); the stored predictions are the input of the scoring re-run below"
  local A=$OUT/alphagenome; mkdir -p "$A/tables" "$A/figures" "$A/a3/tables"
  local TAB=$AG_ROOT/GTEx/alphagenome/MRAS_perbase_analysis/tables A3=$GWAS_ROOT/alphagenome_addon
  cp "$TAB"/*.tsv "$TAB"/*.txt "$A/tables/" 2>/dev/null || true      # copies, not links: the scripts write into this folder
  cp "$A3"/*.tsv "$A3"/*.txt "$A/a3/" 2>/dev/null || true
  sleep 1; touch "$A/.start"                                          # outputs are the files written after this point
  local G=/projects/group/bliu/e1327989/GTEx/alphagenome/MRAS_perbase_analysis
  redirect "$PAPER/07_coloc_gwas_alphagenome/alphagenome/mras_peak_finalize.py" "$A/mras_peak_finalize.py" "\"$G/tables\"" "\"$A/tables\"" "\"$G/figures\"" "\"$A/figures\""
  redirect "$PAPER/07_coloc_gwas_alphagenome/alphagenome/mras_score_variants.py" "$A/mras_score_variants.py" "\"$G/tables\"" "\"$A/tables\""
  redirect "$PAPER/07_coloc_gwas_alphagenome/alphagenome/score_v3_candidates.py" "$A/score_v3_candidates.py" "\"$G/tables\"" "\"$A/tables\"" \
    "\"/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/alphagenome_addon\"" "\"$A/a3\""
  cd "$A"
  for s in mras_peak_finalize mras_score_variants score_v3_candidates; do echo "== $s.py"; "$PYTHON" "$A/$s.py"; done
  "$RS" "$CK/check_05_07.R" alphagenome "$A" "$TAB" "$A3/tables" "$OUT/status/07.3.tsv"
}
substep 07.1 rerun  "puQTL x GWAS colocalization over the LD windows in each tissue (gwas_coloc/14_coloc_tissue.R)" s07_1
substep 07.2 verify "GWAS colocalization: 47,506 tests and 1,125 high-confidence (623 μ, 502 boundary) against Supplementary Table 8" s07_2
substep 07.3 rerun  "MRAS: CAGE and PRO-cap peaks and variant scores from the stored AlphaGenome predictions" s07_3

#### P the package route ################################################################################################
sP() { bash "$PAPER/with_package/run_package_route.sh" "$1"; }
sP_1() { sP cis; }; sP_2() { sP conditional; }; sP_3() { sP trait; }
export OUT TISSUES CHR TRAITS NCORES RSCRIPT LIB PAPER SUPP TRG_ROOT
substep P.1 package "BOHR::bohr_cis() on $CHR in each tissue, against the stored scan and units" sP_1
substep P.2 package "BOHR::bohr_conditional() chained after P.1, against the published signals" sP_2
substep P.3 package "BOHR::bohr_trait() for ${TRT[0]} in each tissue, against the published characteristic tests" sP_3

finish_run
[ "$N_FAIL" -eq 0 ]
