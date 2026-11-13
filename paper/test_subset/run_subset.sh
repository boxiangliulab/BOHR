#!/bin/bash
# Subset test of the analysis code (authored for the repository; not part of the original runs).
#
#   bash paper/test_subset/run_subset.sh paper/test_subset/config_lab.sh
#
# Re-runs the key real-data steps of the paper with the scripts in paper/, as they were run, on a few tissues and one
# chromosome, then compares every output with the stored results of the paper and prints PASS or FAIL per check:
#   1. promoter filter: BOHR::relative_activity() on the promoters of CHR (the package reproduces
#      0103_proActiv_data_process_tissue.R, whose 57 per-tissue proActiv inputs are no longer stored);
#   2. cis mapping: 03_cis_mapping_conditional/scan/run_final2.R for each tissue on CHR;
#   3. independent signals: 03_cis_mapping_conditional/conditional/run_conditional_v3m.R for the units of each tissue and
#      CHR that the final conditional run re-ran, on the scan of step 2;
#   4. the signals of step 3 against the paper's final signal table (14,800 signals).
# Outputs go to OUT only. The stored inputs are read through a shadow folder of symbolic links, so nothing stored is
# written. Takes about 30-60 minutes with 8 cores.
set -eo pipefail
CONFIG=$(readlink -f "${1:-${CONFIG:?usage: run_subset.sh <config>}}")
source "$CONFIG"
HERE=${TEST_DIR:-$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)}   # PBS runs a copy of this script
PAPER=$(dirname "$HERE")

if [ -z "$PBS_JOBID" ] && [ "$SUBMIT" = 1 ] && command -v qsub >/dev/null; then
  mkdir -p "$OUT"
  ## the job reads the values resolved here (including any overrides given on the command line)
  for v in OUT TISSUES CHR NCORES RSCRIPT BOHR_SRC GTEX VCF RUN_ROOT; do printf '%s=%q\n' "$v" "${!v}"; done > "$OUT/config_resolved.sh"
  extra=(); [ -n "$PBS_PROJECT" ] && extra=(-P "$PBS_PROJECT")
  qsub -N bohr_subset -q "$PBS_QUEUE" -l "$PBS_RESOURCES" -l walltime="$PBS_WALLTIME" "${extra[@]}" -j oe \
       -o "$OUT/run_subset.log" -v CONFIG="$OUT/config_resolved.sh",SUBMIT=0,TEST_DIR="$HERE" "$HERE/run_subset.sh"
  echo "submitted; log: $OUT/run_subset.log; summary: $OUT/summary.tsv"
  exit 0
fi

mkdir -p "$OUT"; cd "$OUT"
[ -f "$PAPER/02_bohr_model_and_simulation/engine/bohr2_engine.R" ] || { echo "paper/ not found at $PAPER"; exit 1; }
echo "host $(hostname)  start $(date)  tissues: $TISSUES  $CHR  cores: $NCORES"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1

## code folder: the scripts source each other from one folder (BOHR2_CODE)
mkdir -p code
for f in "$PAPER"/02_bohr_model_and_simulation/engine/bohr2_*.R "$PAPER"/03_cis_mapping_conditional/scan/run_final2.R \
         "$PAPER"/03_cis_mapping_conditional/conditional/run_conditional_v3m.R; do ln -sf "$f" code/; done
export BOHR2_CODE=$OUT/code

## shadow run root: stored inputs linked read-only, new results written here
mkdir -p shadow/results
ln -sfn "$RUN_ROOT/work" shadow/work
ln -sfn "$RUN_ROOT/perm_maps" shadow/perm_maps
ln -sfn "$RUN_ROOT/results/final_nogate_v3m" shadow/results/final_nogate_v3m
mkdir -p shadow/results/final shadow/results/conditional_v3m
export BOHR2_ROOT=$OUT/shadow BOHR2_GTEX=$GTEX BOHR2_VCF=$VCF

## step 1: the BOHR package and the promoter filter
if [ -z "$BOHR_SRC" ]; then BOHR_SRC=$(dirname "$PAPER"); fi
mkdir -p lib
"$(dirname "$RSCRIPT")/R" CMD INSTALL --no-test-load --library=lib "$BOHR_SRC" > step1_install.log 2>&1
echo "== step 1: promoter filter on $CHR $(date)"
R_LIBS=$OUT/lib "$RSCRIPT" "$HERE/step1_filter.R" "$GTEX/v11" "$CHR" step1 > step1_filter.log 2>&1

## step 2: cis mapping (as run)
for T in $TISSUES; do
  echo "== step 2: run_final2.R $T $CHR $(date)"
  "$RSCRIPT" code/run_final2.R freeze "$T" real "$CHR" "$NCORES" > "step2_${T}.log" 2>&1
done

## step 3: independent signals of the re-run units (as run), reading the scan of step 2
for T in $TISSUES; do
  echo "== step 3: run_conditional_v3m.R $T $CHR $(date)"
  "$RSCRIPT" code/run_conditional_v3m.R freeze "$T" real "$CHR" "$NCORES" > "step3_${T}.log" 2>&1
done

## checks
echo "== checks $(date)"
GTEX=$GTEX "$RSCRIPT" "$HERE/check_subset.R" "$OUT" "$RUN_ROOT" "$CHR" $TISSUES | tee check.txt
echo "end $(date)"
