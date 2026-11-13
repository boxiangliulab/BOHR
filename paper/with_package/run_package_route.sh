#!/bin/bash
# Package route of run_paper.sh, step P (authored for the repository; not part of the original runs).
#
# Runs the core chain of the paper with the BOHR package on the subset of run_paper.sh and compares it with the stored
# results of the paper's own scripts: bohr_cis() on one chromosome of each tissue (P.1), bohr_conditional() chained on
# that run (P.2), bohr_trait() for one donor characteristic in each tissue (P.3). It uses the real-data parity scripts of
# the repository (../../workflow/parity/), which rebuild the stored inputs as plain matrices, run the package and print
# the comparison; checks/check_package.R turns their summaries into PASS / FAIL lines.
#
# The package default firth_null_start = TRUE starts each restricted Firth fit of p0 and p1 also from the null fit; the
# paper's engine started it from the full fit only (firth_null_start = FALSE). The two give the same statistics wherever
# the penalised likelihood has one maximum, and the null start finds the higher maximum where there are two (in the paper's
# data, a Brain_Putamen unit, whose spurious second signal it removes; the published total of 98 signals in that tissue
# is what the package gives). P.1 also runs the first tissue with firth_null_start = FALSE, which must equal the paper's
# scan, and counts the tests in which the two settings differ.
#
#   called by run_paper.sh as: run_package_route.sh cis|conditional|trait   (environment: OUT TISSUES CHR TRAITS NCORES RSCRIPT LIB PAPER)
set -eo pipefail
part=$1; W=$(cd "$PAPER/../workflow/parity" && pwd); C=$PAPER/checks; P=$OUT/package; mkdir -p "$P"
read -r -a TIS <<< "$TISSUES"; read -r -a TRT <<< "$TRAITS"
RS=$RSCRIPT; R_() { R_LIBS=$LIB${R_LIBS:+:$R_LIBS} "$RS" "$@"; }
case $part in
  cis)
    for T in "${TIS[@]}"; do
      echo "== $T $CHR: inputs (01_extract_inputs.R)"
      [ -f "$P/in_$T.rds" ] || "$RS" "$W/01_extract_inputs.R" "$T" "$CHR" "$P/in_$T.rds"
      echo "== $T $CHR: bohr_cis(), package default (02_run_bohr_cis.R)"
      R_ "$W/02_run_bohr_cis.R" "$P/in_$T.rds" "$P/cis_$T.rds" "$NCORES" 0 "$LIB"
      "$RS" "$W/03_compare_stored.R" "$P/cis_$T.rds" "$T" "$CHR" "$P/cmp_cis_$T.rds"
      "$RS" "$C/check_package.R" cis "$P/cmp_cis_$T.rds" "$T" "$CHR" "default (firth_null_start = TRUE)" "$OUT/status/P.1.tsv"
    done
    T=${TIS[0]}
    echo "== $T $CHR: bohr_cis(), firth_null_start = FALSE (the paper's engine)"
    BOHR_FIRTH_NULL_START=FALSE R_ "$W/02_run_bohr_cis.R" "$P/in_$T.rds" "$P/cis_${T}_paper_engine.rds" "$NCORES" 0 "$LIB"
    "$RS" "$W/03_compare_stored.R" "$P/cis_${T}_paper_engine.rds" "$T" "$CHR" "$P/cmp_cis_${T}_paper_engine.rds"
    "$RS" "$C/check_package.R" cis_exact "$P/cmp_cis_${T}_paper_engine.rds" "$T" "$CHR" "firth_null_start = FALSE" "$OUT/status/P.1.tsv"
    "$RS" "$C/check_package.R" null_start "$P/cis_${T}_paper_engine.rds" "$P/cis_$T.rds" "$T" "$CHR" "" "$OUT/status/P.1.tsv"
    ;;
  conditional)
    for T in "${TIS[@]}"; do
      echo "== $T $CHR: bohr_conditional() on the bohr_cis() run of P.1 (05_run_conditional.R cis)"
      R_ "$W/05_run_conditional.R" cis "$P/in_$T.rds" "$P/cis_$T.rds" "$P/cond_$T.rds" "$NCORES" "$LIB"
      "$RS" "$W/06_compare_conditional.R" "$P/cond_$T.rds" "$P/cmp_cond_$T.rds"
      "$RS" "$C/check_package.R" conditional "$P/cmp_cond_$T.rds" "$T" "$CHR" "" "$OUT/status/P.2.tsv"
    done
    ;;
  trait)
    for T in "${TIS[@]}"; do
      TR=${TRT[0]}; D=$P/trait_${T}_$TR; mkdir -p "$D"
      echo "== $T x $TR: inputs (07_extract_trait_inputs.R)"
      "$RS" "$W/07_extract_trait_inputs.R" "$D" "$T:$TR"
      IN=$D/trait_${T}_${TR}_inputs.rds
      echo "== $T x $TR: bohr_trait() (08_run_trait.R)"
      R_ "$W/08_run_trait.R" "$IN" "$D/run.rds" "$NCORES" "$LIB"
      "$RS" "$W/09_compare_trait.R" "$D/run.rds" "$D/cmp.rds"
      "$RS" "$C/check_package.R" trait "$D/cmp.rds" "$T" "$TR" "default (firth_null_start = TRUE)" "$OUT/status/P.3.tsv"
    done
    T=${TIS[0]}; TR=${TRT[0]}; D=$P/trait_${T}_$TR
    echo "== $T x $TR: bohr_trait(), firth_null_start = FALSE (the paper's engine)"
    BOHR_FIRTH_NULL_START=FALSE R_ "$W/08_run_trait.R" "$D/trait_${T}_${TR}_inputs.rds" "$D/run_paper_engine.rds" "$NCORES" "$LIB"
    "$RS" "$W/09_compare_trait.R" "$D/run_paper_engine.rds" "$D/cmp_paper_engine.rds"
    "$RS" "$C/check_package.R" trait_exact "$D/cmp_paper_engine.rds" "$T" "$TR" "firth_null_start = FALSE" "$OUT/status/P.3.tsv"
    ;;
  *) echo "usage: run_package_route.sh cis|conditional|trait"; exit 1 ;;
esac
