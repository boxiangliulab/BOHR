#!/bin/bash
# Subset test of the analysis code (authored for the repository; not part of the original runs).
#
#   bash paper/test_subset/run_subset.sh paper/test_subset/config_lab.sh
#
# The short form of ../run_paper.sh: it runs the sub-steps that re-run the core of the paper with the scripts of paper/,
# as they were run, on two tissues and one chromosome, and compares every output with the stored results of the paper:
#   00.1  installs the BOHR package of this repository into OUT/lib;
#   01.2  promoter filter on the promoters of CHR with BOHR::relative_activity() (the package reproduces
#         0103_proActiv_data_process_tissue.R, whose 57 per-tissue proActiv inputs are no longer stored);
#   03.1  per-SNP scan: 03_cis_mapping_conditional/scan/run_final2.R for each tissue on CHR;
#   03.2  independent signals: 03_cis_mapping_conditional/conditional/run_conditional_v3m.R for the units of each tissue
#         and CHR that the final conditional run re-ran, and the signals against the paper's 14,800.
# Outputs go to OUT only; the report is OUT/check.txt, the checks OUT/summary.tsv, the logs OUT/logs/. It is
# run_paper.sh with ONLY="00.1 01.2 03.1 03.2" (one PBS job, about 45 minutes with 8 cores).
set -eo pipefail
CONFIG=$(readlink -f "${1:-${CONFIG:?usage: run_subset.sh <config>}}")
HERE=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)
set -a; source "$CONFIG"; set +a                                   # values of this config win over ../config_lab.sh
STEPS="00 01 03" ONLY="00.1 01.2 03.1 03.2" exec bash "$HERE/../run_paper.sh" "$HERE/../config_lab.sh"
