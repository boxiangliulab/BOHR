# Analysis code of the BOHR paper

This folder holds the analysis code of the manuscript *Genetic regulation of alternative promoter usage with BOHR
reveals distinct modes of transcription initiation*, organised by Methods section. The BOHR R package at the root of
this repository is the reusable form of the method; this folder is the code that produced the paper's results.

## What is here

- **The key analysis scripts, as they were run.** Each script is a copy of the file that ran, taken from the authors'
  storage. Scripts that only draw figures, format Supplementary Tables or Source Data, or check numbers for the text are
  not included, nor are earlier versions whose results the paper does not use, nor sensitivity analyses reported only in
  the Supplementary Information.
- **Cleaning.** Comments were cleaned and nothing else: comments in Chinese were translated into English or removed, and
  comments that referred to internal working notes were removed or corrected. In a few scripts the text of printed
  messages was cleaned the same way. Every edit keeps its line, so line numbers equal those of the original scripts, and
  the code of every edited file was verified identical to its source (R parse data; Python syntax tree and tokens; shell
  lines other than comments). `MANIFEST.tsv` gives, for every file, its source path, the md5 of the source and of the copy,
  the number of edited comment and string lines, and the Methods and Supplementary Methods sections it implements.
- **Paths.** The scripts keep the absolute paths of the original runs (the lab cluster and NSCC). Many read their root
  folders from environment variables (for example `BOHR2_GTEX`, `BOHR2_ROOT`, `BOHR2_VCF` and `BOHR2_CODE`); the folder
  READMEs say which.
- **`run.sh`** in each folder lists the commands in the order they were run, with the arguments of the original jobs.
  These files were written for this repository. The original job scripts are included where they were kept (the NSCC
  `.pbs` jobs of 03, 05 and 07); the others are not.
- **Two helper scripts written for this repository**, for steps whose original code was not kept:
  `01_promoter_quantification/export_txt.R` (the text exports of 0103 and the PEER input) and
  `03_cis_mapping_conditional/inputs/promoter_sets_rule.R` (the per-tissue promoter sets of cis mapping, as a rule). Both
  reproduce the paper's files or sets exactly; they are marked `(authored)` in `MANIFEST.tsv`.
- **`run_paper.sh`** runs the whole chain, folder by folder, as a test: each step's as-run scripts are run again on a
  subset (two tissues, one chromosome, two donor characteristics, one simulated sample size), or, where that is not
  possible, the step's numbers are recomputed from the stored outputs; every result is compared with the stored results
  and the published numbers, with a message, a log and PASS / FAIL / SKIP lines for every sub-step (see "Testing" below).
  `checks/` holds its comparison scripts and `with_package/` its package route. These files were written for this
  repository.
- **`test_subset/`** is the short form of the same test: the promoter filter, the per-SNP scan and the independent
  signals only (sub-steps 01.2, 03.1 and 03.2 of `run_paper.sh`).

## Layout

| folder | Methods (and Supplementary Methods) | what it does |
|---|---|---|
| `01_promoter_quantification/` | §1 | junction reads, proActiv per tissue, promoter filter, relative activity |
| `02_bohr_model_and_simulation/` | §2–3 | the BOHR engine, the check against gamlss, simulations and benchmarking (Fig. 2) |
| `03_cis_mapping_conditional/` | §4, §5, §7 | inputs and covariates, per-SNP scan, units and FDR, independent signals, LD, QTLtools comparison |
| `04_tissue_specificity_traits/` | §6 | tissue-specific promoters (ROKU) and donor-trait tests (Fig. 3) |
| `05_coloc_eqtl/` | §8 | colocalization with GTEx v11 eQTLs (Fig. 4) |
| `06_major_minor_promoters/` | §9 | major and minor promoters, eQTL direction, chromatin marks (Fig. 5) |
| `07_coloc_gwas_alphagenome/` | §10–11 | colocalization with 51 GWAS traits; AlphaGenome at MRAS (Fig. 6) |
| `run_paper.sh`, `config_lab.sh`, `checks/` | | the whole-chain test, with checks |
| `with_package/` | | step P of the test: the same subset with the BOHR package |
| `test_subset/` | | the short test (filter, scan, signals) |

Section numbers are those of the manuscript. Comments in some scripts use the numbering of an earlier draft, in which
GWAS colocalization was §11 and AlphaGenome §12.

## The BOHR package and this folder

The package and this folder do different jobs. This folder is the record of how the paper's numbers were made: the
scripts that ran, unchanged. The package is the method in a form others can run on their own data. The paper's numbers
were not produced with the package, and the scripts here were not rewritten to call it. Instead, the package was
checked against these scripts on the paper's data, value by value, with the scripts in `../workflow/parity/` (inputs
rebuilt from the stored runs, the package run on them, the outputs compared with the stored results). The table gives,
for each step, the as-run scripts, the package function that implements the same step, the parity script and its
result. Numbers are from those parity runs.

| step (Methods) | as-run scripts here | package function | parity script | result on the paper's data |
|---|---|---|---|---|
| promoter filter and relative activity (§1) | `01_promoter_quantification/0103_proActiv_data_process_tissue.R` | `from_proactiv()`, `relative_activity()` | `12_relative_activity_parity.R` | every count of the filter steps of SI §1, the final 10,826 promoters in 4,744 genes and all 10,826 × 22,865 relative activity values identical |
| BOHR likelihood and tests (§2) | `02_bohr_model_and_simulation/engine/bohr2_engine.R`, `bohr2_addon3.R` | `bohr_fit()`, `bohr_test()`, `bohr()`, `bohr_control()` | package tests `test-engine-parity.R`, `test-regression.R` | 540 of 540 engine values bit-identical on synthetic responses covering all four families, separation and the skip rules |
| per-SNP scan, units, ACAT, mirror and estimability rules, BH (§4–5) | `03_cis_mapping_conditional/scan/run_final2.R`, `run_firth_addon.R`, `run_muF.R`; `units/discovery_nogate.R`, `units_v3_estimable.R`, `units_v3m_mirror.R`, `mirror_pairs.R` | `bohr_cis()`, `bohr_combine()`, `acat()`, `mirror_pairs()` | `01_extract_inputs.R`, `02_run_bohr_cis.R`, `03_compare_stored.R` | Kidney_Cortex chr22: all 79,776 tests and 153 units reproduced, unit ACAT P to a relative 7.3 × 10⁻⁸, tissue-wide q to 2.0 × 10⁻¹¹, the same 29 significant units; BA24 chr22: 85,617 tests, the same 97 significant units |
| independent signals (§7) | `03_cis_mapping_conditional/conditional/run_conditional_*.R`, `merge_conditional_*.R` | `bohr_conditional()` | `04_extract_conditional_inputs.R`, `05_run_conditional.R`, `06_compare_conditional.R` | BA24, all 97 significant units: 106 of 106 forward signals and 104 of 104 kept signals with the same lead variant, round and backward call |
| tissue means for tissue specificity (§6) | `04_tissue_specificity_traits/tissue_specificity/gamlss_functions_optimized.R` (`raw_fit()`) | `tissue_means()` | `11_tissue_means_parity.R` | 56 tissues, 430,464 promoter–tissue pairs: the same pairs and rules, means within 0.001 where both fits converged (the paper fitted with gamlss); 5 of 35,870 ROKU calls cross the 0.10 effect filter differently |
| donor-characteristic tests (§6) | `04_tissue_specificity_traits/traits/run_20261005/run_trait.R`, `run_20261104/summarize.R`, `run_20261106_mirror_abs/` | `bohr_trait()` | `07_extract_trait_inputs.R`, `08_run_trait.R`, `09_compare_trait.R` | Breast × SEX and Tibial artery × AGE: all 8,593 and 7,152 tests, P to a relative 1.3 × 10⁻⁹, the same 1,243 and 355 significant tests, the published PDE5A and CSNK1D examples exactly |
| P-value combination (§2, §5) | inside the scripts above | `acat()`, `simes()` | package test `test-cis-parity.R` | identical to the engine's functions, including the underflow and small-P branches |

Steps without a package counterpart, which use the BOHR results but are not part of the method: proActiv quantification
(`0101`–`0102`), the simulations and benchmark of Fig. 2, the QTLtools comparison, colocalization with eQTLs and GWAS
(coloc), the major/minor-promoter and chromatin analysis of Fig. 5, and AlphaGenome at MRAS.

**One deliberate difference.** Each restricted Firth fit of p0 and p1 starts, in the package by default, from the null
fit as well as from the full fit, and keeps the higher maximum (`bohr_control(firth_null_start = TRUE)`). The paper's
engine started from the full fit only; `firth_null_start = FALSE` reproduces it exactly (all engine and fixture tests are
bit-identical with `FALSE`). The two settings differ only where a restricted penalised likelihood has two maxima. In the
parity runs this changed one decision: Brain_Putamen PromoterID_47842 p0, where the full-fit start gives a second signal
and the null start gives one signal, as in the published table, so that the package keeps the published 98 signals in
that tissue. No other unit, signal or characteristic call changed in the parity sets (Kidney_Cortex and BA24 chr22, BA24
in full, the two refitted units, two characteristics). Whether the published total of 14,800 signals is reproduced in
every tissue with the default is known only for the tissues re-run; a re-run of all 50 tissues would answer it.

**Running both.** `run_paper.sh` (below) re-runs the as-run scripts on a subset and compares them with the stored
results. Its step P (`with_package/run_package_route.sh`) runs the same subset with the package through the parity
scripts, with the default and, for the first tissue, with `firth_null_start = FALSE`, and reports how many tests the
two settings change.

## Testing

```sh
bash paper/run_paper.sh paper/config_lab.sh     # from the repository root, on the lab cluster's login node
```

`run_paper.sh` submits itself as one PBS job and works through the folders in order, each split into numbered
sub-steps (00.1 setup; 01.1–07.3; P.1–P.3 for the package route). A sub-step either re-runs as-run scripts on a subset
(*rerun*: Kidney_Cortex and Brain_Substantia_nigra, chr21, SEX and AGE, simulated N = 50) or recomputes the step's numbers
from the stored outputs of all tissues (*verify*), and then compares them with the stored results and the published
Supplementary Tables. Where a re-run is impossible, a SKIP line gives the reason (the per-tissue proActiv results and the
ENCODE archives are no longer stored; ROKU and the GWAS clumping ran on NSCC; AlphaGenome is a remote model). Scripts
that write to fixed folders of the original runs are run as copies in which only those folder strings point into the
output folder (`checks/redirect.py` makes the copy, checks that no other line changes and lists the changed lines in the
log). The inputs are read through symbolic links, never written. Each sub-step has its own log in `OUT/logs/`;
`OUT/check.txt` is the report and `OUT/summary.tsv` has one row per check. Settings are in `config_lab.sh` (paths of the
stored data, subset, PBS) and can be overridden on the command line (`STEPS`, `ONLY`, `KEEP_GOING`, `RESUME`,
`TISSUES`, ...). The run needs
access to the lab's copy of the GTEx data; delete `OUT` afterwards, since it holds per-sample values.

## Data

The analyses use GTEx v11 RNA-seq junction counts, genotypes (GTEx v9 VCF) and donor phenotypes, which are protected
and available through dbGaP (phs000424); GTEx v11 eQTL summary statistics (GTEx portal); ENCODE SCREEN cCREs and
Z-scores; and the GWAS summary statistics listed in Supplementary Table 7. No individual-level data are included here.

## Software

R 4.4–4.5 with data.table, VariantAnnotation, Matrix, coloc 5, TCC (ROKU), gamlss (the gamlss check and the published
scan only) and proActiv 1.16.0; polyester 1.38.0 and STAR for the simulated RNA-seq; Python 3 with pandas, pyarrow and the
AlphaGenome client; PLINK 1.9; QTLtools 1.2.
