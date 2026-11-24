# Analysis code of the BOHR paper

This folder holds the analysis code of the manuscript *Genetic regulation of alternative promoter usage with BOHR
reveals distinct modes of transcription initiation*, organised by Methods section. The BOHR R package at the root of
this repository is the reusable form of the method; this folder is the code that produced the paper's results.

## Version of the analysis

The paper reports the final version of the analysis (November 2026). It ran the scripts of this folder, with the changes
below, on a separate copy of the data tree. Releases v0.1.1 and v0.1.2 of this repository hold the earlier version of the folder.

- **No junction-read floor.** proActiv used GTEx's STAR junction files of every sample with no minimum read count. The
  promoter matrix holds 22,939 promoters in 9,371 genes (earlier version: 10,826 promoters in 4,744 genes).
- **Skip rule of the fits.** A promoter is skipped in a tissue when the standard deviation of its relative activity is
  below 0.005 (earlier: 0.05 or less; `02_bohr_model_and_simulation/engine/bohr2_data.R`). The promoters this adds were
  scanned separately and merged with the first scan (`03/inputs/build_v4.R`, `03/scan/merge_v4.R`). Their tests entered
  the FDR only for the components that passed a permutation calibration fixed beforehand, μ and p0
  (`03/units/units_v3m_mirror.R`).
- **Boundary and read-depth rules.** A promoter is not tested in a tissue when 80% or more of the tissue's analysed
  samples lie at one boundary, or when its gene has a median of fewer than five first-intron junction reads among those
  samples; Benjamini–Hochberg correction is then applied again within each tissue (`03/units/make_version_inputs.R`,
  `build_r5e80_units.R`).
- **Independent signals and eQTL colocalization** were run afresh for all significant units (`03/conditional/`,
  `05_coloc_eqtl/v3m/`).
- **Donor traits.** The trait tests kept the earlier skip rule (standard deviation 0.05 or less) and then apply the
  boundary rule and a gene read-depth rule, with gene medians over all samples of the tissue, and leave out tests of rare
  binary traits (`04/traits/r5e80_20261119/`, `r5e80_ruleb_20261119/`).
- **Tissue specificity** uses `BOHR::tissue_means()` on the GTEx v11 analysis-freeze samples (SMAFRZE = RNASEQ) of the
  genotyped donors in 56 tissues, then ROKU and the 0.10 effect filter (`04/tissue_specificity/roku_freeze51.R`).
- **GWAS clumping** uses an ancestry-matched LD reference: the 799 donors of European ancestry for the 32 European or
  mainly European GWAS, and all 943 donors for the 19 East Asian GWAS (`07/gwas_coloc/`).

Smaller changes that the new data needed (a numeric type for all-missing Firth columns, count assertions turned into
reports, the element set of the chromatin analysis, file names of the new run) are listed in each folder's README.

## What is here

- **The key analysis scripts, as they were run in the final version.** Each script is a copy of the file that ran. Scripts
  that only draw figures, format Supplementary Tables or Source Data, or check numbers for the text are not included, nor
  are side analyses whose outputs the paper does not use, nor sensitivity analyses reported only in the Supplementary
  Information.
- **Cleaning.** Comments were cleaned and nothing else: comments in Chinese were translated into English or removed, and
  comments that referred to internal working notes were removed or corrected. In a few scripts the text of printed
  messages was cleaned the same way. Every edit keeps its line, so line numbers equal those of the original scripts, and
  the code of every edited file was verified identical to its source (R parse data; Python syntax tree and tokens; shell
  lines other than comments). The exceptions are stated in the folder READMEs:
  - In `05_coloc_eqtl/v3m/05_coloc_tissue_v3m.R` and `05_coloc_eqtl/ldwin/clumped_leads/s5_coloc_tissue_v2.R`, one line
    was corrected after the runs. As run, it stored 0 in the `n_int` column (variants in the intersection of the two LD
    windows) of clumped-lead pairs, because data.table resolved `n_int = n_int` to the table's own column. The filter on
    that count used the local value, so no pair or call changes; Supplementary Table 5a carries the recomputed counts.
  - `04_tissue_specificity_traits/tissue_specificity/roku_freeze51.R` calls `TCC::ROKU`; as run, it sourced a copy of
    TCC 1.52.0's `R/ROKU.R`, which gives the same calls.
  - The runs of 2026-11-24 (folders 04 and 07) wrote to a folder whose name is given here as its alias,
    `nofloor_full/consult/freeze_ld_20261124`.
- **`MANIFEST.tsv`** gives, for every file, its source path (the copy that ran), the md5 of the source and of the copy,
  the number of edited comment and string lines, and the Methods and Supplementary Methods sections it implements.
- **Paths.** The scripts keep the absolute paths of the earlier runs (the lab cluster and NSCC). The final run used
  copies in which only folder prefixes and run roots differ, made with `checks/redirect.py`: on the lab cluster,
  `/projects/group/bliu/e1327989/GTEx` and `/home/e1327989/project/GTEx` became `…/GTEx/nofloor_full/GTEx`; on NSCC the
  run roots moved to `/scratch/users/nus/e1327989/nofloor_full/`; the cis-mapping root `QTL_figure/rerun_v2` became
  `rerun_v4` (scans and units of all arms) and `rerun_r5e80` (final units and everything after them), and the GWAS root
  `new_GWAS/coloc_final_nogate_v3` became `coloc_final_nogate_v3_r5e80`. Scripts added for the final run keep their own
  paths. Each folder README gives its mapping. Many scripts read their roots from environment variables (for example
  `BOHR2_GTEX`, `BOHR2_ROOT`, `BOHR2_VCF` and `BOHR2_CODE`).
- **`run.sh`** in each folder lists the commands of the final run in the order they ran, with the arguments of the
  original jobs. These files were written for this repository. The original job scripts are included where they were
  kept (the `.pbs` jobs of 03, 05 and 07); the others are not.

## Layout

| folder | Methods (and Supplementary Methods) | what it does |
|---|---|---|
| `01_promoter_quantification/` | §1 | proActiv per tissue, promoter filter, relative activity; the release comparison (Fig. 1A,B) |
| `02_bohr_model_and_simulation/` | §2–3 | the BOHR engine, the check against gamlss, simulations and benchmarking (Fig. 2) |
| `03_cis_mapping_conditional/` | §4, §5, §7 | inputs and covariates, per-SNP scan, units and FDR, independent signals, LD, QTLtools comparison |
| `04_tissue_specificity_traits/` | §6 | tissue-specific promoters (ROKU) and donor-trait tests (Fig. 3) |
| `05_coloc_eqtl/` | §8 | colocalization with GTEx v11 eQTLs (Fig. 4) |
| `06_major_minor_promoters/` | §9 | major and minor promoters, eQTL direction, chromatin marks (Fig. 5) |
| `07_coloc_gwas_alphagenome/` | §10–11 | colocalization with 51 GWAS traits; AlphaGenome at MRAS (Fig. 6) |
| `checks/redirect.py` | | the folder-prefix redirect that made the run copies |

Section numbers are those of the manuscript. Comments in some scripts use the numbering of an earlier draft, in which
GWAS colocalization was §11 and AlphaGenome §12.

## The BOHR package and this folder

The package and this folder do different jobs. This folder is the record of how the paper's numbers were made: the
scripts that ran. The package is the method in a form others can run on their own data. Apart from the tissue means of
§6, which the final version computes with `BOHR::tissue_means()`, the paper's numbers were not produced with the package,
and the scripts here were not rewritten to call it. Instead, the package was checked against these scripts value by
value, with the scripts in `../workflow/parity/` (inputs rebuilt from the stored runs, the package run on them, the
outputs compared with the stored results). These parity runs used the stored outputs of the earlier version (with the
five-read floor). The table gives, for each step, the scripts here, the package function that implements the same step,
the parity script and its result.

| step (Methods) | scripts here | package function | parity script | result on the earlier version's data |
|---|---|---|---|---|
| promoter filter and relative activity (§1) | `01_promoter_quantification/0103_proActiv_data_process_tissue.R` | `from_proactiv()`, `relative_activity()` | `12_relative_activity_parity.R` | every count of the filter steps, the final 10,826 promoters in 4,744 genes and all 10,826 × 22,865 relative activity values identical |
| BOHR likelihood and tests (§2) | `02_bohr_model_and_simulation/engine/bohr2_engine.R`, `bohr2_addon3.R` | `bohr_fit()`, `bohr_test()`, `bohr()`, `bohr_control()` | package tests `test-engine-parity.R`, `test-regression.R` | 540 of 540 engine values bit-identical on synthetic responses covering all four families, separation and the skip rules |
| per-SNP scan, units, ACAT, mirror and estimability rules, BH (§4–5) | `03_cis_mapping_conditional/scan/run_final2.R`, `run_firth_addon.R`, `run_muF.R`; `units/discovery_nogate.R`, `units_v3_estimable.R`, `units_v3m_mirror.R`, `mirror_pairs.R` | `bohr_cis()`, `bohr_combine()`, `acat()`, `mirror_pairs()` | `01_extract_inputs.R`, `02_run_bohr_cis.R`, `03_compare_stored.R` | Kidney_Cortex chr22: all 79,776 tests and 153 units reproduced, unit ACAT P to a relative 7.3 × 10⁻⁸, tissue-wide q to 2.0 × 10⁻¹¹, the same 29 significant units; BA24 chr22: 85,617 tests, the same 97 significant units |
| independent signals (§7) | `03_cis_mapping_conditional/conditional/run_conditional_v3m.R`, `merge_conditional_*.R` | `bohr_conditional()` | `04_extract_conditional_inputs.R`, `05_run_conditional.R`, `06_compare_conditional.R` | BA24, all 97 significant units: 106 of 106 forward signals and 104 of 104 kept signals with the same lead variant, round and backward call |
| tissue means for tissue specificity (§6) | `04_tissue_specificity_traits/tissue_specificity/roku_freeze51.R` (calls `tissue_means()`) | `tissue_means()` | `11_tissue_means_parity.R` | against the earlier gamlss fits (`raw_fit()`): 56 tissues, 430,464 promoter–tissue pairs, the same pairs and rules, means within 0.001 where both fits converged; 5 of 35,870 ROKU calls cross the 0.10 effect filter differently |
| donor-characteristic tests (§6) | `04_tissue_specificity_traits/traits/run_20261007/run_trait.R`, `run_20261104/summarize.R`, `run_20261106_mirror_abs/` | `bohr_trait()` | `07_extract_trait_inputs.R`, `08_run_trait.R`, `09_compare_trait.R` | Breast × SEX and Tibial artery × AGE: all 8,593 and 7,152 tests, P to a relative 1.3 × 10⁻⁹, the same 1,243 and 355 significant tests |
| P-value combination (§2, §5) | inside the scripts above | `acat()`, `simes()` | package test `test-cis-parity.R` | identical to the engine's functions, including the underflow and small-P branches |

Steps without a package counterpart, which use the BOHR results but are not part of the method: proActiv quantification,
the simulations and benchmark of Fig. 2, the QTLtools comparison, colocalization with eQTLs and GWAS (coloc), the
major/minor-promoter and chromatin analysis of Fig. 5, and AlphaGenome at MRAS.

**Package defaults and the final rules.** The defaults of `bohr_control()` are the earlier skip rules (`min_sd = 0.05`,
`max_boundary = 0.95`). `bohr_control(min_sd = 0.005, max_boundary = 0.80)` applies the final ones; the read-depth rule
is a filter on the input (genes with a median of at least five first-intron junction reads in the tissue).

**One deliberate difference.** Each restricted Firth fit of p0 and p1 starts, in the package by default, from the null
fit as well as from the full fit, and keeps the higher maximum (`bohr_control(firth_null_start = TRUE)`). The paper's
engine started from the full fit only; `firth_null_start = FALSE` reproduces it exactly (all engine and fixture tests are
bit-identical with `FALSE`). The two settings differ only where a restricted penalised likelihood has two maxima. In the
parity runs this changed one decision: Brain_Putamen PromoterID_47842 p0, where the full-fit start gives a second signal
and the null start gives one signal, as in the published table of the earlier version. No other unit, signal or
characteristic call changed in the parity sets (Kidney_Cortex and BA24 chr22, BA24 in full, the two refitted units, two
characteristics).

## Testing

Releases v0.1.1 and v0.1.2 have a whole-chain test (`run_paper.sh`, `checks/`, `with_package/`, `test_subset/`) that re-ran these
scripts on a subset and compared them with the stored outputs of the earlier version. It does not carry over to the final
version (the skip rule and several inputs changed), so it is not included here. The package's own tests are in `../tests/`.

## Data

The analyses use GTEx v11 RNA-seq junction counts, genotypes (GTEx v9 VCF) and donor phenotypes, which are protected
and available through dbGaP (phs000424); GTEx v11 eQTL summary statistics (GTEx portal); ENCODE SCREEN cCREs and
Z-scores; and the GWAS summary statistics listed in Supplementary Table 7. No individual-level data are included here.

## Software

R 4.4–4.5 with data.table, VariantAnnotation, Matrix, coloc 5, TCC (ROKU), gamlss (the gamlss check and the published
scan only) and proActiv 1.16.0; polyester 1.38.0 and STAR for the simulated RNA-seq; Python 3 with pandas, pyarrow and the
AlphaGenome client; PLINK 1.9; QTLtools 1.2.
