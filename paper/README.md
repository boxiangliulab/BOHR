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
  These files were written for this repository; the original job scripts are not included.
- **`test_subset/`** re-runs the key real-data steps on two tissues and one chromosome with one command and compares the
  outputs with the paper's stored results (for collaborators with access to the data).

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
| `test_subset/` | | one-command re-run of a subset, with checks |

Section numbers are those of the manuscript. Comments in some scripts use the numbering of an earlier draft, in which
GWAS colocalization was §11 and AlphaGenome §12.

## The BOHR package and this code

The package reproduces the central steps of this pipeline on the paper's data, checked value by value by the scripts in
`../workflow/parity/` (see its README):

| step | package function | folder here |
|---|---|---|
| promoter filter and relative activity | `relative_activity()`, `from_proactiv()` | `01_promoter_quantification/` |
| per-SNP tests of μ, p0 and p1; units, ACAT and FDR | `bohr_cis()` | `03_cis_mapping_conditional/scan/`, `units/` |
| independent signals | `bohr_conditional()` | `03_cis_mapping_conditional/conditional/` |
| donor-trait tests | `bohr_trait()` | `04_tissue_specificity_traits/traits/` |
| tissue means for tissue specificity | `tissue_means()` | `04_tissue_specificity_traits/tissue_specificity/` |

## Data

The analyses use GTEx v11 RNA-seq junction counts, genotypes (GTEx v9 VCF) and donor phenotypes, which are protected
and available through dbGaP (phs000424); GTEx v11 eQTL summary statistics (GTEx portal); ENCODE SCREEN cCREs and
Z-scores; and the GWAS summary statistics listed in Supplementary Table 7. No individual-level data are included here.

## Software

R 4.4–4.5 with data.table, VariantAnnotation, Matrix, coloc 5, TCC (ROKU), gamlss (the gamlss check and the published
scan only) and proActiv 1.16.0; polyester 1.38.0 and STAR for the simulated RNA-seq; Python 3 with pandas, pyarrow and the
AlphaGenome client; PLINK 1.9; QTLtools 1.2.
