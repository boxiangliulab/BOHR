# Subset test of the analysis code

One command re-runs the key real-data steps of the paper on two tissues and one chromosome, with the scripts of
`paper/` as they were run, and compares every output with the paper's stored results. It is the short form of
`../run_paper.sh`, which tests every step of the paper (see `../README.md`); `run_subset.sh` runs its sub-steps 00.1,
01.2, 03.1 and 03.2.

## Who can run it

The inputs are GTEx v11 individual-level data (dbGaP-protected), kept on the lab cluster in folders readable by the
group `bliu`. You need an account on that cluster in this group. Do not copy the data or the outputs anywhere else; the
outputs (`OUT`) hold per-sample values too, so delete them when the test is done.

## Run

From the root of the repository, on the lab cluster's login node (the script submits one PBS job; nothing heavy runs on
the login node):

```bash
bash paper/test_subset/run_subset.sh paper/test_subset/config_lab.sh
```

It prints the job ID and where the report, log and summary will be. Settings can be changed in `config_lab.sh` or on the command
line, for example `OUT=/scratch/$USER/bohr_test TISSUES="Kidney_Cortex" CHR=chr22 bash paper/test_subset/run_subset.sh
paper/test_subset/config_lab.sh`. The job uses 8 cores and up to 48 GB (queue `long`; the `normal` queue allows only
24 GB) and takes about 30–60 minutes. To run inside an interactive job instead, set `SUBMIT=0`.

## What it does

| sub-step | script (folder in `paper/`) | output in `OUT` |
|---|---|---|
| 00.1 | installs the BOHR package of this repository into `OUT/lib` | `lib/`, `logs/00.1.log` |
| 01.2 | promoter filter of Supplementary Methods §1 on the promoters of `CHR`, with `BOHR::relative_activity()` on the stored outputs of `0103_proActiv_data_process_tissue.R` (`test_subset/step1_filter.R`) | `step1_relative_activity.rds`, `logs/01.2.log` |
| 03.1 | per-SNP BOHR scan of every tested promoter of `CHR` in each tissue (`03_cis_mapping_conditional/scan/run_final2.R`) | `shadow/results/final/freeze/<tissue>/real/<chr>.rds`, `logs/03.1.log` |
| 03.2 | forward conditional selection for the units of each tissue and `CHR` that the final conditional run re-ran (`03_cis_mapping_conditional/conditional/run_conditional_v3m.R`), on the scan of 03.1 | `shadow/results/conditional_v3m/…`, `logs/03.2.log` |
| checks | `test_subset/check_subset.R`, one part after each sub-step | `check.txt` (the report), `summary.tsv` (one row per check) |

Step 01.2 uses the package because the 57 per-tissue proActiv files that `0103` reads are no longer stored; the package
reproduces `0103` exactly (repository `workflow/parity/12_relative_activity_parity.R`). The filter works promoter by
promoter and gene by gene, so one chromosome gives the paper's result for its genes. Steps 03.1 and 03.2 run the original
scripts unchanged. They read the stored inputs through `OUT/shadow`, a folder of symbolic links to the paper's run folder,
and write only under `OUT`.

## Checks

| check | passes when |
|---|---|
| 1a promoter set | the promoters kept on `CHR` equal the paper's final set on `CHR` (`Promoter_info_filter_5`) |
| 1b relative activity values | every value equals the paper's (`filter_5_recal_att_filter`), with the same missing values |
| 2a promoters and rules | per tissue, the same promoters with the same skip reason, family, sample counts, eligibility and SNP counts as the stored scan |
| 2b per-SNP statistics | the same SNP tests, and the log P values and estimates of μ, p0 and p1 (likelihood-ratio and Firth) within 10⁻⁶ of the stored ones; unpenalised p0 and p1 estimates under complete separation diverge and are not compared (those SNPs use the Firth fit, which is) |
| 3 signals of the re-run units | the same signals (promoter, component, order, lead variant), with log P at selection within 10⁻⁶ |
| 4 final signal table | those signals are in the paper's final table of 14,800 independent signals |

`check.txt` lists each sub-step with its log and one line per check (PASS, FAIL or SKIP); after a FAIL it names the log to
read, and the last lines give the totals. Please send `summary.tsv` and `check.txt` (no individual data in either) to the
authors.

## Defaults

Tissues `Kidney_Cortex` and `Brain_Substantia_nigra` (103 and 183 donors), chromosome `chr21`. Any of the 50 tissues
(`work/final/tissues.tsv` of the run folder, rows with `def == freeze`) and any chromosome can be used; larger tissues and
chromosomes take longer, and 03.2 runs only where the final conditional run re-ran a unit (`work/cond_v3m/units_to_run.tsv`).
