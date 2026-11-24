# Colocalization with cis-eQTLs

**Methods §8; Supplementary Methods §8.** Each independent puQTL signal (its leave-one-out scan when the unit has more
than one) is tested by coloc against the GTEx v11 eQTLs of its gene in the same tissue, over the union of the two LD
windows (farthest partner with r² ≥ 0.1, ± 50 kb). For an eGene, the signal is tested against every clumped eQTL lead of
the gene and keeps its best call; for other genes the lead is the GTEx top variant. Each signal is then classed as
shared, distinct, no eQTL or unresolved (Fig. 4E, Supplementary Table 5).

**Run order.** The final run ran the `v3m/` chain once, on the final units and signals (v3c), with every signal new. The
earlier analysis set that `v3m/01_analysis_set_pairs.R` compares with (`results/coloc_eqtl_v3/analysis_set.rds`) was an
empty table, so no coloc result of an earlier run was carried.
1. Set-up (commands in `run.sh`; no script makes these): the empty earlier analysis set, the output folder
   `results/coloc_eqtl_v3m/new_inputs/sig/`, and the link `results/coloc_eqtl_ldwin` to the earlier LD-window run.
2. `v3m/01_analysis_set_pairs.R` (analysis set; pairs of the signals whose LD window is already known) →
   `v3m/02_new_inputs.R` (LD windows, clumpings and eQTL extractions still needed) → `v3m/02b_sig_pairs.py` (GTEx
   significant variants of the eGenes to clump).
3. NSCC job `v3m/03_ld_clump_nscc.pbs`: PLINK LD windows of the new leads, and recursive clumping of the new eGenes with
   `ldwin/clumped_leads/s2_clump_nscc.R`; its outputs are copied to `results/coloc_eqtl_v3m/nscc/`.
4. `v3m/04_pairs_new_all.R` (all pairs and eQTL windows) → `v3m/04b_extract_eqtl.py` (eQTL subsets, top and clump) →
   `v3m/05_coloc_tissue_v3m.R` per tissue → `v3m/06_final_v3m.R` (signal table, gene–tissue classes) →
   `v3m/eqtl_signal_classes.R` (per-signal classes: of 30,043 signals, 7,573 shared, 6,212 distinct, 5,645 no eQTL and
   10,613 unresolved).
5. Links `results/coloc_eqtl_v3c` → `coloc_eqtl_v3m` and `coloc_eqtl_v3m/final_clump_v3c` → `final_clump_v3m`: the
   names that GWAS step 12a (`../07_coloc_gwas_alphagenome/`) and the figure code read. The Figure 5 stage
   (`../06_major_minor_promoters/`) links `final_clump_v3m/coloc_eqtl_signals.tsv.gz` directly.

**Inputs from the earlier run.** LD windows depend only on (tissue, SNP), clumped eQTL leads only on (tissue, gene) and
eQTL subsets only on (tissue, gene, region). The `v3m/` chain therefore reads those of the earlier LD-window run
(`results/coloc_eqtl_ldwin/`, built from the units v3 signals in 2026-10, linked read-only) and computes only what is
missing. In the final run the earlier windows covered 5,134 of the 30,043 signals; 31,929 new LD windows and the
clumping of 8,070 eGenes were computed on NSCC, and eQTL subsets were cut for every gene whose regions the earlier
subsets did not cover. The earlier run is in `ldwin/` and `ldwin/clumped_leads/`, with its analysis set from
`v3/00_build_analysis_set.R`. Its final table (`ldwin/08_final_clump.R`) is read by `v3m/06_final_v3m.R` only for
comparison columns (`fig4e_mirror_v2`, `detail_mirror_v2`, `compare_*.tsv`). `v3/01_extract_eqtl.py` to
`v3/04_classify_all_shared.R` are the earlier signal ± 200 kb run, which `ldwin/07_merge_classify.R` reads.

**Changes for the final run.** Three scripts differ from the earlier run:
- `v3m/01_analysis_set_pairs.R` and `v3m/05_coloc_tissue_v3m.R` read the final units and signals (`final_nogate_v3c`,
  `conditional_final_v3c`) instead of `final_nogate_v3m` and `conditional_final_v3m`.
- `v3c/coloc_final_lib.R` (line 85) makes the main run's Firth columns numeric before they are filled from the add-on:
  when a chromosome has no Firth-eligible promoter the column is all NA, R reads it as logical, and the filled P values
  would become TRUE.
- `v3m/03_ld_clump_nscc.pbs` asks for 24 h, 128 cores and 300 GB (earlier 4 h, 16 cores, 64 GB) and clumps 32 tissues at
  a time (`xargs -P 32`) instead of one after another; each call is unchanged.

**One deviation from the code as run.** `v3m/05_coloc_tissue_v3m.R` line 82 (and line 61 of the earlier
`ldwin/clumped_leads/s5_coloc_tissue_v2.R`) was corrected after the runs: the `n_int` fix described in `../README.md`.
No pair or call changes.

**Paths.** The scripts keep the paths of the earlier runs. In the final run they ran as copies in which only folder
names were changed:

| path in the scripts | folder in the final run |
|---|---|
| `/projects/group/bliu/e1327989/GTEx`, `/home/e1327989/project/GTEx` | `/projects/group/bliu/e1327989/GTEx/nofloor_full/GTEx` (links to the inputs that do not depend on promoter activity, new folders for those that do) |
| `…/GTEx/QTL_figure/rerun_v2` | `…/nofloor_full/GTEx/QTL_figure/rerun_r5e80` |
| `…/GTEx/new_GWAS/coloc_final_nogate_v3/scripts/coloc_final_lib.R` (sourced) | the run copy of `v3c/coloc_final_lib.R` |
| `NDIR` of `coloc_final_lib.R`, `…/GTEx/new_GWAS/coloc_final_nogate_v3` | `…/nofloor_full/GTEx/new_GWAS/coloc_final_nogate_v3_r5e80` |
| NSCC `/scratch/users/nus/e1327989/rerun_v2/ld_ref_v9/v3m` | `/scratch/users/nus/e1327989/nofloor_full/eqtl_ld_r5e80` |

`results/coloc_eqtl_ldwin` of the final run root links to `/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/results/coloc_eqtl_ldwin`.
The library reads `UNITSET` (`v3c` in the final run). coloc is not part of the BOHR package.

**Inputs without a script.** `panel/panel_snps.tsv.gz` (read by `ldwin/03_rarelead_windows.R`, `ldwin/08_final_clump.R`
and `v3m/06_final_v3m.R`) is columns 1, 2 and 4 of the 22 PLINK `.bim` files of the MAF ≥ 0.05 panel
(`LD_ref_v9_plink19`; checked: identical, 6,395,957 variants); `run.sh` gives the command.

**Counts written into the scripts.** `ldwin/08_final_clump.R` asserts 14,783 signals (lines 18 and 52), 50 tissues
(line 21) and 8,283 gene–tissue pairs (line 85), and `ldwin/07_merge_classify.R` 8,283 (lines 55 and 63). These hold for
the earlier run; the final run does not run them. The `v3m/` chain has no counts written in.

`run.sh` lists the commands in the order they were run.


## Files

| file | what it does | inputs → outputs |
|---|---|---|
| `config/traits_51_v3eur.txt` | Lists the 51 GWAS trait files; coloc_final_lib.R reads it when sourced, so the eQTL coloc scripts need it too. | read at load by coloc_final_lib.R |

## v3m/ (the final run)

| file | what it does | inputs → outputs |
|---|---|---|
| `v3m/01_analysis_set_pairs.R` | Builds the analysis set (every backward-kept signal of every significant unit) and the pairs of the signals whose LD window is already known. | units v3c, conditional_final_v3c, v11 eGenes, ld_v3/donors, the earlier analysis set (empty), coloc_eqtl_ldwin windows and pairs → coloc_eqtl_v3m/{analysis_set.rds, missing.tsv, pairs_clump_new.tsv.gz, pairs_top_new.tsv.gz} |
| `v3m/02_new_inputs.R` | Lists the LD windows, clumpings and eQTL extractions still needed. | coloc_eqtl_v3m/{analysis_set.rds, missing.tsv}, coloc_eqtl_ldwin windows and leads, sens_clump pairs, GTEx v11 eGenes → coloc_eqtl_v3m/new_inputs/ |
| `v3m/02b_sig_pairs.py` | Extracts the GTEx v11 significant eQTL variants of the eGenes to clump. | new_inputs/clump_genes.tsv, GTEx v11 eQTLs.signif_pairs → new_inputs/sig/<T>.tsv.gz |
| `v3m/03_ld_clump_nscc.pbs` | NSCC job: LD windows of the new leads (panel; option A' for leads absent from it) and clumping of the new eGenes with s2_clump_nscc.R, 32 tissues at a time. | new_inputs/snps, new_inputs/sig → out/*.win.tsv, clump_out/<T>/{leads.tsv, lead_partners.tsv} |
| `v3m/04_pairs_new_all.R` | Completes the top-variant and clumped-lead pairs of every signal from the step-1 pairs and the NSCC windows and clumps. | analysis_set.rds, pairs_*_new.tsv.gz, nscc/{out, clump_out}, new_inputs → pairs_{top,clump}_new_all.tsv.gz, eqtl_windows_{top,clump}_new.tsv |
| `v3m/04b_extract_eqtl.py` | Cuts the GTEx v11 eQTL all-pairs parquet files to the new windows (top and clumped lead sets). | eqtl_windows_<kind>_new.tsv, GTEx v11 all-pairs parquet → eqtl_subset_<kind>/<T>.tsv.gz |
| `v3m/05_coloc_tissue_v3m.R` | Runs coloc for the signals of one tissue against the GTEx top variant and every clumped eQTL lead. | analysis_set.rds, pairs_*_new_all, eqtl_subset_* and the coloc_eqtl_ldwin subsets, conditional_final_v3c/loo → per_tissue_{top,clump}/<T>.rds |
| `v3m/06_final_v3m.R` | Takes the best call per signal and builds the signal table and the gene–tissue classes, with comparison columns from the earlier table. | analysis_set.rds, per_tissue_*, pairs_clump_new_all, coloc_eqtl_ldwin/final_clump (comparison), panel_snps.tsv.gz, mirror_v3/unit_map_v2_to_v3.tsv.gz → coloc_eqtl_v3m/final_clump_v3m/ |
| `v3m/eqtl_signal_classes.R` | Classes each signal as shared (H4), distinct (H3, different locus or H1), no eQTL (gene not a GTEx v11 eGene in the tissue) or unresolved. | final_clump_v3m/coloc_eqtl_signals.tsv.gz → final_clump_v3m/signal_classes/{signal_classes.tsv.gz, signal_class_counts.tsv} |

## v3c/

| file | what it does | inputs → outputs |
|---|---|---|
| `v3c/coloc_final_lib.R` | Provides the shared puQTL × eQTL and puQTL × GWAS coloc functions (priors p1 = p2 = 10⁻⁴, p12 = 10⁻⁵); `UNITSET` sets the unit set it reads. | sourced by every coloc script of this folder and of ../07_coloc_gwas_alphagenome/; UNITSET=v3c reads final_nogate_v3c and conditional_final_v3c |

## ldwin/ and ldwin/clumped_leads/ (earlier run; inputs of the final run)

| file | what it does | inputs → outputs |
|---|---|---|
| `ldwin/00_leads.R` | Lists the puQTL signal SNPs and GTEx top eQTL variants that need LD windows, per tissue and chromosome. | coloc_eqtl_v3/analysis_set.rds, v11 eGenes → coloc_eqtl_ldwin/{leads.tsv, eqtl_leads.tsv, ldwin_inputs/snps/<T>.<chr>.txt} |
| `ldwin/01_ldwin_nscc.pbs` | NSCC job: PLINK 1.9 r2 >= 0.1 partners within 1 Mb of every lead, in the tissue's donors (MAF >= 0.05 panel). | ldwin/inputs/snps/<tissue>.<chr>.txt → ldwin/out/*.win.tsv, *.absent.txt |
| `ldwin/rarelead.pbs` | NSCC job: bfiles of the unfiltered v9 VCF (MAF >= 0.005) and the LD partners of the eQTL leads absent from the panel. | ldwin/out/*.absent.txt → rarelead/*.ld (read by 03_rarelead_windows.R) |
| `ldwin/02_windows.R` | Assembles LD windows (farthest r2 >= 0.1 partner +/- 50 kb) per lead and per signal-eQTL pair. | leads.tsv, ldwin_out/*.win.tsv, analysis_set.rds, eqtl_leads.tsv → coloc_eqtl_ldwin/{lead_windows.tsv, pair_windows.tsv.gz} |
| `ldwin/03_rarelead_windows.R` | Builds LD windows for eQTL leads absent from the MAF >= 0.05 panel (option A': panel partners only). | lead_windows.tsv, rarelead/*.ld, panel/panel_snps.tsv.gz, pair_windows.tsv.gz → rarelead_windows.tsv, rarelead_pairs.tsv |
| `ldwin/04_pairs.R` | Pairs each analysis signal with its gene's GTEx top eQTL variant and builds the LD-window intersections and unions. | coloc_eqtl_v3/analysis_set.rds, eqtl_leads.tsv, lead_windows.tsv, rarelead_windows.tsv → coloc_eqtl_ldwin/{pairs.tsv.gz, eqtl_windows.tsv} |
| `ldwin/05_extract_eqtl.py` | Cuts the GTEx v11 eQTL all-pairs parquet to the genes and windows the coloc needs, per tissue. | GTEx v11 all-pairs parquet, coloc_eqtl_ldwin/eqtl_windows.tsv → coloc_eqtl_ldwin/eqtl_subset/<tissue>.tsv.gz |
| `ldwin/06_coloc_tissue.R` | Runs puQTL x eQTL coloc for one tissue over each pair's LD-window union (GTEx top-variant leads). | analysis_set.rds, pairs.tsv.gz, eqtl_subset/<T>.tsv.gz, conditional_final_v3/loo → coloc_eqtl_ldwin/per_tissue/<T>.rds |
| `ldwin/07_merge_classify.R` | Merges the top-variant coloc runs and classifies each signal; its non-eGene rows enter the earlier final table. | coloc_eqtl_v3/analysis_set.rds, per_tissue/*.rds, coloc_eqtl_v3/all_shared_rule/ (comparison) → coloc_eqtl_ldwin/{coloc_eqtl_signals*,gene_tissue_classes}.tsv.gz |
| `ldwin/08_final_clump.R` | Builds the earlier final eQTL coloc table (best call per signal over the clumped eQTL leads), read for comparison by v3m/06_final_v3m.R. | sens_clump/{per_tissue_v2, pairs.tsv.gz}, coloc_eqtl_signals.tsv.gz, gene_tissue_classes, panel/panel_snps.tsv.gz → coloc_eqtl_ldwin/final_clump/ |
| `ldwin/clumped_leads/s1_sig_pairs.py` | Extracts the GTEx v11 significant eQTL variants of every eGene in the coloc analysis set, per tissue. | GTEx v11 eQTLs.signif_pairs parquet → coloc_eqtl_ldwin/sens_clump/sig/<tissue>.tsv.gz |
| `ldwin/clumped_leads/s2_clump_nscc.R` | Clumps each eGene's significant eQTL variants recursively (PLINK r2 >= 0.1) and finds each lead's LD partners; the final run calls it from v3m/03_ld_clump_nscc.pbs. | sig/<T>.tsv.gz, LD_ref_v9_plink19 bfiles and donors → <outdir>/<T>/{leads.tsv, lead_partners.tsv} |
| `ldwin/clumped_leads/s2_clump_nscc.pbs` | NSCC array job (one tissue per index) that ran s2_clump_nscc.R in the earlier run. | tissues.txt → s2_clump_nscc.R outputs |
| `ldwin/clumped_leads/s3_pairs.R` | Pairs each eGene signal with every clumped eQTL lead of its gene and builds the lead windows. | pairs.tsv.gz, nscc_out/<T>/{leads,lead_partners}.tsv, tissues.txt → sens_clump/{pairs.tsv.gz, eqtl_windows.tsv} |
| `ldwin/clumped_leads/s4_extract_eqtl.py` | Cuts the GTEx v11 eQTL all-pairs parquet files to the clumped-lead windows, one tissue per call. | GTEx v11 all-pairs parquet, sens_clump/eqtl_windows.tsv → sens_clump/eqtl_subset/<T>.tsv.gz |
| `ldwin/clumped_leads/s5_coloc_tissue_v2.R` | Runs coloc for every eGene signal x clumped eQTL lead pair of one tissue, keeping all coloc columns. | sens_clump/pairs.tsv.gz, sens_clump/eqtl_subset/<T>.tsv.gz, analysis set, LOO scans → sens_clump/per_tissue_v2/<T>.rds |
| `ldwin/clumped_leads/s6_merge_compare.R` | Merges the clumped-lead coloc of all tissues (best call per signal) and compares it with the top-variant and ± 200 kb calls. | sens_clump per-tissue results, coloc_eqtl_signals.tsv.gz → sens_clump/signals_sens.tsv.gz (read by 08_final_clump.R), gene_tissue_sens.tsv.gz, compare_*.tsv |

## v3/ (earlier run)

| file | what it does | inputs → outputs |
|---|---|---|
| `v3/00_build_analysis_set.R` | Builds the analysis set of the earlier run (units v3, 14,783 signals), whose signal SNPs got the LD windows of ldwin/. | final_nogate_v3 (read_units_final, UNITSET v3), conditional_final_v3, GTEx v11 eGenes, ld_v3 donors → results/coloc_eqtl_v3/{analysis_set.rds, eqtl_windows.tsv} |
| `v3/01_extract_eqtl.py` | ± 200 kb run, step 1: cuts the GTEx v11 eQTL all-pairs to the genes and signal ± 200 kb windows of one tissue. | coloc_eqtl_v3/eqtl_windows.tsv, GTEx v11 all-pairs parquet → coloc_eqtl_v3/eqtl_subset/<tissue>.tsv.gz |
| `v3/02_coloc_tissue.R` | ± 200 kb run, step 2: puQTL x eQTL coloc for one tissue. | analysis set, eQTL subset → coloc_eqtl_v3/per_tissue/ |
| `v3/03_merge_classify.R` | ± 200 kb run, step 3: merges the tissues and classifies signals, units and gene-tissue pairs. | per_tissue/ → coloc_eqtl_v3/coloc_eqtl_{signals,units,gene_tissue}.tsv.gz |
| `v3/04_classify_all_shared.R` | ± 200 kb run, step 4: gene-tissue classes under the all-shared rule; its tables are read by ldwin/07_merge_classify.R. | step 3 outputs → coloc_eqtl_v3/all_shared_rule/ |
