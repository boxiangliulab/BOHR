# Colocalization with cis-eQTLs

**Methods §8; Supplementary Methods §8.** Each independent puQTL signal (its leave-one-out scan when the unit has more
than one) is tested by coloc against the GTEx v11 eQTLs of its gene in the same tissue, over the union of the two LD
windows (farthest partner with r² ≥ 0.1, ± 50 kb), against the GTEx top eQTL variant and against every clumped eQTL lead;
the best call per signal gives the gene–tissue classes of Fig. 4E.

**Run order.** The table was built on the units v3 signals (`ldwin/`, 2026-10-12) and updated for the signals that
changed with the mirror rule (`v3m/`) and the two refitted units (`v3c/`); unchanged rows are carried.
1. `v3/00_build_analysis_set.R` (signals and eGenes) → `ldwin/00_leads.R` → PLINK LD windows → `ldwin/02_windows.R`,
   `03_rarelead_windows.R`, `04_pairs.R` → `05_extract_eqtl.py` → `06_coloc_tissue.R` per tissue → `07_merge_classify.R`.
2. Clumped eQTL leads: `ldwin/clumped_leads/s1_sig_pairs.py` → `s2_clump_nscc.R` → `s3_pairs.R` → `s4_extract_eqtl.py` →
   `s5_coloc_tissue_v2.R`, then `ldwin/08_final_clump.R` (best call per signal over all leads).
3. `v3m/`: `01_analysis_set_pairs.R` → `02_new_inputs.R`, `02b_sig_pairs.py` → `04_pairs_new_all.R` → `04b_extract_eqtl.py` →
   `05_coloc_tissue_v3m.R` → `06_final_v3m.R`.
4. `v3c/`: `c1_analysis_set_pairs_v3c.R` → `c5_coloc_tissue_v3c.R` → `c6_final_v3c.R` → `c7_textclean_v3c.R` (the
   paper's table: 1,600 gene–tissue pairs with a distinct signal, 1,663 without a significant eQTL, 1,721 in which every signal colocalized).

`coloc_final_lib.R` holds the shared coloc functions (priors p1 = p2 = 10⁻⁴, p12 = 10⁻⁵); the copy in `v3c/` differs only
in accepting the v3c unit set. coloc is not part of the BOHR package.

`run.sh` lists the commands in the order they were run (paths of the original runs).


## Files

| file | what it does | inputs → outputs |
|---|---|---|
| `coloc_final_lib.R` | Provides the shared coloc functions (load_marginal, coloc_signal, priors) for eQTL and GWAS coloc. | sourced library; reads traits_51_v3eur.txt at load (GWAS sample-size config optional) |
| `config/traits_51_v3eur.txt` | Lists the 51 GWAS trait files; coloc_final_lib.R reads it when sourced, so the eQTL coloc scripts need it too. | read at load by coloc_final_lib.R |

## v3/ and ldwin/

| file | what it does | inputs → outputs |
|---|---|---|
| `v3/00_build_analysis_set.R` | Builds the 2026-10-05 coloc analysis set (units v3, 14,783 signals) behind the 11,983 carried rows of the eQTL coloc table. | final_nogate_v3 (read_units_final, UNITSET v3), conditional_final_v3, GTEx v11 eGenes, ld_v3 donors → results/coloc_eqtl_v3/{analysis_set.rds, eqtl_windows.tsv} |
| `ldwin/00_leads.R` | Lists the puQTL signal SNPs and GTEx top eQTL variants that need LD windows, per tissue and chromosome. | coloc_eqtl_v3/analysis_set.rds, v11 eGenes → coloc_eqtl_ldwin/{leads.tsv, eqtl_leads.tsv, ldwin_inputs/snps/<T>.<chr>.txt} |
| `ldwin/02_windows.R` | Assembles LD windows (farthest r2 >= 0.1 partner +/- 50 kb) per lead and per signal-eQTL pair. | leads.tsv, ldwin_out/*.win.tsv, analysis_set.rds, eqtl_leads.tsv → coloc_eqtl_ldwin/{lead_windows.tsv, pair_windows.tsv.gz} |
| `ldwin/03_rarelead_windows.R` | Builds LD windows for eQTL leads absent from the MAF >= 0.05 panel (option A': panel partners only). | lead_windows.tsv, rarelead/*.ld, panel/panel_snps.tsv.gz, pair_windows.tsv.gz → rarelead_windows.tsv, rarelead_pairs.tsv |
| `ldwin/04_pairs.R` | Pairs each analysis signal with its gene's GTEx top eQTL variant and builds the LD-window intersections and unions. | coloc_eqtl_v3/analysis_set.rds, eqtl_leads.tsv, lead_windows.tsv, rarelead_windows.tsv → coloc_eqtl_ldwin/{pairs.tsv.gz, eqtl_windows.tsv} |
| `ldwin/05_extract_eqtl.py` | Cuts the GTEx v11 eQTL all-pairs parquet to the genes and windows the coloc needs, per tissue. | GTEx v11 all-pairs parquet, coloc_eqtl_ldwin/eqtl_windows.tsv → coloc_eqtl_ldwin/eqtl_subset/<tissue>.tsv.gz |
| `ldwin/06_coloc_tissue.R` | Runs puQTL x eQTL coloc for one tissue over each pair's LD-window union (GTEx top-variant leads). | analysis_set.rds, pairs.tsv.gz, eqtl_subset/<T>.tsv.gz, conditional_final_v3/loo → coloc_eqtl_ldwin/per_tissue/<T>.rds |
| `ldwin/07_merge_classify.R` | Merges the top-variant coloc runs and classifies each signal; only its non-eGene rows reach the final eQTL coloc table. | coloc_eqtl_v3/analysis_set.rds, per_tissue/*.rds, coloc_eqtl_v3/all_shared_rule/ (comparison) → coloc_eqtl_ldwin/{coloc_eqtl_signals*,gene_tissue_classes}.tsv.gz |

## ldwin/clumped_leads/

| file | what it does | inputs → outputs |
|---|---|---|
| `ldwin/clumped_leads/s1_sig_pairs.py` | Extracts the GTEx v11 significant eQTL variants of every eGene in the coloc analysis set, per tissue. | GTEx v11 eQTLs.signif_pairs parquet → coloc_eqtl_ldwin/sens_clump/sig/<tissue>.tsv.gz |
| `ldwin/clumped_leads/s2_clump_nscc.R` | Clumps each eGene's significant eQTL variants recursively (PLINK r2 >= 0.1) and finds each lead's LD partners. | sens_clump/sig/<T>.tsv.gz, LD_ref_v9_plink19 bfiles and donors → <outdir>/<T>/{leads.tsv, lead_partners.tsv} |
| `ldwin/clumped_leads/s3_pairs.R` | Pairs each eGene signal with every clumped eQTL lead of its gene and builds the lead windows. | pairs.tsv.gz, nscc_out/<T>/{leads,lead_partners}.tsv, tissues.txt → sens_clump/{pairs.tsv.gz, eqtl_windows.tsv} |
| `ldwin/clumped_leads/s4_extract_eqtl.py` | Cuts the GTEx v11 eQTL all-pairs parquet files to the clumped-lead windows, one tissue per call. | GTEx v11 all-pairs parquet, sens_clump/eqtl_windows.tsv → sens_clump/eqtl_subset/<T>.tsv.gz |
| `ldwin/clumped_leads/s5_coloc_tissue_v2.R` | Runs coloc for every eGene signal x clumped eQTL lead pair of one tissue, keeping all coloc columns. | sens_clump/pairs.tsv.gz, sens_clump/eqtl_subset/<T>.tsv.gz, analysis set, LOO scans → sens_clump/per_tissue_v2/<T>.rds |
| `ldwin/08_final_clump.R` | Builds the 2026-10-12 final eQTL coloc table (best call per signal over the clumped eQTL leads) whose rows v3m and v3c carry. | sens_clump/{per_tissue_v2, pairs.tsv.gz}, coloc_eqtl_signals.tsv.gz, gene_tissue_classes, panel/panel_snps.tsv.gz → coloc_eqtl_ldwin/final_clump/ |

## v3m/

| file | what it does | inputs → outputs |
|---|---|---|
| `v3m/01_analysis_set_pairs.R` | Builds the mirror-v3 eQTL-coloc analysis set (14,801 signals, carried or new) and the pairs of the new rows. | units v3m, conditional_final_v3m, v11 eGenes, ld_v3/donors, 10-12 analysis set and pairs → coloc_eqtl_v3m/{analysis_set.rds, pairs_*_new.tsv.gz} |
| `v3m/02_new_inputs.R` | Lists the LD windows, clumpings and eQTL extractions still needed for the new rows of the mirror-rule-v3 coloc. | coloc_eqtl_v3m/{analysis_set.rds, missing.tsv}, ldwin windows and leads, sens_clump pairs, GTEx v11 eGenes → coloc_eqtl_v3m/new_inputs/ |
| `v3m/02b_sig_pairs.py` | Extracts the GTEx v11 significant eQTL variants of the eGenes to clump for the new v3m rows. | new_inputs/clump_genes.tsv, GTEx v11 eQTLs.signif_pairs → new_inputs/sig/<T>.tsv.gz |
| `v3m/04_pairs_new_all.R` | Completes the top-variant and clumped-lead pairs of every new v3m row from the step-1 pairs and the NSCC windows and clumps. | analysis_set.rds, pairs_*_new.tsv.gz, nscc/{out,clump_out}, new_inputs → pairs_{top,clump}_new_all.tsv.gz, eqtl_windows_{top,clump}_new.tsv |
| `v3m/04b_extract_eqtl.py` | Cuts the GTEx v11 eQTL all-pairs parquet files to the new v3m windows (top and clumped lead sets). | eqtl_windows_<kind>_new.tsv, GTEx v11 all-pairs parquet → eqtl_subset_<kind>/<T>.tsv.gz |
| `v3m/05_coloc_tissue_v3m.R` | Runs coloc for the new v3m rows of one tissue against the GTEx top variant and every clumped eQTL lead. | analysis_set.rds, pairs_*_new_all, eqtl_subset_*, conditional_final_v3m/loo → per_tissue_{top,clump}/<T>.rds |
| `v3m/06_final_v3m.R` | Merges carried and new coloc rows, takes the best call per signal and builds the Fig. 4E gene-tissue classes. | coloc_eqtl_v3m/{analysis_set.rds, per_tissue_*}, coloc_eqtl_ldwin/final_clump, panel_snps.tsv.gz → coloc_eqtl_v3m/final_clump_v3m/ |

## v3c/

| file | what it does | inputs → outputs |
|---|---|---|
| `v3c/coloc_final_lib.R` | Provides the shared puQTL x GWAS coloc functions; live copy, identical to the cited one except that UNITSET also accepts v3c (Figure 4 conditional refit). | sourced by every coloc script; UNITSET=v3c reads final_nogate_v3c and conditional_final_v3c and writes results_v3c/, coloc_ldunion_20261012/run_v3c/ |
| `v3c/c1_analysis_set_pairs_v3c.R` | Builds the v3c coloc analysis set (14,800 signals, 14,798 carried from v3m) and the pairs of the two refitted units. | final_nogate_v3c, conditional_final_v3c, units_to_run.tsv, v3m and ldwin windows, pairs, leads, eGenes → coloc_eqtl_v3c/{analysis_set.rds, pairs_*_new_all.tsv.gz} |
| `v3c/c5_coloc_tissue_v3c.R` | Runs coloc for the new rows of the two refitted units (Artery_Aorta, Brain_Putamen_basal_ganglia). | coloc_eqtl_v3c analysis set and pairs, conditional_final_v3c/loo, ldwin and v3m eQTL subsets → coloc_eqtl_v3c/per_tissue_{top,clump}/<T>.rds |
| `v3c/c6_final_v3c.R` | Builds the final v3c eQTL coloc table, the gene-tissue classes and the Fig. 4E counts (1,600 / 1,663 / 1,721 / 483 / 2,830). | final_clump_v3m, coloc_eqtl_v3c/per_tissue_*, analysis_set.rds, panel_snps.tsv.gz → coloc_eqtl_v3c/final_clump_v3c/ |
| `v3c/c7_textclean_v3c.R` | Restores the carried rows of the v3c signal-level coloc files verbatim from v3m after c6 rewrote subnormal posteriors. | final_clump_v3m, final_clump_v3c signal files → final_clump_v3c/coloc_eqtl_signals{,_flags}.tsv.gz |
