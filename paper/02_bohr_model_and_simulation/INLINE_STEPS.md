# Figure 2: steps that ran as one-line commands (no script file)

Two steps ran as `Rscript -e` one-liners during the Figure 2 analysis. There is no script file for them, so the code map marks them
NO_CODE. They are copied here verbatim from the record of the commands.

The three commands were later re-run into a scratch folder, reading the same design tables (`sim_50` of each batch).
The output was byte-identical (`cmp`) to the ID files used by every calibration run.

## 1. Simulated activity: 150 calibration genes

Run in `new_simulation/newmethod_bohr2/`. Writes `m1_nothr_g2000/ids_calib_random150.txt`.

```r
set.seed(20261001); writeLines(sort(sprintf("gene_%04d", sample(2000, 150))), "m1_nothr_g2000/ids_calib_random150.txt")
```

## 2. Simulated RNA-seq, batch 1: 150 calibration genes, both promoters

Run in `new_simulation/newmethod_bohr2/`. Writes `m2_nothr/b1/ids_calib_random150genes.txt` and `genes_calib_random150.txt`.

```r
tx <- fread("m2_nothr/b1/sim_<N>/test_output/sim_tx_info_extended.txt")[!is.na(promoter_id)]
g <- unique(tx[, .(gene_id, id = paste0("PromoterID_", promoter_id))]); genes <- sort(unique(g$gene_id))
set.seed(20261001); pick <- sample(genes, 150); ids <- sort(g[gene_id %in% pick, id])
writeLines(ids, "m2_nothr/b1/ids_calib_random150genes.txt"); writeLines(sort(pick), "m2_nothr/b1/genes_calib_random150.txt")
```

The gene list is the same at every N, because genes and transcripts are drawn under `gene_seed`.

## 3. Simulated RNA-seq, batch 2

Batch 2 has a script: `newmethod_bohr2/scripts/runner_m2nothr_local.sh`, lines 10–11. It is the same code with `set.seed(20261002)`.

The analysis reads both promoters of each calibration gene. `prepare_data_v3.R` then drops the minor member of each exact mirror
pair, lines 85–86 and 98.
