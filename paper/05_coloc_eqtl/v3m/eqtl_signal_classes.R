## Per-signal eQTL classes: each conditionally independent puQTL signal is shared / distinct /
## no eQTL / unresolved, without gene x tissue pooling. The gene x tissue rule of 07_merge_classify.R / 08_final_clump.R
## (classify()) applied to one signal; the signal's call is its best call over its eQTL-lead pairs (08_final_clump.R).
##   no eQTL     the gene is not a GTEx v11 eGene in the tissue (qval > 0.05 or not tested)
##   shared      call H4
##   distinct    call H3, different locus (LD windows disjoint), or puQTL-specific (H1)
##   unresolved  ambiguous, eQTL-specific (H2), neither (H0), not assessable, no eQTL rows
## Usage: Rscript eqtl_signal_classes.R <final coloc_eqtl_signals.tsv.gz> <out_dir>
## Output: signal_classes.tsv.gz (one row per analysis signal: keys, call, PP.H4, egene, class), signal_class_counts.tsv
##   (all, by type, by tissue). Aggregate numbers only are printed.
suppressMessages(library(data.table))
a <- commandArgs(TRUE); F <- fread(cmd = paste("zcat", a[1])); OUT <- a[2]; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
F[, class := fcase(egene %in% FALSE | is.na(egene), "no eQTL",
                   call == "shared (H4)", "shared",
                   call %in% c("distinct (H3)", "different locus (LD windows disjoint)", "puQTL-specific (H1)"), "distinct",
                   default = "unresolved")]
lv <- c("shared", "distinct", "no eQTL", "unresolved"); F[, class := factor(class, levels = lv)]
stopifnot(!anyNA(F$class))
fwrite(F[, .(tissue, gene, gene_base, promoter, type, k, signal_snp, stat, egene, call, PP.H4, class)], file.path(OUT, "signal_classes.tsv.gz"), sep = "\t")
cnt <- function(d, by) { x <- d[, .N, by = c(by, "class")]; x <- dcast(x, as.formula(paste(if (length(by)) paste(by, collapse = "+") else ".", "~ class")), value.var = "N", fill = 0)
  for (k in setdiff(lv, names(x))) x[, (k) := 0L]; x[, signals := rowSums(.SD), .SDcols = lv]; x }
C <- rbind(cnt(F, character(0))[, scope := "all"][, level := "all"], cnt(F, "type")[, scope := "type"][, level := type][, type := NULL],
           cnt(F, "tissue")[, scope := "tissue"][, level := tissue][, tissue := NULL], fill = TRUE)
C[, "." := NULL]; setcolorder(C, c("scope", "level", "signals", lv))
fwrite(C, file.path(OUT, "signal_class_counts.tsv"), sep = "\t")
print(C[scope != "tissue"]); cat("calls by class:\n"); print(F[, .N, by = .(class, call)][order(class, -N)])
