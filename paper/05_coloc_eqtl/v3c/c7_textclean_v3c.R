## c7_textclean_v3c.R (2026-11-06): c6 writes carried rows through fread/fwrite, which rewrites posteriors below ~1e-307 (subnormal values,
## e.g. 1.6997e-308 -> 1.9624e-308). For every carried row of final_clump_v3c's signal-level files, restore the v3m line verbatim except
## the `source` field; new rows (conditional refit) keep c6's line. Gene-level files carry no posteriors and are left as written by c6.
suppressPackageStartupMessages(library(data.table))
R2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"
IN3 <- file.path(R2, "results/coloc_eqtl_v3m/final_clump_v3m"); FO <- file.path(R2, "results/coloc_eqtl_v3c/final_clump_v3c")
for (f in c("coloc_eqtl_signals.tsv.gz", "coloc_eqtl_signals_flags.tsv.gz")) {
  a <- readLines(gzfile(file.path(IN3, f))); b <- readLines(gzfile(file.path(FO, f)))
  ha <- strsplit(a[1], "\t")[[1]]; hb <- strsplit(b[1], "\t")[[1]]; stopifnot(identical(ha, hb))
  ka <- c("tissue", "promoter", "type", "k", "signal_snp"); ia <- match(ka, ha); js <- match("source", hb)
  sa <- strsplit(a[-1], "\t", fixed = TRUE); sb <- strsplit(b[-1], "\t", fixed = TRUE)
  keyA <- vapply(sa, function(v) paste(v[ia], collapse = "|"), ""); keyB <- vapply(sb, function(v) paste(v[ia], collapse = "|"), "")
  m <- match(keyB, keyA); carried <- vapply(sb, function(v) v[js] == "carried (final_clump_v3m)", TRUE) & !is.na(m)
  out <- b; n <- 0L
  jg <- match(c("gene_class_fig4e", "gene_class_detail"), hb); jg <- jg[!is.na(jg)]          # gene-level class must not be restored stale
  for (j in which(carried)) { v <- sa[[m[j]]]; stopifnot(identical(v[jg], sb[[j]][jg])); v[js] <- sb[[j]][js]; l <- paste(v, collapse = "\t")
    if (l != b[j + 1]) n <- n + 1L; out[j + 1] <- l }
  con <- gzfile(file.path(FO, f), "w"); writeLines(out, con); close(con)
  cat(f, ": rows", length(b) - 1, "; carried", sum(carried), "(lines restored from v3m:", n, "); new", sum(!carried), "\n") }
