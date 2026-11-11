## Parity of relative_activity() with the paper's promoter filter (0103_proActiv_data_process_tissue.R, SI Section 1).
##
## Inputs (the stored outputs of 0103, GTEx v11 run; one folder):
##   filter_5_rdata.RData            rdata: rowData of every promoter (proActiv annotation + <tissue>.class per tissue)
##   att_sample.info.RData           att: GTEx sample attributes (SMTSD = tissue)
##   trel_att_filter_In.rds          list(meta, rel): proActiv's relative activity, samples x external promoters
##   tabs_att_filter_In.rds          list(meta, abs): proActiv's absolute activity, samples x external promoters
##   Promoter_info_filter_5.RData    info: the paper's final promoter set
##   filter_5_recal_att_filter.RData recal_att_f: the paper's relative activity (external denominator), final set
## The object from_proactiv() would build from the per-tissue proActiv results is assembled from these files with the
## package's own constructor, then relative_activity() runs with its defaults and group = SMTSD.
##
## Usage: Rscript 12_relative_activity_parity.R <input folder> <output folder>
suppressMessages({ library(BOHR); library(Matrix) })
a <- commandArgs(TRUE); inp <- a[1]; out <- a[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time(); say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

load(file.path(inp, "filter_5_rdata.RData"))
load(file.path(inp, "att_sample.info.RData"))
internal <- as.logical(as.character(rdata$internalPromoter))
promoters <- data.frame(id = rownames(rdata), gene = as.character(rdata$geneId), chr = as.character(rdata$seqnames),
                        pos = as.numeric(rdata$start), strand = as.character(rdata$strand), internal = internal,
                        stringsAsFactors = FALSE)
cc <- grep("class", colnames(rdata), value = TRUE)
cls <- as.matrix(as.data.frame(lapply(rdata[, cc], as.character)))
dimnames(cls) <- list(rownames(rdata), sub("^.*\\.rds\\.", "", sub("\\.class$", "", cc)))
ext <- promoters$id[!internal]
say("annotation:", nrow(promoters), "promoters,", length(ext), "external,", ncol(cls), "class columns")

rel <- readRDS(file.path(inp, "trel_att_filter_In.rds"))$rel
abs <- readRDS(file.path(inp, "tabs_att_filter_In.rds"))$abs
samples <- intersect(rownames(att), rownames(rel))                    # 0103: merge(meta_all, rel_all, by = "row.names")
stopifnot(identical(sort(colnames(rel)), sort(ext)), setequal(rownames(abs), rownames(rel)))
samples <- sort(samples)
say("samples:", length(samples), "of", nrow(rel), "in the relative matrix")
to_dense <- function(m) { m <- as.matrix(t(m[samples, ext])); storage.mode(m) <- "double"; m }
R_pa <- to_dense(rel); rm(rel); gc()
A <- to_dense(abs); rm(abs); gc()
say("matrices:", nrow(A), "x", ncol(A))
pa <- BOHR:::new_proactiv(A, R_pa, promoters, cls)
rm(A, R_pa); gc()
group <- as.character(att[samples, "SMTSD"])
say("tissues:", length(unique(group)))

ra <- relative_activity(pa, group = group)
say("relative_activity() done;", nrow(ra$relative), "promoters kept")
print(ra$steps, row.names = FALSE)

## expected (SI Section 1; tables/source_data/figure1/intermediate/filter_reproduction_counts.tsv; C06_promoter_cascade_output.txt)
exp <- data.frame(step = c("annotated", "external", "active", "multi_promoter_gene", "nonzero", "missing", "mean", "min_promoters"),
                  promoters = c(134842, 92557, 41203, 31519, 26508, 25047, 21016, 10826),
                  genes = c(54117, 52037, 28268, NA, NA, NA, NA, 4744),
                  genes_multi = c(24710, 18475, 8284, NA, NA, NA, NA, 4744))
cmp <- merge(exp, ra$steps, by = "step", suffixes = c("_paper", "_package"), sort = FALSE)
cmp$match <- with(cmp, promoters_paper == promoters_package & (is.na(genes_paper) | genes_paper == genes_package) &
                       (is.na(genes_multi_paper) | genes_multi_paper == genes_multi_package))
print(cmp, row.names = FALSE)

load(file.path(inp, "Promoter_info_filter_5.RData"))
set_ok <- setequal(rownames(ra$relative), rownames(info))
say("final set equals Promoter_info_filter_5:", set_ok, "(", length(setdiff(rownames(ra$relative), rownames(info))), "extra,",
    length(setdiff(rownames(info), rownames(ra$relative))), "missing )")

load(file.path(inp, "filter_5_recal_att_filter.RData"))               # recal_att_f: samples x (sample, SMTS, SMTSD, promoters)
P <- as.matrix(recal_att_f[, -(1:3)]); storage.mode(P) <- "double"
ids <- intersect(rownames(ra$relative), colnames(P)); ss <- intersect(samples, rownames(P))
X <- ra$relative[ids, ss]; Y <- t(P[ss, ids])
na_same <- identical(is.na(X), is.na(Y))
d <- abs(X - Y); mx <- max(d, na.rm = TRUE)
Y2 <- Y; Y2[is.nan(Y2)] <- NA
say("values:", length(ids), "promoters x", length(ss), "samples; NA pattern equal:", na_same, "; max |diff|:", signif(mx, 3),
    "; identical (NaN as NA):", identical(unname(X), unname(Y2)))

why <- table(ra$promoters$reason, useNA = "ifany")
print(why)
res <- data.frame(check = c("steps match SI", "final set = paper's", "NA pattern equal", "max abs diff", "identical values",
                            "samples", "tissues", "minutes"),
                  value = c(all(cmp$match), set_ok, na_same, mx, identical(unname(X), unname(Y2)), length(samples),
                            length(unique(group)), round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)))
write.table(cmp, file.path(out, "12_steps_vs_paper.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res, file.path(out, "12_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ra$promoters, gzfile(file.path(out, "12_promoters.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
print(res, row.names = FALSE)
cat("peak memory (gc, MB):", sum(gc()[, 6]), "\n")
