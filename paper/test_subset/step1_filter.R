## Subset test, step 1 (authored for the repository): the promoter filter of SI Section 1 on the promoters of one
## chromosome, with BOHR::relative_activity() on the stored outputs of 0103_proActiv_data_process_tissue.R (proActiv's
## annotation, classes, absolute and relative activity of the external promoters; GTEx v11, 22,865 samples, 57 tissues).
## The filter works promoter by promoter and gene by gene, so one chromosome gives exactly the paper's result for its genes.
## Usage: Rscript step1_filter.R <GTEx v11 folder> <chr> <output prefix>
suppressMessages({ library(BOHR); library(Matrix) })
a <- commandArgs(TRUE); v11 <- a[1]; chr <- a[2]; pre <- a[3]
load(file.path(v11, "filter_5_rdata.RData"))                # rdata: one row per promoter, proActiv classes per tissue
load(file.path(v11, "att_sample.info.RData"))               # att: sample attributes (SMTSD = tissue)
rd <- rdata[as.character(rdata$seqnames) == chr, ]
internal <- as.logical(as.character(rd$internalPromoter))
promoters <- data.frame(id = rownames(rd), gene = as.character(rd$geneId), chr = chr, pos = as.numeric(rd$start),
                        strand = as.character(rd$strand), internal = internal, stringsAsFactors = FALSE)
cc <- grep("class", colnames(rd), value = TRUE)
cls <- as.matrix(as.data.frame(lapply(rd[, cc], as.character))); dimnames(cls) <- list(rownames(rd), cc)
ext <- promoters$id[!internal]
rel <- readRDS(file.path(v11, "trel_att_filter_In.rds"))$rel
samples <- sort(intersect(rownames(att), rownames(rel)))
dense <- function(m) { m <- as.matrix(t(m[samples, ext, drop = FALSE])); storage.mode(m) <- "double"; m }
R_pa <- dense(rel); rm(rel); invisible(gc())
A <- dense(readRDS(file.path(v11, "tabs_att_filter_In.rds"))$abs); invisible(gc())
pa <- BOHR:::new_proactiv(A, R_pa, promoters, cls)           # the object from_proactiv() builds from proActiv results
ra <- relative_activity(pa, group = as.character(att[samples, "SMTSD"]))
saveRDS(list(relative = ra$relative, steps = ra$steps, promoters = ra$promoters, chr = chr), paste0(pre, "_relative_activity.rds"))
print(ra$steps, row.names = FALSE)
cat("kept:", nrow(ra$relative), "promoters in", length(unique(ra$promoters$gene[ra$promoters$kept])), "genes\n")
