## Authored for the final run (not an original script): proActiv per tissue on the GTEx v11 STAR junction files.
## 0102_run_proActiv_tissue.R with every junction (no read floor): the same GENCODE v47 annotation (built with the same
## preparePromoterAnnotation call on the same GTF), the same samples and tissue labels (0102 att_sample.info.RData),
## the same sample order, proActiv() call and complete.cases() step. One tissue per call.
## Usage: Rscript run_proActiv_arm.R <tissue index 1..57> <ncores> <floor|nofloor>; nofloor (every junction) is the final run
suppressMessages({library(proActiv); library(SummarizedExperiment)})
a <- commandArgs(TRUE); k <- as.integer(a[1]); nc <- as.integer(a[2])
arm <- a[3]; SJ <- if (arm == "floor") "/scratch/users/nus/e1327989/junction_floor/filtered_junctions" else "/home/users/nus/e1327989/gtex/gtex_v11/junctions/Star_out_tab"
out_dir <- file.path("/scratch/users/nus/e1327989/junction_floor", arm, "proActiv")
ann <- readRDS("/home/users/nus/e1327989/YuhaoDong/junction_floor/out_test/promoterAnnotation.gencode.v47.rds")
load("/home/users/nus/e1327989/gtex/gtex_v11/att_sample.info.RData")          # att, written by 0102
files <- sort(file.path(SJ, paste0(rownames(att), ".SJ.out.tab.gz")))      # 0102 took Sys.glob order
ids <- sub("\\.SJ\\.out\\.tab\\.gz$", "", basename(files))
stopifnot(all(file.exists(files)), setequal(ids, rownames(att)))
condition <- gsub(" - ", "_", att[ids, "SMTSD"])
tissues <- sort(unique(condition)); tissues <- tissues[vapply(tissues, function(t) sum(condition == t) >= 20, TRUE)]
stopifnot(length(tissues) == 57)
tissue <- tissues[k]; idx <- which(condition == tissue); n_t <- length(idx)
tissue_safe <- gsub("[^A-Za-z0-9_]+", "_", tissue)
out_rds <- file.path(out_dir, paste0("proActiv_", tissue_safe, "_n", n_t, ".rds"))
if (file.exists(out_rds)) { message("[DONE] exists: ", out_rds); quit(save = "no") }
message("[RUN] ", tissue, " (n=", n_t, ") ", Sys.time())
res <- proActiv(files = files[idx], promoterAnnotation = ann, condition = condition[idx], ncores = nc)
if ("promoterCounts" %in% names(assays(res))) res <- res[complete.cases(assays(res)$promoterCounts), ]
saveRDS(res, out_rds)
message("[OK] ", out_rds, " ", Sys.time())
