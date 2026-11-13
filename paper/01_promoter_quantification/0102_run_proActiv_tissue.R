library(proActiv)
library(stringr)
library(SummarizedExperiment)

out_dir <- "~/yuhaodong_colo/GTEx/v11_filter_5_recal/by_SMTSD"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ncores_use <- 16
min_n <- 20

# promoter annotation
promoterAnnotation.gencode.v47 <- preparePromoterAnnotation(
  file = "~/yuhaodong_colo/reference/gencode.v47.annotation.gtf.gz",
  species = "Homo_sapiens"
)
save(promoterAnnotation.gencode.v47,
     file = "~/yuhaodong_colo/GTEx/v11_filter_5_recal/promoterAnnotation.gencode.v47.RData")

# metadata
sample.df  <- read.delim("~/yuhaodong_colo/GTEx/Rawdata/metadata/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt",
                         as.is = TRUE, header = TRUE, row.names = 1)

files <- Sys.glob("~/yuhaodong_colo/GTEx/Rawdata/gtex_v11/filtered_junctions/*.SJ.out.tab.gz")
stopifnot(length(files) > 0)

ids <- str_split_fixed(str_split_fixed(files, "[/]", n = 10)[, 10], "[.]", n = 2)[, 1]
att <- sample.df[ids, , drop = FALSE]
att$SMTSD <- gsub(" - ", "_", att$SMTSD)

save(att, file = "~/yuhaodong_colo/GTEx/v11_filter_5_recal/att_sample.info.RData")

condition <- att$SMTSD
tissues <- sort(unique(condition))

message("Total samples: ", length(files))
message("Total tissues: ", length(tissues))

for (tissue in tissues) {
  idx <- which(condition == tissue)
  n_t <- length(idx)

  if (n_t < min_n) {
    message("[SKIP] ", tissue, " (n=", n_t, " < ", min_n, ")")
    next
  }

  tissue_safe <- gsub("[^A-Za-z0-9_]+", "_", tissue)
  out_rds <- file.path(out_dir, paste0("proActiv_", tissue_safe, "_n", n_t, ".rds"))
  out_log <- file.path(out_dir, paste0("proActiv_", tissue_safe, "_n", n_t, ".log"))
  out_err <- file.path(out_dir, paste0("proActiv_", tissue_safe, "_n", n_t, ".err"))

  if (file.exists(out_rds)) {
    message("[DONE] exists, skip: ", tissue)
    next
  }

  message("[RUN] ", tissue, " (n=", n_t, ")")

  res <- NULL

  # Open the connection
  con_log <- file(out_log, open = "wt")
  con_err <- file(out_err, open = "wt")

  # Make sure the sink/connection is closed whether the run succeeds or fails
  tryCatch({
    sink(file = con_log, split = TRUE)
    sink(file = con_err, type = "message")

    res <- proActiv(
      files = files[idx],
      promoterAnnotation = promoterAnnotation.gencode.v47,
      condition = condition[idx],
      ncores = ncores_use
    )

    if ("promoterCounts" %in% names(assays(res))) {
      res <- res[complete.cases(assays(res)$promoterCounts), ]
    }

    saveRDS(res, out_rds)
    message("[OK] saved: ", out_rds)

  }, error = function(e) {
    message("[FAIL] tissue=", tissue, " error: ", conditionMessage(e))
  }, finally = {
    # Close the sink (mind the order: close the message sink first, then the output sink)
    try(sink(type = "message"), silent = TRUE)
    try(sink(), silent = TRUE)

    # Close the connection
    try(close(con_log), silent = TRUE)
    try(close(con_err), silent = TRUE)
  })
}

message("All done. Results in: ", out_dir)