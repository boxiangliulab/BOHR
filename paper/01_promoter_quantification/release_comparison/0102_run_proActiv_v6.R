library(proActiv)
library(stringr)
library(SummarizedExperiment)

out_dir <- "/home/users/nus/e1327989/scratch/gtex_v6/proActiv_results/by_SMTSD"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ncores_use <- 48
min_n <- 20

# ── Promoter annotation (gencode v19 official, with chr prefix) ──
annot_file <- "/home/users/nus/e1327989/share_data/reference/hg38/gencode.v19.annotation.gtf.gz"
annot_rdata <- "/home/users/nus/e1327989/scratch/gtex_v6/proActiv_results/promoterAnnotation.gencode.v19.RData"

if (file.exists(annot_rdata)) {
  message("Loading cached promoter annotation...")
  load(annot_rdata)
} else {
  message("Preparing promoter annotation from gencode v19...")
  promoterAnnotation.gencode.v19 <- preparePromoterAnnotation(
    file = annot_file,
    species = "Homo_sapiens"
  )
  dir.create(dirname(annot_rdata), recursive = TRUE, showWarnings = FALSE)
  save(promoterAnnotation.gencode.v19, file = annot_rdata)
  message("Saved: ", annot_rdata)
}

# ── Sample metadata (v7 metadata covers v6 samples) ──
sample.df <- read.delim(
  "/data/projects/11003054/share/data/gtex/gtex_v7/GTEx_v7_Annotations_SampleAttributesDS.txt",
  as.is = TRUE, header = TRUE, row.names = 1
)

# ── Junction files ──
sj_dir <- "/home/users/nus/e1327989/scratch/gtex_v6/filtered_junctions"
files <- Sys.glob(file.path(sj_dir, "*.SJ.out.tab"))
stopifnot(length(files) > 0)

# Remove empty files (samples with no junctions passing filter)
file_sizes <- file.info(files)$size
empty_mask <- is.na(file_sizes) | file_sizes == 0
if (any(empty_mask)) {
  message("Removing ", sum(empty_mask), " empty SJ files")
  files <- files[!empty_mask]
}

# Extract sample IDs from filenames
ids <- sub("\\.SJ\\.out\\.tab$", "", basename(files))

# Match to metadata
att <- sample.df[ids, , drop = FALSE]
matched <- !is.na(att$SMTSD)
message("Matched samples: ", sum(matched), " / ", length(ids))

# Keep only matched samples
files <- files[matched]
ids <- ids[matched]
att <- att[matched, , drop = FALSE]
att$SMTSD <- gsub(" - ", "_", att$SMTSD)

save(att, file = file.path(dirname(out_dir), "att_sample_info_v6.RData"))

condition <- att$SMTSD
tissues <- sort(unique(condition))

message("Total samples: ", length(files))
message("Total tissues: ", length(tissues))

# Only run one tissue (set to NULL to run all)
run_tissue <- "Whole Blood"

# ── Run proActiv per tissue ──
for (tissue in tissues) {
  if (!is.null(run_tissue) && tissue != run_tissue) next

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
  con_log <- file(out_log, open = "wt")
  con_err <- file(out_err, open = "wt")

  tryCatch({
    sink(file = con_log, split = TRUE)
    sink(file = con_err, type = "message")

    res <- proActiv(
      files = files[idx],
      promoterAnnotation = promoterAnnotation.gencode.v19,
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
    try(sink(type = "message"), silent = TRUE)
    try(sink(), silent = TRUE)
    try(close(con_log), silent = TRUE)
    try(close(con_err), silent = TRUE)
  })
}

message("All done. Results in: ", out_dir)
