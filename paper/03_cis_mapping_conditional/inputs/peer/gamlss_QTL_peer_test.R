#! /usr/bin/Rscript
#====================================================================
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(source("/home/users/nus/e1327989/YuhaoDong/GTEx/gamlss_functions.R"))
#====================================================================

option_list <- list(
  make_option(c("-n", "--name"), type="character", default= "test",
              help="name of output files [default %default]"),
  make_option(c("-t", "--tissue"), type="character", default= "Muscle_Skeletal",
              help="tissue of output files [default %default]"),
  make_option("--chr", type="character", default= "chr22",
              help="[default %default]"),
  make_option("--flank", type="numeric", default= 500000,
              help="flank region around TSS [default %default]"),
  make_option("--thread", type="numeric", default= 1,
              help="fthread number [default %default]"),
  make_option("--peer", type="numeric", default= 1,
              help="the number of peer factor [default %default]"),
  make_option(c("-p", "--peerfile"), type="character", default= NULL,
              help="peer factor results"),
  make_option("--pca", type="numeric", default= 0,
              help="the number of pca factor [default %default]"),

  make_option("--tissue_col", type="character", default= "SMTSD",
              help="tissue_col [default %default]"),

  make_option(c("-o", "--output"), type="character", default= "./", 
              help="To specify output path [default %default]"))

parser <- OptionParser(usage = "%prog [options]", option_list=option_list,
                       description = "\nThe script of gamlss QTL regression.",
                       epilogue="Yuhao Dong, 2024.\n"
)

arguments <- parse_args(parser)



flank=arguments$flank
peerfile=arguments$peerfile
peer_num=arguments$peer
pca_num=arguments$pca
chr=arguments$chr
thread=arguments$thread
tissue_col=arguments$tissue_col
path=arguments$output
tissue=arguments$tissue
output=paste0(arguments$output, "/", arguments$name, ".", arguments$tissue, ".", arguments$chr)

if (!file.exists(path)) {
  dir.create(path, recursive = TRUE)
  cat("path creat:", path, "\n")
}


#
options(stringsAsFactors = FALSE)

print(paste(arguments$name, arguments$chr, "peer", arguments$peer, "pca", arguments$pca, " start", sep=" "))

### read data================================================================================
mc <- read.table("/home/users/nus/e1327989/YuhaoDong/GTEx/v10_cov_pca.txt", sep="\t", header=T, row.names=1)[20:(19+pca_num)]

if(!is.null(peerfile) & peer_num != 0){
  peer <- read.table(peerfile, header=T,sep="\t")[0:peer_num]
  cov <- merge(mc, peer, by="row.names", all=F)
  colnames(cov)[1] <- "sample"
  factors <- c(colnames(mc)[grep("PC",colnames(mc))], colnames(peer))
}else{
  cov <- mc
  cov$sample <- row.names(cov)
  factors <- c(colnames(mc)[grep("PC",colnames(mc))])
}

recal_att_f <- read.table("/home/users/nus/e1327989/YuhaoDong/GTEx/v11_fillter_recal/filter_5_recal_att_filter.txt", sep="\t", header=T, row.names=1)
recal_att_f$id <- row.names(recal_att_f)

m_data <- merge(cov, recal_att_f, by="sample", all=F)
row.names(m_data) <- m_data$id

# ===========================================================================================
info <- read.table("/home/users/nus/e1327989/YuhaoDong/GTEx/v11_fillter_recal/Promoter_info_filter_5.txt", sep="\t", header=T, row.names=1)
info <- subset(info, seqnames == chr)
m_data <- m_data[colnames(m_data) %in% c(colnames(m_data)[-grep("PromoterID",colnames(m_data))], row.names(info))]
print(paste("Promoters all", ncol(m_data), sep=" "))


m_data <- subset(m_data, m_data[[tissue_col]] %in% tissue)
da <- cbind(m_data['sample'], m_data[,grep("Promoter", colnames(m_data))])

#filter
    s <- as.data.frame(apply(da, 2, FUN=function(x){sd(x[!is.na(x)])} ) )
    colnames(s) <- 'sd'
    ss <- subset(s, sd >= 0.005)
    
    n <- as.data.frame(apply(da, 2, FUN=function(x){sum(is.na(x))}) /nrow(da))
    colnames(n) <- 'N'
    nn <- subset(n, N < 0.4)

    #f_data <- as.data.frame(t(da[,colnames(da) %in% row.names(pp1) & colnames(da) %in% row.names(pp0) & colnames(da) %in% row.names(ss) & colnames(da) %in% row.names(nn)]))
    f_data <- as.data.frame(da[,colnames(da) %in% row.names(ss) & colnames(da) %in% row.names(nn)])
    #f_data <- as.data.frame(da[,colnames(da) %in% row.names(nn)])

data <- cbind(m_data[,-grep("Promoter", colnames(m_data))], f_data)

# Get the TPM data (if available)
if(tissue != "Cells_Leukemia_cell_line_CML"){
  cat("Reading TPM data for", tissue, "...\n")
  tpm <- fread("/home/users/nus/e1327989/gtex/gtex_v11/TPM1_genes.txt")
  tpm <- subset(tpm, SMTSD == tissue)
  data <- cbind(data[,-grep("PromoterID", colnames(data))], 
                       data[, colnames(data) %in% row.names(info[info$geneId %in% tpm$Name,])])
}

# Determine the list of valid promoters in this tissue
tissue_promoters <- colnames(data)[grep("PromoterID_", colnames(data))]
cat("Found", length(tissue_promoters), "promoters for", tissue, "\n")
if(length(tissue_promoters) == 0) {
  stop("Error: No promoters found for the specified tissue.")
}


vcffile <- paste0("/home/users/nus/e1327989/YuhaoDong/GTEx/v9_vcf_MAF_0.05/",chr,".vcf.gz")
#vcf_data <- readVcf(vcffile, "hg38")





# remove objects ============================================================================
objects_to_remove <- c(
  "recal_att_f",     # 1.6GB - already merged into data
  "m_data",          # Already processed into data
  "da",              # Temporary object
  "s", "ss",         # Temporary object used for filtering
  "n", "nn",         # Temporary object used for filtering
  "f_data",          # Temporary object
  "mc",              # Already merged into data
  "peer",            # Already merged into data
  "cov",             # Already merged into data
  "tpm"              # If it exists, it is no longer needed
)

# Remove only the objects that actually exist
existing_objects <- objects_to_remove[objects_to_remove %in% ls()]
if(length(existing_objects) > 0) {
  cat("Removing objects:", paste(existing_objects, collapse=", "), "\n")
  rm(list=existing_objects)
}

# Force garbage collection
gc(verbose = TRUE)

# ===========================================================================


# Set up the parallel environment
num_cores <- thread
#num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Run in parallel
results <- foreach(p = tissue_promoters, .packages = c("gamlss", "GenomicRanges", "VariantAnnotation", "qvalue")) %dopar% {
      tryCatch({
      result <- orig_QTL_regression(
        data=data, info=info,
        tissue_col=tissue_col, flank=flank,
        tissue_name=tissue, activity_col=p, factors=factors, vcffile=vcffile)
      return(list(tissue = tissue, promoter = p, result = result))
    }, error = function(e) {
      return(list(tissue = tissue, promoter = p, error = as.character(e)))
    })

  }

# Shut down the parallel cluster
stopCluster(cl)


processed_results <- list()
for (result in results) {  # Simplified to a single-level loop
    if (is.null(result$error)) {
      processed_results[[paste(result$tissue, result$promoter, sep = "_")]] <- result$result
    } else {
      processed_results[[paste(result$tissue, result$promoter, sep = "_")]] <- result$error
    }
}


saveRDS(processed_results, paste0(output,".factor_peer_", peer_num, "_pca_", pca_num, ".results.rds"))

