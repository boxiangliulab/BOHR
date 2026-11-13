### merge data=================================================
library(SummarizedExperiment)
library(Matrix)
library(data.table)

load("att_sample.info.RData")
att$sample <- vapply(strsplit(rownames(att), "-"), FUN.VALUE = character(1),
                     function(x) paste(x[1], x[2], sep="-"))

rds_files <- list.files("./proActiv/", pattern="\\.rds$", full.names=TRUE)

clean_sample_names <- function(x) {
  x <- gsub("\\.SJ\\.out.*$", "", x)
  x <- gsub("\\.SJ\\.out$", "", x)
  x <- gsub("\\.", "-", x)
  x
}

to_num_mat <- function(x) {
  if (is.data.frame(x)) {
    num <- vapply(x, is.numeric, logical(1))
    as.matrix(x[, num, drop=FALSE])
  } else {
    as.matrix(x)
  }
}

#
rd_list <- list()
rel_list  <- list()
abs_list  <- list()

rd0 <- as.data.frame(rowData(readRDS(rds_files[1])))[,1:8]
pid <- paste0("PromoterID_", rd0$promoterId)

for (f in rds_files) {
  message("Processing: ", basename(f))
  se <- readRDS(f)
  rd <- as.data.frame(rowData(se))[,9:11]
  rownames(rd) <- pid
  rd_list[[basename(f)]] <- rd
}

rdata <- cbind(rd0, do.call(cbind, rd_list))
row.names(rdata) <- paste0("PromoterID_", rdata$promoterId)
save(rdata, file="filter_5_rdata.RData")


info_all <- rdata[1:8]
save(info_all, file=paste0("Promoter_info_all.RData"))

# filter internalPromoter
irdata <- rdata[rdata$internalPromoter == 'FALSE',]

pid_filter <- paste0("PromoterID_", irdata$promoterId)
row.names(irdata) <- pid_filter

# merge data
for (f in rds_files) {
  message("Merging: ", basename(f))
  se <- readRDS(f)

  rel_m <- to_num_mat(assays(se)[[4]])
  abs_m <- to_num_mat(assays(se)[[3]])

  rownames(rel_m) <- pid
  rownames(abs_m) <- pid

  colnames(rel_m) <- clean_sample_names(colnames(rel_m))
  colnames(abs_m) <- clean_sample_names(colnames(abs_m))

  rel_m <- rel_m[pid_filter, , drop=FALSE]
  abs_m <- abs_m[pid_filter, , drop=FALSE]

  trel <- t(Matrix(rel_m, sparse=TRUE))  # samples x promoters
  tabs <- t(Matrix(abs_m, sparse=TRUE))

  rel_list[[basename(f)]]  <- trel
  abs_list[[basename(f)]]  <- tabs
}

meta_all <- att[, c("sample","SMTS","SMTSD"), drop=FALSE]

rel_all <- do.call(rbind, rel_list)
abs_all <- do.call(rbind, abs_list)

base::saveRDS(list(meta = meta_all, rel = rel_all), "trel_att_filter_In.rds")
base::saveRDS(list(meta = meta_all, abs = abs_all), "tabs_att_filter_In.rds")



##absolute activity filter==========================================================
class_col <- na.omit(as.data.frame(rdata[, grep('class', colnames(rdata))]))
class_filter <- class_col[rowSums(class_col != 'Inactive') > 0,]
irdata_filter <- irdata[row.names(class_filter),]

###relative activity filter==========================================================
rel_all_mat <- as.matrix(rel_all)
trel_att_raw <- merge(meta_all, rel_all_mat, by="row.names", all=FALSE)
row.names(trel_att_raw) <- trel_att_raw$Row.names
trel_att_raw <- trel_att_raw[-1]

num <- as.data.frame(table(rdata$geneId))
num <- num[order(num$Freq, decreasing = T),]
x <- as.data.frame(table(num$Freq))
num_1 <- subset(num, Freq == 1)
rdata_num1 <- rdata[rdata$geneId %in% num_1$Var1,]

trel_filter <- trel_att_raw[,row.names(irdata_filter)]
trel_filter <- trel_filter[,!colnames(trel_filter) %in% row.names(rdata_num1)]
trel_filter <- cbind(trel_att_raw[c(1:3)], trel_filter)



###
t <- as.data.frame(table(trel_filter$SMTSD))
t$Var1 <- as.character(t$Var1)
colnames(t) <- c('Tissue', 'Samples')

  for(k in 1:nrow(t)){
      type <- as.character(t[k,1])
      print(paste0(k, '-', type,'-processing'))
      trel_att_t <- subset(trel_filter, SMTSD == type)

      data <- trel_att_t[,-c(1:3)]

      p <- as.data.frame(apply(data > 0, 2, FUN=function(x){sum(x[!is.na(x)])}) /nrow(data))
      colnames(p) <- type

      if(k == 1){
        c <- p
      }else{
        c <- cbind(c, p)
      }

      m <- as.data.frame(apply(data, 2, FUN=function(x){mean(x[!is.na(x)])}) )
      colnames(m) <- type

      if(k == 1){
        cc <- m
      }else{
        cc <- cbind(cc, m)
      }

      s <- as.data.frame(apply(data, 2, FUN=function(x){sd(x[!is.na(x)])}) )
      colnames(s) <- type

      if(k == 1){
        sd <- s
      }else{
        sd <- cbind(sd, s)
      }

      n <- as.data.frame(apply(data, 2, FUN=function(x){sum(is.na(x))}) /nrow(data))
      colnames(n) <- type

      if(k == 1){
        nn <- n
      }else{
        nn <- cbind(nn, n)
      }
  }

# at least one tissue {rel > 0} > 25%
filter2 <- c[rowSums(c > 0.25, na.rm=T) > 0,]
trel_filter2 <- trel_filter[,colnames(trel_filter) %in% row.names(filter2)]
dim(trel_filter2)

# at least one tissue {NA} < 40%
filter3 <- nn[rowSums(nn < 0.4, na.rm=T) > 0,]
trel_filter3 <- trel_filter2[,colnames(trel_filter2) %in% row.names(filter3)]
dim(trel_filter3)

# at least one tissue {mean rel} > 0.1
filter4 <- cc[rowSums(cc > 0.1, na.rm=T) > 0,]
trel_filter4 <- trel_filter3[,colnames(trel_filter3) %in% row.names(filter4)]

trel_filter4 <- cbind(trel_filter[c(1:3)],trel_filter4)
dim(trel_filter4)

# remove all single promoters========================================================
iirdata <- irdata[row.names(irdata) %in% colnames(trel_filter4)[-c(1:3)],]
num <- as.data.frame(table(iirdata$geneId))
num <- num[order(num$Freq, decreasing = T),]
x <- as.data.frame(table(num$Freq))

num_1 <- subset(num, Freq == 1)
iirdata_num1 <- iirdata[!(iirdata$geneId %in% num_1$Var1),]

trel_att <- trel_filter4[, colnames(trel_filter4) %in% row.names(iirdata_num1)]
trel_att <- cbind(trel_filter4[c(1:3)],trel_att)
dim(trel_att)
save(trel_att, file=paste0('filter_5_trel_att.RData'))

abs <- as.matrix(t(abs_all))
abs_att <- abs[row.names(abs) %in% colnames(trel_att),colnames(abs) %in% row.names(trel_att)]
tabs_att <- as.data.frame(t(abs_att))
tabs_att <- cbind(trel_att[1:3], tabs_att)
save(tabs_att, file=paste0('filter_5_tabs_att.RData'))

info <- rdata[row.names(rdata) %in% colnames(trel_att),1:8]
save(info, file=paste0("Promoter_info_filter_5.RData"))








### re-calculate relative activity=====================================================
in_list <- row.names(subset(rdata, internalPromoter == 'FALSE'))
in_info <- info_all[row.names(info_all) %in% in_list,]
in_trel <- trel_att_raw[, colnames(trel_att_raw) %in% c(colnames(trel_att_raw)[1:3], in_list)]


abs <- as.matrix(t(abs_all))

abs_att <- abs[row.names(abs) %in% colnames(in_trel), ]

abs_info <- cbind(in_info[2], abs_att)


library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(Matrix)
library(BiocGenerics) # If a DataFrame needs to be processed

calculate_promoter_ratios_fast <- function(abs_info) {
  # Convert the DataFrame to a data.table
  dt <- as.data.table(as.data.frame(abs_info))
  
  # Extract promoterID
  promoterID <- rownames(abs_info)
  dt[, promoterID := promoterID]
  
  # Move the geneId and promoterID columns to the front
  setcolorder(dt, c("promoterID", "geneId"))
  
  # Convert to a matrix for the calculation
  mat <- as.matrix(dt[, -c(1:2)])  # Exclude the promoterID and geneId columns
  row.names(mat) <- dt$promoterID

  # Create a sparse matrix representing the grouping of each gene
  genes <- as.factor(dt$geneId)
  n <- nrow(mat)
  p <- length(levels(genes))
  i <- seq_len(n)
  j <- as.integer(genes)
  x <- rep(1, n)
  group_mat <- sparseMatrix(i = i, j = j, x = x, dims = c(n, p))
  
  # Compute the total activity of each gene
  gene_totals <- as.matrix(t(group_mat) %*% mat)
  
  # Compute the proportions
  ratios <- mat / gene_totals[genes, ]

  
  # Compute the mean proportion
  average_ratio <- rowMeans(ratios, na.rm = TRUE)
  
  # Combine the results into one data.table
  result <- data.table(
    promoterID = promoterID,
    geneId = dt$geneId,
    as.data.table(ratios),
    average_ratio = average_ratio
  )
  
  # Set the key
  setkey(result, promoterID)
  
  return(result)
}

# Use the function
result <- calculate_promoter_ratios_fast(abs_info)


head(result)[,1:8]

recal <- as.data.frame(result)
row.names(recal) <- recal$promoterID
recal <- recal[-1]
colnames(recal) <- gsub("[.]", "-", colnames(recal))

recal_att <- merge(in_trel[1:3], as.data.frame(t(recal[-1])), by="row.names", all=F)
row.names(recal_att) <- recal_att$Row.names
recal_att <- recal_att[-1]

recal_att$SMTSD <- as.character(recal_att$SMTSD)
recal_att$SMTSD <- gsub("[ ]", "_", recal_att$SMTSD)
recal_att$SMTSD <- gsub("[(]", "", recal_att$SMTSD)
recal_att$SMTSD <- gsub("[)]", "", recal_att$SMTSD)

recal_att$SMTS <- as.character(recal_att$SMTS)
recal_att$SMTS <- gsub("[ ]", "_", recal_att$SMTS)
recal_att$SMTS <- gsub("[(]", "", recal_att$SMTS)
recal_att$SMTS <- gsub("[)]", "", recal_att$SMTS)


save(recal_att, file=paste0('filter_5_recal_att_raw.RData'))


recal_att_f <- cbind(recal_att[1:3], recal_att[,colnames(recal_att) %in% row.names(info)])

save(recal_att_f, file=paste0('filter_5_recal_att_filter.RData'))
