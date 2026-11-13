#conda create -n peer -c conda-forge -c bioconda r-peer
library(peer)
library(parallel)
library(foreach)
library(doParallel)
library(coda)
library(dplyr)

get_factor <- function(data, tissue){
tryCatch({
    print(tissue)
    t_data <- subset(data, SMTSD == tissue)
    row.names(t_data) <- t_data$sample
    
    expr <- t_data[-c(1:3)]

    model = PEER()

    #NULL response means no error here
    PEER_setPhenoMean(model,as.matrix(expr))

    dim(PEER_getPhenoMean(model))

    # Set the number of factors
    PEER_setNk(model,60)

    PEER_getNk(model)
    # Compute the factors
    PEER_update(model)


    factors = PEER_getX(model)

    row.names(factors) <- row.names(expr)
    colnames(factors) <- paste0('peer_', 1:ncol(factors))

    write.table(factors, paste0('v10_filter_5_peer_factors_',tissue,'.txt'), sep='\t', quote=F)

    return(factors)

    }, error = function(e) {
      print(e)
    })
}


# Define the tissue and promoter list to process
#recal_att_f <- read.table('filter_5_recal_att.txt', sep='\t', header=T, row.names=1)
#data <- recal_att_f %>%
#  group_by(sample, SMTSD) %>%
#  summarise(SMTS = first(SMTS),
#            across(starts_with("PromoterID_"), 
#                   .fns = ~mean(., na.rm = TRUE)))

#d <- as.data.frame(data)
#d[is.na(d)] <- 0
#write.table(d, "filter_5_recal_att_filter_SMTSD_mean.txt", sep="\t", row.names=F, col.names=T, quote=F)



d <- read.table("/home/users/nus/e1327989/YuhaoDong/GTEx/peer/v11/filter_5_recal_att_filter_SMTSD_mean.txt", sep='\t', header=T)

tissue_list <- unique(d$SMTSD)
#tissue_list <- c('Cervix_Endocervix',
#'Cervix_Ectocervix',
#'Stomach_Muscularis',
#'Fallopian_Tube',
#'Stomach_Mixed_Cell',
#'Stomach_Mucosa')
#tissue_list <- setdiff(tissue_list, done_list)


# Set up the parallel environment
num_cores <- 50
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Run in parallel
results <- foreach(tissue = tissue_list, .packages = c("peer")) %dopar% {
    factors <- get_factor(d, tissue)
}

# Shut down the parallel cluster
stopCluster(cl)

saveRDS(results, "v10_SMTSD_peer15_factors.rds")

# Process the results
for (i in 1:length(tissue_list)) {
    result <- results[[i]]
    tissue <- gsub(" ","_",tissue_list[i])
    write.table(result, file=paste0(tissue,"_SMTSD_filter_5_recal_peer15.txt"), sep="\t",row.names=T, col.names=T,quote=F)
}


