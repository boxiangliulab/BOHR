library(ACAT)
library(qvalue)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggrepel) 
library(data.table)

tissue_list <- fread("../tissue_list_size20_table.txt")
tissue_list <- as.data.frame(tissue_list)
for(tissue_row in 1:nrow(tissue_list)){
    tissue <- tissue_list[tissue_row,1]
    cat(paste0("Processing ", tissue))

    peer_num <- length(files <- Sys.glob(paste0(tissue, "/sensitivity_test_10k.*rds")))
    if(peer_num==0){next}

    ###
    peer_count <- as.data.frame(1:peer_num)
    colnames(peer_count) <- "peer"
    peer_count$count <- NA

    for(peer in 1:peer_num){
        print(peer)
        t <- readRDS(paste0(tissue, "/sensitivity_test_10k.",tissue,".chr22.factor_peer_",peer,"_pca_5.results.rds"))
        if(length(t) == 0){print("test");next}

        acat_pvalues <- as.data.frame(names(t))
        colnames(acat_pvalues) <- "promoterID"
        acat_pvalues$acat_p <- NA


        for(i in 1:length(t)){
            #print(names(t)[i])
            tt <- t[[i]]
            if(nrow(tt)<=1){next}
            tt[which(tt$lrt_pvalue == 0),] <- .Machine$double.xmin
            acat_p <- ACAT(tt$lrt_pvalue)

            acat_pvalues[i, 2] <- acat_p
        }

        acat_pvalues$acat_q <- qvalue(acat_pvalues$acat_p, lambda=0)$qvalues
        
        # Use lambda=0
        q_values_lambda0 <- qvalue(acat_pvalues$acat_p, lambda=0)$qvalues

        # Try running the default qvalue function
        result <- tryCatch({
        qvalue_result <- qvalue(acat_pvalues$acat_p)
        acat_pvalues$acat_q <- qvalue_result$qvalues
        cat("Default qvalue method ran; estimated π₀ =", qvalue_result$pi0, "\n")
        TRUE  # Indicates success
        }, error = function(e) {
        cat("Default qvalue method failed:", conditionMessage(e), "\n")
        cat("Trying lambda=0 instead...\n")
        FALSE  # Indicates failure
        })

        # If the default method fails, use lambda=0
        if (!result) {
        qvalue_result <- qvalue(acat_pvalues$acat_p, lambda=0)
        acat_pvalues$acat_q <- qvalue_result$qvalues
        }


        sig_num <- nrow(subset(acat_pvalues, acat_q <= 0.05))
        
        peer_count[which(peer_count$peer == peer), 2] <- sig_num

    }
    if(nrow(na.omit(peer_count))==0){next}

    write.table(peer_count, paste0(tissue, "_peer_test_count.txt"), sep="\t", quote=F, row.names=F, col.names=T)



    # Step 1: compute the percentage change in the number of eGenes for each added PEER factor
    # Start from PEER=0 (if present)
    pct_change <- numeric(nrow(peer_count)-1)
    for(i in 1:(nrow(peer_count)-1)) {
    pct_change[i] <- (peer_count$count[i+1] - peer_count$count[i]) / peer_count$count[i] * 100
    }

    # Create a data frame to store the percentage changes
    change_df <- data.frame(
    peer = peer_count$peer[2:nrow(peer_count)],  # Corresponding number of PFs (starting from the second)
    pct_change = pct_change
    )

    # Step 2: smooth the percentage changes with a LOESS fit
    loess_fit <- loess(pct_change ~ peer, data = change_df, span = 0.75)
    smoothed_values <- predict(loess_fit, change_df$peer)

    # Add the fitted values to the data frame
    change_df$smoothed <- smoothed_values

    # Step 3: find the point before the LOESS curve becomes negative
    positive_indices <- which(change_df$smoothed > 0)
    if(length(positive_indices) > 0) {
    optimal_idx <- max(positive_indices)
    optimal_pf <- change_df$peer[optimal_idx]
    } else {
    # If all values are negative, choose the first point
    optimal_pf <- change_df$peer[1]
    }

    # Output the optimal number of PEER factors
    cat("Optimal number of PEER factors (local greedy algorithm):", optimal_pf, "\n")


    # Visualize the results - Plot 1: raw eGene counts with the optimal point marked
    library(ggplot2)

    top5_points <- peer_count %>% 
    arrange(desc(count)) %>% 
    slice(1:5)
    min_count <- min(peer_count$count)
    # Plot 1: line plot of the raw eGene counts
    p1 <- ggplot(peer_count, aes(x = peer, y = count)) +
    geom_line(color = "blue", size = 1) +
    geom_point(color = "blue", size = 3) +
    #geom_vline(xintercept = optimal_pf, linetype = "dashed", color = "red") +
    annotate("text", x = optimal_pf, y = max(peer_count$count)/2, 
            label = paste("best PF =", optimal_pf), color = "red", hjust = -0.1) +
    # Label the five largest points
    geom_text_repel(data = top5_points,
                    aes(x = peer, y = count, label = paste(peer, count, sep="-")),
                    color = "black",
                    size = 3,
                    direction = "x",      # Adjust the position along the x axis only
                    vjust = -0.5,
                    nudge_y = 0.5,       # Fine-tune the y position
                    segment.color = "grey50", # Add guide lines
                    min.segment.length = 0) +
    labs(title = "sig promoter number",
        x = "PEER num",
        y = "sig promoter") +
    theme_minimal()

    # Plot 2: percentage change with the LOESS fit
    p2 <- ggplot(change_df, aes(x = peer)) +
    geom_point(aes(y = pct_change), color = "blue", size = 3) +
    geom_line(aes(y = smoothed), color = "red", size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = optimal_pf, linetype = "dashed", color = "red") +
    annotate("text", x = optimal_pf, y = max(change_df$pct_change)/2, 
            label = paste("best PF =", optimal_pf), color = "red", hjust = -0.1) +
    labs(title = paste0("PEER factor test: LOESS best ", optimal_pf),
        x = "PEER num",
        y = "sig promoter number change (%)") +
    theme_minimal()


    # Show the two plots
    library(gridExtra)

    pdf(paste0(tissue, "_peer_test_count.pdf"))
    grid.arrange(p1, p2, ncol = 1)
    dev.off()

}