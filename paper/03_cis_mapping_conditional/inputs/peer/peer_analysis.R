#!/usr/bin/env Rscript

# peer_analysis.R
# PEER analysis script for a single tissue

library(peer)

# Get the command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop("Usage: Rscript peer_analysis.R <tissue_name> <input_file> <output_dir>")
}

tissue <- args[1]
input_file <- args[2]
output_dir <- args[3]

# Create a safe file name
safe_tissue <- gsub(" ", "_", tissue)

cat(sprintf("Processing tissue: %s\n", tissue))
cat(sprintf("Input file: %s\n", input_file))
cat(sprintf("Output directory: %s\n", output_dir))

# Read the data
d <- read.table(input_file, sep='\t', header=TRUE, stringsAsFactors=FALSE)

# Analysis function
get_factor <- function(data, tissue_name) {
    tryCatch({
        cat(sprintf("Analyzing tissue: %s\n", tissue_name))
        
        # Select the data of the current tissue
        t_data <- subset(data, SMTSD == tissue_name)
        
        if (nrow(t_data) == 0) {
            stop(sprintf("No data found for tissue: %s", tissue_name))
        }
        
        cat(sprintf("Samples found: %d\n", nrow(t_data)))
        
        row.names(t_data) <- t_data$sample
        
        # Extract the expression data
        expr <- t_data[, -c(1:3)]
        
        cat(sprintf("Expression matrix dimension: %d x %d\n", nrow(expr), ncol(expr)))
        
        # Initialize the PEER model
        model <- PEER()
        
        # Set the phenotype data
        PEER_setPhenoMean(model, as.matrix(expr))
        
        cat(sprintf("Phenotype matrix dimension: %s\n", 
                    paste(dim(PEER_getPhenoMean(model)), collapse=" x ")))
        
        # Set the number of factors
        PEER_setNk(model, 60)
        cat(sprintf("Number of factors: %d\n", PEER_getNk(model)))
        
        # Compute the factors
        cat("Running PEER update...\n")
        PEER_update(model)
        cat("PEER update completed!\n")
        
        # Get the factors
        factors <- PEER_getX(model)
        
        row.names(factors) <- row.names(expr)
        colnames(factors) <- paste0('peer_', 1:ncol(factors))
        
        return(factors)
        
    }, error = function(e) {
        cat(sprintf("Error occurred: %s\n", e$message))
        return(NULL)
    })
}

# Run the analysis
factors <- get_factor(d, tissue)

# Save the results
if (!is.null(factors)) {
    output_file <- file.path(output_dir, 
                             paste0('v11_filter_5_peer_factors_', safe_tissue, '.txt'))
    
    write.table(factors, output_file, sep='\t', quote=FALSE)
    
    cat(sprintf("Results saved to: %s\n", output_file))
    cat(sprintf("Factor matrix dimension: %d x %d\n", nrow(factors), ncol(factors)))
    cat("Analysis completed successfully!\n")
} else {
    stop("Analysis failed - no factors generated")
}