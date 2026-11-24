### NO-THRESHOLD copy of v12_simulation_realfloor.R (newmethod_bohr2/README)
### REALFLOOR copy of v12_simulation.R (built by calibration_perm1000/make_realfloor_generators.R)
### library ===============================================================
library(Biostrings)
library(polyester)
library(truncnorm)
library(stringr)
library(dplyr)
library(GenomicRanges)
library(gamlss.dist)

logit <- function(p) log(p/(1-p))
inv_logit <- function(x) exp(x)/(1+exp(x))
### Auxiliary Functions ===============================================================
# Generate random parameters
# V10 changes vs V8:
#   1. middle_value p01: Uniform(0, 0.7) instead of Beta(1.5, 30) (~5%)
#   2. Relaxed p0/p1 type constraints: middle down to 10%, p1/p0 up to 0.6
# V8 changes vs V3:
#   1. Add "one_proportion" effect type
#   2. Balanced baseline p0/p1: Uniform(0.2, 0.8) for boundary types
generate_random_params <- function(effect_type=NULL) {
    if(is.null(effect_type)){
        effect_type <- sample(c("zero_proportion", "one_proportion",
                                "middle_value", "mixed_effect"), 1)
    }
    ## NO-THRESHOLD generation: baseline boundary proportions from ONE common distribution
    ## for every effect type -- p0 + p1 ~ U(0, 0.9), split uniformly. Two runif() draws, as in every realfloor branch,
    ## so all later draws (null flag, effect size, middle mean/sd, genotypes) follow the same random stream.
    p01 <- runif(1, 0, 0.9)
    p0 <- p01*(runif(1, min = 0, max = 1))
    p1 <- p01-p0

    # 10% null, 90% non-zero effect
    if(runif(1) < 0.1) {
        effect_size <- 0
    } else {
        effect_size <- rtruncnorm(1, a=0.001, b=0.5, mean=0.1, sd=0.2)
    }

    params <- list(
        p_zero = p0,
        p_one = p1,
        middle_mean = runif(1, min = 0.1, max = 0.9) ,
        middle_sd = rgamma(1, shape=2, rate=20),
        effect_size = effect_size,
        effect_type =  effect_type,
        maf = rtruncnorm(1, a=0.1, b=0.5, mean=0.25, sd=0.1)
    )

    return(params)
}

# Generate transcript expression ratios
# effect_scale: "prob" (additive on probability) or "link" (additive on logit/log)
generate_transcript_ratios <- function(n_samples, params, effect_scale = "prob") {

    # Generate genotypes
    genotypes <- sample(c(0,1,2), n_samples, replace=TRUE,
                       prob=c((1-params$maf)^2,
                             2*params$maf*(1-params$maf),
                             params$maf^2))

    # Generate transcript 1 proportions
    trans1 <- numeric(n_samples)

    # Pre-determine effect directions ONCE (before genotype loop)
    # Boundary-protecting: push toward 0.5
    if(params$effect_type == "middle_value") {
        dir_mu  <- if(params$middle_mean < 0.5) 1 else -1
        dir_nu  <- 0; dir_tau <- 0
    } else if(params$effect_type == "zero_proportion") {
        dir_mu  <- 0; dir_tau <- 0
        dir_nu  <- if(params$p_zero < 0.5) 1 else -1
    } else if(params$effect_type == "one_proportion") {
        dir_mu  <- 0; dir_nu  <- 0
        dir_tau <- if(params$p_one < 0.5) 1 else -1
    } else if(params$effect_type == "mixed_effect") {
        dir_mu  <- sample(c(-1, 1), 1)
        dir_nu  <- sample(c(-1, 1), 1)
        dir_tau <- 0
    }

    ### Proportional boundary scaling (instead of hard clamping) ================
    et <- params$effect_type
    p0 <- params$p_zero; p1 <- params$p_one
    mid_mean <- params$middle_mean
    p_mid_base <- max(1 - p0 - p1, 0.01)
    nu_base <- p0 / p_mid_base
    tau_base <- p1 / p_mid_base

    sf <- 1.0
    for(iter in 1:20) {
        eff_scaled <- params$effect_size * sf
        valid <- TRUE
        for(gt_val in c(0, 1, 2)) {
            ev <- eff_scaled * gt_val
            if(effect_scale == "prob") {
                if(et == "middle_value") {
                    cm <- mid_mean + dir_mu * ev; cp0 <- p0; cp1 <- p1
                } else if(et == "zero_proportion") {
                    cm <- mid_mean; cp0 <- p0 + dir_nu * ev
                    cp1 <- min(p1, 1 - cp0 - 0.1)
                } else if(et == "one_proportion") {
                    cm <- mid_mean; cp1 <- p1 + dir_tau * ev
                    cp0 <- min(p0, 1 - cp1 - 0.1)
                } else {
                    cm <- mid_mean + dir_mu * ev; cp0 <- p0 + dir_nu * ev
                    cp1 <- min(p1, 1 - cp0 - 0.1)
                }
            } else {
                if(et == "middle_value") {
                    cm <- inv_logit(logit(mid_mean) + dir_mu * ev); cp0 <- p0; cp1 <- p1
                } else if(et == "zero_proportion") {
                    cm <- mid_mean
                    nu_new <- nu_base * exp(dir_nu * ev)
                    cp0 <- nu_new / (1 + nu_new + tau_base)
                    cp1 <- tau_base / (1 + nu_new + tau_base)
                } else if(et == "one_proportion") {
                    cm <- mid_mean
                    tau_new <- tau_base * exp(dir_tau * ev)
                    cp0 <- nu_base / (1 + nu_base + tau_new)
                    cp1 <- tau_new / (1 + nu_base + tau_new)
                } else {
                    cm <- inv_logit(logit(mid_mean) + dir_mu * ev)
                    nu_new <- nu_base * exp(dir_nu * ev)
                    cp0 <- nu_new / (1 + nu_new + tau_base)
                    cp1 <- tau_base / (1 + nu_new + tau_base)
                }
            }
            if(cm < 0.01 || cm > 0.99) { valid <- FALSE; break }
            if(cp0 < 0.01 || cp0 > 0.8) { valid <- FALSE; break }
            if(cp1 < 0 || cp1 > 0.8) { valid <- FALSE; break }
            if(cp0 + cp1 > 0.9) { valid <- FALSE; break }
        }
        if(valid) break
        sf <- sf * 0.8
    }
    eff_size_adj <- params$effect_size * sf

    for(gt in 0:2) {
        idx <- which(genotypes == gt)
        if(length(idx) == 0) next
        n_curr <- length(idx)
        effect_val <- eff_size_adj * gt

        if(effect_scale == "prob") {
            # --- Prob scale: additive on probability, then convert to BEINF ---
            if(et == "middle_value") {
                curr_mid <- mid_mean + dir_mu * effect_val; curr_p0 <- p0; curr_p1 <- p1
            } else if(et == "zero_proportion") {
                curr_mid <- mid_mean; curr_p0 <- p0 + dir_nu * effect_val
                curr_p1 <- max(min(p1, 1 - curr_p0 - 0.1), 0)
            } else if(et == "one_proportion") {
                curr_mid <- mid_mean; curr_p1 <- p1 + dir_tau * effect_val
                curr_p0 <- max(min(p0, 1 - curr_p1 - 0.1), 0)
            } else {
                curr_mid <- mid_mean + dir_mu * effect_val
                curr_p0 <- p0 + dir_nu * effect_val
                curr_p1 <- max(min(p1, 1 - curr_p0 - 0.1), 0)
            }
            curr_mid <- min(max(curr_mid, 0.01), 0.99)
            curr_p0 <- min(max(curr_p0, 0.001), 0.8)
            curr_p1 <- min(max(curr_p1, 0), 0.8)
            # Convert to BEINF parameters
            p_mid <- max(1 - curr_p0 - curr_p1, 0.01)
            mu_beinf <- curr_mid
            nu_beinf <- curr_p0 / p_mid
            tau_beinf <- curr_p1 / p_mid
        } else {
            # --- Link scale: additive on logit(mu)/log(nu)/log(tau) ---
            if(et == "middle_value") {
                mu_beinf <- inv_logit(logit(mid_mean) + dir_mu * effect_val)
                nu_beinf <- nu_base; tau_beinf <- tau_base
            } else if(et == "zero_proportion") {
                mu_beinf <- mid_mean
                nu_beinf <- nu_base * exp(dir_nu * effect_val); tau_beinf <- tau_base
            } else if(et == "one_proportion") {
                mu_beinf <- mid_mean
                nu_beinf <- nu_base; tau_beinf <- tau_base * exp(dir_tau * effect_val)
            } else {
                mu_beinf <- inv_logit(logit(mid_mean) + dir_mu * effect_val)
                nu_beinf <- nu_base * exp(dir_nu * effect_val); tau_beinf <- tau_base
            }
        }

        # Clamp BEINF parameters
        sigma_beinf <- min(params$middle_sd / sqrt(mu_beinf * (1 - mu_beinf) + 1e-10), 0.99)
        sigma_beinf <- max(sigma_beinf, 0.01)
        nu_beinf <- max(nu_beinf, 1e-4)
        tau_beinf <- max(tau_beinf, 1e-4)

        # Generate data with rBEINF
        trans1[idx] <- rBEINF(n_curr, mu = mu_beinf, sigma = sigma_beinf,
                              nu = nu_beinf, tau = tau_beinf)
    }

    trans2 <- 1 - trans1

    # Convert to fold changes
    fold_changes <- matrix(
        c(trans1 * 2, trans2 * 2),
        nrow = 2,
        ncol = n_samples
    )

    # Compute true prob-scale effects (gt=2 vs gt=0) for ground truth
    if(effect_scale == "prob") {
        true_mu  <- dir_mu * eff_size_adj
        true_nu  <- dir_nu * eff_size_adj
        true_tau <- dir_tau * eff_size_adj
    } else {
        true_mu <- inv_logit(logit(mid_mean) + dir_mu * 2 * eff_size_adj) - mid_mean
        nu_at0 <- nu_base; tau_at0 <- tau_base
        nu_at2 <- nu_base * exp(dir_nu * 2 * eff_size_adj)
        tau_at2 <- tau_base * exp(dir_tau * 2 * eff_size_adj)
        true_nu <- nu_at2/(1+nu_at2+tau_at2) - nu_at0/(1+nu_at0+tau_at0)
        true_tau <- tau_at2/(1+nu_at2+tau_at2) - tau_at0/(1+nu_at0+tau_at0)
    }

    return(list(
        genotypes = genotypes,
        transcript1_props = trans1,
        transcript2_props = trans2,
        fold_changes = fold_changes,
        directions = list(mu = dir_mu, nu = dir_nu, tau = dir_tau),
        scale_factor = sf,
        effect_size_adj = eff_size_adj,
        true_effects = list(mu = true_mu, nu = true_nu, tau = true_tau)
    ))
}

### Generate VCF file ===================================================================================
generate_vcf <- function(sim_info, output_file="genotypes.vcf") {
    if(length(sim_info$sim_results) == 0) {
        warning("sim_results is empty, cannot generate VCF")
        return(NULL)
    }

    vcf_header <- c(
        "##fileformat=VCFv4.2",
        paste0("##fileDate=", format(Sys.Date(), "%Y%m%d")),
        "##source=SimulationScript_V10",
        "##reference=GRCh38",
        "##contig=<ID=chr1,length=248956422>",
        "##contig=<ID=chr2,length=242193529>",
        "##contig=<ID=chr3,length=198295559>",
        "##contig=<ID=chr4,length=190214555>",
        "##contig=<ID=chr5,length=181538259>",
        "##contig=<ID=chr6,length=170805979>",
        "##contig=<ID=chr7,length=159345973>",
        "##contig=<ID=chr8,length=145138636>",
        "##contig=<ID=chr9,length=138394717>",
        "##contig=<ID=chr10,length=133797422>",
        "##contig=<ID=chr11,length=135086622>",
        "##contig=<ID=chr12,length=133275309>",
        "##contig=<ID=chr13,length=114364328>",
        "##contig=<ID=chr14,length=107043718>",
        "##contig=<ID=chr15,length=101991189>",
        "##contig=<ID=chr16,length=90338345>",
        "##contig=<ID=chr17,length=83257441>",
        "##contig=<ID=chr18,length=80373285>",
        "##contig=<ID=chr19,length=58617616>",
        "##contig=<ID=chr20,length=64444167>",
        "##contig=<ID=chr21,length=46709983>",
        "##contig=<ID=chr22,length=50818468>",
        "##contig=<ID=chrX,length=156040895>",
        "##contig=<ID=chrY,length=57227415>",
        "##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor Allele Frequency\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
               paste(paste0("sample_", sprintf("%02d", 1:sim_info$simulation_params$n_samples)),
                     collapse="\t"))
    )

    vcf_records <- data.frame(
        CHROM = character(),
        POS = numeric(),
        ID = character(),
        REF = character(),
        ALT = character(),
        QUAL = character(),
        FILTER = character(),
        INFO = character(),
        FORMAT = character(),
        GENOTYPES = character(),
        stringsAsFactors = FALSE
    )

    if(is.null(sim_info$promoter_info)) {
        warning("promoter_info not found, cannot generate VCF")
        return(NULL)
    }

    for(gene_id in names(sim_info$sim_results)) {
        if(!all(c("params", "genotypes") %in% names(sim_info$sim_results[[gene_id]]))) {
            warning("gene_id: ", gene_id, " missing params or genotypes")
            next
        }

        if(!gene_id %in% names(sim_info$transcripts)) {
            warning("gene_id: ", gene_id, " missing transcript info")
            next
        }

        promoter_info <- sim_info$promoter_info %>%
            filter(geneId == gene_id) %>%
            slice(1)

        if(nrow(promoter_info) == 0) {
            warning("No promoter info found for gene: ", gene_id)
            next
        }

        chrom_info <- as.character(promoter_info$seqnames)
        pos_info <- as.integer(promoter_info$start)

        genotypes <- sim_info$sim_results[[gene_id]]$genotypes
        gt_string <- sapply(genotypes, function(g) {
            switch(as.character(g),
                   "0" = "0|0",
                   "1" = "0|1",
                   "2" = "1|1")
        })

        if("maf" %in% names(sim_info$sim_results[[gene_id]]$params)) {
            maf <- format(sim_info$sim_results[[gene_id]]$params$maf,
                         scientific = FALSE,
                         digits = 4)
        } else {
            maf <- "0.5"
        }

        vcf_records <- rbind(vcf_records, data.frame(
            CHROM = chrom_info,
            POS = pos_info,
            ID = paste0(gene_id, "_sQTL_b38"),
            REF = "A",
            ALT = "G",
            QUAL = ".",
            FILTER = "PASS",
            INFO = paste0("MAF=", maf),
            FORMAT = "GT",
            GENOTYPES = paste(gt_string, collapse="\t"),
            stringsAsFactors = FALSE
        ))
    }

    if(nrow(vcf_records) == 0) {
        warning("No valid VCF records")
        return(NULL)
    }

    vcf_records$CHROM_SORT <- factor(vcf_records$CHROM,
                                    levels=paste0("chr", c(1:22, "X", "Y")))
    vcf_records <- vcf_records[order(vcf_records$CHROM_SORT, vcf_records$POS), ]
    vcf_records$CHROM_SORT <- NULL

    vcf_lines <- apply(vcf_records, 1, function(row) {
        paste(row, collapse="\t")
    })

    writeLines(c(vcf_header, vcf_lines), output_file)

    valid_genes <- rownames(vcf_records)

    if(length(valid_genes) > 0) {
        geno_list <- lapply(sim_info$sim_results[valid_genes], function(x) {
            if("genotypes" %in% names(x)) return(x$genotypes) else return(NULL)
        })
        geno_list <- geno_list[!sapply(geno_list, is.null)]

        if(length(geno_list) > 0) {
            geno_matrix <- do.call(rbind, geno_list)
            rownames(geno_matrix) <- names(geno_list)
            colnames(geno_matrix) <- paste0("sample_", sprintf("%02d", 1:ncol(geno_matrix)))
            write.table(geno_matrix,
                       file=sub(".vcf$", ".txt", output_file),
                       sep="\t",
                       quote=FALSE)
        }
    }

    return(list(
        vcf_file = output_file,
        geno_matrix_file = sub(".vcf$", ".txt", output_file)
    ))
}

### Select transcripts from different promoters ======================================
select_transcripts_from_promoters <- function(gene_id, promoter_info, fasta_all) {
    gene_promoters <- promoter_info %>% filter(geneId == gene_id)

    if(nrow(gene_promoters) < 2) {
        warning("Gene ", gene_id, " does not have at least 2 promoters")
        return(NULL)
    }

    selected_promoters <- gene_promoters %>% sample_n(2)

    transcript_ids <- list()
    for(i in 1:2) {
        if(is.list(selected_promoters$txId)) {
            tx_list <- selected_promoters$txId[[i]]
        } else {
            tx_list <- unlist(strsplit(as.character(selected_promoters$txId[i]), ","))
        }

        if(length(tx_list) > 0) {
            transcript_ids[[i]] <- sample(tx_list, 1)
        } else {
            warning("No transcripts found for promoter ", selected_promoters$promoterId[i])
            return(NULL)
        }
    }

    selected_seqs <- DNAStringSet()
    transcript_info <- data.frame(
        transcript_id = character(),
        full_name = character(),
        biotype = character(),
        promoter_id = integer(),
        stringsAsFactors = FALSE
    )

    for(i in 1:2) {
        tx_id <- transcript_ids[[i]]
        tx_idx <- grep(paste0("^", tx_id, " "), names(fasta_all))

        if(length(tx_idx) == 0) {
            warning("Transcript ", tx_id, " not found in FASTA file")
            return(NULL)
        }

        selected_seqs <- c(selected_seqs, fasta_all[tx_idx])
        biotype <- sub(".*transcript_biotype:", "", sub(" gene_symbol.*", "", names(fasta_all[tx_idx])))

        transcript_info <- rbind(transcript_info, data.frame(
            transcript_id = tx_id,
            full_name = names(fasta_all[tx_idx]),
            biotype = biotype,
            promoter_id = selected_promoters$promoterId[i],
            stringsAsFactors = FALSE
        ))
    }

    return(list(
        selected_seqs = selected_seqs,
        transcript_info = transcript_info,
        promoter_info = selected_promoters
    ))
}

### Main function ===============================================================
### V10: separate gene_seed and data_seed for reproducible gene parameters across N
simulate_gene_expression_from_promoters_countmat <- function(
    promoter_info_file,
    n_genes = 10,
    n_samples = NULL,
    fasta_file = "/home/users/nus/e1327989/YuhaoDong/reference/Homo_sapiens.GRCh38.cdna.all.fa",
    out_dir = "test_output",
    mean_expression = 1000,
    effect_type_probs = c(
        zero_proportion = 0.15,
        one_proportion = 0.15,
        middle_value = 0.30,
        mixed_effect = 0.40
    ),
    gene_seed = 42,
    data_seed = 123,
    effect_scale = "prob"
) {

    if(is.null(n_samples)) {
        n_samples <- sample(100:500, 1)
        cat("Randomly selected sample size:", n_samples, "\n")
    }

    cat(sprintf("Gene seed: %d (fixed across all N) | Data seed: %d + %d = %d\n",
                gene_seed, data_seed, n_samples, data_seed + n_samples))

    # ========================================================================
    # PHASE 1: Gene parameters (gene_seed -- FIXED across all sample sizes)
    # ========================================================================
    set.seed(gene_seed)

    # 1. Read promoter info
    cat("Reading promoter information...\n")
    if(is.character(promoter_info_file)) {
        if(grepl("\\.RData$", promoter_info_file)) {
            temp_env <- new.env()
            load(promoter_info_file, envir = temp_env)
            promoter_info <- as.data.frame(temp_env[[ls(temp_env)[1]]])
        } else {
            promoter_info <- read.table(promoter_info_file, header = TRUE, stringsAsFactors = FALSE)
        }
    } else if(is.data.frame(promoter_info_file)) {
        promoter_info <- promoter_info_file
    } else {
        stop("promoter_info_file must be a file path or data.frame")
    }

    required_cols <- c("promoterId", "geneId", "seqnames", "start", "strand", "txId")
    missing_cols <- required_cols[!required_cols %in% colnames(promoter_info)]
    if(length(missing_cols) > 0) {
        stop("Promoter info missing columns: ", paste(missing_cols, collapse=", "))
    }

    # 2. Find eligible genes (>= 2 promoters)
    gene_counts <- table(promoter_info$geneId)
    eligible_genes <- names(gene_counts[gene_counts >= 2])

    if(length(eligible_genes) < n_genes) {
        warning("Only ", length(eligible_genes), " genes have >= 2 promoters. Using all.")
        n_genes <- length(eligible_genes)
    }

    # Select genes (gene_seed)
    selected_genes <- sample(eligible_genes, n_genes)

    # 3. Read FASTA file (deterministic, no RNG)
    cat("Reading FASTA file...\n")
    fasta_all <- readDNAStringSet(fasta_file)

    # 4. Select transcripts from different promoters (gene_seed)
    cat("Selecting transcripts from different promoters...\n")
    selected_seqs <- DNAStringSet()
    transcript_info_list <- list()
    selected_promoters_list <- list()

    for(gene_id in selected_genes) {
        cat("Processing gene:", gene_id, "\n")
        selection <- select_transcripts_from_promoters(gene_id, promoter_info, fasta_all)

        if(is.null(selection)) {
            warning("Cannot select transcripts for gene: ", gene_id)
            next
        }

        selection$transcript_info$gene_id <- gene_id
        selected_seqs <- c(selected_seqs, selection$selected_seqs)
        transcript_info_list[[gene_id]] <- selection$transcript_info
        selected_promoters_list[[gene_id]] <- selection$promoter_info

        cat("\nSelected transcripts for gene", gene_id, ":\n")
        print(selection$transcript_info[, c("transcript_id", "biotype", "promoter_id")])
    }

    n_trans <- length(selected_seqs)
    cat("\nSelected transcripts summary:\n")
    cat("Total genes:", length(transcript_info_list), "\n")
    cat("Total transcripts:", n_trans, "\n")

    writeXStringSet(selected_seqs, "selected_transcripts.fa")

    # 5. Generate effect parameters for each gene (gene_seed)
    cat("\nGenerating expression parameters (gene_seed)...\n")

    gene_transcript_indices <- list()
    current_idx <- 1
    for(gene_id in names(transcript_info_list)) {
        n_gene_trans <- nrow(transcript_info_list[[gene_id]])
        if(n_gene_trans >= 2) {
            gene_transcript_indices[[gene_id]] <- current_idx:(current_idx+1)
            current_idx <- current_idx + n_gene_trans
        } else {
            gene_transcript_indices[[gene_id]] <- current_idx
            current_idx <- current_idx + n_gene_trans
        }
    }

    tx_lengths <- width(selected_seqs)

    # Generate ALL gene parameters under gene_seed
    gene_params_list <- list()
    for(gene_id in names(transcript_info_list)) {
        if(length(gene_transcript_indices[[gene_id]]) < 2) next

        effect_type <- sample(names(effect_type_probs),
                              size = 1,
                              prob = effect_type_probs)
        params <- generate_random_params(effect_type = effect_type)
        gene_params_list[[gene_id]] <- params
    }

    cat(sprintf("Gene parameters generated for %d genes\n", length(gene_params_list)))
    cat(sprintf("Effect type distribution: %s\n",
        paste(names(table(sapply(gene_params_list, function(x) x$effect_type))),
              table(sapply(gene_params_list, function(x) x$effect_type)),
              sep="=", collapse=", ")))

    # ========================================================================
    # PHASE 2: Sample data (data_seed + n_samples -- VARIES per sample size)
    # ========================================================================
    set.seed(data_seed + n_samples)

    cat(sprintf("\nGenerating sample data (data_seed=%d + N=%d)...\n", data_seed, n_samples))

    countmat <- matrix(mean_expression, nrow=n_trans, ncol=n_samples)
    sim_results <- list()

    for(gene_id in names(gene_params_list)) {
        params <- gene_params_list[[gene_id]]

        # Generate transcript ratios (uses data_seed RNG)
        expression_results <- generate_transcript_ratios(n_samples, params, effect_scale = effect_scale)
        genotypes <- expression_results$genotypes
        tx1_props <- expression_results$transcript1_props
        tx2_props <- expression_results$transcript2_props

        idx1 <- gene_transcript_indices[[gene_id]][1]
        idx2 <- gene_transcript_indices[[gene_id]][2]

        # Length-corrected count matrix
        L1 <- tx_lengths[idx1]
        L2 <- tx_lengths[idx2]
        Lmean <- (L1 + L2) / 2
        countmat[idx1, ] <- round(mean_expression * tx1_props * L1 / Lmean)
        countmat[idx2, ] <- round(mean_expression * tx2_props * L2 / Lmean)

        sim_results[[gene_id]] <- list(
            params = params,
            genotypes = genotypes,
            transcript1_props = tx1_props,
            transcript2_props = tx2_props,
            directions = expression_results$directions,
            scale_factor = expression_results$scale_factor,
            effect_size_adj = expression_results$effect_size_adj,
            true_effects = expression_results$true_effects
        )
    }

    # 6. Run polyester simulation
    cat("\nRunning polyester simulation...\n")
    unlink(out_dir, recursive = TRUE)
    dir.create(out_dir, showWarnings = FALSE)

    simulate_experiment_countmat(
        fasta = "selected_transcripts.fa",
        readmat = countmat,
        outdir = out_dir,
        paired = TRUE,
        strand_specific = FALSE,
        error_rate = 0.005,
        seed = data_seed + n_samples  # use data_seed for polyester too
    )

    # 7. Save simulation info
    sim_info <- list(
        genes = selected_genes,
        transcripts = transcript_info_list,
        promoter_info = promoter_info,
        selected_promoters = selected_promoters_list,
        simulation_params = list(
            n_transcripts = n_trans,
            n_samples = n_samples,
            mean_expression = mean_expression,
            gene_seed = gene_seed,
            data_seed = data_seed
        ),
        sim_results = sim_results,
        countmat = countmat,
        simulation_time = Sys.time()
    )
    saveRDS(sim_info, file.path(out_dir, "simulation_info.rds"))

    # Create extended transcript info
    tx_info_extended <- data.frame()

    for(gene_id in names(transcript_info_list)) {
        gene_transcripts <- transcript_info_list[[gene_id]]
        if(gene_id %in% names(sim_results)) {
            gene_params <- sim_results[[gene_id]]$params
            gene_dirs <- sim_results[[gene_id]]$directions
            gene_true_eff <- sim_results[[gene_id]]$true_effects

            for(i in 1:nrow(gene_transcripts)) {
                tx_row <- data.frame(
                    transcript_id = gene_transcripts$transcript_id[i],
                    gene_id = gene_id,
                    promoter_id = gene_transcripts$promoter_id[i],
                    biotype = gene_transcripts$biotype[i],
                    effect_type = gene_params$effect_type,
                    effect_size = sim_results[[gene_id]]$effect_size_adj,
                    true_mu_effect = gene_true_eff$mu,
                    true_nu_effect = gene_true_eff$nu,
                    true_tau_effect = gene_true_eff$tau,
                    maf = gene_params$maf,
                    p_zero = gene_params$p_zero,
                    p_one = gene_params$p_one,
                    middle_mean = gene_params$middle_mean,
                    middle_sd = gene_params$middle_sd,
                    scale_factor = sim_results[[gene_id]]$scale_factor,
                    stringsAsFactors = FALSE
                )
                tx_row$transcript_position <- i
                tx_row$pair_id <- gene_id
                tx_info_extended <- rbind(tx_info_extended, tx_row)
            }
        }
    }

    write.table(tx_info_extended,
                file = file.path(out_dir, "sim_tx_info_extended.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    tx_info_basic <- data.frame(
        transcript_id = tx_info_extended$transcript_id,
        gene = tx_info_extended$gene_id,
        length = sapply(selected_seqs, length),
        effective_length = sapply(selected_seqs, length) - 75 + 1,
        paired_end = TRUE,
        promoter_id = tx_info_extended$promoter_id,
        effect_type = tx_info_extended$effect_type,
        effect_size = tx_info_extended$effect_size,
        maf = tx_info_extended$maf,
        stringsAsFactors = FALSE
    )

    write.table(tx_info_basic,
                file = file.path(out_dir, "sim_tx_info.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    sample_size_info <- data.frame(
        simulation_time = Sys.time(),
        n_samples = n_samples,
        n_genes = length(transcript_info_list),
        n_transcripts = n_trans,
        mean_expression = mean_expression,
        gene_seed = gene_seed,
        data_seed = data_seed
    )

    write.table(sample_size_info,
                file = file.path(out_dir, "sample_size_info.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    # 8. Check output
    out_files <- list.files(out_dir, pattern = "*.fasta")
    cat("\nGenerated output files:\n")
    print(head(out_files))
    cat("Total files:", length(out_files), "\n")

    # 9. Generate VCF file
    cat("\nGenerating VCF file...\n")
    vcf_files <- generate_vcf(sim_info,
                             output_file = file.path(out_dir, "genotypes.vcf"))
    sim_info$genotype_files <- vcf_files

    return(sim_info)
}

### STAR function ===============================================================
generate_star_script <- function(sim_info,
                               star_index_dir,
                               output_dir = "star_output") {
    sim_output_dir <- "test_output"

    if(!is.null(sim_info$genotype_files) && !is.null(sim_info$genotype_files$vcf_file)) {
        sim_output_dir <- dirname(sim_info$genotype_files$vcf_file)
    }

    script <- paste0("#!/bin/bash\n\n",
                    "# Created by R at ", Sys.time(), "\n\n",
                    "# Create output directory\n",
                    "mkdir -p ", output_dir, "\n",
                    "cd ", output_dir, "\n\n",
                    "# Convert FASTA to FASTQ\n",
                    "for f in ../", sim_output_dir, "/sample_*_1.fasta; do\n",
                    "    sample=$(basename $f _1.fasta)\n",
                    "    echo \"Converting $sample to FASTQ format...\"\n",
                    "    seqtk seq -F '#' $f > ${sample}_1.fastq\n",
                    "    seqtk seq -F '#' ${f%_1.fasta}_2.fasta > ${sample}_2.fastq\n",
                    "done\n\n",
                    "# Run STAR - one sample at a time\n",
                    "for sample_prefix in $(ls sample_*_1.fastq | sed 's/_1\\.fastq$//');\n",
                    "do\n",
                    "    echo \"Processing $sample_prefix with STAR...\"\n",
                    "    STAR \\\n",
                    "        --runThreadN 8 \\\n",
                    "        --genomeDir ", star_index_dir, " \\\n",
                    "        --readFilesIn ${sample_prefix}_1.fastq ${sample_prefix}_2.fastq \\\n",
                    "        --outFileNamePrefix ${sample_prefix}_ \\\n",
                    "        --outSAMtype BAM SortedByCoordinate \\\n",
                    "        --quantMode TranscriptomeSAM GeneCounts \\\n",
                    "        --outSAMattributes NH HI AS nM NM MD \\\n",
                    "        --alignSJoverhangMin 8 \\\n",
                    "        --alignSJDBoverhangMin 1\n",
                    "done\n")

    writeLines(script, "run_star.sh")
    system(paste0("chmod +x run_star.sh"))
    cat("Created STAR running script: run_star.sh\n")
}
