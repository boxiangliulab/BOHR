### NO-THRESHOLD copy of v12_sim_run_realfloor_big.R
### BIG copy of v12_sim_run_realfloor.R (built by calibration_perm1000/make_big_M2_runner.R)
### REALFLOOR copy of v12_sim_run.R (built by calibration_perm1000/make_realfloor_generators.R)
#!/usr/bin/env Rscript
### NewSim Method2: Polyester simulation with rBEINF data generation
### Based on V12 proActiv pipeline with two changes:
###   1. rBEINF for generating transcript proportions (instead of rep(0)+rep(1)+rnorm)
###   2. Proportional effect boundary scaling (instead of hard clamping)
###
### Usage: Rscript v12_sim_run.R <n_samples> <sim_dir> [gene_seed] [data_seed] [effect_scale]
### Env: polyester
### =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) stop("Usage: Rscript v12_sim_run.R <n_samples> <sim_dir> [gene_seed] [data_seed] [effect_scale]")

n_samples    <- as.integer(args[1])
sim_dir      <- args[2]
gene_seed    <- ifelse(length(args) >= 3, as.integer(args[3]), 42)
data_seed    <- ifelse(length(args) >= 4, as.integer(args[4]), 123)
effect_scale <- ifelse(length(args) >= 5, args[5], "prob")

cat("=== NewSim Method2 (Polyester + rBEINF): Generate Data ===\n")
cat("Sample size:", n_samples, "\n")
cat("Output directory:", sim_dir, "\n")
cat("Gene seed:", gene_seed, "(fixed across all N)\n")
cat("Data seed:", data_seed, "+ N =", data_seed + n_samples, "\n")
cat("Effect scale:", effect_scale, "\n")

# Create output directory and set working directory
dir.create(sim_dir, recursive = TRUE, showWarnings = FALSE)
setwd(sim_dir)

# Source simulation functions (from method2_polyester directory)
source("/data/projects/11003054/YuhaoDong/GTEx/new_simulation/method2_polyester/v12_simulation_nothr.R")   # NO-THRESHOLD (1200 genes, balanced types as BIG)

# Run simulation
sim_info <- simulate_gene_expression_from_promoters_countmat(
    promoter_info_file = "/home/users/nus/e1327989/YuhaoDong/GTEx/v11_filter_recal/Promoter_info_filter_5.RData",
    n_genes = 1200,                      # BIG: was 300
    n_samples = n_samples,
    fasta_file = "/home/users/nus/e1327989/YuhaoDong/reference/Homo_sapiens.GRCh38.cdna.all.fa",
    out_dir = "test_output",
    mean_expression = 1000,
    effect_type_probs = c(
        zero_proportion = 0.25,          # BIG: was 0.15
        one_proportion = 0.25,           # BIG: was 0.15
        middle_value = 0.25,             # BIG: was 0.30
        mixed_effect = 0.25              # BIG: was 0.40
    ),
    gene_seed = gene_seed,
    data_seed = data_seed,
    effect_scale = effect_scale
)

cat("Step 1 complete: simulation data generated in", sim_dir, "\n")
cat("Number of genes:", length(sim_info$sim_results), "\n")
cat("Number of samples:", n_samples, "\n")
