## Real-data parity, step B (bohrdev R): run BOHR::bohr_cis() on the extracted inputs.
## Usage: Rscript run_bohr_cis.R <inputs.rds> <out.rds> <ncores> [max_promoters] [lib]
## Environment: BOHR_FIRTH_NULL_START=FALSE runs the paper's engine exactly (default TRUE, the package default; stage 3b-1).
a <- commandArgs(TRUE); IN <- a[1]; OUT <- a[2]; NC <- as.integer(a[3]); MAXP <- if (length(a) > 3) as.integer(a[4]) else 0L
LIB <- if (length(a) > 4) a[5] else ""
if (nzchar(LIB)) library(BOHR, lib.loc = LIB) else library(BOHR)
x <- readRDS(IN)
NS <- as.logical(Sys.getenv("BOHR_FIRTH_NULL_START", "TRUE")); ctl <- bohr_control(firth_null_start = NS)
ph <- x$phenotype
if (MAXP > 0) ph <- ph[seq_len(MAXP), , drop = FALSE]
t0 <- proc.time()[["elapsed"]]; c0 <- sum(proc.time()[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE)
## The variant set is the stored run's (VCF pre-filtered at MAF >= 0.05 over all 944 donors). The covariate-count warning (5 covariates > floor(n/30) for small tissues) is expected and recorded.
warn <- character(0)
res <- withCallingHandlers(bohr_cis(ph, x$genotype, x$variants, x$promoters, covariates = x$covariates, major_score = x$major_score,
                                    ncores = NC, all_stats = TRUE, control = ctl),
                           warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") })
res$warnings <- warn
el <- proc.time()[["elapsed"]] - t0
res$cpu <- sum(proc.time()[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE) - c0   # includes the forked workers
res$wall <- el; res$ncores <- NC; res$host <- Sys.info()[["nodename"]]; res$bohr_version <- as.character(utils::packageVersion("BOHR"))
saveRDS(res, OUT)
print(res); cat(sprintf("wall %.1f s on %d cores; CPU %.0f s; %d SNP rows; firth_null_start %s\n", el, NC, res$cpu, nrow(res$snp), NS))
if (length(warn)) cat("warnings:", warn, sep = "\n  ")
