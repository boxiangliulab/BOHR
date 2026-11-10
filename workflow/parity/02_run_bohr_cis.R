## Real-data parity, step B (bohrdev R): run BOHR::bohr_cis() on the extracted inputs.
## Usage: Rscript run_bohr_cis.R <inputs.rds> <out.rds> <ncores> [max_promoters] [lib]
a <- commandArgs(TRUE); IN <- a[1]; OUT <- a[2]; NC <- as.integer(a[3]); MAXP <- if (length(a) > 3) as.integer(a[4]) else 0L
LIB <- if (length(a) > 4) a[5] else ""
if (nzchar(LIB)) library(BOHR, lib.loc = LIB) else library(BOHR)
x <- readRDS(IN)
ph <- x$phenotype
if (MAXP > 0) ph <- ph[seq_len(MAXP), , drop = FALSE]
t0 <- proc.time()[["elapsed"]]; c0 <- sum(proc.time()[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE)
res <- bohr_cis(ph, x$genotype, x$variants, x$promoters, covariates = x$covariates, major_score = x$major_score,
                ncores = NC, all_stats = TRUE)
el <- proc.time()[["elapsed"]] - t0
res$cpu <- sum(proc.time()[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE) - c0   # includes the forked workers
res$wall <- el; res$ncores <- NC; res$host <- Sys.info()[["nodename"]]; res$bohr_version <- as.character(utils::packageVersion("BOHR"))
saveRDS(res, OUT)
print(res); cat(sprintf("wall %.1f s on %d cores; CPU %.0f s; %d SNP rows\n", el, NC, res$cpu, nrow(res$snp)))
