## Real-data parity of bohr_trait(), step B (bohrdev R with BOHR installed): one tissue x trait from 07_extract_trait_inputs.R,
## with the absolute activity (option B labels) and the paper's canonical mirror members.
## Usage: Rscript 08_run_trait.R <inputs.rds> <out.rds> <ncores> [lib]
a <- commandArgs(TRUE); IN <- a[1]; OUT <- a[2]; NC <- as.integer(a[3]); LIB <- if (length(a) > 3) a[4] else ""
if (nzchar(LIB)) library(BOHR, lib.loc = LIB) else library(BOHR)
x <- readRDS(IN)
t0 <- proc.time()[["elapsed"]]; c0 <- sum(proc.time()[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE)
res <- bohr_trait(x$phenotype, x$traits, x$covariates, promoters = x$promoters, major_score = x$major_score, absolute = x$absolute, ncores = NC)
res$wall <- proc.time()[["elapsed"]] - t0
res$cpu <- sum(proc.time()[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE) - c0
res$tissue <- x$tissue; res$trait <- x$trait
saveRDS(res, OUT)
print(res); cat(sprintf("wall %.1f s on %d cores; CPU %.0f s\n", res$wall, NC, res$cpu))
