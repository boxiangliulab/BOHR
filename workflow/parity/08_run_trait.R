## Real-data parity of bohr_trait(), step B (bohrdev R with BOHR installed): one tissue x trait from 07_extract_trait_inputs.R,
## with the absolute activity (option B labels) and the paper's canonical mirror members.
## Usage: Rscript 08_run_trait.R <inputs.rds> <out.rds> <ncores> [lib]
## Environment: BOHR_FIRTH_NULL_START=FALSE runs the paper's engine exactly (default TRUE, the package default; stage 3b-1).
a <- commandArgs(TRUE); IN <- a[1]; OUT <- a[2]; NC <- as.integer(a[3]); LIB <- if (length(a) > 3) a[4] else ""
if (nzchar(LIB)) library(BOHR, lib.loc = LIB) else library(BOHR)
x <- readRDS(IN)
NS <- as.logical(Sys.getenv("BOHR_FIRTH_NULL_START", "TRUE"))
t0 <- proc.time()[["elapsed"]]; c0 <- sum(proc.time()[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE)
res <- bohr_trait(x$phenotype, x$traits, x$covariates, promoters = x$promoters, major_score = x$major_score, absolute = x$absolute, ncores = NC,
                  control = bohr_control(firth_null_start = NS))
res$wall <- proc.time()[["elapsed"]] - t0
res$cpu <- sum(proc.time()[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE) - c0
res$tissue <- x$tissue; res$trait <- x$trait
saveRDS(res, OUT)
print(res); cat(sprintf("wall %.1f s on %d cores; CPU %.0f s; firth_null_start %s\n", res$wall, NC, res$cpu, NS))
