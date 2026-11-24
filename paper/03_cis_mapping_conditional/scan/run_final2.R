## Full rerun (configuration chosen by the pilot: C3 + Firth for p0/p1). One (definition, tissue, arm, chr).
## Usage: Rscript run_final.R <def: freeze|pilot> <tissue> <arm: real|perm0> <chr> <ncores>
## Engine and scheduling as run_pilot.R (split likelihood, sigma ~ 1, two-phase promoter/SNP-chunk
## parallelism). For promoters eligible for p0 (f0 >= 0.35, BEINF/BEZI) or p1 (f1 >= 0.35, BEINF/BEOI) the
## Firth penalised LRT for nu/tau is added (bohr2_addon3.R: alternating full/restricted search,
## multi-start); these are the p-values
## used for p0/p1 discovery. Output (atomic): results/final/<def>/<tissue>/<arm>/<chr>.rds
suppressPackageStartupMessages(library(parallel))
CODE <- Sys.getenv("BOHR2_CODE", "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/code")
source(file.path(CODE, "bohr2_data.R")); source(file.path(CODE, "bohr2_engine.R")); source(file.path(CODE, "bohr2_addon3.R"))
a <- commandArgs(TRUE); def <- a[1]; tis <- a[2]; arm <- a[3]; chr <- a[4]; NC <- as.integer(a[5])
out_dir <- file.path(Sys.getenv("FINAL_ROOT", file.path(V2, "results", "final")), def, tis, arm)
out <- file.path(out_dir, paste0(chr, ".rds"))
if (file.exists(out)) { cat("exists:", out, "\n"); quit(save = "no") }
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
t_start <- Sys.time()

td <- readRDS(file.path(V2, "work", "final", def, paste0(tis, ".rds"))); K <- td$K
if (arm == "real") { A <- td$act; P <- td$peer } else {
  m <- fread(file.path(V2, "perm_maps", "final", def, paste0(tis, ".", arm, ".perm_map.tsv")))
  stopifnot(setequal(m$slot, td$U), setequal(m$donor, td$U), !anyDuplicated(m$donor))
  d <- setNames(m$donor, m$slot)[td$U]; A <- td$act[d, , drop = FALSE]; P <- td$peer[d, , drop = FALSE]; rownames(A) <- rownames(P) <- td$U }
C_all <- cov_matrix(td, P, K)
prs <- colnames(td$act)[td$info[colnames(td$act), "seqnames"] == chr]
prs <- prs[order(td$info[prs, "start"])]
if (as.integer(Sys.getenv("FINAL_MAXPROM", "0")) > 0) {                         # smoke tests: never into the production tree
  if (Sys.getenv("FINAL_ROOT") == "") stop("FINAL_MAXPROM requires FINAL_ROOT"); prs <- head(prs, as.integer(Sys.getenv("FINAL_MAXPROM"))) }

SNP_COLS <- c("mu_stat", "nu_stat", "tau_stat", "cat_joint_stat", "joint_stat", "joint_df",
              "mu_lrt_logp", "nu_lrt_logp", "tau_lrt_logp", "joint_logp", "mu_wald_logp", "nu_wald_logp", "tau_wald_logp",
              "mu_est", "mu_se", "nu_est", "nu_se", "tau_est", "tau_se", "mu_pe", "p0_pe", "p1_pe",
              "mu_conv", "nu_conv", "tau_conv", "mu_sep", "nu_sep", "tau_sep", "mu_iter", "cat_iter",
              "ll_beta_full", "ll_cat_full", "cat_complete_sep")
FIRTH_COLS <- c("nu_firth_stat", "nu_firth_logp", "nu_firth_est", "nu_firth_se", "nu_firth_conv",
                "tau_firth_stat", "tau_firth_logp", "tau_firth_est", "tau_firth_se", "tau_firth_conv", "firth_rounds")
CHUNK <- 100L
prep_promoter <- function(p) {
  t0 <- proc.time()[["elapsed"]]
  base <- list(promoter = p, pos = td$info[p, "start"], gene = td$info[p, "geneId"], family = NA_character_, skip = NA_character_,
               n = NA_integer_, n0 = NA_integer_, n1 = NA_integer_, n_cov = NA_integer_, n_cov_dropped = NA_integer_,
               elig_mu = FALSE, elig_p0 = FALSE, elig_p1 = FALSE, null_beta_conv = NA, null_cat_conv = NA,
               null_ll_beta = NA_real_, null_ll_cat = NA_real_, n_snp_window = NA_integer_, n_snp_tested = 0L,
               n_snp_error = 0L, n_firth_error = 0L, geno_missing = NA_integer_, sec = NA_real_)
  z <- prepare_y(A[, p], ncol(C_all) - 1L)
  if (!is.na(z$skip)) { base$skip <- z$skip; base$sec <- proc.time()[["elapsed"]] - t0; return(list(prom = base)) }
  if (!any(z$y > 0 & z$y < 1)) { base$skip <- "no_interior"; base$sec <- proc.time()[["elapsed"]] - t0; return(list(prom = base)) }
  s <- names(z$y); y <- unname(z$y); Cs <- orth_design(C_all[s, , drop = FALSE])
  n <- length(y); n0 <- sum(y == 0); n1 <- sum(y == 1)
  base$family <- z$family; base$n <- n; base$n0 <- n0; base$n1 <- n1
  base$n_cov <- ncol(Cs) - 1L; base$n_cov_dropped <- length(attr(Cs, "dropped"))
  base$elig_mu <- (n - n0 - n1) / n >= 0.70
  base$elig_p0 <- n0 / n >= 0.35 && z$family %in% c("BEINF", "BEZI")
  base$elig_p1 <- n1 / n >= 0.35 && z$family %in% c("BEINF", "BEOI")
  nul <- promoter_null(y, Cs, z$family)
  base$null_beta_conv <- if (!is.null(nul$beta)) nul$beta$converged else NA
  base$null_cat_conv  <- if (!is.null(nul$cat)) nul$cat$converged else NA
  base$null_ll_beta   <- if (!is.null(nul$beta)) nul$beta$ll else NA_real_
  base$null_ll_cat    <- if (!is.null(nul$cat)) nul$cat$ll else NA_real_
  pos <- td$info[p, "start"]; geno <- read_geno(chr, pos - FLANK, pos + FLANK - 1L, td$U)
  if (is.null(geno)) { base$skip <- "no_snp"; base$sec <- proc.time()[["elapsed"]] - t0; return(list(prom = base)) }
  idx <- cis_index(geno, pos); base$n_snp_window <- length(idx); base$geno_missing <- geno$n_missing
  G <- geno$G[s, idx, drop = FALSE]
  cnt <- apply(G, 2, function(x) tabulate(x + 1L, nbins = 3L)); keep <- which(apply(cnt, 2, min) >= MIN_GENO)
  if (!length(keep)) { base$skip <- "no_snp"; base$sec <- proc.time()[["elapsed"]] - t0; return(list(prom = base)) }
  base$n_snp_pass <- length(keep); base$sec <- proc.time()[["elapsed"]] - t0
  list(prom = base, y = y, Cs = Cs, nul = nul, G = G[, keep, drop = FALSE], cnt = cnt[, keep, drop = FALSE],
       snp_pos = geno$start[idx[keep]], firth = base$elig_p0 || base$elig_p1)
}
prep_one <- function(p) tryCatch(prep_promoter(p), error = function(e) list(prom = list(promoter = p, skip = paste("ERROR:", conditionMessage(e)))))
prep <- mclapply(prs, prep_one, mc.cores = NC, mc.preschedule = FALSE)
dead <- which(!vapply(prep, is.list, TRUE))                   # a worker died (e.g. killed) -> retry serially
for (i in dead) prep[[i]] <- tryCatch(prep_one(prs[i]), error = function(e) NULL)
incomplete <- sum(!vapply(prep, is.list, TRUE))
for (i in which(!vapply(prep, is.list, TRUE))) prep[[i]] <- list(prom = list(promoter = prs[i], skip = "ERROR: worker died twice"))
tasks <- list()
for (i in seq_along(prep)) if (!is.null(prep[[i]]$G)) {
  m <- ncol(prep[[i]]$G); w <- c(BEINF = 3, BEZI = 1, BEOI = 1, BE = 0.7)[[prep[[i]]$prom$family]] * prep[[i]]$prom$n *
    ncol(prep[[i]]$Cs)^1.5 * (if (isTRUE(prep[[i]]$firth)) 4 else 1)
  for (ch in split(seq_len(m), ceiling(seq_len(m) / CHUNK))) tasks[[length(tasks) + 1]] <- list(i = i, j = ch, w = w * length(ch))
}
tasks <- tasks[order(-vapply(tasks, `[[`, 0, "w"))]
run_chunk <- function(t) {
  P <- prep[[t$i]]; t0 <- proc.time()[["elapsed"]]; nerr <- 0L; fe <- 0L
  rows <- lapply(t$j, function(j) {
    r <- tryCatch(snp_tests(P$y, P$Cs, P$G[, j], P$nul), error = function(e) NULL)
    if (is.null(r)) { nerr <<- nerr + 1L; return(NULL) }
    f <- setNames(rep(list(NA), length(FIRTH_COLS)), FIRTH_COLS); ferr <- FALSE
    if (isTRUE(P$firth)) { ff <- tryCatch(firth_snp_tests(P$y, P$Cs, P$G[, j], P$nul), error = function(e) NULL)
      if (is.null(ff)) { ferr <- TRUE; fe <<- fe + 1L; f$nu_firth_conv <- FALSE; f$tau_firth_conv <- FALSE } else f <- ff[FIRTH_COLS] }
    af <- (P$cnt[2, j] + 2 * P$cnt[3, j]) / (2 * length(P$y))
    c(list(snp = colnames(P$G)[j], snp_pos = P$snp_pos[j], alt_af = af, maf = min(af, 1 - af), min_geno = min(P$cnt[, j]), firth_error = ferr),
      r[SNP_COLS], f)
  })
  list(i = t$i, nerr = nerr, ferr = fe, sec = proc.time()[["elapsed"]] - t0, snp = rbindlist(rows[!vapply(rows, is.null, TRUE)], fill = TRUE))
}
chunk_one <- function(t) tryCatch(run_chunk(t), error = function(e) list(i = t$i, nerr = length(t$j), ferr = 0L, sec = NA_real_, snp = NULL))
res <- mclapply(tasks, chunk_one, mc.cores = NC, mc.preschedule = FALSE)
dead <- which(!vapply(res, is.list, TRUE))                    # a worker died -> retry those chunks once
if (length(dead)) res[dead] <- mclapply(tasks[dead], chunk_one, mc.cores = min(NC, length(dead)), mc.preschedule = FALSE)
incomplete <- incomplete + sum(!vapply(res, is.list, TRUE))
res <- res[vapply(res, is.list, TRUE)]
for (r in res) { b <- prep[[r$i]]$prom
  b$n_snp_error <- b$n_snp_error + r$nerr; b$n_firth_error <- b$n_firth_error + r$ferr
  b$n_snp_tested <- b$n_snp_tested + (if (is.null(r$snp)) 0L else nrow(r$snp))
  b$sec <- b$sec + (if (is.finite(r$sec)) r$sec else 0); prep[[r$i]]$prom <- b }
prom <- rbindlist(lapply(prep, `[[`, "prom"), fill = TRUE)
## bookkeeping: every SNP that passed the genotype filter was either tested or counted as an error
if ("n_snp_pass" %in% names(prom)) incomplete <- incomplete + prom[!is.na(n_snp_pass) & n_snp_tested + n_snp_error != n_snp_pass, .N]
snp <- rbindlist(lapply(res, function(r) { x <- r$snp; if (is.null(x) || !nrow(x)) NULL else cbind(promoter = prep[[r$i]]$prom$promoter, x) }), fill = TRUE)
meta <- list(def = def, tissue = tis, K = K, arm = arm, chr = chr, n_samples = length(td$U), n_promoters = length(prs), cores = NC,
             wall_sec = as.numeric(difftime(Sys.time(), t_start, units = "secs")), host = Sys.info()[["nodename"]], date = format(Sys.time()),
             engine = "bohr2_engine.R (split likelihood, sigma ~ 1) + bohr2_addon3.R alternating-search Firth for p0/p1-eligible promoters")
meta$incomplete <- incomplete
if (incomplete > 0) {                                          # never publish a chromosome with lost work
  saveRDS(list(snp = snp, prom = prom, meta = meta), sub("\\.rds$", ".incomplete.rds", out))
  cat("INCOMPLETE:", incomplete, "promoters/chunks lost after retry; not published\n"); quit(save = "no", status = 0) }
tmp <- paste0(out, ".tmp"); saveRDS(list(snp = snp, prom = prom, meta = meta), tmp); file.rename(tmp, out)
cat(sprintf("%s %s K%d %s %s: %d promoters, %d SNP rows, %.1f min\n", def, tis, K, arm, chr, length(prs), nrow(snp), meta$wall_sec / 60))
