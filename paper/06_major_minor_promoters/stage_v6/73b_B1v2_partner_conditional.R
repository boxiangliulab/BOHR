## 73b_B1v2_partner_conditional.R -- B1 refit. Same donors, covariates and genotype rule as 73.
## Changes from 73:
##  1. Designated partner. In "two-promoter" genes (n_prom_gene = 2 among the tested promoters) signals can have extra
##     non-internal promoters in the absolute matrix (untested, low activity), so 73's "others" sum (A_O) is not the
##     designated partner. Here the partner P (mirror table) is fitted on its own and the extra promoters are summed
##     separately (X). A_O (all others) is kept: it is what usage = 1 (p1) is defined on.
##  2. Conditional fits. The genotypes of the other retained signals (backward_keep) of the same unit
##     (tissue x promoter x type) enter as covariates (multi-signal units).
##  3. Level part (OLS on donors with A > 0) for own and partner whenever >= 30 donors are positive, with no zero-fraction
##     gate (= the full OLS when the promoter is always on); on-state (logistic) when 2-98% of donors are zero.
##  4. Effects kept as b, se and SD of the response, for equivalence tests ("flat" = CI inside +-delta SD, not P > 0.05).
##  5. Usage of the own and of the partner promoter refitted, so the partner's usage direction is measured, not assumed.
## Output: out/B1v2_multipromoter/B1v2_fits.tsv.gz (one row per signal)
suppressPackageStartupMessages({ library(data.table); library(VariantAnnotation); library(parallel) })
source(file.path(Sys.getenv("F5S", "/home/e1327989/f5_stage_20261106_v6"), "code/00_helpers.R"))
OUT <- file.path(OUTROOT, "B1v2_multipromoter"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
cat("start", format(Sys.time()), "\n")

sig <- fr(file.path(IN, "delivery/signals.tsv.gz"), select = c("tissue", "promoter", "type", "k", "snp", "chr", "backward_keep", "loo_est"))
fl  <- fr(file.path(IN, "mirror/promoter_flags.tsv.gz"))[, .(tissue, promoter, gene, n_prom_gene, mirror_partner, exact_mirror)]
sig <- merge(sig, fl, by = c("tissue", "promoter")); stopifnot(nrow(sig) > 0)   # the per-signal table of delivery v6 must not be empty
sig[, pos := as.integer(sub("^chr[^_]+_([0-9]+)_.*$", "\\1", snp))]
sig[, unit := paste(tissue, promoter, type, sep = "|")]
sig[, bk := as.logical(backward_keep)]
sig[, cond_snps := { s <- snp; b <- bk; vapply(seq_along(s), function(i) paste(setdiff(unique(s[-i][b[-i] %in% TRUE]), s[i]), collapse = ";"), "") }, by = unit]
cat("signals with a conditioning set:", sig[cond_snps != "", .N], "\n")

e <- new.env(); load(file.path(IN, "v11/Promoter_info_all.RData"), envir = e); ia <- as.data.table(e$info_all, keep.rownames = "pid"); rm(e)
x <- readRDS(file.path(IN, "v11/tabs_att_filter_In.rds")); abs <- x$abs; rm(x); invisible(gc())
ia <- ia[internalPromoter == FALSE & pid %in% colnames(abs)]
genes <- unique(sig$gene); keep <- ia[geneId %in% genes, pid]
abs <- abs[, keep, drop = FALSE]; invisible(gc())
gp_of <- split(ia[geneId %in% genes, pid], ia[geneId %in% genes, geneId])
TT <- fr(file.path(IN, "freeze/tissues.tsv"))[def == "freeze"]

gt_dosage <- function(gt) { D <- matrix(NA_integer_, nrow(gt), ncol(gt), dimnames = dimnames(gt))
  D[gt %in% c("0|0", "0/0")] <- 0L; D[gt %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L; D[gt %in% c("1|1", "1/1")] <- 2L; D }
impute_modal <- function(D) { for (i in which(rowSums(is.na(D)) > 0)) { tb <- tabulate(D[i, ] + 1L, nbins = 3L); D[i, is.na(D[i, ])] <- which.max(tb) - 1L }; D }
## OLS of y on X (last column = dosage) over finite y: b, se, z of the dosage, SD of y, n
fit1 <- function(X, y, nmin = 30L) { ok <- is.finite(y); n <- sum(ok); na <- c(b = NA_real_, se = NA_real_, z = NA_real_, sdy = NA_real_, n = n)
  if (n < nmin) return(na); Xo <- X[ok, , drop = FALSE]; yo <- y[ok]; p <- ncol(Xo); if (n <= p + 2) return(na)
  f <- tryCatch(lm.fit(Xo, yo), error = function(e) NULL); if (is.null(f) || f$rank < p) return(na)
  s2 <- sum(f$residuals^2) / (n - p); Rinv <- backsolve(qr.R(f$qr), diag(p)); j <- which(f$qr$pivot == p)
  b <- unname(f$coefficients[p]); se <- sqrt(sum(Rinv[j, ]^2) * s2); c(b = b, se = se, z = b / se, sdy = sd(yo), n = n) }
glm_z <- function(X, y) { if (mean(y) < 0.02 || mean(y) > 0.98) return(NA_real_)
  f <- tryCatch(suppressWarnings(glm.fit(X, y, family = binomial())), error = function(e) NULL); if (is.null(f) || !f$converged) return(NA_real_)
  s <- tryCatch(summary.glm(structure(f, class = c("glm", "lm")))$coefficients, error = function(e) NULL)
  if (is.null(s) || nrow(s) < ncol(X)) return(NA_real_); s[nrow(s), 3] }
put <- function(v, nm) setNames(as.list(v), paste0(c("b_", "se_", "z_", "sd_", "n_"), nm))

run_tissue <- function(t) {
  td <- readRDS(file.path(IN, "freeze", paste0(t, ".rds"))); K <- TT[tissue == t, K_C3]
  st <- sig[tissue == t]; if (!nrow(st)) return(NULL)
  Z <- scale(cbind(td$pcs, td$peer[, seq_len(K), drop = FALSE])); C <- cbind(1, Z)
  G <- list()
  for (ch in unique(st$chr)) {
    s <- st[chr == ch]; prm <- ScanVcfParam(which = GRanges(ch, IRanges(s$pos, s$pos)), samples = td$U, fixed = NA, info = NA, geno = "GT")
    v <- readVcf(TabixFile(file.path(IN, "vcf", paste0(ch, ".vcf.gz"))), "hg38", prm)
    gt <- geno(v)$GT[, td$U, drop = FALSE]; gt <- gt[!duplicated(rownames(gt)), , drop = FALSE]
    D <- impute_modal(gt_dosage(gt)); for (id in intersect(rownames(D), s$snp)) G[[id]] <- D[id, ]
  }
  Asub <- abs[td$row_id, , drop = FALSE]
  out <- vector("list", nrow(st))
  for (i in seq_len(nrow(st))) {
    r <- st[i]; key <- list(tissue = t, promoter = r$promoter, type = r$type, k = r$k)
    g <- G[[r$snp]]; if (is.null(g)) { out[[i]] <- c(key, status = "no genotype"); next }
    pr <- gp_of[[r$gene]]; if (!(r$promoter %in% pr)) { out[[i]] <- c(key, status = "promoter not in abs"); next }
    A <- as.matrix(Asub[, pr, drop = FALSE]); At <- A[, r$promoter]
    P <- if (r$n_prom_gene == 2 && !is.na(r$mirror_partner) && r$mirror_partner %in% pr) r$mirror_partner else NA_character_
    oth <- setdiff(pr, r$promoter); xs <- setdiff(oth, P)
    AO <- if (length(oth)) rowSums(A[, oth, drop = FALSE]) else rep(0, nrow(A))
    AP <- if (!is.na(P)) A[, P] else rep(NA_real_, nrow(A))
    AX <- if (length(xs)) rowSums(A[, xs, drop = FALSE]) else rep(0, nrow(A))
    AT <- log2(rowSums(2^A - 1) + 1); tot <- rowSums(A)
    ut <- ifelse(tot > 0, At / tot, NA_real_); uP <- if (!is.na(P)) ifelse(tot > 0, AP / tot, NA_real_) else rep(NA_real_, nrow(A))
    ## conditioning genotypes (other retained signals of the unit); dropped if they make the design rank deficient
    cs <- if (nzchar(r$cond_snps)) intersect(strsplit(r$cond_snps, ";", fixed = TRUE)[[1]], names(G)) else character(0)
    Gc <- if (length(cs)) do.call(cbind, lapply(cs, function(s) as.numeric(G[[s]]))) else NULL
    X <- cbind(C, Gc, dose = as.numeric(g)); cond_dropped <- FALSE
    if (qr(X)$rank < ncol(X)) { X <- cbind(C, dose = as.numeric(g)); cond_dropped <- length(cs) > 0 }
    if (qr(X)$rank < ncol(X)) { out[[i]] <- c(key, status = "rank deficient"); next }
    gperm <- { set.seed(1e5 + i); sample(as.numeric(g)) }; Xp <- X; Xp[, ncol(Xp)] <- gperm
    res <- c(key, status = "ok", n = nrow(A), n_prom_abs = length(pr), partner = P, n_extra = length(xs), n_cond = if (cond_dropped) 0L else length(cs),
             cond_dropped = cond_dropped, zero_t = mean(At == 0), zero_P = if (!is.na(P)) mean(AP == 0) else NA_real_, extra_on = mean(AX > 0),
             put(fit1(X, At), "t"), put(fit1(X, AP), "P"), put(fit1(X, AX), "X"), put(fit1(X, AO), "O"), put(fit1(X, AT), "T"),
             put(fit1(X, ut), "ut"), put(fit1(X, uP), "uP"),
             put(fit1(X, ifelse(At > 0, At, NA_real_)), "levt"), put(fit1(X, ifelse(AP > 0, AP, NA_real_)), "levP"),
             zon_t = glm_z(X, as.numeric(At > 0)), zon_P = if (!is.na(P)) glm_z(X, as.numeric(AP > 0)) else NA_real_,
             zperm_t = fit1(Xp, At)[["z"]], zperm_P = fit1(Xp, AP)[["z"]], zperm_O = fit1(Xp, AO)[["z"]], zperm_T = fit1(Xp, AT)[["z"]],
             zperm_levP = fit1(Xp, ifelse(AP > 0, AP, NA_real_))[["z"]], zperm_levt = fit1(Xp, ifelse(At > 0, At, NA_real_))[["z"]])
    out[[i]] <- res
  }
  cat(t, nrow(st), "signals done", format(Sys.time()), "\n")
  rbindlist(out, fill = TRUE)
}
RES <- rbindlist(mclapply(sort(unique(sig$tissue)), function(t) tryCatch(run_tissue(t), error = function(e) { cat("ERROR", t, conditionMessage(e), "\n"); NULL }),
                          mc.cores = 2, mc.preschedule = FALSE), fill = TRUE)
fwrite(RES, file.path(OUT, "B1v2_fits.tsv.gz"), sep = "\t")
cat("rows", nrow(RES), "; status:\n"); print(RES[, .N, by = status])
cat("two-promoter signals with a partner fit:", RES[!is.na(partner), .N], "; with extra promoters:", RES[!is.na(partner) & n_extra > 0, .N],
    "; conditional fits:", RES[n_cond > 0, .N], "; conditioning dropped (rank):", RES[cond_dropped == TRUE, .N], "\n")
cat("end", format(Sys.time()), "\n")
