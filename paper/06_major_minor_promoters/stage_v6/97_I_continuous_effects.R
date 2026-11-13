## 97_I_continuous_effects.R -- follow-ups to a review (2026-11-03 evening).
## I1  F2 within-promoter pairs on continuous effects: standardised per-allele effects at S, N and gene output (b / SD of the
##     response, signed by the S-usage direction) for the boundary and the mu signal of each pair; paired mean / median
##     differences with a gene bootstrap; and "N opposite" split into strict reciprocal (S significant) vs N only.
## I2  Alternative phenotype: linear share n_S / (n_S + n_N) of normalised counts (not log2), and log ratio; does each class's
##     lead move the linear share in the S-usage direction? (two-promoter signals; same donors / covariates / conditioning as 73b)
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(VariantAnnotation); library(parallel) })
source(file.path(Sys.getenv("F5S", "/home/e1327989/f5_stage_20261106_v6"), "code/00_helpers.R"))
OUT <- file.path(OUTROOT, "I_continuous_effects"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE); B <- 2000L; ZC <- 1.96
F <- fr(file.path(OUTROOT, "B1v2_multipromoter/B1v2_fits.tsv.gz")); F[, sid := paste(tissue, promoter, type, k, sep = "|")]
P <- fr(file.path(OUTROOT, "B1v2_multipromoter/B1v2_signal_patterns.tsv.gz"),
        select = c("sid", "gene", "two", "route", "promoter", "partner", "S_is_own", "sS", "patternS", "zS", "zN", "zTS", "dN", "dT"))[two == TRUE]
P[, S_prom := fifelse(S_is_own, promoter, partner)]   # v4: compare the switching promoter itself (cross-member exact pairs: mu on canonical, p0 on partner)
X <- merge(P, F[, .(sid, b_t, sd_t, b_P, sd_P, se_t, se_P)], by = "sid")
X[, dS := fifelse(S_is_own, b_t / sd_t, b_P / sd_P) * sS]
pr <- fread(file.path(OUTROOT, "A1_location/A1_within_unit_pairs.tsv"), select = c("sid_b", "sid_mu", "gene", "n_prom_gene"))[n_prom_gene == 2]
v <- function(s, col) X[match(s, X$sid)][[col]]
for (cc in c("dS", "dN", "dT", "zS", "zN", "S_is_own", "S_prom", "patternS")) { pr[, paste0(cc, "_b") := v(sid_b, cc)]; pr[, paste0(cc, "_mu") := v(sid_mu, cc)] }
Q <- pr[S_prom_b == S_prom_mu & !is.na(dN_b) & !is.na(dN_mu)]
cat("F2 pairs with continuous effects:", nrow(Q), "\n")
pstat <- function(d, a, b, lab, seed, fun = mean) { d <- d[is.finite(get(a)) & is.finite(get(b))]; dd <- d[, .(gene, x = get(a) - get(b))]
  ug <- unique(dd$gene); set.seed(seed); e <- fun(dd$x)
  bt <- replicate(B, { g <- sample(ug, length(ug), TRUE); fun(dd[data.table(gene = g), on = "gene", allow.cartesian = TRUE]$x) })
  data.table(comparison = lab, stat = deparse(substitute(fun)), n_pairs = nrow(d), boundary = fun(d[[a]]), mu = fun(d[[b]]), diff = e,
             lo = quantile(bt, 0.025), hi = quantile(bt, 0.975), p_boot = min(1, max(2 * min(mean(bt <= 0), mean(bt >= 0)), 1 / B))) }
Q[, `:=`(recip_b = as.numeric(zS_b >= ZC & zN_b <= -ZC), recip_mu = as.numeric(zS_mu >= ZC & zN_mu <= -ZC),
         nonly_b = as.numeric(abs(zS_b) < ZC & zN_b <= -ZC), nonly_mu = as.numeric(abs(zS_mu) < ZC & zN_mu <= -ZC),
         coact_b = as.numeric(zS_b >= ZC & zN_b >= ZC), coact_mu = as.numeric(zS_mu >= ZC & zN_mu >= ZC))]
R1 <- rbind(pstat(Q, "dS_b", "dS_mu", "S effect (SD per allele, S-usage signed)", 1L), pstat(Q, "dN_b", "dN_mu", "N effect", 2L),
            pstat(Q, "dT_b", "dT_mu", "gene output effect", 3L), pstat(Q, "dN_b", "dN_mu", "N effect", 4L, median),
            pstat(Q, "dT_b", "dT_mu", "gene output effect", 5L, median),
            pstat(Q, "recip_b", "recip_mu", "strict reciprocal (S up, N down, both |z| >= 1.96)", 6L),
            pstat(Q, "nonly_b", "nonly_mu", "N only (N down, S not significant)", 7L),
            pstat(Q, "coact_b", "coact_mu", "coactivation (S up, N up)", 8L))
print(R1); writeT(R1, file.path(OUT, "I1_F2_continuous.tsv")); fwrite(Q, file.path(OUT, "I1_F2_pairs.tsv.gz"), sep = "\t")
## all two-promoter signals: N opposite split
A <- X[route %in% c("onoff_S_minor", "mu")]
A[, Ncls := fifelse(zN <= -ZC & zS >= ZC, "reciprocal", fifelse(zN <= -ZC, "N_only_or_S_against", fifelse(zN >= ZC, "N_same", "N_not_sig")))]
splitN <- A[, .N, by = .(route, Ncls)][, share := N / sum(N), by = route]; print(splitN); writeT(splitN, file.path(OUT, "I1_N_opposite_split_all.tsv"))
fwrite(A[, .(sid, gene, route, dS, dN, dT, zS, zN, zTS)], file.path(OUT, "I1_aligned_effects.tsv.gz"), sep = "\t")
## ---------------- I2 linear share and log ratio
sig <- fr(file.path(IN, "delivery/signals.tsv.gz"), select = c("tissue", "promoter", "type", "k", "snp", "chr", "backward_keep"))
sig[, sid := paste(tissue, promoter, type, k, sep = "|")][, unit := paste(tissue, promoter, type, sep = "|")][, bk := as.logical(backward_keep)]
sig[, cond_snps := { s <- snp; b <- bk; vapply(seq_along(s), function(i) paste(setdiff(unique(s[-i][b[-i] %in% TRUE]), s[i]), collapse = ";"), "") }, by = unit]
Z <- merge(X[, .(sid, gene, route, S_is_own, sS)], sig[, .(sid, tissue, promoter, type, snp, chr, cond_snps)], by = "sid")
Z <- merge(Z, F[, .(sid, partner)], by = "sid")
Z[, `:=`(Sprom = fifelse(S_is_own, promoter, partner), Nprom = fifelse(S_is_own, partner, promoter))]
x <- readRDS(file.path(IN, "v11/tabs_att_filter_In.rds")); abs <- x$abs; rm(x); invisible(gc())
abs <- abs[, intersect(colnames(abs), unique(c(Z$Sprom, Z$Nprom))), drop = FALSE]; invisible(gc())
TT <- fr(file.path(IN, "freeze/tissues.tsv"))[def == "freeze"]
gt_dosage <- function(gt) { D <- matrix(NA_integer_, nrow(gt), ncol(gt), dimnames = dimnames(gt))
  D[gt %in% c("0|0", "0/0")] <- 0L; D[gt %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L; D[gt %in% c("1|1", "1/1")] <- 2L; D }
impute_modal <- function(D) { for (i in which(rowSums(is.na(D)) > 0)) { tb <- tabulate(D[i, ] + 1L, nbins = 3L); D[i, is.na(D[i, ])] <- which.max(tb) - 1L }; D }
lmz <- function(X, y) { ok <- is.finite(y); if (sum(ok) < 30 || var(y[ok]) == 0) return(NA_real_); Xo <- X[ok, , drop = FALSE]; yo <- y[ok]; p <- ncol(Xo)
  f <- lm.fit(Xo, yo); if (f$rank < p) return(NA_real_); s2 <- sum(f$residuals^2) / (length(yo) - p); if (s2 <= 0) return(NA_real_)
  Rinv <- backsolve(qr.R(f$qr), diag(p)); j <- which(f$qr$pivot == p); unname(f$coefficients[p]) / sqrt(sum(Rinv[j, ]^2) * s2) }
run_tissue <- function(t) { td <- readRDS(file.path(IN, "freeze", paste0(t, ".rds"))); K <- TT[tissue == t, K_C3]; st <- Z[tissue == t]; if (!nrow(st)) return(NULL)
  C <- cbind(1, scale(cbind(td$pcs, td$peer[, seq_len(K), drop = FALSE])))
  need <- unique(c(st$snp, unlist(strsplit(st$cond_snps[nzchar(st$cond_snps)], ";", fixed = TRUE))))
  ss <- data.table(snp = need)[, `:=`(chr = sub("_.*$", "", snp), pos = as.integer(sub("^chr[^_]+_([0-9]+)_.*$", "\\1", snp)))]; G <- list()
  for (ch in unique(ss$chr)) { s <- ss[chr == ch]
    vv <- readVcf(TabixFile(file.path(IN, "vcf", paste0(ch, ".vcf.gz"))), "hg38", ScanVcfParam(which = GRanges(ch, IRanges(s$pos, s$pos)), samples = td$U, fixed = NA, info = NA, geno = "GT"))
    gt <- geno(vv)$GT[, td$U, drop = FALSE]; gt <- gt[!duplicated(rownames(gt)), , drop = FALSE]; D <- impute_modal(gt_dosage(gt)); for (id in intersect(rownames(D), s$snp)) G[[id]] <- D[id, ] }
  Asub <- as.matrix(abs[td$row_id, intersect(colnames(abs), unique(c(st$Sprom, st$Nprom))), drop = FALSE])
  rbindlist(lapply(seq_len(nrow(st)), function(i) { r <- st[i]; g <- G[[r$snp]]; if (is.null(g) || !all(c(r$Sprom, r$Nprom) %in% colnames(Asub))) return(NULL)
    cs <- if (nzchar(r$cond_snps)) intersect(strsplit(r$cond_snps, ";", fixed = TRUE)[[1]], names(G)) else character(0)
    Xd <- cbind(C, if (length(cs)) do.call(cbind, lapply(cs, function(s) as.numeric(G[[s]]))) else NULL, dose = as.numeric(g)); if (qr(Xd)$rank < ncol(Xd)) Xd <- cbind(C, dose = as.numeric(g))
    nS <- 2^Asub[, r$Sprom] - 1; nN <- 2^Asub[, r$Nprom] - 1; tot <- nS + nN
    data.table(sid = r$sid, z_lin = lmz(Xd, ifelse(tot > 0, nS / tot, NA_real_)) * r$sS, z_lr = lmz(Xd, Asub[, r$Sprom] - Asub[, r$Nprom]) * r$sS,
               z_log2u = lmz(Xd, ifelse(Asub[, r$Sprom] + Asub[, r$Nprom] > 0, Asub[, r$Sprom] / (Asub[, r$Sprom] + Asub[, r$Nprom]), NA_real_)) * r$sS) })) }
RES <- rbindlist(mclapply(sort(unique(Z$tissue)), function(t) tryCatch(run_tissue(t), error = function(e) { cat("ERROR", t, conditionMessage(e), "\n"); NULL }), mc.cores = 2, mc.preschedule = FALSE))
Y <- merge(Z[, .(sid, gene, route)], RES, by = "sid"); fwrite(Y, file.path(OUT, "I2_linear_share_fits.tsv.gz"), sep = "\t")
tabs <- list()
add <- function(d, groups, lab, seed) { groups <- intersect(groups, unique(d$group)); r <- ratio_table(d, groups, B, seed, lab)
  tabs[[length(tabs) + 1]] <<- rbind(r$groups[, .(label, group, n_signals, estimate, lo, hi, contrast = NA_character_, p_boot = NA_real_)],
    if (!is.null(r$contrasts)) r$contrasts[, .(label, group = NA_character_, n_signals = NA_real_, estimate = difference, lo, hi, contrast, p_boot)]) }
rts <- c("onoff_S_minor", "onoff_S_major", "mu")
for (m in c("z_log2u", "z_lin", "z_lr")) {
  add(Y[!is.na(get(m)), .(gene, group = route, num = as.numeric(get(m) >= ZC), den = 1)], rts, paste(m, ">= 1.96 in the S-usage direction"), 1L)
  add(Y[!is.na(get(m)), .(gene, group = route, num = as.numeric(get(m) <= -ZC), den = 1)], rts, paste(m, "<= -1.96 (against the S-usage direction)"), 2L) }
add(Y[!is.na(z_log2u) & z_log2u >= ZC & !is.na(z_lin), .(gene, group = route, num = as.numeric(z_lin >= ZC), den = 1)], rts, "given a log2-usage effect: the linear share moves too", 3L)
TAB <- rbindlist(tabs); writeT(TAB, file.path(OUT, "I2_summary.tsv"))
fmtrow <- function(x) if (is.na(x$contrast)) sprintf("  - %s: %s [%s, %s] (n = %s)", x$group, pct(x$estimate), pct(x$lo), pct(x$hi), format(round(x$n_signals, 1))) else
  sprintf("  - %s: %+.1f pp [%+.1f, %+.1f], P = %s", x$contrast, 100 * x$estimate, 100 * x$lo, 100 * x$hi, formatC(x$p_boot, format = "g", digits = 2))
L <- c("# I -- continuous effects in F2 pairs; N-opposite split; linear-share phenotype", "", "## I1 F2 pairs (boundary vs mu signal of the same promoter)",
       paste0("- ", R1[, sprintf("%s (%s): boundary %.3f vs mu %.3f, diff %+.3f [%+.3f, %+.3f], P = %s (n = %d)", comparison, stat, boundary, mu, diff, lo, hi, formatC(p_boot, format = "g", digits = 2), n_pairs)]),
       "", "## All aligned signals: N response classes", paste0("- ", splitN[order(route, Ncls), sprintf("%s %s: %.1f%% (n = %d)", route, Ncls, 100 * share, N)]), "", "## I2 phenotype scales")
for (lab in unique(TAB$label)) { L <- c(L, paste0("- ", lab)); x <- TAB[label == lab]; for (i in seq_len(nrow(x))) L <- c(L, fmtrow(x[i])) }
writeLines(L, file.path(OUT, "I_report.md")); cat(L, sep = "\n")
