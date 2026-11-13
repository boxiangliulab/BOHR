## crosstype_conditional_v6.R (2026-11-06): copy of crosstype_conditional_v5.R for Figure 5 delivery v6 (conditional refit v3c). Only the delivery path differs.
## crosstype_conditional_v5.R (2026-11-06): copy of crosstype_conditional_v4.R for Figure 5 delivery v5 (abs-activity mirror labels, units v3a). Only the delivery path differs.
## crosstype_conditional_v4.R (2026-11-04): copy of crosstype_conditional.R writing Figure 5 delivery v4 (mirror rule v3: exact switches = p0 on the silent TSS;
## in boundary_only_partner_pheno, columns canonical / partner mean the promoter holding the boundary unit / the other promoter).
## CONDITIONAL cross-type tests for pairs of mu and boundary signals of the same promoter, units v3.
## Pair set = Figure 5's 738 pairs:
##   (i)  within-promoter: crosstype_ld mu-p0 / mu-p1 signal pairs with r2 < 0.2 (718 pairs, 574 promoter x tissue);
##   (ii) rule-2 (boundary-only exact pairs): canonical boundary signal x partner mu signal with r2 < 0.2
##        (boundary_only_partner_pheno; 20 pairs, 17 units).
## For every SNP of a pair (the mu signal s_mu and the boundary signal s_b), ONE conditioning set is used for its own-type
## and cross-type tests: the pair partner plus all OTHER retained signals (backward_keep == TRUE, all three components) of
## the tested promoter. The partner is itself retained in 732 of 738 pairs; in the other 6 one SNP of the pair failed
## backward, and it is added explicitly (s_*_backward_keep flags), so that SNP has a pair-specific set.
##   (i)  tested promoter = the unit's promoter;
##   (ii) tested promoter = the PARTNER B. Its boundary indicators are the canonical A's with the labels swapped
##        (p0(B) = p1(A), same indicator, same estimate). The conditioning set is B's retained signals plus A's retained
##        boundary signals. b_type keeps the canonical label; b_type_tested is the label on B.
## A SNP that is a signal of two components of the same promoter has one conditioning set and is fitted once
## (role = "boundary_signal;mu_signal").
## Tests (plan B: conditioning SNPs appended to the covariates of both likelihood parts, as in the conditional analysis):
##   mu = F-corrected LRT (F = (n_mid - p) expm1(W / n_mid), p = interior rank + 1) + the chi2 version
##   p0 = nu Firth PLR, p1 = tau Firth PLR (bohr2_addon3.R firth_snp_tests, multi-start)
## Testability: p0 needs zeros, p1 needs ones (family after the <= 5% endpoint rule); mu needs n_mid > p + 1. A part is also
## not testable when the tested SNP lies in the span of the conditioning design (all donors for p0/p1, interior donors for
## mu): its effect is not identifiable given the conditioning SNPs. NA means "not testable" (flag *_testable = FALSE, text
## in *_na_reason, "none" when tested); a fitted-but-unconverged test is flagged *_conv = FALSE, not set to NA.
## *_ml_sep = separation flag of the unpenalised fit (Firth estimates stay finite). interior_rank_deficient = the
## conditioning design is rank-deficient on the interior donors (the mu LRT stays valid; mu_se can be NA).
## Diagnostics: max_r2_cond = the tested SNP's highest r2 with any conditioning SNP (the shared LD standard of
## code/ld_tool.R: the tissue's freeze donors, modal-imputed dosages); cond_ld08_snps /
## cond_ld08_types = conditioning SNPs with r2 >= 0.8 to the tested SNP and the components they are signals of
## (tested promoter's labels; "_failed" = failed backward). Conditional tests there are underpowered by construction.
## Probability-scale effects (dosage 0 -> 2) with delta-method SEs from the same fits: *_pe, *_pe_se (code/pe_delta.R).
## Output: results/fig5_delivery_v6/crosstype_conditional_{snps,pairs}.tsv.gz
suppressPackageStartupMessages({ library(data.table); library(parallel) })
CODE <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/code"
source(file.path(CODE, "bohr2_data.R")); source(file.path(CODE, "bohr2_engine.R")); source(file.path(CODE, "bohr2_addon3.R")); source(file.path(CODE, "pe_delta.R"))
source(file.path(CODE, "ld_tool.R"))
OUT <- file.path(V2, "results/fig5_delivery_v6"); NC <- as.integer(Sys.getenv("NC", "8"))
rd <- function(f) fread(cmd = paste("zcat", file.path(OUT, f)))
S <- rd("signals.tsv.gz"); LD <- rd("crosstype_ld.tsv.gz"); BP <- rd("boundary_only_partner_pheno.tsv.gz")
KEEP <- S[backward_keep == TRUE, .(tissue, promoter, type, k, snp)]
SW <- c(mu = "mu", p0 = "p1", p1 = "p0")
## pair set
P1 <- LD[type_a == "mu" & type_b %in% c("p0", "p1") & r2 < 0.2, .(pair_kind = "within", tissue, tested_promoter = promoter, mu_promoter = promoter, b_promoter = promoter,
        b_type = type_b, mu_k = k_a, s_mu = snp_a, b_k = k_b, s_b = snp_b, r2)]
P2 <- unique(BP[r2 < 0.2, .(pair_kind = "rule2", tissue, tested_promoter = partner, mu_promoter = partner, b_promoter = canonical,
        b_type = canonical_type, mu_k, s_mu = mu_snp, b_k = boundary_k, s_b = boundary_snp, r2)])
PR <- rbind(P1, P2)[, b_type_tested := fifelse(pair_kind == "rule2", unname(SW[b_type]), b_type)]
## stable pair IDs: pairs of the first delivery keep their f-number (work/crosstype_pair_ids_v1.tsv, 738 pairs);
## pairs that enter later (e.g. after the PLINK LD standard) continue the series
IDM <- fread(file.path(V2, "work/crosstype_pair_ids_v1.tsv"))
PR <- merge(PR, IDM, by = c("pair_kind", "tissue", "tested_promoter", "s_mu", "s_b"), all.x = TRUE, sort = FALSE)
nmax <- max(as.integer(sub("f", "", IDM$pair_id))); nn <- PR[is.na(pair_id), .N]
if (nn) PR[is.na(pair_id), pair_id := sprintf("f%04d", nmax + seq_len(nn))]
cat("pair IDs: kept", PR[pair_id %in% IDM$pair_id, .N], "; new", nn, "; first-delivery pairs no longer in the set:", paste(setdiff(IDM$pair_id, PR$pair_id), collapse = ","), "\n")
cat("pairs:", nrow(PR), "(within", nrow(P1), ", rule2", nrow(P2), "); both k = 1:", PR[mu_k == 1 & b_k == 1, .N], "\n")
## conditioning set of a SNP on its tested promoter, and the components (tested promoter's labels) it is itself a signal of
cond_of <- function(t, p, pk, bprom) { x <- KEEP[tissue == t & promoter == p]$snp
  if (pk == "rule2") x <- c(x, KEEP[tissue == t & promoter == bprom & type %in% c("p0", "p1")]$snp); unique(x) }
csid <- function(t, p, pk, bprom, s, partner) paste(sort(setdiff(union(cond_of(t, p, pk, bprom), partner), s)), collapse = ";")
types_of <- function(t, p, pk, bprom, s) { x <- KEEP[tissue == t & promoter == p & snp == s]$type
  if (pk == "rule2") x <- c(x, SW[KEEP[tissue == t & promoter == bprom & snp == s & type %in% c("p0", "p1")]$type]); paste(sort(unique(x)), collapse = ";") }
PR[, `:=`(cs_mu = mapply(csid, tissue, tested_promoter, pair_kind, b_promoter, s_mu, s_b, USE.NAMES = FALSE),
          cs_b = mapply(csid, tissue, tested_promoter, pair_kind, b_promoter, s_b, s_mu, USE.NAMES = FALSE),
          s_mu_signal_types = mapply(types_of, tissue, tested_promoter, pair_kind, b_promoter, s_mu, USE.NAMES = FALSE),
          s_b_signal_types = mapply(types_of, tissue, tested_promoter, pair_kind, b_promoter, s_b, USE.NAMES = FALSE))]
BK <- S[, .(tissue, promoter, type, k, snp, backward_keep)]
PR <- merge(PR, BK[type == "mu", .(tissue, mu_promoter = promoter, mu_k = k, s_mu = snp, s_mu_backward_keep = backward_keep)], by = c("tissue", "mu_promoter", "mu_k", "s_mu"), all.x = TRUE)
PR <- merge(PR, BK[, .(tissue, b_promoter = promoter, b_type = type, b_k = k, s_b = snp, s_b_backward_keep = backward_keep)], by = c("tissue", "b_promoter", "b_type", "b_k", "s_b"), all.x = TRUE)
cat("pairs with a SNP that failed backward:", PR[!(s_mu_backward_keep & s_b_backward_keep), .N], "; backward flag missing:", PR[is.na(s_mu_backward_keep) | is.na(s_b_backward_keep), .N], "\n")
in_set <- function(s, id) s %in% strsplit(id, ";", fixed = TRUE)[[1]]
cat("pairs whose partner SNP is in the conditioning set (both sides):", sum(mapply(in_set, PR$s_b, PR$cs_mu) & mapply(in_set, PR$s_mu, PR$cs_b)), "of", nrow(PR), "\n")
## SNP tasks: one per (tissue, tested promoter, snp, conditioning set)
TK <- rbind(PR[, .(tissue, tested_promoter, snp = s_mu, cond_set_id = cs_mu, role = "mu_signal", pair_kind, b_promoter, st = s_mu_signal_types)],
            PR[, .(tissue, tested_promoter, snp = s_b, cond_set_id = cs_b, role = "boundary_signal", pair_kind, b_promoter, st = s_b_signal_types)])
TK <- TK[, .(role = paste(sort(unique(role)), collapse = ";"), pair_kind = paste(sort(unique(pair_kind)), collapse = ";"), b_promoter = b_promoter[1],
             snp_signal_types = paste(sort(unique(unlist(strsplit(st, ";", fixed = TRUE)))), collapse = ";")), by = .(tissue, tested_promoter, snp, cond_set_id)]
cat("SNP tasks:", nrow(TK), "; SNPs with two roles:", TK[role == "boundary_signal;mu_signal", .N], "\n")
pos_of <- function(id) as.integer(vapply(strsplit(id, "_", fixed = TRUE), `[`, "", 2))
## component labels of every forward signal (tested promoter's labels; rule 2 adds the canonical's boundary signals, swapped)
LAB <- S[, .(tissue, promoter, snp, type, fail = fifelse(backward_keep, "", "_failed"))]
labels_of <- function(t, p, pk, bprom) { x <- LAB[tissue == t & promoter == p, .(snp, lab = paste0(type, fail))]
  if (grepl("rule2", pk)) x <- rbind(x, LAB[tissue == t & promoter == bprom & type %in% c("p0", "p1"), .(snp, lab = paste0(SW[type], fail))]); x }
one_tissue <- function(t) {
  td <- readRDS(file.path(V2, "work/final/freeze", paste0(t, ".rds"))); C_all <- cov_matrix(td, td$peer, td$K); x <- TK[tissue == t]; out <- list()
  for (p in unique(x$tested_promoter)) {
    z <- prepare_y(td$act[, p], ncol(C_all) - 1L); if (!is.na(z$skip)) { message(t, " ", p, " skip ", z$skip); next }
    smp <- names(z$y); y <- unname(z$y); ch <- as.character(td$info[p, "seqnames"]); fam <- z$family; m <- y > 0 & y < 1; n_mid <- sum(m)
    xp <- x[tested_promoter == p]
    for (i in seq_len(nrow(xp))) { sn <- xp$snp[i]
      cs <- if (nzchar(xp$cond_set_id[i])) strsplit(xp$cond_set_id[i], ";", fixed = TRUE)[[1]] else character(0); ids <- c(cs, sn); ps <- pos_of(ids)
      g <- read_geno(ch, min(ps), max(ps), td$U); ok <- !is.null(g) && all(ids %in% colnames(g$G))
      row <- data.table(xp[i], family = fam, n = length(y), n_mid = n_mid, n_cond = length(cs), geno_ok = ok)
      if (ok) { G <- g$G[smp, ids, drop = FALSE]
        gr <- tryCatch(ld_raw_region(ch, min(ps), max(ps), td$U), error = function(e) NULL)   # LD = PLINK rule (ld_tool.R): raw calls, all freeze donors
        r2c <- if (length(cs) && !is.null(gr) && all(ids %in% colnames(gr))) ld_r(gr[, sn], gr[, cs, drop = FALSE])^2 else rep(NA_real_, length(cs))
        Cc <- orth_design(cbind(C_all[smp, , drop = FALSE], G[, cs, drop = FALSE])); dropped <- length(attr(Cc, "dropped"))
        nul <- promoter_null(y, Cc, fam); rk_m <- qr(Cc[m, , drop = FALSE])$rank; pp <- rk_m + 1L
        mu_id <- qr(cbind(Cc[m, , drop = FALSE], G[m, sn]))$rank > rk_m; b_id <- qr(cbind(Cc, G[, sn]))$rank > qr(Cc)$rank
        r <- tryCatch(snp_tests(y, Cc, G[, sn], nul), error = function(e) NULL)
        f <- tryCatch(firth_snp_tests_ff(y, Cc, G[, sn], nul), error = function(e) NULL)
        mp <- tryCatch(pe_mu_fit(y, Cc, G[, sn], nul), error = function(e) NULL)
        bp <- if (!is.null(f)) tryCatch(pe_firth(f$ff, Cc, fam), error = function(e) NULL) else NULL
        mu_n <- n_mid > pp + 1L; p0_f <- fam %in% c("BEZI", "BEINF"); p1_f <- fam %in% c("BEOI", "BEINF")
        mu_ok <- mu_n && mu_id; p0_ok <- p0_f && b_id; p1_ok <- p1_f && b_id
        W <- if (!is.null(r)) max(unname(r$mu_stat), 0) else NA_real_
        hi <- cs[which(r2c >= 0.8)]; lb <- labels_of(t, p, xp$pair_kind[i], xp$b_promoter[i])
        row[, `:=`(cond_cols_dropped = dropped, p_interior = pp, interior_rank_deficient = rk_m < ncol(Cc),
          max_r2_cond = if (any(is.finite(r2c))) max(r2c, na.rm = TRUE) else NA_real_,
          max_r2_cond_snp = if (any(is.finite(r2c))) cs[which.max(r2c)] else NA_character_,
          cond_ld08_snps = paste(hi, collapse = ";"), cond_ld08_types = paste(sort(unique(lb[snp %in% hi]$lab)), collapse = ";"),
          mu_testable = mu_ok, p0_testable = p0_ok, p1_testable = p1_ok,
          mu_na_reason = if (!mu_n) "too few interior donors (n_mid <= p_interior + 1)" else if (!mu_id) "not identifiable: tested SNP in the span of the conditioning design (interior donors)" else "none",
          p0_na_reason = if (!p0_f) paste0("no zeros after the 5% endpoint rule (family ", fam, ")") else if (!b_id) "not identifiable: tested SNP in the span of the conditioning design" else "none",
          p1_na_reason = if (!p1_f) paste0("no ones after the 5% endpoint rule (family ", fam, ")") else if (!b_id) "not identifiable: tested SNP in the span of the conditioning design" else "none",
          mu_logp = if (mu_ok && is.finite(W)) pf((n_mid - pp) * expm1(W / n_mid), 1, n_mid - pp, lower.tail = FALSE, log.p = TRUE) else NA_real_,
          mu_logp_chi = if (mu_ok && !is.null(r)) unname(r$mu_lrt_logp) else NA_real_,
          mu_est = if (mu_ok && !is.null(r)) unname(r$mu_est) else NA_real_, mu_se = if (mu_ok && !is.null(r)) unname(r$mu_se) else NA_real_,
          mu_conv = if (mu_ok && !is.null(r)) isTRUE(r$mu_conv) else NA, mu_ml_sep = if (mu_ok && !is.null(r)) isTRUE(r$mu_sep) else NA,
          mu_pe = if (mu_ok && !is.null(mp)) mp$pe else NA_real_, mu_pe_se = if (mu_ok && !is.null(mp)) mp$pe_se else NA_real_,
          p0_pe = if (p0_ok && !is.null(bp)) bp$p0_pe else NA_real_, p0_pe_se = if (p0_ok && !is.null(bp)) bp$p0_pe_se else NA_real_,
          p1_pe = if (p1_ok && !is.null(bp)) bp$p1_pe else NA_real_, p1_pe_se = if (p1_ok && !is.null(bp)) bp$p1_pe_se else NA_real_,
          mu_est_check = if (mu_ok && !is.null(mp)) mp$est - unname(r$mu_est) else NA_real_,
          p0_logp = if (p0_ok && !is.null(f)) unname(f$nu_firth_logp) else NA_real_, p0_est = if (p0_ok && !is.null(f)) unname(f$nu_firth_est) else NA_real_,
          p0_se = if (p0_ok && !is.null(f)) unname(f$nu_firth_se) else NA_real_, p0_conv = if (p0_ok && !is.null(f)) isTRUE(f$nu_firth_conv) else NA,
          p0_ml_sep = if (p0_ok && !is.null(r)) isTRUE(r$nu_sep) else NA,
          p1_logp = if (p1_ok && !is.null(f)) unname(f$tau_firth_logp) else NA_real_, p1_est = if (p1_ok && !is.null(f)) unname(f$tau_firth_est) else NA_real_,
          p1_se = if (p1_ok && !is.null(f)) unname(f$tau_firth_se) else NA_real_, p1_conv = if (p1_ok && !is.null(f)) isTRUE(f$tau_firth_conv) else NA,
          p1_ml_sep = if (p1_ok && !is.null(r)) isTRUE(r$tau_sep) else NA,
          cat_complete_sep = if (!is.null(r) && !is.null(r$cat_complete_sep)) isTRUE(r$cat_complete_sep) else NA,
          fit_error = is.null(r) || is.null(f) || is.null(mp) || is.null(bp))] }
      out[[length(out) + 1]] <- row } }
  rbindlist(out, fill = TRUE) }
X <- rbindlist(mclapply(unique(TK$tissue), function(t) tryCatch(one_tissue(t), error = function(e) { message(t, ": ", conditionMessage(e)); NULL }), mc.cores = NC, mc.preschedule = FALSE), fill = TRUE)
for (v in c("mu", "p0", "p1")) { rs <- paste0(v, "_na_reason"); lp <- X[[paste0(v, "_logp")]]
  X[geno_ok %in% FALSE, (rs) := "genotypes missing"]; X[fit_error %in% TRUE, (rs) := "fit error"]
  X[get(rs) == "none" & is.na(lp), (rs) := "test returned NA"] }
setorder(X, tissue, tested_promoter, snp); fwrite(X, file.path(OUT, "crosstype_conditional_snps.tsv.gz"), sep = "\t")
## pair view: both SNPs, joined on their own conditioning set; headline own/cross columns on the tested promoter's labels
VAL <- c("role", "snp_signal_types", "n_cond", "max_r2_cond", "max_r2_cond_snp", "cond_ld08_snps", "cond_ld08_types", "cond_cols_dropped", "interior_rank_deficient",
         "mu_na_reason", "p0_na_reason", "p1_na_reason",
         "mu_logp", "mu_logp_chi", "mu_est", "mu_se", "mu_pe", "mu_pe_se", "mu_conv", "mu_testable",
         "p0_logp", "p0_est", "p0_se", "p0_pe", "p0_pe_se", "p0_conv", "p0_ml_sep", "p0_testable",
         "p1_logp", "p1_est", "p1_se", "p1_pe", "p1_pe_se", "p1_conv", "p1_ml_sep", "p1_testable", "cat_complete_sep")
side <- function(s, cs, pre) setnames(X[, c("tissue", "tested_promoter", "snp", "cond_set_id", VAL), with = FALSE], c("snp", "cond_set_id", VAL), c(s, cs, paste0(pre, VAL)))
PP <- merge(PR, side("s_mu", "cs_mu", "at_mu_"), by = c("tissue", "tested_promoter", "s_mu", "cs_mu"), all.x = TRUE)
PP <- merge(PP, side("s_b", "cs_b", "at_b_"), by = c("tissue", "tested_promoter", "s_b", "cs_b"), all.x = TRUE)
pick <- function(b, a0, a1) fifelse(b == "p0", a0, a1)
PP[, `:=`(own_mu_logp = at_mu_mu_logp, own_mu_est = at_mu_mu_est, own_mu_se = at_mu_mu_se, own_mu_pe = at_mu_mu_pe, own_mu_pe_se = at_mu_mu_pe_se, own_mu_conv = at_mu_mu_conv,
          cross_b_at_mu_logp = pick(b_type_tested, at_mu_p0_logp, at_mu_p1_logp), cross_b_at_mu_est = pick(b_type_tested, at_mu_p0_est, at_mu_p1_est),
          cross_b_at_mu_se = pick(b_type_tested, at_mu_p0_se, at_mu_p1_se), cross_b_at_mu_conv = pick(b_type_tested, at_mu_p0_conv, at_mu_p1_conv),
          cross_b_at_mu_pe = pick(b_type_tested, at_mu_p0_pe, at_mu_p1_pe), cross_b_at_mu_pe_se = pick(b_type_tested, at_mu_p0_pe_se, at_mu_p1_pe_se),
          own_b_logp = pick(b_type_tested, at_b_p0_logp, at_b_p1_logp), own_b_est = pick(b_type_tested, at_b_p0_est, at_b_p1_est),
          own_b_se = pick(b_type_tested, at_b_p0_se, at_b_p1_se), own_b_conv = pick(b_type_tested, at_b_p0_conv, at_b_p1_conv),
          own_b_pe = pick(b_type_tested, at_b_p0_pe, at_b_p1_pe), own_b_pe_se = pick(b_type_tested, at_b_p0_pe_se, at_b_p1_pe_se),
          cross_mu_at_b_logp = at_b_mu_logp, cross_mu_at_b_est = at_b_mu_est, cross_mu_at_b_se = at_b_mu_se, cross_mu_at_b_conv = at_b_mu_conv,
          cross_mu_at_b_pe = at_b_mu_pe, cross_mu_at_b_pe_se = at_b_mu_pe_se,
          s_mu_is_b_signal = mapply(function(b, s) b %in% strsplit(s, ";", fixed = TRUE)[[1]], b_type_tested, s_mu_signal_types, USE.NAMES = FALSE),
          s_b_is_mu_signal = mapply(function(s) "mu" %in% strsplit(s, ";", fixed = TRUE)[[1]], s_b_signal_types, USE.NAMES = FALSE))]
## locus flags: the SNP itself, or a retained signal in r2 >= 0.8 with it, is a signal of the other arm's component. Then the
## conditional cross test is absorbed by that signal: a large p there means "not separable", not "no cross effect".
has <- function(lbl, types) mapply(function(l, s) l %in% strsplit(s, ";", fixed = TRUE)[[1]], lbl, fifelse(is.na(types), "", types), USE.NAMES = FALSE)
PP[, `:=`(s_mu_locus_b_signal = s_mu_is_b_signal | has(b_type_tested, at_mu_cond_ld08_types),
          s_b_locus_mu_signal = s_b_is_mu_signal | has(rep("mu", .N), at_b_cond_ld08_types))]
setcolorder(PP, c("pair_id", "pair_kind", "tissue", "tested_promoter", "mu_promoter", "b_promoter", "b_type", "b_type_tested", "mu_k", "s_mu", "b_k", "s_b", "r2",
                  grep("^own_|^cross_", names(PP), value = TRUE), "s_mu_backward_keep", "s_b_backward_keep", "s_mu_locus_b_signal", "s_b_locus_mu_signal",
                  "s_mu_is_b_signal", "s_b_is_mu_signal", "s_mu_signal_types", "s_b_signal_types", "at_mu_max_r2_cond", "at_b_max_r2_cond", "at_mu_cond_ld08_types", "at_b_cond_ld08_types"))
setorder(PP, pair_id); fwrite(PP, file.path(OUT, "crosstype_conditional_pairs.tsv.gz"), sep = "\t")
hd <- c("own_mu_logp", "cross_b_at_mu_logp", "own_b_logp", "cross_mu_at_b_logp")
cat("SNP rows", nrow(X), "; geno ok", X[geno_ok == TRUE, .N], "; fit errors", X[fit_error == TRUE, .N], "; unconverged mu/p0/p1:", X[mu_conv == FALSE, .N], X[p0_conv == FALSE, .N], X[p1_conv == FALSE, .N],
    "; pairs", nrow(PP), "(unique", uniqueN(PP$pair_id), "); pairs with an NA headline test", PP[!complete.cases(PP[, ..hd]), .N], "\n")
cat("pe check: max |mu est (pe refit) - mu est (test)|", suppressWarnings(max(abs(X$mu_est_check), na.rm = TRUE)), "; pe present but SE NA: mu", X[!is.na(mu_pe) & is.na(mu_pe_se), .N],
    "p0", X[!is.na(p0_pe) & is.na(p0_pe_se), .N], "p1", X[!is.na(p1_pe) & is.na(p1_pe_se), .N], "; logp present but pe NA:", X[(!is.na(mu_logp) & is.na(mu_pe)) | (!is.na(p0_logp) & is.na(p0_pe)) | (!is.na(p1_logp) & is.na(p1_pe)), .N], "\n")
