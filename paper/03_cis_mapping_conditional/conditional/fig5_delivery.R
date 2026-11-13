## Figure 5 delivery from the FINAL conditional analysis (results/conditional_final, built by merge_conditional_ng.R).
## Spec: Figure 5 (2026-10-01 "shape we need"). Freeze definition, real arm, unit level,
## ungated + mirror table v2 de-duplicated + mu estimable (n_mid >= 80).
## Output (results/fig5_delivery/):
##   units.tsv.gz        units with gate flags, mirror flags (merged_component, twin, eq_class, ...), minor_count, signal counts
##   signals.tsv.gz      one row per signal: rank k, SNP, p at selection, backward columns, MARGINAL probability-scale effect
##                       of the unit's own type at that SNP (mu_pe / p0_pe / p1_pe, dosage 0 -> 2) and marginal link-scale est
##   genotype_pheno.tsv.gz   per signal SNP: n0 n1 n2, m0 m1 m2, sd0 sd1 sd2 of the promoter's RAW relative usage
##                       (freeze donors, before the <= 5% endpoint deletion), genotype = alt-allele dosage
##   crosstype_ld.tsv.gz  r2 between every signal of one type and every signal of another type on the same promoter x tissue
##   boundary_only_partner_pheno.tsv.gz  boundary-only exact pairs where mu(partner) and the canonical's p0/p1 unit are
##                       both significant: partner's own raw-usage genotype stats at each canonical boundary-signal SNP,
##                       and r2 to each mu(partner) signal
## Usage: NC=8 Rscript code/fig5_delivery.R [tissue ...]   (default: all tissues)
suppressPackageStartupMessages({ library(data.table); library(parallel) })
CODE <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/code"; source(file.path(CODE, "bohr2_data.R")); source(file.path(CODE, "ld_tool.R"))
FIN <- Sys.getenv("COND_FINAL", file.path(V2, "results/conditional_final")); OUT <- Sys.getenv("FIG5_OUT", file.path(V2, "results/fig5_delivery"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
NC <- as.integer(Sys.getenv("NC", "8")); args <- commandArgs(TRUE)
UF <- fread(cmd = paste("zcat", file.path(FIN, "freeze_real_units.tsv.gz"))); SF <- fread(cmd = paste("zcat", file.path(FIN, "freeze_real_signals.tsv.gz")))
TIS <- if (length(args)) args else sort(unique(UF$tissue))
pos_of <- function(id) as.integer(vapply(strsplit(id, "_", fixed = TRUE), `[`, "", 2))
gstat <- function(y, g) { ok <- !is.na(y) & !is.na(g); y <- y[ok]; g <- g[ok]
  f <- function(d, fun) { v <- y[g == d]; if (length(v)) fun(v) else NA_real_ }
  c(n0 = sum(g == 0), n1 = sum(g == 1), n2 = sum(g == 2), m0 = f(0, mean), m1 = f(1, mean), m2 = f(2, mean),
    sd0 = f(0, sd), sd1 = f(1, sd), sd2 = f(2, sd)) }
one_tissue <- function(t) {
  td <- readRDS(file.path(V2, "work/final/freeze", paste0(t, ".rds")))
  V <- as.data.table(readRDS(file.path(V2, Sys.getenv("UNITS_DIR", "results/final_nogate_v2"), "freeze", paste0(t, ".real.units.rds"))))
  u <- merge(UF[tissue == t], V[, .(promoter, type, gene, n_mid, mid_share, f0, f1, elig_mu, elig_p0, elig_p1, old_gate, estimable, minor_count,
                                   boundary_share = if ("boundary_share" %in% names(V)) boundary_share else NA_real_,
                                   n_prom_gene, mirror_pair_id, boundary_exact, interior_exact, mirror_role, mirror_partner, merged_component, twin, eq_class)],
             by = c("promoter", "type"), all.x = TRUE)
  s <- SF[tissue == t]
  ## marginal effects of each signal SNP from the main run (own type)
  chrs <- unique(s$chr)
  M <- rbindlist(lapply(chrs, function(ch) { x <- readRDS(file.path(V2, "results/final/freeze", t, "real", paste0(ch, ".rds")))$snp
    x <- as.data.table(x)[paste(promoter, snp) %in% paste(s$promoter, s$snp), .(promoter, snp, mu_pe, p0_pe, p1_pe, mu_est, nu_est, tau_est, alt_af, maf)] }))
  s <- merge(s, M, by = c("promoter", "snp"), all.x = TRUE)
  s[, `:=`(p_at_selection = exp(logp_at_selection), pe_marginal = fifelse(type == "mu", mu_pe, fifelse(type == "p0", p0_pe, p1_pe)),
           est_marginal_link = fifelse(type == "mu", mu_est, fifelse(type == "p0", nu_est, tau_est)))]
  ## genotypes: one read per promoter, spanning its signals (and its mirror partner's signals, for the partner table)
  prom_snps <- s[, .(snps = list(unique(snp))), by = promoter]
  G_of <- function(p, extra = character(0)) { ids <- unique(c(prom_snps[promoter == p]$snps[[1]], extra)); ch <- as.character(td$info[p, "seqnames"])
    ps <- pos_of(ids); g <- read_geno(ch, min(ps), max(ps), td$U); if (is.null(g)) return(NULL); g$G[, intersect(ids, g$id), drop = FALSE] }
  ## LD (2026-10-03, author: PLINK as the standard): raw calls, donors missing either SNP dropped per pair (ld_tool.R)
  Graw_of <- function(p, extra = character(0)) { ids <- unique(c(prom_snps[promoter == p]$snps[[1]], extra)); ch <- as.character(td$info[p, "seqnames"])
    ps <- pos_of(ids); D <- ld_raw_region(ch, min(ps), max(ps), td$U); if (is.null(D)) return(NULL); D[, intersect(ids, colnames(D)), drop = FALSE] }
  GP <- list(); LD <- list(); BP <- list()
  for (p in unique(s$promoter)) { G <- tryCatch(G_of(p), error = function(e) NULL); if (is.null(G)) next
    y <- td$act[rownames(G), p]
    sp <- s[promoter == p][order(match(type, c("mu", "p0", "p1")), k)]      # pairs come out as mu-p0, mu-p1, p0-p1
    GP[[p]] <- rbindlist(lapply(seq_len(nrow(sp)), function(i) if (sp$snp[i] %in% colnames(G)) as.data.table(as.list(c(gstat(y, G[, sp$snp[i]]))))[, `:=`(promoter = p, type = sp$type[i], k = sp$k[i], snp = sp$snp[i])]))
    if (uniqueN(sp$type) >= 2) { pr <- CJ(a = seq_len(nrow(sp)), b = seq_len(nrow(sp)))[a < b][sp$type[a] != sp$type[b]]
      if (nrow(pr)) { Gr <- Graw_of(p)
        LD[[p]] <- pr[, .(promoter = p, type_a = sp$type[a], k_a = sp$k[a], snp_a = sp$snp[a], type_b = sp$type[b], k_b = sp$k[b], snp_b = sp$snp[b],
          r = mapply(function(x, z) if (!is.null(Gr) && all(c(x, z) %in% colnames(Gr))) ld_r(Gr[, x], Gr[, z]) else NA_real_, sp$snp[a], sp$snp[b]))][, r2 := r^2] } } }
  ## boundary-only partner table
  bo <- u[mirror_role == "canonical" & boundary_exact == TRUE & interior_exact == FALSE & type %in% c("p0", "p1")]
  for (i in seq_len(nrow(bo))) { A <- bo$promoter[i]; B <- bo$mirror_partner[i]
    if (!nrow(u[promoter == B & type == "mu"])) next
    sa <- s[promoter == A & type == bo$type[i]]; sb <- s[promoter == B & type == "mu"]
    G <- tryCatch(G_of(A, extra = sb$snp), error = function(e) NULL); if (is.null(G)) next
    Gr <- tryCatch(Graw_of(A, extra = sb$snp), error = function(e) NULL)
    yB <- td$act[rownames(G), B]
    BP[[paste(A, bo$type[i])]] <- rbindlist(lapply(seq_len(nrow(sa)), function(j) rbindlist(lapply(seq_len(nrow(sb)), function(m) {
      st <- if (sa$snp[j] %in% colnames(G)) gstat(yB, G[, sa$snp[j]]) else setNames(rep(NA_real_, 9), c("n0","n1","n2","m0","m1","m2","sd0","sd1","sd2"))
      data.table(gene = bo$gene[i], eq_class = bo$eq_class[i], canonical = A, partner = B, canonical_type = bo$type[i],
                 partner_type = c(p0 = "p1", p1 = "p0")[[bo$type[i]]], boundary_k = sa$k[j], boundary_snp = sa$snp[j], t(st),
                 mu_k = sb$k[m], mu_snp = sb$snp[m],
                 r2 = if (!is.null(Gr) && all(c(sa$snp[j], sb$snp[m]) %in% colnames(Gr))) ld_r(Gr[, sa$snp[j]], Gr[, sb$snp[m]])^2 else NA_real_) })))) }
  list(units = u, signals = s, gp = rbindlist(GP, fill = TRUE)[, tissue := t], ld = rbindlist(LD, fill = TRUE), bp = rbindlist(BP, fill = TRUE))
}
RES <- mclapply(TIS, function(t) tryCatch(one_tissue(t), error = function(e) { message(t, ": ", conditionMessage(e)); NULL }), mc.cores = NC, mc.preschedule = FALSE)
names(RES) <- TIS; ok <- !vapply(RES, is.null, TRUE); if (any(!ok)) cat("FAILED tissues:", TIS[!ok], "\n")
pick <- function(k) rbindlist(lapply(names(RES)[ok], function(t) { d <- RES[[t]][[k]]; if (NROW(d) && !"tissue" %in% names(d)) d[, tissue := t]; d }), fill = TRUE)
U <- pick("units"); S <- pick("signals"); GP <- pick("gp"); LD <- pick("ld"); BP <- pick("bp")
fwrite(U, file.path(OUT, "units.tsv.gz"), sep = "\t"); fwrite(S, file.path(OUT, "signals.tsv.gz"), sep = "\t")
fwrite(GP, file.path(OUT, "genotype_pheno.tsv.gz"), sep = "\t"); fwrite(LD, file.path(OUT, "crosstype_ld.tsv.gz"), sep = "\t")
fwrite(BP, file.path(OUT, "boundary_only_partner_pheno.tsv.gz"), sep = "\t")
cat("tissues", sum(ok), "; units", nrow(U), "; signals", nrow(S), "; genotype rows", nrow(GP), "; cross-type pairs", nrow(LD), "; boundary-only partner rows", nrow(BP), "\n")
