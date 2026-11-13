## Shared LD resource, units v3 (author 2026-10-03: Figures 4, 5 and 6 use ONE LD standard = PLINK 1.9; see code/ld_tool.R).
## Missing genotype calls stay missing and each pair uses the donors called at both SNPs (PLINK's rule), so the proxy
## table equals the PLINK run of work/plink_v3 (checked below). Version 1 (modal-filled genotypes) is archived in
## results/ld_v3/archive_modal_20261003/.
## Output: results/ld_v3/
##   donors/<tissue>.txt, donors/ALL.txt     the donor set of each tissue (td$U = the model's donors) and their union
##   signal_proxies.tsv.gz                  every signal x every SNP the main run tested for that promoter with r2 >= 0.8
##                                          (signal itself included, r = 1): signed r and r2. Same rows as
##                                          fig5_delivery_v3/genotype_pheno_proxies (checked).
##   signal_pairs.tsv.gz                    every pair of signals in the same tissue within 1 Mb (any promoter, any type):
##                                          signed r and r2. Includes all crosstype_ld pairs (checked).
## Usage: NC=16 Rscript code/ld_v3.R
suppressPackageStartupMessages({ library(data.table); library(parallel) })
CODE <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/code"; source(file.path(CODE, "bohr2_data.R")); source(file.path(CODE, "ld_tool.R"))
NC <- as.integer(Sys.getenv("NC", "16")); R2MIN <- 0.8; PAIR_MAX <- 1e6
D3 <- file.path(V2, "results/fig5_delivery_v3"); dir.create(file.path(LD_DIR, "donors"), recursive = TRUE, showWarnings = FALSE)
SF <- fread(cmd = paste("zcat", file.path(D3, "signals.tsv.gz")))[, .(tissue, promoter, type, k, snp, chr, backward_keep)]
TIS <- sort(unique(SF$tissue))
one_tissue <- function(t) {
  td <- readRDS(file.path(V2, "work/final/freeze", paste0(t, ".rds"))); s <- SF[tissue == t]; U <- td$U
  writeLines(U, file.path(LD_DIR, "donors", paste0(t, ".txt")))
  ## (1) proxies within the promoter's tested SNPs (as fig5_proxies.R, plus the sign of r)
  PX <- list()
  for (ch in unique(s$chr)) {
    M <- as.data.table(readRDS(file.path(V2, "results/final/freeze", t, "real", paste0(ch, ".rds")))$snp)[promoter %in% s$promoter, .(promoter, snp)]
    for (p in unique(s[chr == ch]$promoter)) { pos <- td$info[p, "start"]
      G <- tryCatch(ld_raw_region(ch, pos - FLANK, pos + FLANK - 1L, U), error = function(e) NULL); if (is.null(G)) next
      G <- G[, colnames(G) %in% M[promoter == p]$snp, drop = FALSE]   # the tested SNPs (all inside the cis window)
      sp <- s[promoter == p]
      PX[[p]] <- rbindlist(lapply(seq_len(nrow(sp)), function(i) { sn <- sp$snp[i]; if (!sn %in% colnames(G)) return(NULL)
        r <- ld_r(G[, sn], G); px <- names(r)[is.finite(r) & r^2 >= R2MIN]; px <- unique(c(sn, px))
        data.table(tissue = t, promoter = p, type = sp$type[i], k = sp$k[i], signal_snp = sn, proxy_snp = px, r = ifelse(px == sn, 1, unname(r[px]))) })) } }
  ## (2) signal pairs within 1 Mb, same tissue
  u <- unique(s[, .(snp, chr)])[, pos := ld_pos(snp)]; Dm <- ld_raw_dosage(u$snp, U)
  pr <- u[u, on = .(chr), allow.cartesian = TRUE, nomatch = 0L][snp < i.snp & abs(pos - i.pos) <= PAIR_MAX, .(snp_a = snp, snp_b = i.snp, dist = abs(pos - i.pos))]
  pr[, r := as.numeric(mapply(function(a, b) if (a %in% colnames(Dm) && b %in% colnames(Dm)) ld_r(Dm[, a], Dm[, b]) else NA_real_, snp_a, snp_b))]   # as.numeric: 0 pairs
  list(px = rbindlist(PX)[, r2 := r^2], pr = pr[, tissue := t][, r2 := r^2], n_donors = length(U), missing_dosage = setdiff(u$snp, colnames(Dm))) }
RES <- mclapply(TIS, function(t) tryCatch(one_tissue(t), error = function(e) { message(t, ": ", conditionMessage(e)); NULL }), mc.cores = NC, mc.preschedule = FALSE)
names(RES) <- TIS; ok <- !vapply(RES, is.null, TRUE); if (any(!ok)) cat("FAILED tissues:", TIS[!ok], "\n")
writeLines(sort(unique(unlist(lapply(TIS[ok], ld_donors)))), file.path(LD_DIR, "donors", "ALL.txt"))
PX <- rbindlist(lapply(RES[ok], `[[`, "px")); setorder(PX, tissue, promoter, type, k, -r2, proxy_snp)
fwrite(PX, file.path(LD_DIR, "signal_proxies.tsv.gz"), sep = "\t")
## pair table: expand SNP pairs to every (signal_a, signal_b) combination of the same tissue (a SNP can carry several signals)
PR <- rbindlist(lapply(RES[ok], `[[`, "pr"))
SA <- SF[, .(tissue, snp_a = snp, promoter_a = promoter, type_a = type, k_a = k, keep_a = backward_keep)]
SB <- SF[, .(tissue, snp_b = snp, promoter_b = promoter, type_b = type, k_b = k, keep_b = backward_keep)]
PP <- merge(merge(PR, SA, by = c("tissue", "snp_a"), allow.cartesian = TRUE), SB, by = c("tissue", "snp_b"), allow.cartesian = TRUE)
## same-SNP signal pairs (one SNP that is a signal of two units) have r = 1 by definition
SS <- merge(SA, SB, by.x = c("tissue", "snp_a"), by.y = c("tissue", "snp_b"), allow.cartesian = TRUE)[paste(promoter_a, type_a, k_a) < paste(promoter_b, type_b, k_b)][, `:=`(snp_b = snp_a, dist = 0L, r = 1, r2 = 1)]
PP <- rbind(PP, SS, use.names = TRUE)[, same_promoter := promoter_a == promoter_b]
setcolorder(PP, c("tissue", "promoter_a", "type_a", "k_a", "snp_a", "promoter_b", "type_b", "k_b", "snp_b", "dist", "r", "r2", "same_promoter", "keep_a", "keep_b"))
setorder(PP, tissue, promoter_a, type_a, k_a, promoter_b, type_b, k_b); fwrite(PP, file.path(LD_DIR, "signal_pairs.tsv.gz"), sep = "\t")
## check against the PLINK run (results/ld_v3/plink/signal_proxies_plink.tsv.gz, from work/plink_v3 out2)
PL <- fread(cmd = paste("zcat", file.path(LD_DIR, "plink/signal_proxies_plink.tsv.gz")))[, .(tissue, promoter, type, k, signal_snp, proxy_snp, r_plink = r)]
cx <- merge(PX, PL, by = c("tissue", "promoter", "type", "k", "signal_snp", "proxy_snp"), all = TRUE)
cat("proxies: rows", nrow(PX), "; PLINK rows", nrow(PL), "; only here", cx[is.na(r_plink), .N], "; only PLINK", cx[is.na(r), .N],
    "; max |r diff|", signif(max(abs(cx$r - cx$r_plink), na.rm = TRUE), 3), "; sign mismatches", cx[sign(r) != sign(r_plink), .N], "\n")
cat("signal pairs", nrow(PP), "(same promoter", PP[same_promoter == TRUE, .N], "); r2 >= 0.8:", PP[r2 >= 0.8, .N], "; tissues", sum(ok),
    "; SNPs without dosage:", sum(lengths(lapply(RES[ok], `[[`, "missing_dosage"))), "\n")
