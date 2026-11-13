## signals_partner.tsv.gz of the Figure 5 delivery. The v3 table was written ad hoc; this script
## reproduces it (checked identical on delivery v3) so it can be rebuilt for delivery v4 (mirror rule v3).
## One row per signal: snp position, the unit promoter's TSS / strand (v11/Promoter_info_filter_5.txt, TSS = start),
## dist_to_tss = snp_pos - tss; mirror pair columns from the promoter flags (every two-promoter gene; exactness flags).
## Usage: FIG5_OUT=<delivery dir> FLAGS=<promoter_flags.tsv.gz> Rscript code/fig5_signals_partner.R
suppressPackageStartupMessages(library(data.table))
V2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"; OUT <- Sys.getenv("FIG5_OUT")
FLAGS <- Sys.getenv("FLAGS", file.path(V2, "results/mirror/promoter_flags.tsv.gz"))
S <- fread(cmd = paste("zcat", file.path(OUT, "signals.tsv.gz")))[, .(tissue, promoter, type, k, snp)]
info <- fread("/projects/group/bliu/e1327989/GTEx/v11/Promoter_info_filter_5.txt")[, .(promoter = V1, tss = start, strand)]
FL <- fread(cmd = paste("zcat", FLAGS))[, .(tissue, promoter, mirror_pair_id, mirror_role, mirror_partner, boundary_exact, interior_exact, exact_mirror)]
S[, snp_pos := as.integer(sub("^chr[^_]+_([0-9]+)_.*$", "\\1", snp))]
X <- merge(S, info, by = "promoter", all.x = TRUE); X[, dist_to_tss := snp_pos - tss]
X <- merge(X, FL, by = c("tissue", "promoter"), all.x = TRUE)
X <- merge(X, info[, .(mirror_partner = promoter, partner_tss = tss, partner_strand = strand)], by = "mirror_partner", all.x = TRUE)
X[, dist_to_partner_tss := snp_pos - partner_tss]
setcolorder(X, c("tissue", "promoter", "type", "k", "snp", "snp_pos", "tss", "strand", "dist_to_tss", "mirror_pair_id", "mirror_role", "mirror_partner",
                 "partner_tss", "partner_strand", "dist_to_partner_tss", "boundary_exact", "interior_exact", "exact_mirror"))
setorder(X, tissue, promoter, type, k)
fwrite(X, file.path(OUT, "signals_partner.tsv.gz"), sep = "\t"); cat("rows", nrow(X), "\n")
