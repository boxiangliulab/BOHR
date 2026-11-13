## Figure 4 panel e, LD-window coloc, step 2: assemble the LD window of every lead (from 01_ldwin_nscc.pbs) and the
## pair-level window table (puQTL signal x its gene's eQTL lead). Window = [min(lead, partners) - 50 kb,
## max(lead, partners) + 50 kb]; a lead with no partner at r2 >= 0.1 gets lead +/- 50 kb. Leads absent from the v9
## reference (not in the MAF >= 0.05 panel) are listed; they have no LD window.
suppressPackageStartupMessages(library(data.table))
V2  <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"; OUT <- file.path(V2, "results/coloc_eqtl_ldwin")
PAD <- 50000L
pos_of <- function(s) as.integer(sub("^chr[^_]+_([0-9]+)_.*$", "\\1", s))
L <- fread(file.path(OUT, "leads.tsv"))
rd <- function(pat, f) rbindlist(lapply(list.files(file.path(OUT, "ldwin_out"), pat, full.names = TRUE), function(x) {
  d <- f(x); if (nrow(d)) d[, tissue := sub("\\.chr[0-9]+\\..*$", "", basename(x))]; d }), fill = TRUE)
W  <- rd("\\.win\\.tsv$", function(x) fread(x, colClasses = list(character = "snp")))
AB <- rd("\\.absent\\.txt$", function(x) { r <- readLines(x); data.table(snp = r) })
L <- merge(L, W[, .(tissue, snp, n_partners, min_partner_pos, max_partner_pos)], by = c("tissue", "snp"), all.x = TRUE)
L[, absent := paste(tissue, snp) %in% AB[, paste(tissue, snp)]]
L[, pos := pos_of(snp)]
L[absent == FALSE & is.na(n_partners), `:=`(n_partners = 0L, min_partner_pos = pos, max_partner_pos = pos)]
L[absent == FALSE, `:=`(win_start = pmin(pos, min_partner_pos) - PAD, win_end = pmax(pos, max_partner_pos) + PAD)]
fwrite(L, file.path(OUT, "lead_windows.tsv"), sep = "\t")
cat("leads:", nrow(L), "; absent from the v9 panel:", L[absent == TRUE, .N], "\n"); print(L[, .N, by = .(role, absent)])
cat("window width (kb), present leads:\n"); print(L[absent == FALSE, quantile((win_end - win_start) / 1000, c(0, .1, .25, .5, .75, .9, .99, 1))])
## pairs: each analysis signal x its gene's eQTL lead
A <- readRDS(file.path(V2, "results/coloc_eqtl_v3/analysis_set.rds"))[, .(tissue, gene_base, promoter, type, k, signal_snp, egene)]
G <- fread(file.path(OUT, "eqtl_leads.tsv"))[, .(tissue, gene_base, eqtl_lead)]
P <- merge(A, G, by = c("tissue", "gene_base"))
P <- merge(P, L[, .(tissue, signal_snp = snp, p_absent = absent, p_start = win_start, p_end = win_end)], by = c("tissue", "signal_snp"))
P <- merge(P, L[, .(tissue, eqtl_lead = snp, e_absent = absent, e_start = win_start, e_end = win_end)], by = c("tissue", "eqtl_lead"))
P[, `:=`(int_start = pmax(p_start, e_start), int_end = pmin(p_end, e_end), uni_start = pmin(p_start, e_start), uni_end = pmax(p_end, e_end))]
P[, overlap_bp := pmax(0L, int_end - int_start)]
fwrite(P, file.path(OUT, "pair_windows.tsv.gz"), sep = "\t")
cat("\npairs:", nrow(P), "; eQTL lead absent:", P[e_absent == TRUE, .N], "; positional overlap > 0:", P[e_absent == FALSE & overlap_bp > 0, .N],
    "; no overlap:", P[e_absent == FALSE & overlap_bp == 0, .N], "\n")
print(P[egene == TRUE, .(n = .N, e_absent = sum(e_absent), no_overlap = sum(!e_absent & overlap_bp == 0))])
cat("union width (kb):\n"); print(P[e_absent == FALSE, quantile((uni_end - uni_start) / 1000, c(.1, .5, .9, .99, 1))])
