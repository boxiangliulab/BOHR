## 132_marks_by_usage_role.R -- author 2026-11-04: show the three ENCODE signals separately (DNase, H3K27ac, H3K4me3) instead of the
## mutually exclusive states, which hide H3K27ac (a TSS that is H3K27ac+ and H3K4me3+ is called promoter-like, so "enhancer-like" =
## H3K27ac+ without H3K4me3 is rare by construction). Same element, samples and threshold as 122a / 129: the ENCODE element covering
## the TSS, the QTL tissue's matched biosamples (exact + approximate), mark called when z > 1.64. Per (tissue, promoter) unit: the
## fraction of its matched biosamples in which each mark is called; then the mean over units with a gene-clustered bootstrap CI.
## Role = usage-defined major / minor promoter of two-promoter genes (mirror_pairs tissue means), as in 129.
F5S <- Sys.getenv("F5S", "/home/e1327989/f5_stage_20261106_v6")
src <- readLines(file.path(F5S, "code/122a_tissue_matched_chromatin_v4.R"))
eval(parse(text = src[1:(grep("^summ <- function", src) - 1)]))
mp <- MAPb[match %in% c("exact", "approximate")]
MK <- rbindlist(lapply(split(TP, TP$tissue), function(x) { bs <- mp[gtex == x$tissue[1]]$biosample_value
  if (!length(bs)) return(NULL); has <- !is.na(x$acc) & x$acc %in% rownames(A); if (!any(has)) return(NULL); x <- x[has]
  x[, `:=`(n_bios = length(bs), DNase = rowMeans(D[acc, bs, drop = FALSE] > TH), H3K27ac = rowMeans(K[acc, bs, drop = FALSE] > TH),
           H3K4me3 = rowMeans(A[acc, bs, drop = FALSE] > TH))] }))
MK <- grp(MK)
MP <- fread(cmd = sprintf("zcat %s", file.path(F5S, "in/mirror/mirror_pairs.tsv.gz")), select = c("tissue", "canonical", "partner"))   # major = canonical = higher usage pooled over all tissues (author 2026-11-04: one definition, as Figure 4)
RO <- rbind(MP[, .(tissue, promoter = canonical, role = "major")], MP[, .(tissue, promoter = partner, role = "minor")])
MK <- merge(MK, RO, by = c("tissue", "promoter"))
res <- rbindlist(lapply(c("major", "minor"), function(r) rbindlist(lapply(GL, function(g) { d <- MK[role == r & group == g]
  rbindlist(lapply(c("DNase", "H3K27ac", "H3K4me3"), function(m) { f <- function(x) mean(x[[m]]); b <- gboot(d, f)
    data.table(role = r, group = g, mark = m, units = nrow(d), genes = uniqueN(d$gene0), mean_frac = f(d), lo = quantile(b, .025), hi = quantile(b, .975)) })) }))))
fwrite(res, file.path(OUTD, "132_marks_by_usage_role.tsv"), sep = "\t")
## H3K27ac among open TSSs, split by H3K4me3 (shows where H3K27ac goes)
options(width = 200); print(dcast(res[, v := sprintf("%.1f%% [%.1f, %.1f]", 100 * mean_frac, 100 * lo, 100 * hi)], role + group + units ~ mark, value.var = "v"))
