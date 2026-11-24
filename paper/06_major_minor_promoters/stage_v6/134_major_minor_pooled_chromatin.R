## 134_major_minor_pooled_chromatin.R -- the main panel compares only the major and the minor promoter (no QTL-class split;
## the class-split tables 129 / 132 are kept as tables). Every tested (tissue, promoter) of a two-promoter gene counted once; role =
## usage-defined major / minor (mirror_pairs tissue means). Read-outs, gene-clustered bootstrap B = 2000:
##   marks: mean share of the tissue's matched ENCODE samples with DNase / H3K27ac / H3K4me3 called at the TSS element (z > 1.64);
##   states: closed / enhancer-like (H3K27ac+ H3K4me3-) / promoter-like / open-unmarked, as 122a;
##   H3K27ac-to-H3K4me3 balance: share of samples with H3K27ac called among samples with either mark (how "enhancer-leaning" the TSS is).
F5S <- Sys.getenv("F5S", "/home/e1327989/f5_stage_20261106_v6")
src <- readLines(file.path(F5S, "code/122a_tissue_matched_chromatin_v4.R"))
eval(parse(text = src[1:(grep("^summ <- function", src) - 1)]))
mp <- MAPb[match %in% c("exact", "approximate")]
U1 <- unique(TP[, .(tissue, promoter, gene0, acc)])
ST <- state_for(U1, c("exact", "approximate"))
MK <- rbindlist(lapply(split(U1, U1$tissue), function(x) { bs <- mp[gtex == x$tissue[1]]$biosample_value
  if (!length(bs)) return(NULL); has <- !is.na(x$acc) & x$acc %in% rownames(A); if (!any(has)) return(NULL); x <- x[has]
  k27 <- K[x$acc, bs, drop = FALSE] > TH; k4 <- A[x$acc, bs, drop = FALSE] > TH
  x[, `:=`(DNase = rowMeans(D[acc, bs, drop = FALSE] > TH), H3K27ac = rowMeans(k27), H3K4me3 = rowMeans(k4),
           k27_balance = ifelse(rowSums(k27 | k4) > 0, rowSums(k27 & !k4) / rowSums(k27 | k4), NA_real_))] }))
MP <- fread(cmd = sprintf("zcat %s", file.path(F5S, "in/mirror/mirror_pairs.tsv.gz")), select = c("tissue", "canonical", "partner"))   # major = canonical = higher usage pooled over all tissues (one definition, as Figure 4)
RO <- rbind(MP[, .(tissue, promoter = canonical, role = "major")], MP[, .(tissue, promoter = partner, role = "minor")])
ST <- merge(ST, RO, by = c("tissue", "promoter")); MK <- merge(MK, RO, by = c("tissue", "promoter"))
bt <- function(d, f) { b <- gboot(d, f); c(est = f(d), lo = unname(quantile(b, .025, na.rm = TRUE)), hi = unname(quantile(b, .975, na.rm = TRUE))) }
R1 <- rbindlist(lapply(c("major", "minor"), function(r) rbindlist(lapply(c("DNase", "H3K27ac", "H3K4me3", "k27_balance"), function(m) { d <- MK[role == r & is.finite(get(m))]
  v <- bt(d, function(x) mean(x[[m]])); data.table(role = r, what = m, kind = "mark", units = nrow(d), est = v[["est"]], lo = v[["lo"]], hi = v[["hi"]]) }))))
R2 <- rbindlist(lapply(c("major", "minor"), function(r) rbindlist(lapply(c("no accessible element", "enhancer-like", "promoter-like", "accessible, unmarked"), function(st) { d <- ST[role == r]
  v <- bt(d, function(x) mean(x$state == st)); data.table(role = r, what = st, kind = "state", units = nrow(d), est = v[["est"]], lo = v[["lo"]], hi = v[["hi"]]) }))))
R <- rbind(R1, R2); fwrite(R, file.path(OUTD, "134_major_minor_pooled_chromatin.tsv"), sep = "\t")
options(width = 200); print(R[, .(role, what, kind, units, est = round(100 * est, 2), lo = round(100 * lo, 2), hi = round(100 * hi, 2))])
