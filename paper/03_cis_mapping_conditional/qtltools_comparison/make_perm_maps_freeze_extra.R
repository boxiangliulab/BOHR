## Extra permutation maps for the freeze definition (author 2026-10-02, plan B: 50 tissues x 2 more permutations, BOHR and
## QTLtools). Same construction as prep_all.R's freeze perm0 (slot = genotype owner, donor = whose phenotype + PEER row the
## slot receives; sample(U) without replacement), seed 20261002 + 1000 * k + i (k = permutation index, i = tissue index in
## the sorted tissue list). Existing map files are NEVER overwritten: Adipose_Subcutaneous and Whole_Blood already have
## perm1-3 maps (the 2026-10-01 extra mu permutations) and are reused as they are.
## Output: perm_maps/final/freeze/<T>.perm{1,2}.perm_map.tsv; rows appended to perm_maps/MANIFEST.tsv
suppressPackageStartupMessages(library(data.table))
V2 <- "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2"; OUTM <- file.path(V2, "perm_maps/final/freeze")
TIS <- sort(sub("[.]rds$", "", list.files(file.path(V2, "work/final/freeze"), pattern = "[.]rds$"))); man <- list()
for (k in 1:2) for (i in seq_along(TIS)) { t <- TIS[i]; U <- readRDS(file.path(V2, "work/final/freeze", paste0(t, ".rds")))$U
  f <- file.path(OUTM, sprintf("%s.perm%d.perm_map.tsv", t, k))
  if (file.exists(f)) { m <- fread(f); how <- "existing (reused)" } else {
    set.seed(20261002L + 1000L * k + i); m <- data.table(slot = U, donor = sample(U)); fwrite(m, f, sep = "\t"); how <- sprintf("generated seed=%d", 20261002L + 1000L * k + i) }
  stopifnot(setequal(m$slot, U), setequal(m$donor, U), !anyDuplicated(m$donor), !anyDuplicated(m$slot))
  man[[length(man) + 1]] <- data.table(tissue = t, arm = paste0("perm", k), n_slots = nrow(m), n_fixed_points = sum(m$slot == m$donor), how = how,
                                       source = "code/make_perm_maps_freeze_extra.R (freeze)", md5 = unname(tools::md5sum(f))) }
M <- rbindlist(man); mf <- file.path(V2, "perm_maps/MANIFEST.tsv")
old <- fread(mf); new <- M[!paste(tissue, arm, source) %in% paste(old$tissue, old$arm, old$source)]
if (nrow(new)) fwrite(new, mf, sep = "\t", append = TRUE)
print(M[, .N, by = .(arm, how = sub(" seed=.*", "", how))]); cat("max fixed points:", max(M$n_fixed_points), "; manifest rows appended:", nrow(new), "\n")
