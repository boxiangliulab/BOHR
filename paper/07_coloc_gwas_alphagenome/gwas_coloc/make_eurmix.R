## Ancestry-matched LD for the GWAS loci: European and mostly European GWAS (Supplementary Table 7a, ancestry_group
## European / Multi-ancestry) take the European-donor clumping (clump_EUR, 799 donors); East Asian GWAS keep the 943-donor
## clumping (clump_ALL, byte-identical to the NSCC outputs of the earlier run). A lead monomorphic in the 799 European donors has no PLINK r2 row:
## it is kept as a lead with no partners (step 12 then uses lead +- 50 kb), recorded as n_partners = 0 and pos from its ID.
suppressMessages(library(data.table))
K <- "/projects/group/bliu/e1327989/GTEx/nofloor_full/consult/freeze_ld_20261124"; OUT <- file.path(K, "gwas/clump_EURmix")
S7 <- "/projects/group/bliu/e1327989/GTEx/manuscript/BOHR_v7/build/tables/supplementary/Supplementary_Table_S7a_gwas_provenance.tsv"
s7 <- fread(S7)[, .(trait = sub("\\.txt$", "", gwas_file), ancestry_group)]
tr <- readLines(file.path(K, "gwas/traits.txt")); stopifnot(setequal(tr, s7$trait))
dir.create(OUT, showWarnings = FALSE)
res <- rbindlist(lapply(tr, function(t) {
  g <- s7[trait == t, ancestry_group]; d <- file.path(OUT, t)
  if (g == "East Asian") { if (!file.exists(d)) file.symlink(file.path(K, "gwas/clump_ALL", t), d); return(data.table(trait = t, ancestry_group = g, source = "ALL", filled = 0L)) }
  dir.create(d, showWarnings = FALSE)
  file.copy(file.path(K, "gwas/clump_EUR", t, c("leads.tsv", "DONE")), d, overwrite = TRUE)
  q <- fread(file.path(K, "gwas/clump_EUR", t, "lead_partners.tsv"), colClasses = list(character = "lead"))
  n <- q[is.na(pos), .N]
  q[is.na(pos), `:=`(pos = as.integer(sub("^chr[0-9]+_([0-9]+)_.*$", "\\1", lead)), n_partners = 0L)]
  fwrite(q, file.path(d, "lead_partners.tsv"), sep = "\t")
  data.table(trait = t, ancestry_group = g, source = "EUR", filled = n) }))
fwrite(res, file.path(K, "gwas/clump_EURmix_sources.tsv"), sep = "\t")
print(res[, .(traits = .N, leads_filled = sum(filled)), by = .(ancestry_group, source)])
