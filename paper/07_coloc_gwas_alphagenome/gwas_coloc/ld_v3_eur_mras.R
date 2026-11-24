## Follow-up to ld_v3.R (run after it):
##   (1) EUR donor lists: donors/<tissue>.EUR.txt and ALL.EUR.txt (799 donors, the LD reference of the European GWAS) = the donor
##       set restricted to self-reported RACE == 3 (metadata/v10_cov_pca_full_GRU.txt).
##   (2) The MRAS check under the shared standard: v784 (chr3_138348784_C_T_b38) against the hypertension GWAS lead
##       (chr3_138336294_C_T_b38), in every tissue, in ALL, and in the EUR versions. Output: results/ld_v3/mras_v784_gwas_lead.tsv
suppressPackageStartupMessages(library(data.table))
source("/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/code/bohr2_data.R"); source(file.path(V2, "code/ld_tool.R"))
PH <- fread(file.path(GTEX, "metadata/v10_cov_pca_full_GRU.txt"), select = c("SUBJID", "RACE"))
eur <- PH[RACE == 3]$SUBJID
sets <- sub("[.]txt$", "", basename(Sys.glob(file.path(LD_DIR, "donors", "*.txt")))); sets <- sets[!grepl("[.]EUR$", sets)]
for (s in sets) writeLines(intersect(ld_donors(s), eur), file.path(LD_DIR, "donors", paste0(s, ".EUR.txt")))
cat("EUR lists:", length(sets), "; ALL", length(ld_donors("ALL")), "donors, EUR", length(ld_donors("ALL", TRUE)), "\n")
P <- data.table(tissue = sets, snp_a = "chr3_138348784_C_T_b38", snp_b = "chr3_138336294_C_T_b38")
R <- merge(ld_pairs(P)[, .(tissue, n_donors, r, r2)], ld_pairs(P, eur = TRUE)[, .(tissue, n_donors_eur = n_donors, r_eur = r, r2_eur = r2)], by = "tissue")
fwrite(R[order(tissue != "ALL", tissue)], file.path(LD_DIR, "mras_v784_gwas_lead.tsv"), sep = "\t")
print(R[tissue == "ALL"]); print(R[tissue != "ALL", .(tissues = .N, r2_min = round(min(r2), 3), r2_median = round(median(r2), 3), r2_max = round(max(r2), 3),
                                                     r2_eur_min = round(min(r2_eur), 3), r2_eur_max = round(max(r2_eur), 3))])
