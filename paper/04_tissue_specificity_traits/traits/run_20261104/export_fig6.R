## =============================================================================
## export_fig6.R -- export of the run 20261104 trait regression for Figure 6 (Supplementary Tables S9a/S9b: direction consistency, "double significant") (mirror rule v3)
## Input: results/summary/units_real.rds (all test units, including partner units removed by mirror de-duplication)
## Output: results/export/trait_regression_units.tsv.gz + README_trait_regression.txt
## =============================================================================
source("/projects/group/bliu/e1327989/GTEx/tissue_regression/trait_regression_20261104/code/data.R")
OUT <- file.path(PATHS$res, "export"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
U <- readRDS(file.path(PATHS$res, "summary", "units_real.rds"))
E <- U[testable == TRUE, .(tissue, trait, promoter, gene, type, p, logp, est, se, pe,
                           direction = sign(pe), final, q_final, sig_final, primary, q_all = q, sig_all = sig,
                           dedup_drop, duplicate_of, dedup_drop_v2, mirror_role, mirror_partner, boundary_exact, interior_exact,
                           n, n0, n1, n_mid, share, n_donors, trait_class, family)]
setorder(E, tissue, trait, promoter, type)
f <- file.path(OUT, "trait_regression_units.tsv.gz"); fwrite(E, f, sep = "\t")
cat("exported", nrow(E), "units to", f, "\n")
writeLines(c(
"Figure 3 trait regression, run 20261104 = the fits of run 20261007 with mirror rule v3 -- export for Figure 6 (Supplementary Tables S9a/S9b)",
paste0("Built ", format(Sys.time()), " from results/summary/units_real.rds. One row per (tissue, trait, promoter, type); every testable unit."),
"",
"WHAT CHANGED FROM v13 (tissue_regression/v13/output/v13_regression_<tissue>.regression_results_<trait>.{rds,txt})",
"  * parameter names: v13 mu / nu / tau  ->  run 20261007 type = mu / p0 / p1  (p0 = zero-inflation, p1 = one-inflation component)",
"  * there is NO joint likelihood-ratio test in run 20261007 (v13 lr_pvalue has no counterpart). Each unit is a 1-df test:",
"      mu : beta-part LRT with the small-sample F correction ; p0/p1 : Firth penalised LRT",
"  * p is the per-unit nominal p-value; logp = natural log(p).",
"  * samples: GTEx v11 freeze (SMAFRZE = RNASEQ), one sample per donor (highest RIN); v13 used every sample row.",
"  * covariates: SEX, AGE, PC1-5, DTHHRDY (factor; levels < 10 donors merged into the reference), MHABNWBC;",
"    the tested trait is dropped from the covariates. Donor phenotypes from the dbGaP GRU table (2026-10-03).",
"",
"EFFECT COLUMNS (sign is comparable across types: positive = the trait increases that component)",
"  est, se : link-scale coefficient for the contrast below (mu: logit of the interior mean; p0/p1: log-odds vs interior)",
"  pe      : probability-scale effect for the same contrast (mu: change in interior mean; p0/p1: change in P(y=0) / P(y=1))",
"  direction = sign(pe)",
"  contrast: continuous traits (AGE, BMI) = mean -> mean + 1 SD ;",
"            binary traits = unexposed (0) -> exposed (1) ; SEX = male -> female",
"",
"SIGNIFICANCE (Figure 3 final rules, trait_regression_20261007/DECISIONS.md)",
"  final     = primary & mu n_mid >= 80 & p0/p1 boundary share 0.10-0.90 & p0/p1 not completely separated & tissue >= 50 donors",
"  q_final   = BH within tissue x trait x type on 'final' ; sig_final = q_final <= 0.05",
"  q_all/sig_all = BH on all primary units (no estimability rules), for reference only",
"  For a nominal 'trait-regression P < 0.05' use p directly (all testable units are included, final or not).",
"",
"MIRROR PAIRS, rule v3 (author 2026-11-04; Figure 4 results/mirror_v3):",
"  boundary_exact two-promoter pairs: p0 is kept on BOTH promoters and p1 is dropped on both (p1(A) == p0(B)), so every exact",
"  on/off switch is counted once, as p0 on the promoter that loses all reads. v2 p1(major) -> v3 p0(minor); v2 p0(major) unchanged (major / minor = higher / lower mean usage pooled over tissues = mirror_role canonical / partner).",
"  mu (interior_exact): unchanged, kept on the major promoter, mu(minor) dropped. Non-exact pairs and genes with >= 3 promoters: unchanged.",
"  Dropped units stay here with dedup_drop = TRUE and duplicate_of = <twin promoter>:<type>; their p equals the twin's p (p0/p1 twins: same 0/1 response, same effect sign; mu twins: opposite sign).",
"  dedup_drop_v2 = the flag under the old rule v2. The v2 -> v3 unit map is results/summary/unit_map_v2_to_v3.tsv.gz.",
"  BH (q_final) was re-run within tissue x trait x type as the labels now fall. Run 20261007 (rule v2) is kept unchanged.",
"",
"TRAITS (20): SEX AGE BMI MHDRNKSTS MHSMKSTS MHASTHMA MHCOPD MHCVD MHHRTDIS MHHTN MHT1D MHT2D MHDPRSSN MHHRTATT",
"             MHRNLFLR MHARTHTS MHPNMNIA MHCANCERNM MHLVRDIS MHCOCAINE5",
"  Drinking and smoking are MHDRNKSTS / MHSMKSTS (yes/no; unexposed = non-drinker / non-smoker) since run 20261007;",
"  they replaced MHDRNKNMB / MHSMKNMB (counts per self-chosen period, recorded for drinkers / smokers only)."
), file.path(OUT, "README_trait_regression.txt"))
cat("README written\n")
