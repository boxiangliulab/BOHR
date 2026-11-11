## Real-data parity of bohr_trait(), step C (original R + data.table): compare one tissue x trait with the published export
## (trait_regression_20261106_mirror_abs/results/export_B/trait_regression_units.tsv.gz: P, estimates, q, significance and labels)
## and with the stored per-promoter rows (trait_regression_20261104/results/real/<tissue>/<trait>.rds: skip reasons).
## Prints counts, largest differences and agreement rates only (no sample or donor IDs).
## Usage: Rscript 09_compare_trait.R <run.rds> [summary.rds]
options(width = 200)
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); R <- readRDS(a[1]); OUTS <- if (length(a) > 1) a[2] else ""
TIS <- R$tissue; TR <- R$trait
EXP <- "/projects/group/bliu/e1327989/GTEx/tissue_regression/trait_regression_20261106_mirror_abs/results/export_B/trait_regression_units.tsv.gz"
RAW <- file.path("/projects/group/bliu/e1327989/GTEx/tissue_regression/trait_regression_20261104/results/real", TIS, paste0(TR, ".rds"))
E <- fread(cmd = sprintf("zcat %s | awk -F'\\t' 'NR == 1 || ($1 == \"%s\" && $2 == \"%s\")'", EXP, TIS, TR))
P <- as.data.table(R$results); setnames(P, "component", "type")
cat(sprintf("==== %s x %s: export_B rows %d, package rows %d\n", TIS, TR, nrow(E), nrow(P)))
M <- merge(E, P, by = c("promoter", "type"), suffixes = c(".st", ".pk"), all = TRUE)
cat(sprintf("units matched %d; only in export %d; only in package %d\n", M[!is.na(p.st) & !is.na(p.pk), .N], M[is.na(p.pk), .N], M[is.na(p.st), .N]))
K <- M[!is.na(p.st) & !is.na(p.pk)]
dd <- function(a, b, lab) { ok <- is.finite(a) & is.finite(b); d <- abs(a[ok] - b[ok])
  data.table(stat = lab, both = sum(ok), na_mismatch = sum(is.finite(a) != is.finite(b)), exact = sum(d == 0), max_abs_diff = if (any(ok)) max(d) else NA_real_,
             max_rel_diff = if (any(ok)) max(d / pmax(abs(a[ok]), 1e-300)) else NA_real_, agree_1e6 = if (any(ok)) mean(d <= 1e-6 * pmax(1, abs(a[ok]))) else NA_real_) }
print(rbind(dd(K$logp, K$log_p, "log P"), dd(K$p.st, K$p.pk, "P"), dd(K$est.st, K$est.pk, "estimate (per contrast)"), dd(K$se.st, K$se.pk, "SE"),
            dd(K$pe, K$effect, "effect"), dd(K$q_final, K$q, "q (within tissue x trait x component)")), digits = 3)
eqv <- function(a, b) sum((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b))
## mirror role: the export names it for boundary_exact pairs only; the package (stage 3b-1) for every pair, as bohr_cis()
KB <- K[boundary_exact.st %in% TRUE]
cat(sprintf("agreement over %d units: eligible (final) %d; significant %d; mirror drop %d; BH family %d; relabelled %d; n/n0/n1/n_mid %d; n_donors %d; boundary_exact %d; interior_exact %d; mirror role (boundary_exact pairs) %d of %d; roles of other pairs named by the package only: %d\n",
            nrow(K), eqv(K$final, K$eligible), eqv(K$sig_final, K$significant), eqv(K$dedup_drop.st, K$dedup_drop.pk), eqv(K$bh_family.st, K$bh_family.pk),
            eqv(K$abs_relabel.st, K$abs_relabel.pk), K[n.st == n.pk & n0.st == n0.pk & n1.st == n1.pk & n_mid.st == n_mid.pk, .N], eqv(K$n_donors.st, K$n_donors.pk),
            eqv(K$boundary_exact.st, K$boundary_exact.pk), eqv(K$interior_exact.st, K$interior_exact.pk),
            eqv(fifelse(KB$mirror_role.st == "", NA_character_, KB$mirror_role.st), KB$mirror_role.pk), nrow(KB),
            K[!(boundary_exact.st %in% TRUE) & !is.na(mirror_role.pk), .N]))
cat(sprintf("significant: export %d (mu %d, p0 %d, p1 %d), package %d (mu %d, p0 %d, p1 %d); eligible: export %d, package %d\n",
            E[sig_final == TRUE, .N], E[sig_final == TRUE & type == "mu", .N], E[sig_final == TRUE & type == "p0", .N], E[sig_final == TRUE & type == "p1", .N],
            P[significant == TRUE, .N], P[significant == TRUE & type == "mu", .N], P[significant == TRUE & type == "p0", .N], P[significant == TRUE & type == "p1", .N],
            E[final == TRUE, .N], P[eligible == TRUE, .N]))
qs <- K[final == TRUE & eligible == TRUE]
cat(sprintf("smallest |log(q/0.05)| among eligible units: %.3g\n", qs[, min(abs(log(q_final / 0.05)))]))
L <- R$labels
cat(sprintf("absolute-activity labels: %d significant exact-pair events judged; reported as p1 at the partner: %d (export: %d rows with abs_relabel, %d events)\n",
            if (is.null(L)) 0L else nrow(L), if (is.null(L)) 0L else sum(L$changed), E[abs_relabel == TRUE, .N], E[abs_relabel == TRUE, .N] / 2))
## per-promoter skips against the stored raw rows
X <- as.data.table(readRDS(RAW)$res); Fp <- as.data.table(R$features)
F2 <- merge(X[, .(promoter, skip_st = sub(":.*", "", skip), family_st = family)], Fp[, .(promoter, skip_pk = sub(":.*", "", skip), family_pk = family)], by = "promoter", all = TRUE)
F2[skip_pk == "trait_rare_in_feature", skip_pk := "trait_rare_in_promoter"]
cat(sprintf("features: stored %d, package %d, matched %d; same skip reason %d; same family %d\n", nrow(X), nrow(Fp), F2[!is.na(skip_st) | !is.na(family_st)][!is.na(skip_pk) | !is.na(family_pk), .N],
            eqv(F2$skip_st, F2$skip_pk), eqv(F2$family_st, F2$family_pk)))
print(F2[, .N, by = .(skip_st, skip_pk)][order(-N)])
## published examples
ex <- rbind(data.table(promoter = "PromoterID_34648", type = "p0"), data.table(promoter = "PromoterID_34646", type = "p1"),
            data.table(promoter = "PromoterID_36592", type = "mu"), data.table(promoter = "PromoterID_36602", type = "mu"))
print(merge(ex, K[, .(promoter, type, p_export = signif(p.st, 4), p_package = signif(p.pk, 4), q_export = signif(q_final, 4), q_package = signif(q, 4),
                      est_export = signif(est.st, 5), est_package = signif(est.pk, 5), sig_export = sig_final, sig_package = significant)], by = c("promoter", "type")))
if (nzchar(OUTS)) saveRDS(list(units = K, features = F2), OUTS)
cat(sprintf("package run: wall %.0f s, CPU %.0f s\n", R$wall, R$cpu))
