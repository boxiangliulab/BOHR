## run_paper.sh step P (authored for the repository): the summaries of the parity scripts (workflow/parity/03, 06, 09) as
## PASS / FAIL lines. Thresholds for significance come from Supplementary Table 4a (environment SUPP).
## Usage: Rscript check_package.R cis|cis_exact <cmp.rds> <tissue> <chr> <label> <status.tsv>
##        Rscript check_package.R null_start <run FALSE.rds> <run TRUE.rds> <tissue> <chr> "" <status.tsv>
##        Rscript check_package.R conditional <cmp.rds> <tissue> <chr> "" <status.tsv>
##        Rscript check_package.R trait <cmp.rds> <tissue> <trait> "" <status.tsv>
a <- commandArgs(TRUE)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
mode <- a[1]; st <- a[length(a)]
thr_of <- function(t) supp(Sys.getenv("SUPP"), "S4a")[SMTSD == t, unit_p_threshold]
pv <- c("mu_F_logp", "mu_lrt_logp", "nu_firth_logp", "tau_firth_logp", "nu_lrt_logp", "tau_lrt_logp")
if (mode %in% c("cis", "cis_exact")) {
  S <- readRDS(a[2]); t <- a[3]; ch <- a[4]; lab <- a[5]; TB <- as.data.table(S$snp); U <- as.data.table(S$units)
  g <- function(cl, f) TB[column == cl][[f]]
  mu <- c("mu_F_logp", "mu_lrt_logp", "mu_est", "mu_se")
  add(sprintf("%s %s, %s: μ tests of every SNP", t, ch, lab), all(vapply(mu, g, 0, "agree_1e6") == 1) && all(vapply(mu, g, 0, "na_mismatch") == 0),
      sprintf("%d tests; agree within 1e-6: %s", g("mu_F_logp", "both_finite"), paste(sprintf("%s %.4f", mu, vapply(mu, g, 0, "agree_1e6")), collapse = ", ")))
  fi <- c("nu_firth_logp", "tau_firth_logp")
  ag <- vapply(fi, g, 0, "agree_1e6"); nm <- vapply(fi, g, 0, "na_mismatch")
  add(sprintf("%s %s, %s: p0 and p1 Firth tests of every SNP", t, ch, lab), all(nm == 0) && (mode == "cis" || all(ag == 1)),
      sprintf("p0 %d tests, agree within 1e-6 %.4f; p1 %d tests, %.4f%s", g("nu_firth_logp", "both_finite"), ag[1], g("tau_firth_logp", "both_finite"), ag[2],
              if (mode == "cis" && any(ag < 1)) " (the null start of the package default changes these; see the paper-engine run and the comparison below)" else ""))
  M <- U[!is.na(p) & !is.na(p_pk)]; thr <- thr_of(t)
  M[, tie := abs(leader_logp - llp_pk) <= 1e-8 * pmax(1, abs(leader_logp))]   # leaders with the same P: the choice between them is arbitrary
  rd <- M[, max(abs(p - p_pk) / p)]
  e <- M[estimable == TRUE & dedup_drop == FALSE]
  dec <- e[, sum((p <= thr) != (p_pk <= thr))]
  add(sprintf("%s %s, %s: units (ACAT P, leader, rules, significance at the tissue's threshold)", t, ch, lab),
      nrow(M) == nrow(U) && M[, all(estimable == est_pk & dedup_drop == dd_pk)] && dec == 0 &&
        (mode == "cis" || (rd <= 1e-6 && M[, all(leader == leader_pk | tie)])),   # with the default, P values and leaders may move where the null start finds a higher maximum
      sprintf("%d units, %d matched; max relative difference of the unit P %.3g; same leader %d%s; significant at P <= %.3g: stored %d, package %d, %d differ",
              nrow(U), nrow(M), rd, M[leader == leader_pk, .N],
              if (M[leader != leader_pk & tie, .N]) sprintf(" (other %d: a tie, two variants with the same P, e.g. in perfect LD)", M[leader != leader_pk & tie, .N]) else "",
              thr, e[p <= thr, .N], e[p_pk <= thr, .N], dec))
} else if (mode == "null_start") {
  A <- readRDS(a[2]); B <- readRDS(a[3]); t <- a[4]; ch <- a[5]
  sa <- as.data.table(A$snp); sb <- as.data.table(B$snp); k <- c("promoter", "variant"); k <- intersect(k, names(sa)); if (!length(k)) k <- c("promoter", "snp")
  m <- merge(sa, sb, by = k, suffixes = c(".f", ".t"))
  cols <- intersect(c("nu_firth_log_p", "tau_firth_log_p", "nu_firth_logp", "tau_firth_logp"), names(sa))
  dif <- vapply(cols, function(cl) { x <- m[[paste0(cl, ".f")]]; y <- m[[paste0(cl, ".t")]]; sum(abs(x - y) > 1e-6 * pmax(1, abs(x)), na.rm = TRUE) }, 0)
  ua <- as.data.table(A$units); ub <- as.data.table(B$units); ku <- intersect(c("promoter", "component", "type"), names(ua))
  mu <- merge(ua, ub, by = ku, suffixes = c(".f", ".t")); pc <- intersect(c("p", "acat_p"), names(ua))[1]
  thr <- thr_of(t); dec <- mu[, sum((get(paste0(pc, ".f")) <= thr) != (get(paste0(pc, ".t")) <= thr), na.rm = TRUE)]
  add(sprintf("%s %s: package default against the paper's engine (firth_null_start TRUE vs FALSE)", t, ch), dec == 0,
      sprintf("%d SNP rows; Firth log P differing by more than 1e-6: %s; units whose significance changes: %d of %d",
              nrow(m), paste(sprintf("%s %d", cols, dif), collapse = ", "), dec, nrow(mu)))
} else if (mode == "conditional") {
  S <- readRDS(a[2]); t <- a[3]; ch <- a[4]; K <- as.data.table(S$signals); M <- as.data.table(S$units)
  ks <- K[bk %in% TRUE, paste(promoter, type, k, snp)]; kp <- K[bk_pk %in% TRUE, paste(promoter, type, k, snp_pk)]
  add(sprintf("%s %s: independent signals kept after the backward step", t, ch), setequal(ks, kp) && M[, all(nk == nk_pk, na.rm = TRUE)],
      sprintf("%d units; signals kept: stored (v3c, published) %d, package %d, the same %d", nrow(M), length(ks), length(kp), length(intersect(ks, kp))))
} else if (mode %in% c("trait", "trait_exact")) {
  S <- readRDS(a[2]); t <- a[3]; tr <- a[4]; lab <- a[5]; K <- as.data.table(S$units)
  eqv <- function(x, y) sum((is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y))
  d <- K[is.finite(logp) & is.finite(log_p), abs(logp - log_p) <= 1e-6 * pmax(1, abs(logp))]
  ## with the default, the null start may move Firth P values (and with them the testability of a unit); the decisive
  ## comparison is significance; with firth_null_start = FALSE everything must agree
  ## eligibility can differ where the unpenalised p0/p1 fit of the promoter is separated: the stored run (NSCC) returned no
  ## statistic there and the unit was not counted as testable, the same script on this platform returns one (as in 04.2)
  de <- K[!((is.na(final) & is.na(eligible)) | (!is.na(final) & !is.na(eligible) & final == eligible))]
  raw <- file.path(Sys.getenv("TRG_ROOT", "/projects/group/bliu/e1327989/GTEx/tissue_regression"), "trait_regression_20261104/results/real", t, paste0(tr, ".rds"))
  sep_ok <- if (nrow(de) && file.exists(raw)) { X <- as.data.table(readRDS(raw)$res)
    all(mapply(function(pr, ty) { r <- X[promoter == pr]; nrow(r) == 1 && isTRUE(if (ty == "p0") r$nu_sep else if (ty == "p1") r$tau_sep else FALSE) || isTRUE(r$cat_complete_sep) },
               de$promoter, de$type)) } else nrow(de) == 0
  add(sprintf("%s x %s, %s: characteristic tests (eligibility, significance, P)", t, tr, lab),
      eqv(K$sig_final, K$significant) == nrow(K) && (mode == "trait" || (sep_ok && all(d))),
      paste0(sprintf("%d units; same eligibility %d, same significance %d (significant: published %d, package %d); log P within 1e-6: %d of %d",
              nrow(K), eqv(K$final, K$eligible), eqv(K$sig_final, K$significant), K[sig_final == TRUE, .N], K[significant == TRUE, .N], sum(d), length(d)),
      if (nrow(de)) sprintf("; eligibility differs for %d units (%s), all with a separated unpenalised fit in the stored run: %s", nrow(de),
                            paste(unique(de$promoter), collapse = ", "), sep_ok) else ""))
}
finish(st)
