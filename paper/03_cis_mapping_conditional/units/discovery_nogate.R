## Ungated discovery: no proportion gate in the primary analysis (gates are display filters only).
## Rule firthF as in final_discovery.R, but units = EVERY (promoter, type) whose part exists after the standard <= 5%
## endpoint rule: mu for every testable promoter (finite mu_F); p0 for BEZI/BEINF; p1 for BEOI/BEINF. Firth p for
## p0/p1 comes from the main run (old-gate promoters) or from the Firth add-on (all other promoters with a 0/1 part).
## Per unit: ACAT over the window of the own-type p (cap 1 - 1/M), BH over all units of the tissue, q <= 0.05.
## Old gate flags are kept as columns (elig_mu: interior share >= 0.70; elig_p0 / elig_p1: share >= 0.35).
## Also stored per unit: the leader's p and effect for the OTHER two types (same SNP), to study multi-type effects.
## Mirror rule (results/mirror/README.txt, table v2 with per-component exactness): a partner's mu unit
## duplicates the canonical's when interior_exact; its p0/p1 units duplicate the canonical's p1/p0 when boundary_exact.
## The PRIMARY BH (q) runs on the de-duplicated set. Dropped units carry duplicate_of and q = NA; surviving merged
## units carry merged_component (interior / boundary), twin and eq_class. q_full = BH over all units (reference only).
## Usage: NC=4 Rscript code/discovery_nogate.R <def> <arm> <tissue> [tissue ...]
## Output: results/final_nogate/<def>/<tissue>.<arm>.units.rds (only when main, muF and add-on are complete for the arm)
suppressPackageStartupMessages({ library(data.table); library(parallel) })
V2 <- Sys.getenv("BOHR2_ROOT", "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2")
a <- commandArgs(TRUE); DEF <- a[1]; ARM <- a[2]; TIS <- a[-(1:2)]; NC <- as.integer(Sys.getenv("NC", "4"))
OUT <- file.path(V2, "results", "final_nogate", DEF); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
acat_log <- function(logp) { p <- exp(logp[is.finite(logp)]); if (!length(p)) return(NA_real_); if (any(p == 0)) return(0)   # = final_discovery.R
  if (length(p) > 1) p <- pmin(p, 1 - 1 / length(p))
  sm <- p < 1e-15; t <- numeric(length(p)); t[!sm] <- tan((0.5 - p[!sm]) * pi); t[sm] <- 1 / (p[sm] * pi); st <- mean(t)
  if (st > 1e15) return(1 / (st * pi)); pcauchy(st, lower.tail = FALSE) }
FC <- c("nu_firth_logp", "nu_firth_est", "tau_firth_logp", "tau_firth_est", "nu_firth_se", "tau_firth_se")   # + the Firth standard errors
one_chr <- function(tis, i) {
  f <- function(root) file.path(V2, "results", root, DEF, tis, ARM, sprintf("chr%d.rds", i))
  M <- readRDS(f("final")); F <- readRDS(f("final_muF")); A <- readRDS(f("final_firth_addon"))
  s <- as.data.table(M$snp)[, c("promoter", "snp", "mu_est", "mu_se", "mu_pe", "p0_pe", "p1_pe", FC), with = FALSE]   # + mu_se
  s <- merge(s, as.data.table(F$snp)[, .(promoter, snp, mu_F_logp)], by = c("promoter", "snp"), all.x = TRUE)
  if (NROW(A$snp)) { ad <- as.data.table(A$snp)[, c("promoter", "snp", FC), with = FALSE]; setnames(ad, FC, paste0("a_", FC))
    s <- merge(s, ad, by = c("promoter", "snp"), all.x = TRUE)
    for (k in FC) set(s, j = k, value = as.numeric(s[[k]]))   # numeric first: an all-NA column is logical, and := would store the add-on's values as TRUE
    for (k in FC) s[is.na(get(k)) & !is.na(get(paste0("a_", k))), (k) := get(paste0("a_", k))]
    s[, paste0("a_", FC) := NULL] }
  list(snp = s, prom = as.data.table(M$prom)[, .(promoter, gene, family, skip, n, n0, n1, elig_mu, elig_p0, elig_p1)])
}
PCOL <- c(mu = "mu_F_logp", p0 = "nu_firth_logp", p1 = "tau_firth_logp")
EST <- c(mu = "mu_est", p0 = "nu_firth_est", p1 = "tau_firth_est"); PE <- c(mu = "mu_pe", p0 = "p0_pe", p1 = "p1_pe")
SE <- c(mu = "mu_se", p0 = "nu_firth_se", p1 = "tau_firth_se")   # standard error column of each type's estimate
for (tis in TIS) {
  out <- file.path(OUT, paste(tis, ARM, "units.rds", sep = "."))
  roots <- c("final", "final_muF", "final_firth_addon")
  ok <- all(vapply(roots, function(r) all(file.exists(file.path(V2, "results", r, DEF, tis, ARM, sprintf("chr%d.rds", 1:22)))), TRUE))
  if (!ok) { cat(tis, ARM, ": inputs incomplete, skipped\n"); next }
  X <- mclapply(1:22, function(i) one_chr(tis, i), mc.cores = NC)
  S <- rbindlist(lapply(X, `[[`, "snp")); P <- rbindlist(lapply(X, `[[`, "prom"))[is.na(skip) & !is.na(family)]
  types <- rbind(P[, .(promoter, type = "mu")], P[family %in% c("BEZI", "BEINF"), .(promoter, type = "p0")], P[family %in% c("BEOI", "BEINF"), .(promoter, type = "p1")])
  U <- rbindlist(lapply(names(PCOL), function(ty) { col <- PCOL[[ty]]; s <- S[promoter %in% types[type == ty, promoter]]
    s[, { j <- which.min(replace(get(col), !is.finite(get(col)), Inf))
          .(p = acat_log(get(col)), n_snp = sum(is.finite(get(col))), leader = snp[j], leader_logp = get(col)[j], leader_est = get(EST[[ty]])[j], leader_se = get(SE[[ty]])[j],   # leader_se: standard error of the leader's estimate (for odds-ratio intervals)
            leader_pe = get(PE[[ty]])[j], mu_logp_at_leader = mu_F_logp[j], p0_logp_at_leader = nu_firth_logp[j], p1_logp_at_leader = tau_firth_logp[j],
            mu_est_at_leader = mu_est[j], p0_est_at_leader = nu_firth_est[j], p1_est_at_leader = tau_firth_est[j]) }, by = promoter][, type := ty] }))
  U <- U[is.finite(p)]
  FL <- fread(cmd = paste("zcat", file.path(V2, "results/mirror/promoter_flags.tsv.gz")))[tissue == tis,
          .(promoter, n_prom_gene, mirror_pair_id, boundary_exact, interior_exact, mirror_role, mirror_partner, drop_mu, drop_boundary)]
  U <- merge(U, FL, by = "promoter", all.x = TRUE)
  for (k in c("boundary_exact", "interior_exact", "drop_mu", "drop_boundary")) U[is.na(get(k)), (k) := FALSE]
  ## per-component de-duplication (mirror table v2): mu(partner) == mu(canonical) when interior_exact;
  ## p0(partner) == p1(canonical) and p1(partner) == p0(canonical) when boundary_exact
  U[, dedup_drop := fifelse(type == "mu", drop_mu, drop_boundary)]
  U[, duplicate_of := fifelse(dedup_drop, paste0(mirror_partner, ":", c(mu = "mu", p0 = "p1", p1 = "p0")[type]), NA_character_)]
  U[, merged_component := fifelse(mirror_role %in% "canonical" & type == "mu" & interior_exact, "interior",
                         fifelse(mirror_role %in% "canonical" & type != "mu" & boundary_exact, "boundary", NA_character_))]
  U[, twin := fifelse(!is.na(merged_component), paste0(mirror_partner, ":", c(mu = "mu", p0 = "p1", p1 = "p0")[type]), NA_character_)]
  U[, eq_class := fifelse(!is.na(merged_component) | dedup_drop, paste0(mirror_pair_id, "|", fifelse(type == "mu", "interior",
                   fifelse((mirror_role == "canonical" & type == "p0") | (mirror_role == "partner" & type == "p1"), "canonical_off", "canonical_full"))), NA_character_)]
  U[, q_full := p.adjust(p, "BH")]
  U[, q := NA_real_]; U[dedup_drop == FALSE, q := p.adjust(p, "BH")]
  U <- merge(U, P[, .(promoter, gene, family, n, n0, n1, elig_mu, elig_p0, elig_p1)], by = "promoter")
  U[, `:=`(n_mid = n - n0 - n1, mid_share = (n - n0 - n1) / n, f0 = n0 / n, f1 = n1 / n,
           old_gate = fifelse(type == "mu", elig_mu, fifelse(type == "p0", elig_p0, elig_p1)), def = DEF, tissue = tis, arm = ARM)]
  tmp <- paste0(out, ".tmp"); saveRDS(U, tmp); file.rename(tmp, out)
  cat(sprintf("%s %s: %d ungated units (%d after mirror de-dup; old-gate %d); significant %d (old-gate subset %d); thr %.3g\n", tis, ARM, nrow(U),
              sum(!U$dedup_drop), sum(U$old_gate), sum(U$q <= 0.05, na.rm = TRUE), sum(U$q <= 0.05 & U$old_gate, na.rm = TRUE), U[q <= 0.05, if (.N) max(p) else NA_real_]))
}
