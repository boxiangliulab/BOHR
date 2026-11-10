#' Simulate a complete synthetic input set
#'
#' Generates every input of the genomics layer, for demonstrations and tests:
#' individuals with genotypes and two traits, several tissues with
#' per-tissue covariates, genes with two or three promoters, their absolute
#' and relative activity, a TSS table and the true effects. With `out_dir`,
#' everything is also written to files in the formats that [bohr_cis()]
#' reads, including a small VCF file. No real data ship with the package; this
#' is not the simulation benchmark of the BOHR paper.
#'
#' **Genes.** Genes are placed on `n_chromosomes` chromosomes, 400 kb apart.
#' Their types cycle through four layouts of external promoters (A is the
#' primary promoter):
#' \describe{
#'   \item{`exact`}{two promoters, A and B; the relative activities sum to 1,
#'     so the pair is an exact mirror pair (`boundary_exact` and
#'     `interior_exact` in [mirror_pairs()]).}
#'   \item{`three`}{three promoters, all tested: B and C share what A leaves,
#'     each silent in some samples.}
#'   \item{`boundary_only`}{A and B tested; a third promoter C, active only
#'     when A is neither silent nor alone, is not tested but stays in the
#'     denominator: the pair is `boundary_exact` but not `interior_exact`.}
#'   \item{`nonexact`}{A and B tested; C, active in about 30\% of samples, is
#'     not tested but stays in the denominator: not a mirror pair.}
#' }
#' Every third gene also has an internal promoter (`internal = TRUE` in the
#' TSS table), which has an absolute activity but is not in the denominator.
#' The untested and internal promoters are in `absolute` and `tss` but not in
#' `relative`, the matrix that is tested.
#'
#' **Activities.** In each sample, A's relative activity is 0 with
#' probability \eqn{p_0}, 1 with probability \eqn{p_1}, and otherwise a beta
#' value with mean \eqn{\mu} (dispersion 0.3, within 0.005 and 0.995), the
#' model of [simulate_response()]. \eqn{p_0}, \eqn{p_1} and \eqn{\mu} change
#' linearly with the dosage of the gene's causal variant (`effects`, per
#' allele, probability scale); \eqn{\mu} also depends on the covariates `PC1`
#' and `F1` (logit scale). B and C take shares of the rest. Absolute
#' activities are the shares times a gene-level total (rounded to 4
#' decimals); a gene is not expressed (all promoters 0, relative activity
#' `NA`) in about 2\% of samples. The relative activity is each external
#' promoter's absolute activity over the sum of its gene's external
#' promoters, rounded to 12 decimals (so that it is the same after a round
#' trip through a file).
#'
#' **Effects.** The first `round(prop_causal * n_genes)` genes (so every gene
#' type is among them, and among the others) have a causal variant (allele
#' frequency 0.3 to 0.5) within 50 kb of A, acting on p0, mu, p1 or p0 and mu
#' (in turn); in every tissue for alternate causal genes, in the first tissue
#' only for the others. The two variants on each side of a causal variant are
#' in linkage with it (copies with 10\% of the genotypes redrawn), so a unit's
#' leader need not be the causal variant. The binary trait raises p1 of A in
#' the second gene in the first tissue (by 0.15), and the continuous trait
#' raises mu of A in the third gene in every tissue (by 0.08 per standard
#' deviation). The truth table gives, for every tested
#' promoter affected (A and the promoters that share its gene), the true
#' change of each component on the scale of [bohr_test()]'s effects: dosage 2
#' versus 0, binary trait 1 versus 0, continuous trait one standard deviation
#' above the mean versus the mean, averaged over the tissue's samples.
#'
#' Random numbers come from the current R session; call [set.seed()] first for
#' reproducible data.
#'
#' @param n_individuals Number of individuals (donors).
#' @param n_tissues Number of tissues (`T1`, `T2`, ...); each individual
#'   contributes a sample to each tissue with probability 0.85.
#' @param n_genes Number of genes.
#' @param n_chromosomes Number of chromosomes (`chr1`, ...).
#' @param variant_spacing Distance between neighbouring variants, in base
#'   pairs.
#' @param prop_causal Share of genes with a causal variant.
#' @param effects Per-allele effects of a causal variant on mu, p0 and p1
#'   (probability scale), a named vector with values between -0.2 and 0.2.
#'   The defaults are large, so that the demonstrations find the effects with
#'   a few hundred samples.
#' @param out_dir Optional folder (created if needed) to which the input set
#'   is written; see Value.
#' @return A list of class `"bohr_sim"`:
#'   \describe{
#'     \item{`relative`}{Relative activity of the tested promoters (rows) in
#'       every sample (columns) of every tissue: the `phenotype` of
#'       [bohr_cis()].}
#'     \item{`absolute`}{Absolute activity of every promoter (tested,
#'       untested and internal).}
#'     \item{`samples`}{The sample table: `sample`, `individual` and `tissue`.}
#'     \item{`tss`}{The TSS table, one row per promoter: `id`, `gene`, `chr`,
#'       `pos`, `strand`, `internal` and `tested`.}
#'     \item{`genotype`, `variants`}{Dosages (variants x individuals, about
#'       1\% missing) and their positions (`id`, `chr`, `pos`, `end`, `ref`,
#'       `alt`; about 10\% are indels).}
#'     \item{`covariates`}{A list by tissue of covariate matrices (samples x
#'       `PC1`, `PC2`, `F1`).}
#'     \item{`traits`}{One row per individual: `binary_trait` (0/1) and
#'       `continuous_trait`.}
#'     \item{`truth`}{The true effects: `tissue`, `gene`, `gene_type`,
#'       `promoter`, `role` (`"primary"`, `"mirror"` or `"shared"`),
#'       `predictor` (variant ID or trait), `predictor_type`, `effect_mu`,
#'       `effect_p0` and `effect_p1`. Promoters without a row have no
#'       effect.}
#'     \item{`files`}{With `out_dir`, the paths written: `relative_activity.tsv`,
#'       `absolute_activity.tsv`, `samples.tsv`, `tss.tsv`, `genotypes.tsv`,
#'       `genotypes.vcf`, `variants.tsv`, `traits.tsv`, `truth.tsv` and the
#'       folder `covariates` with one `<tissue>.tsv` per tissue. Matrices have
#'       their row IDs in the first column.}
#'   }
#' @seealso [bohr_cis()], [simulate_response()] for one response.
#' @examples
#' set.seed(1)
#' sim <- simulate_bohr(n_individuals = 100, n_tissues = 2, n_genes = 4)
#' sim
#' head(sim$truth)
#' ## one gene's promoters in the first samples
#' sim$relative[1:2, 1:6]
#' ## written to files
#' d <- file.path(tempdir(), "bohr_sim")
#' sim <- simulate_bohr(n_individuals = 60, n_tissues = 2, n_genes = 2, out_dir = d)
#' basename(sim$files)
#' unlink(d, recursive = TRUE)
#' @export
simulate_bohr <- function(n_individuals = 200, n_tissues = 3, n_genes = 12, n_chromosomes = 2, variant_spacing = 50000,
                          prop_causal = 0.5, effects = c(mu = 0.1, p0 = 0.2, p1 = 0.2), out_dir = NULL) {
  n_individuals <- check_count(n_individuals, "n_individuals", min = 20)
  n_tissues <- check_count(n_tissues, "n_tissues", min = 1); n_genes <- check_count(n_genes, "n_genes", min = 1)
  n_chromosomes <- check_count(n_chromosomes, "n_chromosomes", min = 1)
  if (!is.numeric(variant_spacing) || length(variant_spacing) != 1L || !is.finite(variant_spacing) || variant_spacing < 100)
    stop("`variant_spacing` must be a single number >= 100.", call. = FALSE)
  if (!is.numeric(prop_causal) || length(prop_causal) != 1L || !(prop_causal >= 0 && prop_causal <= 1)) stop("`prop_causal` must lie in [0, 1].", call. = FALSE)
  if (!is.numeric(effects) || !setequal(names(effects), c("mu", "p0", "p1")) || any(!is.finite(effects)) || any(abs(effects) > 0.2))
    stop("`effects` must be a vector c(mu = , p0 = , p1 = ) of per-allele effects between -0.2 and 0.2.", call. = FALSE)
  if (!is.null(out_dir) && !is_path(out_dir)) stop("`out_dir` must be the path of one folder.", call. = FALSE)
  w <- nchar(n_individuals)
  ind <- sprintf(paste0("IND%0", max(3L, w), "d"), seq_len(n_individuals)); tis <- paste0("T", seq_len(n_tissues))

  ## ---- genes, promoters and variants ----
  types <- rep(c("exact", "three", "boundary_only", "nonexact"), length.out = n_genes)
  gchr <- paste0("chr", rep(seq_len(n_chromosomes), length.out = n_genes))
  gk <- stats::ave(seq_len(n_genes), gchr, FUN = seq_along)                 # rank of the gene on its chromosome
  gstart <- as.integer(3e5 + (gk - 1) * 4e5); gstrand <- ifelse(seq_len(n_genes) %% 2 == 1, "+", "-")
  gid <- sprintf("G%03d", seq_len(n_genes))
  P <- list(); pid <- 0L
  for (j in seq_len(n_genes)) {
    roles <- switch(types[j], exact = c("A", "B"), c("A", "B", "C"))
    if (j %% 3 == 0) roles <- c(roles, "I")
    off <- c(A = 0L, B = 15000L, C = 30000L, I = 60000L)[roles] * (if (gstrand[j] == "+") 1L else -1L)
    ids <- sprintf("P%04d", pid + seq_along(roles)); pid <- pid + length(roles)
    P[[j]] <- data.frame(id = ids, gene = gid[j], chr = gchr[j], pos = gstart[j] + off, strand = gstrand[j], internal = roles == "I",
                         tested = roles %in% c("A", "B") | (roles == "C" & types[j] == "three"), role = roles, stringsAsFactors = FALSE)
  }
  P <- do.call(rbind, P); rownames(P) <- NULL
  chr_len <- vapply(paste0("chr", seq_len(n_chromosomes)), function(ch) as.integer(max(c(0, gstart[gchr == ch])) + 7e5), 0L)
  V <- do.call(rbind, lapply(names(chr_len), function(ch) data.frame(chr = ch, pos = as.integer(seq(variant_spacing / 2, chr_len[[ch]], by = variant_spacing)),
                                                                      stringsAsFactors = FALSE)))
  nv <- nrow(V); bases <- c("A", "C", "G", "T")
  V$ref <- sample(bases, nv, replace = TRUE)
  V$alt <- vapply(V$ref, function(r) sample(setdiff(bases, r), 1L), "")
  indel <- stats::runif(nv) < 0.1
  V$ref[indel] <- vapply(which(indel), function(i) paste0(V$alt[i], paste(sample(bases, sample(1:3, 1L), replace = TRUE), collapse = "")), "")
  V$end <- V$pos + nchar(V$ref) - 1L
  V$id <- paste(V$chr, V$pos, V$ref, V$alt, sep = "_")
  V <- V[, c("id", "chr", "pos", "end", "ref", "alt")]
  freq <- stats::runif(nv, 0.15, 0.5)
  ## causal variants: one per causal gene, within 50 kb of its primary promoter
  causal <- seq_len(n_genes) <= max(round(prop_causal * n_genes), as.integer(prop_causal > 0))
  cvar <- rep(NA_integer_, n_genes)
  for (j in which(causal)) {
    on <- which(V$chr == gchr[j]); target <- gstart[j] + stats::runif(1, -5e4, 5e4)
    cvar[j] <- on[which.min(abs(V$pos[on] - target))]; freq[cvar[j]] <- stats::runif(1, 0.3, 0.5)
  }
  G <- t(vapply(freq, function(q) as.integer(stats::rbinom(n_individuals, 2L, q)), integer(n_individuals)))
  dimnames(G) <- list(V$id, ind)
  ## linkage: the two variants on each side of a causal variant copy its genotypes, with 10% of them redrawn
  for (cv in stats::na.omit(cvar)) for (nb in cv + c(-2L, -1L, 1L, 2L)) if (nb >= 1L && nb <= nv && V$chr[nb] == V$chr[cv] && !nb %in% cvar) {
    G[nb, ] <- G[cv, ]; re <- stats::runif(n_individuals) < 0.1
    G[nb, re] <- as.integer(stats::rbinom(sum(re), 2L, freq[cv])); freq[nb] <- freq[cv]
  }
  ccomp <- rep(list(character(0)), n_genes); cwhere <- rep("all", n_genes); k <- 0L
  for (j in which(causal)) { k <- k + 1L; ccomp[[j]] <- list("p0", "mu", "p1", c("p0", "mu"))[[(k - 1L) %% 4L + 1L]]
    cwhere[j] <- if (k %% 2L == 1L || n_tissues == 1L) "all" else tis[1] }

  ## ---- individuals: traits and genotype PCs; samples and tissue factors ----
  traits <- data.frame(individual = ind, binary_trait = stats::rbinom(n_individuals, 1L, 0.5),
                       continuous_trait = dec(stats::rnorm(n_individuals, 50, 10), 2), stringsAsFactors = FALSE)
  pcs <- matrix(dec(stats::rnorm(2 * n_individuals), 6), n_individuals, 2, dimnames = list(ind, c("PC1", "PC2")))
  S <- do.call(rbind, lapply(tis, function(tt) {
    inn <- ind[stats::runif(n_individuals) < 0.85]
    if (length(inn) < 20L) inn <- ind[sort(sample.int(n_individuals, min(n_individuals, 20L)))]
    data.frame(sample = paste0(inn, "-", tt), individual = inn, tissue = tt, stringsAsFactors = FALSE)
  }))
  rownames(S) <- NULL; ns <- nrow(S)
  f1 <- dec(stats::rnorm(ns), 6)
  cov_list <- lapply(stats::setNames(tis, tis), function(tt) {
    i <- which(S$tissue == tt)
    m <- cbind(pcs[S$individual[i], , drop = FALSE], F1 = f1[i]); rownames(m) <- S$sample[i]
    m
  })

  ## ---- activities ----
  gi <- match(S$individual, ind)
  bt <- traits$binary_trait[gi]; ct <- (traits$continuous_trait[gi] - 50) / 10
  shift <- 0.25 * pcs[S$individual, "PC1"] + 0.3 * f1                            # covariate effect on logit(mu) of A
  share <- matrix(0, nrow(P), ns, dimnames = list(P$id, S$sample))
  truth <- list()
  for (j in seq_len(n_genes)) {
    pj <- P[P$gene == gid[j], , drop = FALSE]
    p0b <- stats::runif(1, 0.1, 0.2); p1b <- if (types[j] == "exact" && j %% 2 == 1 && !"p1" %in% ccomp[[j]]) 0 else stats::runif(1, 0.06, 0.12)
    mub <- stats::runif(1, 0.3, 0.6) + stats::rnorm(n_tissues, 0, 0.03); names(mub) <- tis
    g <- if (causal[j]) G[cvar[j], gi] else rep(0L, ns)
    act <- if (causal[j]) S$tissue %in% (if (cwhere[j] == "all") tis else cwhere[j]) else rep(FALSE, ns)
    eff <- c(mu = 0, p0 = 0, p1 = 0); eff[ccomp[[j]]] <- effects[ccomp[[j]]]
    tau1 <- if (j == 2L) 0.15 * (S$tissue == tis[1]) else rep(0, ns)          # binary trait on p1 of A (gene 2, first tissue)
    kap <- if (j == 3L) 0.08 else 0                                           # continuous trait on mu of A (gene 3)
    comp_A <- function(g, b, c) {
      p0 <- pmin(pmax(p0b + eff[["p0"]] * g * act, 0), 0.9); p1 <- pmin(pmax(p1b + eff[["p1"]] * g * act + tau1 * b, 0), 0.9)
      tot <- p0 + p1; sc <- ifelse(tot > 0.95, 0.95 / tot, 1); p0 <- p0 * sc; p1 <- p1 * sc
      m <- pmin(pmax(mub[S$tissue] + eff[["mu"]] * g * act + kap * c, 0.05), 0.95)
      list(p0 = p0, p1 = p1, mu = stats::plogis(stats::qlogis(m) + shift))
    }
    A <- comp_A(g, bt, ct)
    u <- stats::runif(ns); phi <- (1 - 0.3^2) / 0.3^2
    a <- pmin(pmax(stats::rbeta(ns, A$mu * phi, (1 - A$mu) * phi), 0.005), 0.995)
    a[u < A$p0] <- 0; a[u >= A$p0 & u < A$p0 + A$p1] <- 1
    ## the other external promoters share what A leaves (B gets the share s of 1 - a)
    pr <- split_par(types[j])
    s <- switch(types[j], exact = rep(1, ns),
                boundary_only = ifelse(a %in% c(0, 1), 1, 1 - stats::runif(ns, 0.05, 0.25)),
                nonexact = 1 - stats::runif(ns, 0.05, 0.25) * (stats::runif(ns) < pr$pi),
                three = { v <- pmin(pmax(stats::rbeta(ns, 0.5 * 6, 0.5 * 6), 0.02), 0.98); r <- stats::runif(ns)
                          v[r < pr$r] <- 0; v[r >= pr$r & r < pr$r + pr$q] <- 1; v })
    rid <- stats::setNames(pj$id, pj$role)
    share[rid[["A"]], ] <- a; share[rid[["B"]], ] <- (1 - a) * s
    if ("C" %in% names(rid)) share[rid[["C"]], ] <- (1 - a) * (1 - s)
    ## truth: the change of each tested promoter's components, averaged over the tissue's samples
    preds <- list()
    if (causal[j]) preds[[V$id[cvar[j]]]] <- list(type = "dosage", lo = function() comp_A(0, bt, ct), hi = function() comp_A(2, bt, ct))
    if (j == 2L) preds$binary_trait <- list(type = "binary", lo = function() comp_A(g, 0, ct), hi = function() comp_A(g, 1, ct))
    if (j == 3L) preds$continuous_trait <- list(type = "continuous", lo = function() comp_A(g, bt, 0), hi = function() comp_A(g, bt, 1))
    for (pn in names(preds)) {
      lo <- preds[[pn]]$lo(); hi <- preds[[pn]]$hi()
      for (ro in intersect(c("A", "B", "C"), pj$role[pj$tested])) {
        clo <- derived_comp(lo, types[j], ro, pr); chi <- derived_comp(hi, types[j], ro, pr)
        for (tt in tis) {
          i <- S$tissue == tt; d <- vapply(c("mu", "p0", "p1"), function(cc) mean(chi[[cc]][i] - clo[[cc]][i]), 0)
          if (all(abs(d) < 1e-12)) next
          truth[[length(truth) + 1L]] <- data.frame(tissue = tt, gene = gid[j], gene_type = types[j], promoter = rid[[ro]],
                                                    role = c(A = "primary", B = if (types[j] == "exact") "mirror" else "shared", C = "shared")[[ro]],
                                                    predictor = pn, predictor_type = preds[[pn]]$type,
                                                    effect_mu = d[["mu"]], effect_p0 = d[["p0"]], effect_p1 = d[["p1"]], stringsAsFactors = FALSE)
        }
      }
    }
  }
  ## absolute activity: shares times a gene total; a gene is silent in ~2% of samples; internal promoters on their own
  absm <- matrix(0, nrow(P), ns, dimnames = dimnames(share))
  for (j in seq_len(n_genes)) {
    r <- which(P$gene == gid[j] & !P$internal)
    tot <- exp(stats::rnorm(ns, log(stats::runif(1, 3, 9)), 0.4)) * (stats::runif(ns) >= 0.02)
    absm[r, ] <- dec(sweep(share[r, , drop = FALSE], 2, tot, "*"), 4)
    ri <- which(P$gene == gid[j] & P$internal)
    if (length(ri)) absm[ri, ] <- dec(stats::runif(ns, 0.5, 3) * (stats::runif(ns) < 0.6), 4)
  }
  ext <- !P$internal
  den <- rowsum(absm[ext, , drop = FALSE], P$gene[ext], reorder = FALSE)
  rel <- absm[P$tested, , drop = FALSE] / den[P$gene[P$tested], , drop = FALSE]
  rel[!is.finite(rel)] <- NA_real_
  rel <- dec(rel, 12); dimnames(rel) <- list(P$id[P$tested], S$sample)
  G_obs <- G; G_obs[sample.int(length(G_obs), round(0.01 * length(G_obs)))] <- NA_integer_
  truth <- if (length(truth)) do.call(rbind, truth) else
    data.frame(tissue = character(0), gene = character(0), gene_type = character(0), promoter = character(0), role = character(0),
               predictor = character(0), predictor_type = character(0), effect_mu = numeric(0), effect_p0 = numeric(0), effect_p1 = numeric(0))
  rownames(truth) <- NULL
  out <- list(relative = rel, absolute = absm, samples = S, tss = P[, c("id", "gene", "chr", "pos", "strand", "internal", "tested")],
              genotype = G_obs, variants = V, covariates = cov_list, traits = traits, truth = truth, files = NULL)
  if (!is.null(out_dir)) out$files <- write_sim(out, out_dir, chr_len)
  structure(out, class = "bohr_sim")
}

## x rounded to d decimals through its decimal string: the double a file reader returns for the written value
dec <- function(x, d) { i <- !is.na(x); x[i] <- as.numeric(sprintf(paste0("%.", d, "f"), x[i])); x }
## parameters of the split of 1 - A between B and C (see simulate_bohr Details)
split_par <- function(type) list(pi = 0.3, ec = 0.15, r = 0.15, q = 0.15, m = 0.5)
## components (p0, p1, mu per sample) of promoter `role` given those of A, from the generative model of each gene type
derived_comp <- function(A, type, role, pr) {
  p0 <- A$p0; p1 <- A$p1; mu <- A$mu; mid <- 1 - p0 - p1
  if (role == "A") return(A)
  if (type == "exact") return(list(p0 = p1, p1 = p0, mu = 1 - mu))
  if (type == "boundary_only") return(list(p0 = p1, p1 = p0, mu = (1 - mu) * (1 - pr$ec)))
  if (type == "nonexact") {
    num <- mid * (1 - mu) * (1 - pr$pi * pr$ec) + p0 * pr$pi * (1 - pr$ec); den <- mid + p0 * pr$pi
    return(list(p0 = p1, p1 = p0 * (1 - pr$pi), mu = num / den))
  }
  ## three: B = (1 - A) s, C = (1 - A)(1 - s); P(s = 0) = r, P(s = 1) = q, else mean m
  r <- pr$r; q <- pr$q; m <- if (role == "B") pr$m else 1 - pr$m
  if (role == "C") { tmp <- r; r <- q; q <- tmp }                            # C = 0 when s = 1, C = 1 when s = 0
  num <- mid * (1 - mu) * (q + (1 - q - r) * m) + p0 * (1 - q - r) * m
  den <- mid * (1 - r) + p0 * (1 - q - r)
  list(p0 = p1 + (1 - p1) * r, p1 = p0 * q, mu = num / den)
}

## write the input set in the formats bohr_cis() reads; returns the paths
write_sim <- function(x, out_dir, chr_len) {
  dir.create(file.path(out_dir, "covariates"), recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out_dir)) stop("Cannot create `out_dir`: ", out_dir, call. = FALSE)
  f <- function(name) file.path(out_dir, name)
  mat <- function(m, first) { d <- data.frame(rownames(m), m, check.names = FALSE, stringsAsFactors = FALSE); names(d)[1] <- first; d }
  write_tsv(mat(x$relative, "id"), f("relative_activity.tsv"))
  write_tsv(mat(x$absolute, "id"), f("absolute_activity.tsv"))
  write_tsv(x$samples, f("samples.tsv"))
  write_tsv(x$tss, f("tss.tsv"))
  write_tsv(mat(x$genotype, "id"), f("genotypes.tsv"))
  write_tsv(x$variants, f("variants.tsv"))
  write_tsv(x$traits, f("traits.tsv"))
  write_tsv(x$truth, f("truth.tsv"))
  cv <- vapply(names(x$covariates), function(tt) { p <- file.path(out_dir, "covariates", paste0(tt, ".tsv")); write_tsv(mat(x$covariates[[tt]], "sample"), p); p }, "")
  write_vcf(x$genotype, x$variants, f("genotypes.vcf"), chr_len)
  c(relative = f("relative_activity.tsv"), absolute = f("absolute_activity.tsv"), samples = f("samples.tsv"), tss = f("tss.tsv"),
    genotype = f("genotypes.tsv"), vcf = f("genotypes.vcf"), variants = f("variants.tsv"), traits = f("traits.tsv"),
    truth = f("truth.tsv"), covariates = file.path(out_dir, "covariates"), stats::setNames(cv, paste0("covariates_", names(cv))))
}
## a minimal VCF 4.2 with unphased GT calls (dosage = number of ALT alleles; NA -> ./.)
write_vcf <- function(G, V, path, chr_len) {
  gt <- matrix(c("0/0", "0/1", "1/1")[G + 1L], nrow(G)); gt[is.na(gt)] <- "./."
  hdr <- c("##fileformat=VCFv4.2", "##source=BOHR::simulate_bohr (synthetic data)",
           sprintf("##contig=<ID=%s,length=%d>", names(chr_len), as.integer(chr_len)),
           "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
           paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", colnames(G)), collapse = "\t"))
  o <- order(match(V$chr, names(chr_len)), V$pos)
  body <- paste(V$chr[o], as.integer(V$pos[o]), V$id[o], V$ref[o], V$alt[o], ".", "PASS", ".", "GT",
                do.call(paste, c(as.data.frame(gt[o, , drop = FALSE], stringsAsFactors = FALSE), sep = "\t")), sep = "\t")
  writeLines(c(hdr, body), path)
  invisible(path)
}

#' @export
print.bohr_sim <- function(x, ...) {
  cat("BOHR synthetic input set\n")
  cat(sprintf("  %d individuals; %d tissues (%s); %d samples\n", nrow(x$traits), length(x$covariates),
              paste(sprintf("%s: %d", names(x$covariates), vapply(x$covariates, nrow, 0L)), collapse = ", "), nrow(x$samples)))
  ty <- tapply(x$tss$gene, x$tss$gene, length)
  cat(sprintf("  %d genes; %d promoters (%d tested, %d internal)\n", length(ty), nrow(x$tss), sum(x$tss$tested), sum(x$tss$internal)))
  cat(sprintf("  %d variants on %d chromosome(s); covariates %s\n", nrow(x$variants), length(unique(x$variants$chr)),
              paste(colnames(x$covariates[[1]]), collapse = ", ")))
  cat(sprintf("  true effects: %d rows (%d promoters)\n", nrow(x$truth), length(unique(x$truth$promoter))))
  if (!is.null(x$files)) cat("  written to", dirname(x$files[["relative"]]), "\n")
  invisible(x)
}
