## =============================================================================
## LDscnR-multi/R/08_structure_alignment.R
##
## PRE-SCAN DIAGNOSTIC: how much of the tested phenotype is explained by
## genetic/sampling structure, before any association test runs. Generalises
## LDscnR-paper/module_9sp/R/02b_env_structure_alignment.R (confirmed, via
## LDscnR_manuscript/generate_empirical_values.R, as the actual source of the
## manuscript's empirical alignment numbers) from that script's hardcoded
## 3sp/9sp bundle layout into this repo's dataset-folder convention.
##
## Per materials_and_methods.tex's "Environment--structure alignment" and
## Discussion.tex: "Environment-structure alignment is therefore useful before
## a genome scan begins... it reveals when phenotype and ancestry provide
## nearly interchangeable explanations." A WARNING DIAGNOSTIC, not a gate --
## this function never blocks or fails a dataset over its own result, it only
## reports.
##
## Three measures, all following the reference script exactly:
##   grm_axes_r2_population           R^2 of population-mean phenotype ~ the
##                                     top `n_axes` eigenvectors of the
##                                     population-level GRM
##   structure_r2_population          R^2 of population-mean phenotype ~
##                                     structure_group, population-weighted
##   structure_r2_individual_weighted same, individual-weighted
##
## against two nulls:
##   "within-group permutation"  population-mean phenotype shuffled within
##                                structure_group (matches the manuscript's
##                                primary null exactly)
##   "GRM-matched continuous"    a continuous surrogate drawn from the GRM's
##                                own eigenbasis (s = V %*% (sqrt(lambda) *
##                                rnorm(n))), residualised against the tested
##                                phenotype -- "does alignment this strong
##                                arise just from sharing the GRM's own
##                                covariance structure?" Computed
##                                unconditionally (cheap, needs only GRM +
##                                phenotype), not gated to any one panel the
##                                way the reference script's 9sp-only version
##                                is.
##
## NO COVARIATE ADJUSTMENT: LDscnR-multi has no Covar/residualisation
## mechanism (see R/03_stage_A_emmax.R), so this reports alignment for the
## phenotype actually tested (input/phenotype.rds) only -- the reference
## script's "raw habitat" line, not its 9sp-specific "lineage-adjusted"
## variant.
## =============================================================================

.r2 <- function(y, X) { if (!length(y) || stats::var(y) == 0) return(NA_real_)
  summary(stats::lm(y ~ X))$r.squared }
.adj_r2 <- function(y, X) { if (!length(y) || stats::var(y) == 0) return(NA_real_)
  summary(stats::lm(y ~ X))$adj.r.squared }
.group_r2 <- function(y, group) {
  if (!length(y) || stats::var(y) == 0 || length(unique(group)) < 2L) return(NA_real_)
  summary(stats::lm(y ~ factor(group)))$r.squared
}
.pop_mean <- function(y, pop, pops) as.numeric(tapply(y, pop, mean)[pops])

## Read an individual-level named input (same contract as phenotype.rds --
## see R/01_dataset_io.R), realigned to `ind_ids` order, erroring if any
## individual is missing.
.read_individual_input <- function(f, ind_ids) {
  x <- readRDS(f)
  if (is.null(names(x))) stop(f, ": must be a NAMED vector (names = genotypes' rownames).")
  if (!identical(names(x), ind_ids)) x <- x[ind_ids]
  if (anyNA(names(x))) stop(f, ": does not cover every individual in genotypes.")
  x
}

## Population-level GRM projection: block-average the individual-level GRM
## into a population x population matrix (Ind = 1/n_p indicator), eigen-
## decompose, return the top `n_axes` eigenvectors (population order = `pops`).
.grm_population_axes <- function(GRM, pop, pops, n_axes) {
  n_pop <- length(pops)
  pop_idx <- match(pop, pops)
  Ind <- matrix(0, length(pop), n_pop)
  Ind[cbind(seq_along(pop), pop_idx)] <- 1
  Ind <- sweep(Ind, 2, colSums(Ind), "/")
  K_pop <- crossprod(Ind, GRM %*% Ind)
  K_pop <- (K_pop + t(K_pop)) / 2
  eg <- eigen(K_pop, symmetric = TRUE)
  eg$vectors[, seq_len(min(n_axes, n_pop - 1L)), drop = FALSE]
}

.observed_row <- function(y_ind, pop, pops, group, axes) {
  y_pop <- .pop_mean(y_ind, pop, pops)
  group_pop <- group[match(pops, pop)]
  data.table::data.table(
    n_individuals = length(y_ind), n_populations = length(pops),
    n_structure_groups = data.table::uniqueN(group_pop),
    grm_axes_r2_population = .r2(y_pop, axes),
    grm_axes_adjusted_r2_population = .adj_r2(y_pop, axes),
    structure_r2_population = .group_r2(y_pop, group_pop),
    structure_r2_individual_weighted = .group_r2(y_ind, group))
}

.null_row <- function(scheme, draw, y_ind, pop, pops, group, axes) {
  y_pop <- .pop_mean(y_ind, pop, pops)
  group_pop <- group[match(pops, pop)]
  data.table::data.table(
    scheme = scheme, draw = draw,
    grm_axes_r2_population = .r2(y_pop, axes),
    structure_r2_population = .group_r2(y_pop, group_pop),
    structure_r2_individual_weighted = .group_r2(y_ind, group))
}

.summarise_metric <- function(null_draws, observed, metric) {
  null_draws[, {
    obs <- observed[[metric]]
    x <- get(metric)
    q <- stats::quantile(x, c(0.025, 0.975), na.rm = TRUE)
    list(observed = obs, null_mean = mean(x, na.rm = TRUE), null_median = stats::median(x, na.rm = TRUE),
        null_sd = stats::sd(x, na.rm = TRUE), null_q025 = q[1], null_q975 = q[2],
        observed_within_null95 = obs >= q[1] - 1e-12 && obs <= q[2] + 1e-12,
        proportion_null_ge_observed = mean(x >= obs - 1e-12, na.rm = TRUE))
  }, by = scheme][, metric := metric]
}

## @param dataset_dir Path to one dataset folder. Needs input/population.rds
##   (required) and input/phenotype.rds; input/structure_group.rds is
##   optional (defaults to population.rds itself -- the within-group null
##   then has nothing to shuffle and degenerates to the observed value,
##   which is reported, not an error).
## @param cfg Resolved config; reads grm_method (via build_stage1()), seed.
## @param n_draws Null draws per scheme (default 1000L, matching the
##   reference script).
## @param n_axes Population-level GRM eigenvectors to use (default 5L,
##   capped at n_populations - 1).
## @param force Rebuild even if the receipt says nothing changed.
## @return list(observed, null_draws, null_summary, out_dir).
check_structure_alignment <- function(dataset_dir, cfg = resolve_config(dataset_dir),
                                      n_draws = 1000L, n_axes = 5L, force = FALSE) {
  pop_f <- file.path(dataset_dir, "input", "population.rds")
  if (!file.exists(pop_f))
    stop("check_structure_alignment() needs ", pop_f, " (named vector, individual -> population ID).")
  grp_f <- file.path(dataset_dir, "input", "structure_group.rds")

  s1 <- build_stage1(dataset_dir, cfg)
  d  <- read_dataset(dataset_dir, need_phenotype = TRUE)
  ind_ids <- rownames(d$genotypes)

  population <- as.character(.read_individual_input(pop_f, ind_ids))
  structure_group <- if (file.exists(grp_f)) as.character(.read_individual_input(grp_f, ind_ids))
                     else population
  bad <- tapply(structure_group, population, function(g) length(unique(g)) > 1L)
  if (any(bad)) stop("structure_group.rds is not constant within population for: ",
                     paste(names(bad)[bad], collapse = ", "))

  out_dir <- file.path(dataset_dir, "output", "alignment")
  inputs <- c(file.path(dataset_dir, "input", c("genotypes.rds", "map.rds", "phenotype.rds", "population.rds")),
             if (file.exists(grp_f)) grp_f, receipt_path(file.path(dataset_dir, "cache")))
  params <- list(n_draws = n_draws, n_axes = n_axes, seed = cfg$seed)

  if (!force && !stage_stale(out_dir, inputs, params, label = "alignment")) {
    return(list(observed = data.table::fread(file.path(out_dir, "observed.csv")),
               null_draws = data.table::fread(file.path(out_dir, "null_draws.csv")),
               null_summary = data.table::fread(file.path(out_dir, "null_summary.csv")),
               out_dir = out_dir))
  }

  pops <- unique(population)
  say("    [alignment] %d individuals, %d populations, %d structure group(s)\n",
      length(ind_ids), length(pops), length(unique(structure_group)))
  axes <- .grm_population_axes(s1$GRM, population, pops, n_axes)
  y <- d$phenotype

  observed <- .observed_row(y, population, pops, structure_group, axes)

  say("    [alignment] null 1/2: %d within-group population-mean permutations\n", n_draws)
  pop_group <- structure_group[match(pops, population)]
  y_pop_obs <- .pop_mean(y, population, pops)
  null_perm <- data.table::rbindlist(lapply(seq_len(n_draws), function(b) {
    set.seed(b)
    y_pop_shuf <- stats::ave(y_pop_obs, pop_group, FUN = sample)
    y_ind_shuf <- y_pop_shuf[match(population, pops)]
    .null_row("within-group permutation", b, y_ind_shuf, population, pops, structure_group, axes)
  }))

  say("    [alignment] null 2/2: %d GRM-matched continuous surrogates\n", n_draws)
  eg <- eigen(s1$GRM, symmetric = TRUE)
  values <- pmax(eg$values, 0); vectors <- eg$vectors
  null_grm <- data.table::rbindlist(lapply(seq_len(n_draws), function(b) {
    set.seed(b)
    s <- as.numeric(vectors %*% (sqrt(values) * stats::rnorm(length(y))))
    y_null <- as.numeric(stats::resid(stats::lm(s ~ y)))
    .null_row("GRM-matched continuous", b, y_null, population, pops, structure_group, axes)
  }))

  null_draws <- data.table::rbindlist(list(null_perm, null_grm))
  null_summary <- data.table::rbindlist(lapply(
    c("grm_axes_r2_population", "structure_r2_population", "structure_r2_individual_weighted"),
    function(m) .summarise_metric(null_draws, observed, m)))
  data.table::setcolorder(null_summary, c("scheme", "metric", "observed", "null_mean", "null_median",
                                         "null_sd", "null_q025", "null_q975",
                                         "observed_within_null95", "proportion_null_ge_observed"))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(observed, file.path(out_dir, "observed.csv"))
  data.table::fwrite(null_draws, file.path(out_dir, "null_draws.csv"))
  data.table::fwrite(null_summary, file.path(out_dir, "null_summary.csv"))
  write_receipt(out_dir, inputs = inputs, params = params,
               outputs = file.path(out_dir, c("observed.csv", "null_draws.csv", "null_summary.csv")))

  list(observed = observed, null_draws = null_draws, null_summary = null_summary, out_dir = out_dir)
}
