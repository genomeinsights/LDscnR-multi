## =============================================================================
## LDscnR-multi/R/03_stage_A_emmax.R
##
## STAGE A: GTs + map + phenotype -> EMMAX observed and permuted-null p-values.
##
## Two arms, following module_3sp's validated pattern (LDscnR-paper/module_3sp/
## R/03_EMMAX.R) rather than picking one: "unit" tests one association per
## stage-1 cluster (ld_unit_matrix()'s consensus-dosage variable, the cheaper
## and more powerful arm); "simes" tests every marker directly and aggregates
## per cluster via Simes' method inside ld_outlier_test(). Both write a
## p_obs/p_perm pair under output/pvalues/<engine>/, ready for Stage B.
##
## "unit" ALSO writes a marker-level companion, p_obs_marker.rds -- the same
## per-marker EMMAX scan (same GRM, same phenotype) "simes" computes for its
## own test, but here purely for DISPLAY resolution. The "unit" TEST itself
## has no per-marker p (p_obs is one value per cluster-summary-variable), so
## without this a per-SNP plot of the unit arm can only show every member of a
## cluster at its cluster's flat, broadcast value. module_3sp/R_figures/
## figure_manhattan.R makes exactly this choice for the same reason (its
## `pm_emmax` is plotted regardless of which statistic was actually tested).
## Cheap when "simes" is also requested: the scan is identical, so it is
## computed once and reused, not run twice.
##
## The GRM (built from stage1$pruned -- the stage-1 cluster representatives,
## not a separate ld_w threshold: module_3sp/R/00_config.R's own measurement
## found this "the easy call", see that file's GRM_BASIS section) lives in
## build_stage1() now, not here -- it's phenotype-free and
## check_structure_alignment() (R/08_structure_alignment.R) needs the same
## one. This file just reads `s1$GRM`.
## =============================================================================

## Run permutation index `b` through `perm_fun` and `emmax_fast`, seeded the
## same way module_3sp's p_perm closures are (`set.seed(bb)` per draw) so a
## rerun with the same B reproduces the same surrogate matrix.
.emmax_perm_matrix <- function(prep, y, perm_fun, B, cores = 1L) {
  one <- function(b) { set.seed(b); emmax_fast(prep, perm_fun(b, y)) }
  cols <- if (cores > 1L) parallel::mclapply(seq_len(B), one, mc.cores = cores) else lapply(seq_len(B), one)
  do.call(cbind, cols)
}

## @param dataset_dir Path to one dataset folder (must have input/phenotype.rds).
## @param cfg Resolved config; reads statistics, unit_repr, size_floor,
##   grm_method, b_unit, b_simes, seed, cores.
## @param force Rebuild every arm even if its receipt is current.
## @param engine_suffix Appended to the "unit" arm's output directory name
##   only (output/pvalues/emmax_unit<engine_suffix>/) -- used by
##   run_floor_sweep() (R/09_floor_sweep.R) to keep each swept size_floor's
##   unit-arm output distinct, since that arm's p-values genuinely depend on
##   size_floor (ld_unit_matrix()'s columns change with it). The "simes" arm
##   never depends on size_floor at this stage (only Stage B's aggregation
##   does), so it is deliberately NOT suffixed -- one shared computation is
##   reused, receipt-gated, across every floor in a sweep.
## @return Named list, one element per computed statistic ("unit"/"simes"),
##   each list(p_obs, p_perm, out_dir).
run_stage_A <- function(dataset_dir, cfg = resolve_config(dataset_dir), force = FALSE, engine_suffix = "") {
  s1 <- build_stage1(dataset_dir, cfg)
  d  <- read_dataset(dataset_dir, need_phenotype = TRUE)
  perm_fun <- resolve_perm_fun(dataset_dir)
  GRM <- s1$GRM

  out <- list()
  ## Depends on cache/'s own receipt, not just the raw input files: without
  ## this, changing a decay_args/cr_rho setting rebuilds stage1 (build_stage1()
  ## catches that on its own params) but this stage's receipt saw no change in
  ## genotypes/map/phenotype and reported "up to date", silently keeping
  ## association p-values computed against the OLD stage1/GRM. Found by an
  ## external audit; run_stage_B() already had this dependency, this stage
  ## didn't.
  stage_inputs <- c(file.path(dataset_dir, "input", "genotypes.rds"),
                    file.path(dataset_dir, "input", "map.rds"),
                    file.path(dataset_dir, "input", "phenotype.rds"),
                    receipt_path(file.path(dataset_dir, "cache")))
  perm_fun_file <- file.path(dataset_dir, "input", "perm_fun.R")
  if (file.exists(perm_fun_file)) stage_inputs <- c(stage_inputs, perm_fun_file)

  ## Memoised marker-level scan: built the first time either arm needs it,
  ## reused by the second so "unit" and "simes" never pay for it twice.
  .marker_scan <- local({
    cache <- NULL
    function() {
      if (is.null(cache)) {
        prep <- emmax_setup(s1$genotypes, GRM)
        p_obs <- emmax_fast(prep, d$phenotype)
        names(p_obs) <- colnames(s1$genotypes)
        cache <<- list(prep = prep, p_obs = p_obs)
      }
      cache
    }
  })

  if ("unit" %in% cfg$statistics) {
    out_dir <- file.path(dataset_dir, "output", "pvalues", paste0("emmax_unit", engine_suffix))
    params <- list(cr_rho = cfg$cr_rho, size_floor = cfg$size_floor, unit_repr = cfg$unit_repr,
                   grm_method = cfg$grm_method, b = cfg$b_unit, seed = cfg$seed)
    if (force || stage_stale(out_dir, stage_inputs, params, label = paste0("emmax_unit", engine_suffix))) {
      say("    [emmax_unit%s] ld_unit_matrix(repr = \"%s\", size_floor = %d) + emmax_fast, %d permutations\n",
          engine_suffix, cfg$unit_repr, cfg$size_floor, cfg$b_unit)
      um <- ld_unit_matrix(s1$genotypes, s1$stage1, s1$map, size_floor = cfg$size_floor,
                           repr = cfg$unit_repr)
      prep <- emmax_setup(um, GRM)
      p_obs <- emmax_fast(prep, d$phenotype)
      p_perm <- .emmax_perm_matrix(prep, d$phenotype, perm_fun, cfg$b_unit, cfg$cores)
      p_obs_marker <- .marker_scan()$p_obs   ## display-only companion; see file header
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(p_obs, file.path(out_dir, "p_obs.rds"))
      saveRDS(p_perm, file.path(out_dir, "p_perm.rds"))
      saveRDS(p_obs_marker, file.path(out_dir, "p_obs_marker.rds"))
      saveRDS(params, file.path(out_dir, "params.rds"))
      write_receipt(out_dir, inputs = stage_inputs, params = params,
                    outputs = file.path(out_dir, c("p_obs.rds", "p_perm.rds", "p_obs_marker.rds")))
    }
    out$unit <- list(p_obs = readRDS(file.path(out_dir, "p_obs.rds")),
                     p_perm = readRDS(file.path(out_dir, "p_perm.rds")),
                     p_obs_marker = readRDS(file.path(out_dir, "p_obs_marker.rds")), out_dir = out_dir)
  }

  if ("simes" %in% cfg$statistics) {
    out_dir <- file.path(dataset_dir, "output", "pvalues", "emmax_simes")
    params <- list(grm_method = cfg$grm_method, b = cfg$b_simes, seed = cfg$seed)
    if (force || stage_stale(out_dir, stage_inputs, params, label = "emmax_simes")) {
      say("    [emmax_simes] marker-level emmax_fast, %d permutations (rescans every marker per draw)\n",
          cfg$b_simes)
      ms <- .marker_scan()
      p_obs <- ms$p_obs
      p_perm <- .emmax_perm_matrix(ms$prep, d$phenotype, perm_fun, cfg$b_simes, cfg$cores)
      rownames(p_perm) <- colnames(s1$genotypes)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(p_obs, file.path(out_dir, "p_obs.rds"))
      saveRDS(p_perm, file.path(out_dir, "p_perm.rds"))
      saveRDS(params, file.path(out_dir, "params.rds"))
      write_receipt(out_dir, inputs = stage_inputs, params = params,
                    outputs = file.path(out_dir, c("p_obs.rds", "p_perm.rds")))
    }
    out$simes <- list(p_obs = readRDS(file.path(out_dir, "p_obs.rds")),
                      p_perm = readRDS(file.path(out_dir, "p_perm.rds")), out_dir = out_dir)
  }

  out
}
