## =============================================================================
## LDscnR-multi/R/09_floor_sweep.R
##
## SIZE_FLOOR SENSITIVITY SWEEP. Whether this is cheap depends on which arm:
##
##   "simes": CHEAP. Stage A's marker-level EMMAX scan (R/03_stage_A_emmax.R)
##   never depends on size_floor -- only Stage B's per-cluster aggregation
##   does. run_stage_A() writes that arm to a size_floor-independent,
##   UNSUFFIXED path ("emmax_simes"), so calling it again at a different
##   floor just reports "up to date" and costs nothing; only Stage B
##   (ld_outlier_test's aggregation) reruns per floor.
##
##   "unit": NOT CHEAP. ld_unit_matrix()'s columns are themselves determined
##   by which Stage-1 clusters clear size_floor, so a different floor is a
##   genuinely different design matrix -- emmax_setup()/emmax_fast() and every
##   permutation must be redone. run_stage_A()'s `engine_suffix` (added for
##   this) keeps each floor's unit-arm output distinct
##   (output/pvalues/emmax_unit_floor<N>/) rather than overwriting the default
##   run, at the cost of B_UNIT permutations PER SWEPT FLOOR.
##
## Floors are deduplicated by their RESOLVED INTEGER VALUE, not by the
## requested factor: on a small panel (e.g. this repo's own worked examples,
## where cfg$size_floor derives to 1 -- see R/00_config.R's DERIVED
## SIZE_FLOOR) a 0.5x factor floors to the same value as 1x, and the two
## naturally collapse onto one shared computation rather than a wasted repeat.
## =============================================================================

## @param dataset_dir Path to one dataset folder.
## @param factors Multipliers applied to `cfg$size_floor` (default: half,
##   default, double). Each resolves to `max(1L, round(base_floor * factor))`.
## @param cfg Resolved config (resolve_config(dataset_dir) by default); its
##   OWN `size_floor` is the sweep's base/"default" point.
## @param annotation,chrom_lengths,force Passed through to every
##   run_stage_B() call, same meaning as run_dataset()'s.
## @param plot Write a combined output/figures/manhattan_floor_sweep.{png,pdf}
##   (one panel per (floor, statistic), via ldm_manhattan()).
## @return list(summary (data.table, one row per (floor, statistic): floor,
##   factor, n_units_tested, n_significant, n_regions), floors (the resolved
##   integer floors actually run), figure).
run_floor_sweep <- function(dataset_dir, factors = c(0.5, 1, 2), cfg = resolve_config(dataset_dir),
                            annotation = NULL, chrom_lengths = NULL, plot = TRUE, force = FALSE) {
  if (!file.exists(file.path(dataset_dir, "input", "phenotype.rds")))
    stop("run_floor_sweep() needs input/phenotype.rds (Stage A must be able to run).")

  base_floor <- cfg$size_floor
  floor_of <- vapply(factors, function(f) as.integer(max(1L, round(base_floor * f))), 0L)
  by_floor <- split(factors, floor_of)   # dedupe: factors that resolve to the same floor share one run
  floors <- as.integer(names(by_floor))
  say("=== %s: floor sweep -- base %d, factors %s -> floors %s ===\n",
      basename(dataset_dir), base_floor, paste(factors, collapse = ","), paste(floors, collapse = ","))

  s1 <- build_stage1(dataset_dir, cfg, force = force)   # shared across every floor
  all_engines <- character()
  rows <- vector("list", length(floors))

  for (i in seq_along(floors)) {
    fl <- floors[i]
    say("\n--- floor = %d (factor(s) %s) ---\n", fl, paste(by_floor[[i]], collapse = ", "))
    cfg_fl <- utils::modifyList(cfg, list(size_floor = fl))
    suffix <- sprintf("_floor%d", fl)

    stageA <- run_stage_A(dataset_dir, cfg_fl, force = force, engine_suffix = suffix)

    ## statistic = "unit" explicit: run_stage_B()'s own default guesses from
    ## `engine == "emmax_unit"` exactly, which this floor-suffixed name never
    ## matches.
    unit_engine <- paste0("emmax_unit", suffix)
    res_unit <- run_stage_B(dataset_dir, unit_engine, statistic = "unit", cfg = cfg_fl,
                            annotation = annotation, chrom_lengths = chrom_lengths, force = force)

    simes_engine <- paste0("emmax_simes", suffix)
    simes_dir <- file.path(dataset_dir, "output", "pvalues", "emmax_simes")
    res_simes <- run_stage_B(dataset_dir, simes_engine,
                             p_obs = readRDS(file.path(simes_dir, "p_obs.rds")),
                             p_perm = readRDS(file.path(simes_dir, "p_perm.rds")),
                             statistic = "simes", cfg = cfg_fl, annotation = annotation,
                             chrom_lengths = chrom_lengths, force = force)

    all_engines <- c(all_engines, unit_engine, simes_engine)
    rows[[i]] <- data.table::rbindlist(list(
      data.table::data.table(floor = fl, statistic = "unit", n_units_tested = nrow(res_unit$test$units),
                             n_significant = sum(res_unit$test$units$significant),
                             n_regions = nrow(res_unit$test$regions)),
      data.table::data.table(floor = fl, statistic = "simes", n_units_tested = nrow(res_simes$test$units),
                             n_significant = sum(res_simes$test$units$significant),
                             n_regions = nrow(res_simes$test$regions))))
  }

  summary <- data.table::rbindlist(rows)
  out_dir <- file.path(dataset_dir, "output")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(summary, file.path(out_dir, "floor_sweep_summary.csv"))
  say("\nwrote %s\n", file.path(out_dir, "floor_sweep_summary.csv"))

  fig <- NULL
  if (isTRUE(plot)) {
    fig <- ldm_manhattan(dataset_dir, engines = all_engines, alpha = cfg$alpha, ncol = 2)
    fig_dir <- file.path(out_dir, "figures")
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    h <- max(3, 3 * ceiling(length(all_engines) / 2))
    ggplot2::ggsave(file.path(fig_dir, "manhattan_floor_sweep.png"), fig, width = 20, height = h,
                    dpi = 150, limitsize = FALSE)
    ggplot2::ggsave(file.path(fig_dir, "manhattan_floor_sweep.pdf"), fig, width = 20, height = h,
                    limitsize = FALSE)
    say("wrote %s\n", file.path(fig_dir, "manhattan_floor_sweep.png"))
  }

  list(summary = summary, floors = floors, figure = fig)
}
