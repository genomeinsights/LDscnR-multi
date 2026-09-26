## =============================================================================
## LDscnR-multi/R/06_run_dataset.R
##
## ORCHESTRATE ONE DATASET FOLDER: stage1 clustering -> the structure-
## alignment diagnostic (only if input/population.rds exists) -> Stage A
## (only if input/phenotype.rds exists) -> Stage B for every available
## p-value engine (Stage A's own arms plus every external_pvalues/*.rds) ->
## two manhattan figures (q-value, and the same panels with ld_w_095 --
## local-LD support, independent of any test -- on the y-axis instead).
## Every step checks its own receipt first (see R/00_config.R), so calling
## this again on an unchanged dataset is a no-op walk over receipts.
## =============================================================================

## @param dataset_dir Path to one dataset folder.
## @param engines Optional character vector restricting which engines Stage B
##   runs on (default: every engine available -- Stage A's arms plus every
##   external_pvalues/*.rds).
## @param cfg Resolved config (resolve_config(dataset_dir) by default).
## @param annotation,chrom_lengths Optional, passed through to every
##   run_stage_B() call (ld_region_rotation()).
## @param plot Write output/figures/manhattan.{png,pdf} (q-value) and
##   output/figures/manhattan_ld_w.{png,pdf} (same panels/faceting/colouring,
##   ld_w_095 on the y-axis instead) via ldm_manhattan().
## @param force Rebuild every stage even if its receipt is current.
## @return list(stage1, alignment, stageA, results (named by engine), figure,
##   figure_ld_w).
run_dataset <- function(dataset_dir, engines = NULL, cfg = resolve_config(dataset_dir),
                        annotation = NULL, chrom_lengths = NULL, plot = TRUE, force = FALSE) {
  say("=== %s ===\n", basename(dataset_dir))
  ## check_ldscnr()'s own default (stop_on_fail = TRUE unless LDSCNR_LAX is
  ## set) is deliberately NOT overridden here: a batch of hundreds of
  ## datasets should stop on the FIRST call against a mismatched/dirty
  ## LDscnR install, not print a warning and keep running through all of
  ## them against results that may not be reproducible.
  invisible(check_ldscnr())

  s1 <- build_stage1(dataset_dir, cfg, force = force)

  alignment <- NULL
  has_pop <- file.exists(file.path(dataset_dir, "input", "population.rds"))
  has_pheno <- file.exists(file.path(dataset_dir, "input", "phenotype.rds"))
  if (has_pop && has_pheno) {
    alignment <- check_structure_alignment(dataset_dir, cfg, force = force)
  } else if (has_pop) {
    say("    [alignment] skipped -- population.rds present but no input/phenotype.rds\n")
  }

  stageA <- NULL
  if (has_pheno) {
    stageA <- run_stage_A(dataset_dir, cfg, force = force)
  } else {
    say("    [stage A] skipped -- no input/phenotype.rds\n")
  }

  ## Exactly the two canonical Stage A engines, not every subdirectory under
  ## output/pvalues/ -- a run_floor_sweep() call (R/09_floor_sweep.R) leaves
  ## floor-suffixed siblings there (e.g. "emmax_unit_floor4") that are a
  ## separate, purpose-built analysis with their own explicit `statistic`;
  ## picking them up here would run them through run_stage_B()'s name-based
  ## default guess ("emmax_unit_floor4" != "emmax_unit" -> guessed "simes"),
  ## which is wrong for a unit-level p_obs and errors on the length mismatch.
  pv_dir <- file.path(dataset_dir, "output", "pvalues")
  emmax_engines <- intersect(c("emmax_unit", "emmax_simes"),
                             if (dir.exists(pv_dir)) list.dirs(pv_dir, recursive = FALSE, full.names = FALSE)
                             else character())
  ext_engines <- names(list_external_pvalues(dataset_dir, s1$map))
  all_engines <- unique(c(emmax_engines, ext_engines))
  if (!is.null(engines)) all_engines <- intersect(all_engines, engines)

  if (!length(all_engines)) {
    warning("run_dataset(", basename(dataset_dir), "): no p-value engine available -- ",
           "no input/phenotype.rds for Stage A, and no external_pvalues/*.rds.")
    return(list(stage1 = s1, alignment = alignment, stageA = stageA, results = list(), figure = NULL))
  }

  results <- stats::setNames(lapply(all_engines, function(eng)
    run_stage_B(dataset_dir, eng, cfg = cfg, annotation = annotation,
               chrom_lengths = chrom_lengths, force = force)), all_engines)

  fig <- fig_ld_w <- NULL
  if (isTRUE(plot)) {
    fig_dir <- file.path(dataset_dir, "output", "figures")
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    h <- max(3, 3 * length(all_engines))

    fig <- ldm_manhattan(dataset_dir, engines = all_engines, alpha = cfg$alpha)
    ggplot2::ggsave(file.path(fig_dir, "manhattan.png"), fig, width = 12, height = h,
                    dpi = 150, limitsize = FALSE)
    ggplot2::ggsave(file.path(fig_dir, "manhattan.pdf"), fig, width = 12, height = h, limitsize = FALSE)
    say("    [figure] %s\n", file.path(fig_dir, "manhattan.png"))

    ## Same panels/faceting/colouring as the q-value plot above, just
    ## ld_w_095 (local-LD support, no association test involved) on the
    ## y-axis -- lets a significant region be compared directly against
    ## local LD, independent of any engine's own significance.
    fig_ld_w <- ldm_manhattan(dataset_dir, engines = all_engines, value = "ld_w_095")
    ggplot2::ggsave(file.path(fig_dir, "manhattan_ld_w.png"), fig_ld_w, width = 12, height = h,
                    dpi = 150, limitsize = FALSE)
    ggplot2::ggsave(file.path(fig_dir, "manhattan_ld_w.pdf"), fig_ld_w, width = 12, height = h, limitsize = FALSE)
    say("    [figure] %s\n", file.path(fig_dir, "manhattan_ld_w.png"))
  }

  list(stage1 = s1, alignment = alignment, stageA = stageA, results = results,
      figure = fig, figure_ld_w = fig_ld_w)
}
