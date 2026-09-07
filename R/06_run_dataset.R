## =============================================================================
## LDscnR-multi/R/06_run_dataset.R
##
## ORCHESTRATE ONE DATASET FOLDER: stage1 clustering -> Stage A (only if
## input/phenotype.rds exists) -> Stage B for every available p-value engine
## (Stage A's own arms plus every external_pvalues/*.rds) -> the manhattan
## figure. Every step checks its own receipt first (see R/00_config.R), so
## calling this again on an unchanged dataset is a no-op walk over receipts.
## =============================================================================

## @param dataset_dir Path to one dataset folder.
## @param engines Optional character vector restricting which engines Stage B
##   runs on (default: every engine available -- Stage A's arms plus every
##   external_pvalues/*.rds).
## @param cfg Resolved config (resolve_config(dataset_dir) by default).
## @param annotation,chrom_lengths Optional, passed through to every
##   run_stage_B() call (ld_region_rotation()).
## @param plot Write output/figures/manhattan.{png,pdf} via ldm_manhattan().
## @param force Rebuild every stage even if its receipt is current.
## @return list(stage1, stageA, results (named by engine), figure).
run_dataset <- function(dataset_dir, engines = NULL, cfg = resolve_config(dataset_dir),
                        annotation = NULL, chrom_lengths = NULL, plot = TRUE, force = FALSE) {
  say("=== %s ===\n", basename(dataset_dir))
  invisible(check_ldscnr(stop_on_fail = FALSE))

  s1 <- build_stage1(dataset_dir, cfg, force = force)

  stageA <- NULL
  if (file.exists(file.path(dataset_dir, "input", "phenotype.rds"))) {
    stageA <- run_stage_A(dataset_dir, cfg, force = force)
  } else {
    say("    [stage A] skipped -- no input/phenotype.rds\n")
  }

  pv_dir <- file.path(dataset_dir, "output", "pvalues")
  emmax_engines <- if (dir.exists(pv_dir)) list.dirs(pv_dir, recursive = FALSE, full.names = FALSE) else character()
  ext_engines <- names(list_external_pvalues(dataset_dir, s1$map))
  all_engines <- unique(c(emmax_engines, ext_engines))
  if (!is.null(engines)) all_engines <- intersect(all_engines, engines)

  if (!length(all_engines)) {
    warning("run_dataset(", basename(dataset_dir), "): no p-value engine available -- ",
           "no input/phenotype.rds for Stage A, and no external_pvalues/*.rds.")
    return(list(stage1 = s1, stageA = stageA, results = list(), figure = NULL))
  }

  results <- stats::setNames(lapply(all_engines, function(eng)
    run_stage_B(dataset_dir, eng, cfg = cfg, annotation = annotation,
               chrom_lengths = chrom_lengths, force = force)), all_engines)

  fig <- NULL
  if (isTRUE(plot)) {
    fig <- ldm_manhattan(dataset_dir, engines = all_engines, alpha = cfg$alpha)
    fig_dir <- file.path(dataset_dir, "output", "figures")
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    h <- max(3, 3 * length(all_engines))
    ggplot2::ggsave(file.path(fig_dir, "manhattan.png"), fig, width = 12, height = h,
                    dpi = 150, limitsize = FALSE)
    ggplot2::ggsave(file.path(fig_dir, "manhattan.pdf"), fig, width = 12, height = h, limitsize = FALSE)
    say("    [figure] %s\n", file.path(fig_dir, "manhattan.png"))
  }

  list(stage1 = s1, stageA = stageA, results = results, figure = fig)
}
