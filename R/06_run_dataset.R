## =============================================================================
## LDscnR-multi/R/06_run_dataset.R
##
## ORCHESTRATE ONE DATASET FOLDER: stage1 clustering -> the LD-decay report
## (decay_summary.csv/decay_print.txt/decay_curves.pdf, see
## R/10_ld_decay_report.R -- always runs, genotype-only like stage1 itself)
## -> the structure-alignment diagnostic (only if input/population.rds
## exists) -> the floor profile (R/11_floor_profile.R -- only if
## input/phenotype.rds exists; see @param floor_profile to disable) -> Stage
## A/B for every available p-value engine (the floor profile's own canonical
## floor for "emmax_unit"/"emmax_simes", plus every external_pvalues/*.rds at
## that same floor) -> two manhattan figures (q-value, and the same panels
## with ld_w_095 -- local-LD support, independent of any test -- on the
## y-axis instead), labelled with each region's floor-stability tier where
## the floor profile found one.
## Every step checks its own receipt first (see R/00_config.R), so calling
## this again on an unchanged dataset is a no-op walk over receipts.
##
## WHY THE FLOOR PROFILE'S CANONICAL FLOOR IS NOW THIS FUNCTION'S OWN
## `size_floor`, NOT `cfg`'s own derivation: run_floor_profile()
## (R/11_floor_profile.R) picks size_floor by an actual, verified property --
## the test-count reduction its own Stage-1 clustering achieves -- while
## resolve_config()'s own default is a marker-count heuristic with no such
## grounding. Once the floor profile runs anyway (cheap: the two side floors
## are observed-only, no permutations -- see that file's header), reporting
## results from a DIFFERENT, unrelated floor alongside it would be confusing,
## not merely redundant. `cfg$size_floor` (whatever resolve_config() gave) is
## overridden to the floor profile's own canonical floor before Stage A/B run
## -- resolve_config() and run_floor_sweep() are both untouched; only THIS
## function's own default reporting floor changes when floor_profile = TRUE.
## =============================================================================

## @param dataset_dir Path to one dataset folder.
## @param engines Optional character vector restricting which engines Stage B
##   runs on (default: every engine available -- Stage A's arms plus every
##   external_pvalues/*.rds).
## @param cfg Resolved config (resolve_config(dataset_dir) by default).
## @param annotation,chrom_lengths Optional, passed through to every
##   run_stage_B() call (ld_region_rotation()) and to the floor profile's own
##   canonical-floor call.
## @param floor_profile TRUE (default) runs run_floor_profile() and adopts
##   its canonical (99.7%-target) floor as `cfg$size_floor` for this call's
##   own Stage A/B and manhattan figures -- see the file header. FALSE keeps
##   `cfg$size_floor` exactly as resolve_config() gave it (today's
##   behaviour), and skips the floor profile entirely.
## @param floor_profile_targets,floor_profile_canonical_target Passed through
##   to run_floor_profile() when `floor_profile = TRUE`.
## @param plot Write output/figures/manhattan.{png,pdf} (q-value) and
##   output/figures/manhattan_ld_w.{png,pdf} (same panels/layout/colouring,
##   ld_w_095 on the y-axis instead) via ldm_manhattan().
## @param force Rebuild every stage even if its receipt is current.
## @return list(stage1, ld_decay, alignment, floor_profile, stageA, results
##   (named by engine), figure, figure_ld_w).
run_dataset <- function(dataset_dir, engines = NULL, cfg = resolve_config(dataset_dir),
                        annotation = NULL, chrom_lengths = NULL,
                        floor_profile = TRUE, floor_profile_targets = c(0.995, 0.997, 0.999),
                        floor_profile_canonical_target = 0.997,
                        plot = TRUE, force = FALSE) {
  say("=== %s ===\n", basename(dataset_dir))
  ## check_ldscnr()'s own default (stop_on_fail = TRUE unless LDSCNR_LAX is
  ## set) is deliberately NOT overridden here: a batch of hundreds of
  ## datasets should stop on the FIRST call against a mismatched/dirty
  ## LDscnR install, not print a warning and keep running through all of
  ## them against results that may not be reproducible.
  invisible(check_ldscnr())

  s1 <- build_stage1(dataset_dir, cfg, force = force)
  ## force = FALSE, always, not `force`: report_ld_decay() calls build_stage1()
  ## again internally (same as check_structure_alignment()/run_stage_A() do --
  ## see R/02_stage1_cluster.R's header on why every stage rederives it rather
  ## than threading `s1` through call signatures), but stage1's receipt was
  ## JUST written by the call above, so a second force = TRUE here would
  ## force a full, genuinely expensive stage1 rebuild a SECOND time for
  ## nothing -- caught via full-genome timing (this step read as 386s,
  ## against structure alignment's 12s, on a run where force = TRUE).
  ## force = FALSE always hits the fast cache-read path instead, correctly,
  ## since nothing has changed since the line above just wrote it.
  ld_decay <- report_ld_decay(dataset_dir, cfg, force = FALSE)

  alignment <- NULL
  has_pop <- file.exists(file.path(dataset_dir, "input", "population.rds"))
  has_pheno <- file.exists(file.path(dataset_dir, "input", "phenotype.rds"))
  if (has_pop && has_pheno) {
    alignment <- check_structure_alignment(dataset_dir, cfg, force = force)
  } else if (has_pop) {
    say("    [alignment] skipped -- population.rds present but no input/phenotype.rds\n")
  }

  fp <- NULL
  stageA <- NULL
  if (has_pheno) {
    if (isTRUE(floor_profile)) {
      fp <- run_floor_profile(dataset_dir, targets = floor_profile_targets,
                              canonical_target = floor_profile_canonical_target, cfg = cfg,
                              annotation = annotation, chrom_lengths = chrom_lengths,
                              canonical_unsuffixed = TRUE, force = force)
      ## Every call below in this function now sees the floor profile's own
      ## canonical floor, not resolve_config()'s -- see file header. The
      ## run_stage_A()/run_stage_B() calls that follow just hit the receipt
      ## the floor profile's own canonical-floor pass (canonical_unsuffixed =
      ## TRUE, above) already wrote at the standard unsuffixed paths; they are
      ## not redundant work, just this function's own return value/API shape.
      cfg <- utils::modifyList(cfg, list(size_floor = fp$canonical_floor))
      say("    [floor profile] canonical size_floor = %d (%.1f%% target)\n",
          fp$canonical_floor, 100 * floor_profile_canonical_target)
    }
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
    return(list(stage1 = s1, ld_decay = ld_decay, alignment = alignment, floor_profile = fp,
               stageA = stageA, results = list(), figure = NULL))
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

    ## Same panels/layout/colouring as the q-value plot above, just
    ## ld_w_095 (local-LD support, no association test involved) on the
    ## y-axis -- lets a significant region be compared directly against
    ## local LD, independent of any engine's own significance.
    fig_ld_w <- ldm_manhattan(dataset_dir, engines = all_engines, value = "ld_w_095")
    ggplot2::ggsave(file.path(fig_dir, "manhattan_ld_w.png"), fig_ld_w, width = 12, height = h,
                    dpi = 150, limitsize = FALSE)
    ggplot2::ggsave(file.path(fig_dir, "manhattan_ld_w.pdf"), fig_ld_w, width = 12, height = h, limitsize = FALSE)
    say("    [figure] %s\n", file.path(fig_dir, "manhattan_ld_w.png"))
  }

  list(stage1 = s1, ld_decay = ld_decay, alignment = alignment, floor_profile = fp,
      stageA = stageA, results = results, figure = fig, figure_ld_w = fig_ld_w)
}
