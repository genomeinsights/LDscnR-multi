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
##   run_floor_sweep() (R/09_floor_sweep.R) and run_floor_profile()
##   (R/11_floor_profile.R) to keep each floor's unit-arm output distinct,
##   since that arm's p-values genuinely depend on size_floor
##   (ld_unit_matrix()'s columns change with it).
## @param permute Only affects the "unit" arm. FALSE skips
##   .emmax_perm_matrix() entirely (no p_perm.rds is written or listed in
##   this call's receipt outputs) -- still builds p_obs/p_obs_marker via a
##   single emmax_fast() call, i.e. a fast, floor-specific OBSERVED scan with
##   no permutation cost. Used by run_floor_profile()'s non-canonical
##   sensitivity floors, which need a per-floor observed test but never a
##   permutation-based one. The "simes" arm is deliberately NOT touched by
##   this flag -- see run_floor_profile()'s own header for why its shared,
##   floor-independent observed scan is read directly rather than routed
##   through a "permute" mode here (avoids entangling that arm's receipt
##   params, which are keyed on `simes_floor`, with an unrelated concept).
## @param simes_floor size_floor used ONLY to decide which markers are
##   "eligible" for the "simes" arm's PERMUTATION scan (see that arm's own
##   comment below) -- defaults to cfg$size_floor, the common case.
##   run_floor_sweep() overrides this to the LOWEST floor in its sweep: since
##   a lower floor only ever makes MORE clusters eligible, one compact scan
##   at the sweep's minimum floor is a superset covering every higher floor
##   tested afterward too, preserving the "simes is cheap to reuse across a
##   floor sweep" design this arm has always had. Part of the "simes" arm's
##   own receipt params (below), deliberately: a later call at a LOWER
##   simes_floor needs a STRICTLY LARGER eligible set than an earlier scan
##   covered, and without this the receipt would falsely read "up to date"
##   against the old, now-too-narrow compact scan -- silently starving
##   ld_outlier_perm()'s surrogate count for any newly-eligible unit's member
##   markers the old scan never covered.
## @return Named list, one element per computed statistic ("unit"/"simes"),
##   each list(p_obs, p_perm/p_perm_compact, out_dir) -- "unit"'s `p_perm` is
##   NULL when `permute = FALSE`.
run_stage_A <- function(dataset_dir, cfg = resolve_config(dataset_dir), force = FALSE,
                        engine_suffix = "", simes_floor = cfg$size_floor, permute = TRUE) {
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
                   grm_method = cfg$grm_method, b = cfg$b_unit, seed = cfg$seed, permute = permute)
    if (force || stage_stale(out_dir, stage_inputs, params, label = paste0("emmax_unit", engine_suffix))) {
      say("    [emmax_unit%s] ld_unit_matrix(repr = \"%s\", size_floor = %d) + emmax_fast%s\n",
          engine_suffix, cfg$unit_repr, cfg$size_floor,
          if (permute) sprintf(", %d permutations", cfg$b_unit) else " (observed only, no permutations)")
      um <- ld_unit_matrix(s1$genotypes, s1$stage1, s1$map, size_floor = cfg$size_floor,
                           repr = cfg$unit_repr)
      prep <- emmax_setup(um, GRM)
      p_obs <- emmax_fast(prep, d$phenotype)
      p_perm <- if (permute) .emmax_perm_matrix(prep, d$phenotype, perm_fun, cfg$b_unit, cfg$cores) else NULL
      p_obs_marker <- .marker_scan()$p_obs   ## display-only companion; see file header
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(p_obs, file.path(out_dir, "p_obs.rds"))
      perm_f <- file.path(out_dir, "p_perm.rds")
      ## Stale permutation file from an earlier permute = TRUE call at this
      ## SAME (floor-suffixed) out_dir: removed, not left orphaned next to a
      ## receipt that no longer promises it -- same defensive cleanup
      ## pattern as the "simes" arm's old_full removal just below.
      if (permute) saveRDS(p_perm, perm_f) else if (file.exists(perm_f)) file.remove(perm_f)
      saveRDS(p_obs_marker, file.path(out_dir, "p_obs_marker.rds"))
      saveRDS(params, file.path(out_dir, "params.rds"))
      write_receipt(out_dir, inputs = stage_inputs, params = params,
                    outputs = file.path(out_dir, c("p_obs.rds", "p_obs_marker.rds", if (permute) "p_perm.rds")))
    }
    perm_f <- file.path(out_dir, "p_perm.rds")
    out$unit <- list(p_obs = readRDS(file.path(out_dir, "p_obs.rds")),
                     p_perm = if (file.exists(perm_f)) readRDS(perm_f) else NULL,
                     p_obs_marker = readRDS(file.path(out_dir, "p_obs_marker.rds")), out_dir = out_dir)
  }

  if ("simes" %in% cfg$statistics) {
    out_dir <- file.path(dataset_dir, "output", "pvalues", "emmax_simes")
    ## The OBSERVED scan (p_obs) always covers every marker -- it is also
    ## marker_p/p_display's source (R/04_stage_B_outlier.R), used for
    ## genome-wide manhattan resolution regardless of which markers are ever
    ## actually tested. Only the PERMUTATION scan below is restricted: a
    ## unit's Simes p-value (LDscnR's .ld_outlier_tested_units()) is built
    ## strictly from its own member markers' p-values, and a marker only
    ## ever becomes a member once its cluster clears size_floor
    ## (.ld_outlier_units()) -- so a marker outside every size_floor-clearing
    ## cluster can NEVER contribute to any tested unit's surrogate p-value,
    ## observed or permuted. Rescanning it 200 times over was pure waste, and
    ## at full-genome scale it was the dominant cost: on the full 3sp panel
    ## (LDscnR-paper/module_3sp), only 16,890 of 790,578 markers (2.1%) sit
    ## in a cluster clearing size_floor = 8, yet emmax_setup()'s working
    ## object was being built over all 790,578 -- multiple GiB just for that
    ## one object, which a forked mclapply worker (cores > 1) then multiplies
    ## by copy-on-write, once per worker. Restricting emmax_setup() itself to
    ## the eligible markers shrinks that object by the same ~47x, not just
    ## the final result.
    cl <- data.table::as.data.table(s1$stage1$clusters)
    nl <- if ("n_loci" %in% names(cl)) cl$n_loci else cl$n_snps
    eligible <- unique(unlist(cl$members[nl >= simes_floor], use.names = FALSE))
    params <- list(grm_method = cfg$grm_method, b = cfg$b_simes, seed = cfg$seed,
                   simes_floor = simes_floor, n_eligible = length(eligible))
    if (force || stage_stale(out_dir, stage_inputs, params, label = "emmax_simes")) {
      say("    [emmax_simes] marker-level emmax_fast (observed, all %s markers); %d permutations restricted to %s markers (%.1f%%) in clusters clearing simes_floor = %d\n",
          format(ncol(s1$genotypes), big.mark = ","), cfg$b_simes, format(length(eligible), big.mark = ","),
          100 * length(eligible) / ncol(s1$genotypes), simes_floor)
      ms <- .marker_scan()
      p_obs <- ms$p_obs
      prep_compact <- emmax_setup(s1$genotypes[, eligible, drop = FALSE], GRM)
      p_perm_compact <- .emmax_perm_matrix(prep_compact, d$phenotype, perm_fun, cfg$b_simes, cfg$cores)
      rownames(p_perm_compact) <- eligible
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(p_obs, file.path(out_dir, "p_obs.rds"))
      saveRDS(p_perm_compact, file.path(out_dir, "p_perm_compact.rds"))
      ## Stale full-length format from a pre-compact-scan run: removed, not
      ## left behind -- run_stage_B() prefers p_perm_compact.rds when both
      ## exist, but a leftover full p_perm.rds is still bytes on disk this
      ## arm no longer needs and could confuse a reader inspecting the folder.
      old_full <- file.path(out_dir, "p_perm.rds")
      if (file.exists(old_full)) file.remove(old_full)
      saveRDS(params, file.path(out_dir, "params.rds"))
      write_receipt(out_dir, inputs = stage_inputs, params = params,
                    outputs = file.path(out_dir, c("p_obs.rds", "p_perm_compact.rds")))
    }
    out$simes <- list(p_obs = readRDS(file.path(out_dir, "p_obs.rds")),
                      p_perm_compact = readRDS(file.path(out_dir, "p_perm_compact.rds")), out_dir = out_dir)
  }

  out
}
