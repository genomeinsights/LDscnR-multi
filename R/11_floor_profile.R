## =============================================================================
## LDscnR-multi/R/11_floor_profile.R
##
## FLOOR-STABILITY LABELS FOR THE REPORTED REGIONS -- deliberately separate
## from run_floor_sweep() (R/09_floor_sweep.R), which reruns full
## permutations at every floor and is not suitable for hundreds of datasets.
## run_floor_profile() runs the pipeline's own, unmodified Stage A/B once with
## permutations at the REPORTED floor (cfg$size_floor), and observed-only
## passes (no permutations, no rotation) at two comparison floors. It never
## changes the reported floor.
##
## THE THREE FLOORS (PK, 2026-10-07): the reported floor is the dataset's own
## size_floor -- by default the marker-density rule, round(markers per Mb /
## 250), see resolve_config() -- and it is compared with HALF and TWICE that
## floor (`factors`, default c(0.5, 2)). The low side may be 1 (singletons
## included): it is a comparison floor only, never the reported one. This
## replaces the earlier 99.5/99.7/99.9% test-reduction targets; on the LDscnR
## paper's simulations the 0.5x/2x design gave reported-region precision of
## ~0.2 (1/3), ~0.5 (2/3) and ~0.8 (3/3), and recovery at the LOWER floor was
## the more informative check. grid.csv still records each floor's achieved
## test-count reduction and marker coverage, as information only.
##
## COST: the reported floor runs Stage A/B exactly as run_dataset() always
## has (permutations, rotation). The consensus ("unit") arm is fitted ONCE
## more, observed only, at the low floor; the high floor's units are a
## subset of the low floor's and each unit's EMMAX p-value does not depend on
## which other units are tested, so the high floor reuses those p-values
## (selected by core_snp) and only BH + Stage-2 assembly are redone. The
## "simes" arm reuses the shared, floor-independent marker scan at every
## floor. Neither side floor gets a permutation null.
##
## RECOVERY RULE: a reported region is recovered at a comparison floor when
## at least one of its significant Stage-1 units is also significant there
## AND the region overlaps an outlier region assembled at that floor (the
## manuscript's rule, as in LDscnR-paper/module_3sp/R/helpers_floor_tiers.R).
## Units are matched by `core_snp`, each unit's representative marker, which
## does not depend on the floor. `unit_id` cannot be used, because it is a
## position index into one floor's own size-filtered unit table.
## =============================================================================

## Achieved test-count reduction for one candidate integer floor.
.floor_reduction <- function(n_loci, floor, n_markers) 1 - sum(n_loci >= floor) / n_markers

## The three profile floors for reported floor `canon`: low = max(1,
## round(factors[1] * canon)), high = round(factors[2] * canon).
## @return data.table(role, floor, achieved_reduction, n_units_tested,
##   marker_coverage), one row per role (low/reported/high).
.profile_floors <- function(clusters, n_markers, canon, factors = c(0.5, 2)) {
  clusters <- data.table::as.data.table(clusters)
  nl <- if ("n_loci" %in% names(clusters)) clusters$n_loci else clusters$n_snps
  fl <- c(max(1L, as.integer(round(factors[1] * canon))), as.integer(canon),
          as.integer(round(factors[2] * canon)))
  if (anyDuplicated(fl)) stop("floor profile: low/reported/high floors must differ, got ", paste(fl, collapse = "/"))
  data.table::data.table(role = c("low", "reported", "high"), floor = fl,
                         achieved_reduction = vapply(fl, function(f) .floor_reduction(nl, f, n_markers), 0),
                         n_units_tested = vapply(fl, function(f) sum(nl >= f), 0L),
                         marker_coverage = vapply(fl, function(f) sum(nl[nl >= f]) / n_markers, 0))
}

## Significant units at one floor, with their floor-independent core_snp.
## test$units follows .ld_outlier_units() row order (checked), so core_snp is
## attached by position.
.sig_units_with_core <- function(test, stage1, map, floor) {
  uf <- LDscnR:::.ld_outlier_units(stage1, map, floor)
  u <- data.table::copy(data.table::as.data.table(test$units))
  stopifnot("outlier test units must follow .ld_outlier_units() order" = identical(u$unit_id, uf$unit_id))
  u[, core_snp := uf$core_snp]
  u[significant == TRUE]
}

## Reported-region stability table for one statistic ("unit" or "simes").
## `test_by_floor` is a named list (names = floor, as character) of that
## floor's own run_stage_B()$test. Returns one row per region at the REPORTED
## floor. A region is recovered at a comparison floor when at least one of
## the significant Stage-1 units inside its span is also significant there
## (matched by core_snp) AND the region overlaps an outlier region assembled
## at that floor. This is the manuscript's rule and the one used by its
## modules (LDscnR-paper/module_3sp/R/helpers_floor_tiers.R) and simulations.
.floor_stability_table <- function(test_by_floor, canon_floor, side_floors, stage1, map) {
  canon <- test_by_floor[[as.character(canon_floor)]]
  if (is.null(canon) || !nrow(canon$regions)) return(data.table::data.table())
  reg <- data.table::copy(data.table::as.data.table(canon$regions))
  canon_sig <- .sig_units_with_core(canon, stage1, map, canon_floor)
  side <- stats::setNames(lapply(side_floors, function(fl) {
    t <- test_by_floor[[as.character(fl)]]
    if (is.null(t)) return(NULL)
    list(sig = .sig_units_with_core(t, stage1, map, fl), regions = data.table::as.data.table(t$regions))
  }), as.character(side_floors))

  recovered <- function(fl) vapply(seq_len(nrow(reg)), function(i) {
    r <- reg[i]; s <- side[[as.character(fl)]]
    if (is.null(s)) return(FALSE)
    cores <- canon_sig[Chr == r$Chr & to >= r$from & from <= r$to, core_snp]
    any(cores %in% s$sig$core_snp) && nrow(s$regions[Chr == r$Chr & to >= r$from & from <= r$to]) > 0
  }, logical(1))

  out <- reg[, .(region_id = sprintf("%s:%d-%d", Chr, as.integer(from), as.integer(to)), Chr, from, to)]
  for (fl in side_floors) data.table::set(out, j = sprintf("tick_%03d", round(fl)), value = recovered(fl))
  tick_cols <- grep("^tick_", names(out), value = TRUE)
  out[, n_ticks := 1L + rowSums(.SD), .SDcols = tick_cols]
  out[, `:=`(stability_tier = sprintf("%d/%d", n_ticks, length(side_floors) + 1L), span = to - from)]
  data.table::setorder(out, -n_ticks, -span)
  out[]
}

## @param dataset_dir Path to one dataset folder (must have input/phenotype.rds).
## @param cfg Resolved config; `cfg$size_floor` is the REPORTED floor (the
##   marker-density rule unless set explicitly). Other fields -- statistics,
##   grm_method, b_unit, b_simes, alpha, assembly, score_threshold,
##   distance_threshold, gap, n_rotations, rotation_scheme, seed, cores --
##   are used as given.
## @param factors Multipliers for the low and high comparison floors
##   (default half and twice the reported floor).
## @param annotation,chrom_lengths Passed to the REPORTED floor's
##   run_stage_B() calls only (ld_region_rotation() is a reported-floor-only,
##   structure-aware analysis, same as the permutation null).
## @param canonical_unsuffixed FALSE (default) keeps EVERY floor -- the
##   reported one included -- at its own `_floor<N>`-suffixed output paths,
##   never touching a dataset's existing results; a standalone
##   run_floor_profile() call is always safe to run alongside whatever else
##   the dataset folder already has. TRUE writes the REPORTED floor's Stage
##   A/B to the traditional UNSUFFIXED paths (output/pvalues/emmax_unit/,
##   output/stageB_emmax_unit/, ...) instead -- used by run_dataset()'s own
##   internal call. The two comparison floors stay `_floor<N>`-suffixed.
## @param force Rebuild every floor's Stage A/B even if their receipts are current.
## @return list(grid, canonical_floor, floors, regions = list(unit=, simes=),
##   summary, out_dir).
run_floor_profile <- function(dataset_dir, cfg = resolve_config(dataset_dir), factors = c(0.5, 2),
                              annotation = NULL, chrom_lengths = NULL,
                              canonical_unsuffixed = FALSE, force = FALSE) {
  if (!file.exists(file.path(dataset_dir, "input", "phenotype.rds")))
    stop("run_floor_profile() needs input/phenotype.rds (the association scan must be able to run).")

  ## force NOT forwarded here -- same reason run_stage_A()/run_stage_B()/
  ## check_structure_alignment() never forward theirs either: when called
  ## from run_dataset(), stage1's receipt was JUST written by that function's
  ## own build_stage1(force = force) call moments earlier.
  s1 <- build_stage1(dataset_dir, cfg)
  canon_floor <- as.integer(cfg$size_floor)
  grid <- .profile_floors(s1$stage1$clusters, nrow(s1$map), canon_floor, factors)
  low_floor <- grid$floor[grid$role == "low"]; high_floor <- grid$floor[grid$role == "high"]
  side_floors <- c(low_floor, high_floor)

  out_dir <- file.path(dataset_dir, "output", "floor_profile")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(grid, file.path(out_dir, "grid.csv"))
  say("=== %s: floor profile -- floors %d / %d (reported) / %d ===\n",
      basename(dataset_dir), low_floor, canon_floor, high_floor)

  test_by_floor <- list(unit = list(), simes = list())
  canon_suffix <- if (canonical_unsuffixed) "" else sprintf("_floor%d", canon_floor)

  ## ---- reported floor: the full Stage A/B, unchanged -----------------------
  say("\n--- floor = %d (reported) ---\n", canon_floor)
  cfg_c <- utils::modifyList(cfg, list(size_floor = canon_floor))
  if ("unit" %in% cfg$statistics) {
    run_stage_A(dataset_dir, utils::modifyList(cfg_c, list(statistics = "unit")),
               engine_suffix = canon_suffix, permute = TRUE, force = force)
    r <- run_stage_B(dataset_dir, paste0("emmax_unit", canon_suffix), statistic = "unit", cfg = cfg_c,
                     annotation = annotation, chrom_lengths = chrom_lengths, force = force)
    test_by_floor$unit[[as.character(canon_floor)]] <- r$test
  }
  simes_dir <- file.path(dataset_dir, "output", "pvalues", "emmax_simes")
  simes_p_obs <- NULL
  if ("simes" %in% cfg$statistics) {
    ## The shared marker scan's permutations are restricted to markers in
    ## clusters clearing simes_floor; the reported floor's permutation null
    ## is the only one needed, so simes_floor stays at the reported floor.
    run_stage_A(dataset_dir, utils::modifyList(cfg_c, list(statistics = "simes")),
               simes_floor = canon_floor, force = force)
    simes_p_obs <- readRDS(file.path(simes_dir, "p_obs.rds"))
    ## p_obs/p_perm_compact passed explicitly: run_stage_A()'s "simes" arm
    ## always writes to the UNSUFFIXED emmax_simes/ directory, so a
    ## floor-suffixed Stage B engine name would never find it on its own.
    r <- run_stage_B(dataset_dir, paste0("emmax_simes", canon_suffix), p_obs = simes_p_obs,
                     p_perm_compact = readRDS(file.path(simes_dir, "p_perm_compact.rds")),
                     statistic = "simes", cfg = cfg_c,
                     annotation = annotation, chrom_lengths = chrom_lengths, force = force)
    test_by_floor$simes[[as.character(canon_floor)]] <- r$test
  }

  ## ---- comparison floors: observed only ------------------------------------
  if ("unit" %in% cfg$statistics) {
    ## One observed consensus fit at the LOW floor; the high floor is a
    ## subset of its units (matched by core_snp, which is floor-independent).
    low_suffix <- sprintf("_floor%d", low_floor)
    cfg_l <- utils::modifyList(cfg, list(size_floor = low_floor))
    say("\n--- floor = %d (comparison, observed only; also supplies floor %d) ---\n", low_floor, high_floor)
    a_low <- run_stage_A(dataset_dir, utils::modifyList(cfg_l, list(statistics = "unit")),
                         engine_suffix = low_suffix, permute = FALSE, force = force)$unit
    ## p-value vectors follow .ld_outlier_units() ROW order (= ld_unit_matrix
    ## column order), not unit_id values, so map positions through it.
    core_low <- LDscnR:::.ld_outlier_units(s1$stage1, s1$map, low_floor)$core_snp
    stopifnot(length(core_low) == length(a_low$p_obs))
    p_by_core <- stats::setNames(as.numeric(a_low$p_obs), core_low)
    for (fl in side_floors) {
      cfg_fl <- utils::modifyList(cfg, list(size_floor = fl))
      core_fl <- LDscnR:::.ld_outlier_units(s1$stage1, s1$map, fl)$core_snp
      stopifnot("comparison floor units must be a subset of the low floor's" = all(core_fl %in% names(p_by_core)))
      r <- run_stage_B(dataset_dir, sprintf("emmax_unit_floor%d", fl), p_obs = unname(p_by_core[core_fl]),
                       p_display = as.numeric(a_low$p_obs_marker[s1$map$marker]),
                       statistic = "unit", cfg = cfg_fl, observed_only = TRUE, force = force)
      test_by_floor$unit[[as.character(fl)]] <- r$test
    }
  }
  if ("simes" %in% cfg$statistics) {
    for (fl in side_floors) {
      say("\n--- floor = %d (comparison, observed only, simes) ---\n", fl)
      r <- run_stage_B(dataset_dir, sprintf("emmax_simes_floor%d", fl), p_obs = simes_p_obs,
                       statistic = "simes", cfg = utils::modifyList(cfg, list(size_floor = fl)),
                       observed_only = TRUE, force = force)
      test_by_floor$simes[[as.character(fl)]] <- r$test
    }
  }

  ## ---- reported-region stability tables -------------------------------------
  regions <- list()
  summary_rows <- list()
  for (stat in intersect(c("unit", "simes"), cfg$statistics)) {
    tab <- .floor_stability_table(test_by_floor[[stat]], canon_floor, side_floors, s1$stage1, s1$map)
    regions[[stat]] <- tab
    data.table::fwrite(tab, file.path(out_dir, sprintf("canonical_regions_%s.csv", stat)))
    canon_test <- test_by_floor[[stat]][[as.character(canon_floor)]]
    tiers <- if (nrow(tab)) table(factor(tab$n_ticks, levels = 1:3)) else c(0L, 0L, 0L)
    summary_rows[[stat]] <- data.table::data.table(
      method = stat, floor_low = low_floor, floor = canon_floor, floor_high = high_floor,
      reduction_low = grid$achieved_reduction[1], reduction = grid$achieved_reduction[2],
      reduction_high = grid$achieved_reduction[3],
      n_units_tested = nrow(canon_test$units), n_significant = sum(canon_test$units$significant),
      n_regions = nrow(canon_test$regions),
      n_tier_1 = as.integer(tiers[1]), n_tier_2 = as.integer(tiers[2]), n_tier_3 = as.integer(tiers[3]))
  }
  summary <- data.table::rbindlist(summary_rows, fill = TRUE)
  data.table::fwrite(summary, file.path(out_dir, "summary.csv"))
  say("\nwrote %s\n", out_dir)

  list(grid = grid, canonical_floor = canon_floor, floors = grid$floor,
      regions = regions, summary = summary, out_dir = out_dir)
}
