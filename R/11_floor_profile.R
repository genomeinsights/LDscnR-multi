## =============================================================================
## LDscnR-multi/R/11_floor_profile.R
##
## A LIGHTWEIGHT, OBSERVED-ONLY size_floor SENSITIVITY CHECK -- deliberately
## separate from run_floor_sweep() (R/09_floor_sweep.R), which reruns full
## permutations at every floor and is not suitable for hundreds of datasets.
## run_floor_profile() runs the pipeline's OWN, UNMODIFIED Stage A/B functions
## exactly once with permutations (the "canonical" floor), and twice more
## observed-only (no permutations, no rotation) at two side floors, all three
## floors chosen from the dataset's ACTUAL Stage-1 cluster-size distribution
## rather than from marker count alone.
##
## This does NOT change what `size_floor` means anywhere else in the repo:
## resolve_config()'s own derivation (R/00_config.R), run_dataset(), and
## run_floor_sweep() are all untouched. This is an additional, opt-in report,
## not a redefinition -- "Keep the existing explicit size_floor option for
## reproducibility and older analyses."
##
## THE THREE FLOORS are chosen by ACHIEVED TEST-COUNT REDUCTION, not a fixed
## rule: for candidate integer floor f, reduction(f) = 1 - (number of Stage-1
## clusters with n_loci >= f) / (number of assayed markers) -- the fraction of
## per-marker tests that clustering-plus-floor avoids running at all. Targets
## default to 99.5% / 99.7% / 99.9%; the middle one (99.7%) is canonical --
## the ONLY floor that gets the structure-aware permutation null and Stage-2
## rotation test. Verified directly against real cluster-size data: 99.7%
## resolves to floor 8 on examples/3sp_chr1_chr4, floor 12 on
## examples/9sp_chr19_chr20, and floor 7 on the full 3sp panel
## (LDscnR-paper/module_3sp) -- the last of which is NOT the manuscript's own
## floor of 8 (closer to a 99.8% target there). This function reports that
## kind of discrepancy; it never silently revises anything upstream of it.
##
## REGIONS ARE MATCHED ACROSS FLOORS BY STAGE-1 CLUSTER IDENTITY, NOT SPAN
## OVERLAP: a unit's own `unit_id` is a position index into that floor's own
## size_floor-filtered cluster subset (.ld_outlier_units(), LDscnR) and is
## therefore NOT comparable across floors -- two different clusters can be
## "unit_id 5" at two different floors. `core_snp` (each cluster's designated
## representative marker, assigned once during clustering, independent of any
## floor) IS stable, and is what "is this canonical region's cluster also
## significant at the other floor" is actually asked about. Physical span
## overlap is used only as a secondary sanity check on a recovered hit, per
## the user's own instruction -- not as the primary matching key.
## =============================================================================

## Achieved test-count reduction for one candidate integer floor.
.floor_reduction <- function(n_loci, floor, n_markers) 1 - sum(n_loci >= floor) / n_markers

## For each target reduction, the integer floor whose ACHIEVED reduction is
## closest. Ties (only possible when no cluster has the exact size that would
## separate two candidate floors) are broken toward the LARGER floor --
## consistent with this pipeline's existing "exclude the marginal case"
## philosophy (e.g. the size_floor >= 2 hard minimum, R/00_config.R).
## @return data.table(target, floor, achieved_reduction, n_units_tested,
##   marker_coverage), one row per target -- floors are NOT deduplicated here;
##   run_floor_profile() only runs Stage A/B once per DISTINCT floor value,
##   but every target's own row is preserved so a repeat is visible, not hidden.
.resolve_floor_targets <- function(clusters, n_markers, targets) {
  clusters <- data.table::as.data.table(clusters)
  nl <- if ("n_loci" %in% names(clusters)) clusters$n_loci else clusters$n_snps
  max_floor <- max(nl)
  floors <- seq_len(max_floor)
  n_units <- vapply(floors, function(f) sum(nl >= f), 0L)
  reduction <- 1 - n_units / n_markers
  coverage <- vapply(floors, function(f) sum(nl[nl >= f]) / n_markers, 0)

  data.table::rbindlist(lapply(targets, function(t) {
    d <- abs(reduction - t)
    floor <- max(floors[d == min(d)])   # tie-break: larger floor
    data.table::data.table(target = t, floor = floor,
                           achieved_reduction = reduction[floor],
                           n_units_tested = n_units[floor],
                           marker_coverage = coverage[floor])
  }))
}

## Stable Stage-1 cluster identity for a set of unit_ids AT A GIVEN FLOOR --
## see file header for why unit_id itself isn't comparable across floors.
.unit_core_snp <- function(stage1_clusters, floor, unit_ids) {
  stage1_clusters <- data.table::as.data.table(stage1_clusters)
  nl <- if ("n_loci" %in% names(stage1_clusters)) stage1_clusters$n_loci else stage1_clusters$n_snps
  cl_f <- stage1_clusters[nl >= floor]
  cl_f$core_snp[unit_ids]
}

## Canonical-region stability table for one statistic ("unit"/"simes").
## `snp_results_by_floor` is a named list (names = floor, as character) of
## that floor's own `run_stage_B()$snp_results`. Returns one row per
## CANONICAL region (never per side-floor-only region -- this profile is
## about how robust the PRIMARY, reported findings are, not a census of every
## floor's own discoveries).
.floor_stability_table <- function(snp_results_by_floor, canon_floor, side_floors, stage1_clusters) {
  canon_sr <- snp_results_by_floor[[as.character(canon_floor)]]
  if (is.null(canon_sr) || !any(!is.na(canon_sr$region_id))) return(data.table::data.table())

  canon_units <- unique(canon_sr[!is.na(region_id), .(unit_id, region_id, Chr)])
  canon_units[, core_snp := .unit_core_snp(stage1_clusters, canon_floor, unit_id)]
  spans <- canon_sr[!is.na(region_id), .(from = min(Pos), to = max(Pos)), by = .(region_id, Chr)]

  ## Per side floor: significant units' stable cluster identity, joined to
  ## THAT FLOOR'S OWN region span (for the physical-overlap check) -- a
  ## cluster's own significance and its Stage-2 grouping are BOTH floor
  ## specific, even though the cluster's identity (core_snp) is not.
  side_lookup <- stats::setNames(lapply(side_floors, function(fl) {
    sr <- snp_results_by_floor[[as.character(fl)]]
    if (is.null(sr)) return(NULL)
    sig <- unique(sr[!is.na(unit_id) & significant == TRUE, .(unit_id, region_id)])
    if (!nrow(sig)) return(NULL)
    sig[, core_snp := .unit_core_snp(stage1_clusters, fl, unit_id)]
    reg_spans <- sr[!is.na(region_id), .(from = min(Pos), to = max(Pos)), by = region_id]
    sig[reg_spans, on = "region_id", `:=`(side_from = i.from, side_to = i.to)]
    sig
  }), as.character(side_floors))

  tick_col  <- function(fl) sprintf("tick_%03d", round(fl))
  check_col <- function(fl) sprintf("overlap_check_%03d", round(fl))

  rows <- lapply(unique(canon_units$region_id), function(rid) {
    core_snps <- canon_units[region_id == rid, core_snp]
    sp <- spans[region_id == rid]
    row <- data.table::data.table(region_id = rid, Chr = sp$Chr, from = sp$from, to = sp$to)
    n_ticks <- 1L   # the canonical floor itself always counts
    for (fl in side_floors) {
      su <- side_lookup[[as.character(fl)]]
      hit <- if (!is.null(su)) su[core_snp %in% core_snps][1] else NULL
      recovered <- !is.null(hit) && nrow(hit) > 0 && !is.na(hit$core_snp[1])
      overlap_ok <- NA
      if (isTRUE(recovered)) {
        overlap_ok <- hit$side_from[1] <= sp$to && hit$side_to[1] >= sp$from
        n_ticks <- n_ticks + 1L
      }
      ## data.table::set(), not row[[col]] <- value: the latter is a base-R
      ## list assignment that data.table only tolerates via a defensive
      ## shallow copy (and warns about it, harmlessly but noisily, every
      ## single call) -- set() adds the column by reference, as intended.
      data.table::set(row, j = tick_col(fl), value = isTRUE(recovered))
      data.table::set(row, j = check_col(fl), value = overlap_ok)
    }
    row[, `:=`(n_ticks = n_ticks,
              stability_tier = sprintf("%d/%d", n_ticks, length(side_floors) + 1L))]
    row
  })

  out <- data.table::rbindlist(rows, fill = TRUE)
  ## -n_ticks/-span, not -stability_tier: that column is a display STRING
  ## ("3/3") -- setorder() needs the underlying numeric count to sort
  ## descending, and setorder() takes column names, not arbitrary expressions
  ## like -(to - from), hence the explicit `span` column (kept in the output,
  ## not just a sort key -- "physical span" was asked for as its own column).
  out[, span := to - from]
  data.table::setorder(out, -n_ticks, -span)
  out[]
}

## @param dataset_dir Path to one dataset folder (must have input/phenotype.rds).
## @param targets Reduction targets; the pipeline's own Stage A/B (with
##   permutations and, if given, ld_region_rotation()) runs ONLY at
##   `canonical_target`. The other two run observed-only: a fast, per-floor
##   EMMAX-consensus scan and Stage-2 assembly, but no permutation null, no
##   rotation test -- see run_stage_A(permute=)/run_stage_B(observed_only=).
## @param canonical_target Must be one of `targets`.
## @param cfg Resolved config (its OWN `size_floor` is ignored here --
##   every floor this function tests comes from `targets`, not from
##   resolve_config()'s derivation or an explicit override in it. Other
##   fields -- statistics, grm_method, b_unit, b_simes, alpha, assembly,
##   score_threshold, distance_threshold, gap, n_rotations, rotation_scheme,
##   seed, cores -- are used as given).
## @param annotation,chrom_lengths Passed to the CANONICAL floor's
##   run_stage_B() calls only (ld_region_rotation() is a canonical-only,
##   structure-aware analysis, same as the permutation null).
## @param canonical_unsuffixed FALSE (default) keeps EVERY floor -- canonical
##   included -- at its own `_floor<N>`-suffixed output paths, never touching
##   a dataset's existing default-floor results; a standalone
##   run_floor_profile() call is always safe to run alongside whatever else
##   the dataset folder already has. TRUE writes the CANONICAL floor's Stage
##   A/B to the traditional UNSUFFIXED paths (output/pvalues/emmax_unit/,
##   output/stageB_emmax_unit/, ...) instead -- used by run_dataset()'s own
##   internal call, where the canonical floor IS now that dataset's primary,
##   reported result (see run_dataset()'s own header for why). The two side
##   floors stay `_floor<N>`-suffixed either way.
## @param force Rebuild every floor's Stage A/B even if their receipts are current.
## @return list(grid, canonical_floor, floors, regions = list(unit=, simes=),
##   summary, out_dir).
run_floor_profile <- function(dataset_dir, targets = c(0.995, 0.997, 0.999),
                              canonical_target = 0.997, cfg = resolve_config(dataset_dir),
                              annotation = NULL, chrom_lengths = NULL,
                              canonical_unsuffixed = FALSE, force = FALSE) {
  if (!file.exists(file.path(dataset_dir, "input", "phenotype.rds")))
    stop("run_floor_profile() needs input/phenotype.rds (Stage A must be able to run).")
  if (!(canonical_target %in% targets))
    stop("`canonical_target` (", canonical_target, ") must be one of `targets` (",
         paste(targets, collapse = ", "), ").")

  ## force NOT forwarded here -- same reason run_stage_A()/run_stage_B()/
  ## check_structure_alignment() never forward theirs either (and the exact
  ## bug report_ld_decay() had until it was fixed, R/06_run_dataset.R's own
  ## header): when called from run_dataset(), stage1's receipt was JUST
  ## written by that function's own build_stage1(force = force) call moments
  ## earlier, so forwarding `force` here would force a full, expensive stage1
  ## rebuild a SECOND time for nothing. build_stage1()'s own staleness check
  ## still rebuilds correctly if anything upstream genuinely changed.
  s1 <- build_stage1(dataset_dir, cfg)
  grid <- .resolve_floor_targets(s1$stage1$clusters, nrow(s1$map), targets)
  grid[, is_canonical := target == canonical_target]

  out_dir <- file.path(dataset_dir, "output", "floor_profile")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(grid, file.path(out_dir, "grid.csv"))

  canon_floor <- grid$floor[grid$is_canonical]
  distinct_floors <- unique(grid$floor)
  side_floors <- setdiff(distinct_floors, canon_floor)
  say("=== %s: floor profile -- targets %s -> floors %s (canonical = %d) ===\n",
      basename(dataset_dir), paste(sprintf("%.1f%%", 100 * targets), collapse = ", "),
      paste(distinct_floors, collapse = ","), canon_floor)
  if (length(distinct_floors) < length(targets))
    say("    note: %d target(s) share a floor with another -- reported once, not rerun.\n",
        length(targets) - length(distinct_floors))

  snp_results_by_floor <- list(unit = list(), simes = list())
  test_by_floor <- list(unit = list(), simes = list())
  simes_p_obs <- NULL

  ## Canonical first, always -- guarantees the "simes" arm's shared observed
  ## scan (output/pvalues/emmax_simes/p_obs.rds) exists before any side floor
  ## needs to reuse it directly.
  ordered_floors <- c(canon_floor, side_floors)

  for (fl in ordered_floors) {
    is_canon <- identical(fl, canon_floor)
    suffix <- if (is_canon && canonical_unsuffixed) "" else sprintf("_floor%d", fl)
    cfg_fl <- utils::modifyList(cfg, list(size_floor = fl))
    say("\n--- floor = %d (%s) ---\n", fl, if (is_canon) "canonical" else "observed-only")

    if ("unit" %in% cfg$statistics) {
      run_stage_A(dataset_dir, utils::modifyList(cfg_fl, list(statistics = "unit")),
                 engine_suffix = suffix, permute = is_canon, force = force)
      r <- run_stage_B(dataset_dir, paste0("emmax_unit", suffix), statistic = "unit", cfg = cfg_fl,
                       annotation = if (is_canon) annotation else NULL,
                       chrom_lengths = chrom_lengths, observed_only = !is_canon, force = force)
      test_by_floor$unit[[as.character(fl)]] <- r$test
      snp_results_by_floor$unit[[as.character(fl)]] <- r$snp_results
    }

    if ("simes" %in% cfg$statistics) {
      simes_dir <- file.path(dataset_dir, "output", "pvalues", "emmax_simes")
      if (is_canon) {
        run_stage_A(dataset_dir, utils::modifyList(cfg_fl, list(statistics = "simes")),
                   simes_floor = fl, force = force)
        simes_p_obs <- readRDS(file.path(simes_dir, "p_obs.rds"))
        ## p_obs/p_perm_compact passed explicitly, not left to run_stage_B()'s
        ## own auto-resolution: that logic looks for output/pvalues/<engine>/,
        ## using the STAGE B engine NAME given below -- but run_stage_A()'s
        ## "simes" arm always writes to the UNSUFFIXED emmax_simes/ directory
        ## regardless of engine_suffix (its observed scan is floor-independent
        ## and deliberately shared), so a floor-suffixed Stage B engine name
        ## here would never find it (exactly the mismatch this fixed after
        ## first hitting it directly on this example).
        r <- run_stage_B(dataset_dir, paste0("emmax_simes", suffix), p_obs = simes_p_obs,
                         p_perm_compact = readRDS(file.path(simes_dir, "p_perm_compact.rds")),
                         statistic = "simes", cfg = cfg_fl,
                         annotation = annotation, chrom_lengths = chrom_lengths, force = force)
      } else {
        ## Reused directly, not routed back through run_stage_A(): that
        ## arm's observed scan is floor-independent and already on disk from
        ## the canonical pass above -- see file header for why this sidesteps
        ## run_stage_A()'s simes receipt (keyed on simes_floor) entirely
        ## rather than entangling it with an unrelated "observed_only" idea.
        r <- run_stage_B(dataset_dir, paste0("emmax_simes", suffix), p_obs = simes_p_obs,
                         statistic = "simes", cfg = cfg_fl, observed_only = TRUE, force = force)
      }
      test_by_floor$simes[[as.character(fl)]] <- r$test
      snp_results_by_floor$simes[[as.character(fl)]] <- r$snp_results
    }
  }

  ## ---- canonical-region stability tables --------------------------------
  regions <- list()
  summary_rows <- list()
  for (stat in intersect(c("unit", "simes"), cfg$statistics)) {
    tab <- .floor_stability_table(snp_results_by_floor[[stat]], canon_floor, side_floors, s1$stage1$clusters)
    regions[[stat]] <- tab
    data.table::fwrite(tab, file.path(out_dir, sprintf("canonical_regions_%s.csv", stat)))

    canon_test <- test_by_floor[[stat]][[as.character(canon_floor)]]
    row <- data.table::data.table(method = stat)
    for (i in seq_len(nrow(grid))) {
      tgt <- round(grid$target[i] * 1000)
      data.table::set(row, j = sprintf("floor_%03d", tgt), value = grid$floor[i])
      data.table::set(row, j = sprintf("reduction_%03d", tgt), value = grid$achieved_reduction[i])
    }
    row[, `:=`(n_units_tested_canonical = nrow(canon_test$units),
              n_significant_canonical = sum(canon_test$units$significant),
              n_regions_canonical = nrow(canon_test$regions))]
    summary_rows[[stat]] <- row
  }
  summary <- data.table::rbindlist(summary_rows, fill = TRUE)
  data.table::fwrite(summary, file.path(out_dir, "summary.csv"))
  say("\nwrote %s\n", out_dir)

  list(grid = grid, canonical_floor = canon_floor, floors = distinct_floors,
      regions = regions, summary = summary, out_dir = out_dir)
}
