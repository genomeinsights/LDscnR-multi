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

## Stable Stage-1 cluster identity for a set of unit_ids AT A GIVEN FLOOR --
## see file header for why unit_id itself isn't comparable across floors.
## unit_id VALUES (as carried by test$units and snp_results) are row indices
## into the size-filtered cluster table: .ld_outlier_units() assigns unit_id =
## cluster-table row, then sorts its rows by position. So index the filtered
## cluster table by the unit_id value -- never by row position in
## .ld_outlier_units() output (that is the order of p-value vectors and unit
## matrices, used separately below).
.unit_core_snp <- function(stage1, map, floor, unit_ids) {
  cl <- data.table::as.data.table(stage1$clusters)
  nl <- if ("n_loci" %in% names(cl)) cl$n_loci else cl$n_snps
  cl[nl >= floor]$core_snp[unit_ids]
}

## Canonical-region stability table for one statistic ("unit"/"simes").
## `snp_results_by_floor` is a named list (names = floor, as character) of
## that floor's own `run_stage_B()$snp_results`. Returns one row per
## CANONICAL region (never per side-floor-only region -- this profile is
## about how robust the PRIMARY, reported findings are, not a census of every
## floor's own discoveries).
.floor_stability_table <- function(snp_results_by_floor, canon_floor, side_floors, stage1, map) {
  canon_sr <- snp_results_by_floor[[as.character(canon_floor)]]
  if (is.null(canon_sr) || !any(!is.na(canon_sr$region_id))) return(data.table::data.table())

  canon_units <- unique(canon_sr[!is.na(region_id), .(unit_id, region_id, Chr)])
  canon_units[, core_snp := .unit_core_snp(stage1, map, canon_floor, unit_id)]
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
    sig[, core_snp := .unit_core_snp(stage1, map, fl, unit_id)]
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
    stop("run_floor_profile() needs input/phenotype.rds (Stage A must be able to run).")

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

  snp_results_by_floor <- list(unit = list(), simes = list())
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
    snp_results_by_floor$unit[[as.character(canon_floor)]] <- r$snp_results
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
    snp_results_by_floor$simes[[as.character(canon_floor)]] <- r$snp_results
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
      snp_results_by_floor$unit[[as.character(fl)]] <- r$snp_results
    }
  }
  if ("simes" %in% cfg$statistics) {
    for (fl in side_floors) {
      say("\n--- floor = %d (comparison, observed only, simes) ---\n", fl)
      r <- run_stage_B(dataset_dir, sprintf("emmax_simes_floor%d", fl), p_obs = simes_p_obs,
                       statistic = "simes", cfg = utils::modifyList(cfg, list(size_floor = fl)),
                       observed_only = TRUE, force = force)
      test_by_floor$simes[[as.character(fl)]] <- r$test
      snp_results_by_floor$simes[[as.character(fl)]] <- r$snp_results
    }
  }

  ## ---- reported-region stability tables -------------------------------------
  regions <- list()
  summary_rows <- list()
  for (stat in intersect(c("unit", "simes"), cfg$statistics)) {
    tab <- .floor_stability_table(snp_results_by_floor[[stat]], canon_floor, side_floors, s1$stage1, s1$map)
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
