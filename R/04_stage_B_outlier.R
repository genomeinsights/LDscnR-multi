## =============================================================================
## LDscnR-multi/R/04_stage_B_outlier.R
##
## STAGE B: any p-value source -> ld_outlier_test() -> per-marker results +
## region table. This is what makes Stage 2 engine-agnostic in practice, not
## just in theory: an "emmax_unit"/"emmax_simes" engine from Stage A and an
## external engine dropped in external_pvalues/*.rds (LFMM, or anything else
## that produces one p-value per marker) go through the exact same call, the
## same way module_3sp/R/04_lfmm.R feeds LFMM's precomputed p-values through
## the identical ld_outlier_test() call 03_EMMAX.R uses for EMMAX.
##
## PER-MARKER OUTPUT. ld_outlier_test()$units is one row per TESTED CLUSTER,
## not per marker -- turning that back into a per-marker table
## (snp_results.csv) uses EXACT CLUSTER MEMBERSHIP, not a spatial join on
## coordinates, at BOTH levels this function assigns:
##
## - unit_id: replicated from ld_outlier_test()'s own internal filter
##   (LDscnR:::.ld_outlier_units(), R/ld_outlier_internal.R) -- the SAME
##   size_floor filter over stage1$clusters, in ORIGINAL row order, assigns
##   unit_id = seq_len(nrow(filtered)) before that internal table is resorted
##   by Chr/from, so row i of the SAME filter here is unit_id i regardless of
##   what order test$units itself ended up in.
## - region_id (assembly = "stage2_discovered" only): replicated from
##   ld_outlier_test()'s own internal ld_prune_and_eMLG() call
##   (.stage2_unit_to_region(), below) -- NOT a span-containment join. An
##   external audit found that even after the unit_id fix, joining a unit's
##   OWN span against a region's span is still ambiguous whenever regions
##   overlap physically (quantified: 6 significant units in the bundled 3sp
##   LFMM output, 4-6 in each bundled 9sp EMMAX output, fit inside more than
##   one region's bounds) -- Stage-2 groups are membership-based too, so their
##   bounding spans can overlap the same way Stage-1 clusters' can.
##   assembly = "physical" doesn't need this: its regions come from a plain
##   gap-merge of already-sorted, already-disjoint unit spans
##   (.physical_merge()), which cannot produce overlapping output by
##   construction, so the existing span join stays correct there.
##
## POSITION-INDEPENDENT ASSIGNMENT. Every column below is assigned either (a)
## directly onto a freshly-copied `map`-ordered table, before anything reorders
## it, when the source is a plain vector aligned to map$marker by CONTRACT
## (marker_p/marker_q/p_display when marker-aligned), or (b) via a data.table
## update-join keyed on `marker` (`out[lookup, on = "marker", (cols) := ...]`),
## which matches by NAME and leaves `out`'s row order untouched, never via
## merge() + a manual setorder() meant to restore it. An earlier version did
## the latter: it built `out` via merge() (whose row order is undefined, not
## `map`'s), called setorder(out, Chr, Pos) to put it back, and only THEN
## assigned `marker_p := as.numeric(p_obs)` positionally -- correct only when
## `map` happened to already be sorted by (Chr, Pos), which the dataset
## contract never required. An external audit reproduced the swap directly
## on a two-marker case.
.build_snp_results <- function(map, test, p_obs, statistic, stage1, size_floor,
                               GTs = NULL, LD_decay = NULL, p_display = NULL) {
  out <- data.table::as.data.table(map)[, .(marker, Chr, Pos)]   # map's own row order, untouched from here on

  ## marker_p/marker_q/p_display: assigned NOW, positionally, while `out`'s
  ## row order is still provably identical to `map`'s (the contract p_obs/
  ## p_display are aligned to) -- before any join below can disturb it.
  out[, statistic := statistic]
  ## Carried straight from `map` (already in this same row order, positional
  ## and safe for the same reason marker_p/p_display are, above) -- lets
  ## ldm_manhattan(value = "ld_w_095") plot local-LD support without a
  ## separate per-dataset data source; NA if `map` predates this column
  ## (an older cached stage1, or a caller supplying its own `map`).
  out[, ld_w_095 := if ("ld_w_095" %in% names(map)) as.numeric(map$ld_w_095) else NA_real_]
  if (statistic == "simes") {
    if (length(p_obs) != nrow(map))
      stop("statistic = \"simes\": p_obs must have one value per marker, aligned to map (",
           length(p_obs), " vs ", nrow(map), ").")
    out[, marker_p := as.numeric(p_obs)]
    out[, marker_q := stats::p.adjust(marker_p, method = "BH")]
    out[, `:=`(p_display = marker_p, q_display = marker_q)]
  } else {
    out[, `:=`(marker_p = NA_real_, marker_q = NA_real_)]
    if (!is.null(p_display)) {
      if (length(p_display) != nrow(map))
        stop("p_display must have one value per marker, aligned to map (",
             length(p_display), " vs ", nrow(map), ").")
      out[, p_display := as.numeric(p_display)]
      out[, q_display := stats::p.adjust(p_display, method = "BH")]
    }   # else: p_display/q_display assigned below, from unit_p/unit_q (also a join, not a position)
  }

  ## unit_id/unit_p/unit_q/significant: a NAME-KEYED update-join (`on =
  ## "marker"`), so it cannot depend on -- or disturb -- `out`'s row order.
  cl <- data.table::as.data.table(stage1$clusters)
  nl <- if ("n_loci" %in% names(cl)) cl$n_loci else cl$n_snps
  cl_f <- cl[nl >= size_floor]
  stopifnot("stage1/size_floor do not reproduce test$units' own filter -- pass the SAME stage1 and size_floor ld_outlier_test() was called with" =
              nrow(cl_f) == nrow(test$units))
  out[, `:=`(unit_id = NA_integer_, unit_p = NA_real_, unit_q = NA_real_, significant = NA)]
  if (nrow(cl_f)) {
    units_pq <- test$units[, .(unit_id, unit_p = p, unit_q = q, significant)]
    member_map <- data.table::data.table(
      marker = unlist(cl_f$members, use.names = FALSE),
      unit_id = rep.int(seq_len(nrow(cl_f)), lengths(cl_f$members)))
    member_map <- merge(member_map, units_pq, by = "unit_id", all.x = TRUE)
    out[member_map, on = "marker",
       `:=`(unit_id = i.unit_id, unit_p = i.unit_p, unit_q = i.unit_q, significant = i.significant)]
  }
  if (statistic == "unit" && is.null(p_display)) out[, `:=`(p_display = unit_p, q_display = unit_q)]

  ## region_id: exact Stage-2 membership for "stage2_discovered" (needs
  ## GTs/LD_decay to reconstruct -- see .stage2_unit_to_region()); the
  ## existing span join for "physical", where it's already exact (see header).
  out[, region_id := NA_character_]
  if (nrow(test$regions)) {
    if (identical(test$params$assembly, "stage2_discovered")) {
      u2r <- .stage2_unit_to_region(test, stage1, map, GTs, LD_decay)
      if (!is.null(u2r)) out[u2r, on = "unit_id", region_id := i.region_id]
    } else {
      regions <- data.table::copy(test$regions)
      regions[, region_id := sprintf("%s:%.0f-%.0f", Chr, from, to)]
      data.table::setkey(regions, Chr, from, to)
      uspan <- out[!is.na(unit_id), .(marker, Chr, from = Pos, to = Pos)]
      ## "physical" assembly's regions are disjoint gap-merges of the units'
      ## own spans, so a marker's own position (equally, its unit's span) is
      ## an exact, unambiguous join here -- see the file header.
      rj <- data.table::foverlaps(uspan, regions[, .(Chr, from, to, region_id)],
                                  by.x = c("Chr", "from", "to"), type = "within",
                                  mult = "first", nomatch = NA)
      out[rj, on = "marker", region_id := i.region_id]
    }
  }

  out[, tested := !is.na(unit_id)]
  data.table::setcolorder(out, c("marker", "Chr", "Pos", "statistic",
                                 "marker_p", "marker_q", "unit_p", "unit_q", "p_display", "q_display",
                                 "unit_id", "region_id", "significant", "tested", "ld_w_095"))
  out[]
}

## Exact Stage-2 group membership for assembly = "stage2_discovered" --
## ld_outlier_test()'s own return exposes only the merged region SPANS, not
## which significant units were merged into which; this reruns its internal
## ld_prune_and_eMLG() call to recover that. Mirrors that call's OWN
## parameters exactly (LDscnR::ld_outlier_test(), the assembly =
## "stage2_discovered" branch) -- flagged, not hidden, as fragile to a future
## change in that specific internal call. Returns a data.table (unit_id,
## region_id), one row per significant unit, or NULL when there is nothing to
## assemble.
.stage2_unit_to_region <- function(test, stage1, map, GTs, LD_decay) {
  p <- test$params
  sig <- test$units[significant == TRUE]
  if (!nrow(sig) || !nrow(test$regions)) return(NULL)
  if (is.null(GTs) || is.null(LD_decay))
    stop(".stage2_unit_to_region() needs GTs and LD_decay to reconstruct assembly = \"stage2_discovered\"'s groups.")

  cl <- data.table::as.data.table(stage1$clusters)
  nl <- if ("n_loci" %in% names(cl)) cl$n_loci else cl$n_snps
  cl_sig <- cl[nl >= p$size_floor][sig$unit_id]
  mk_sig <- unlist(cl_sig$members, use.names = FALSE)
  ms_sig <- data.table::as.data.table(stage1$map_snp)[marker %chin% mk_sig]
  sub <- structure(list(map_snp = ms_sig, clusters = cl_sig, pruned = cl_sig$core_snp),
                   class = "ld_complexity_reduction")
  pr <- ld_prune_and_eMLG(GTs = GTs[, mk_sig, drop = FALSE], stage1 = sub,
                          ld_w_col = "n_loci", ld_w_threshold = 0,
                          LD_decay = LD_decay, min_r2_rho = stage1$params$rho,
                          score_threshold = p$score_threshold,
                          distance_threshold = p$distance_threshold,
                          compute_unflagged_eMLG = FALSE, min_n_loci_eMLG = 1,
                          min_n_loci_flag = 1, cores = 1)
  g <- data.table::as.data.table(pr$groups)

  ## Each significant unit's own representative marker (core_snp -- always
  ## one of that cluster's own members) identifies which reconstructed group
  ## it landed in: this assembly consolidates whole clusters, never fragments
  ## one across two groups, so any single member marker resolves it exactly.
  marker_to_grp <- data.table::data.table(
    marker = unlist(g$members, use.names = FALSE),
    grp_id = rep.int(seq_len(nrow(g)), lengths(g$members)))
  rep_marker <- data.table::data.table(unit_id = sig$unit_id, marker = cl_sig$core_snp)
  rep_marker <- merge(rep_marker, marker_to_grp, by = "marker", all.x = TRUE)
  if (anyNA(rep_marker$grp_id))
    stop(".stage2_unit_to_region(): ", sum(is.na(rep_marker$grp_id)),
        " significant unit(s)' representative marker was not found in any reconstructed ",
        "Stage-2 group -- this function's ld_prune_and_eMLG() call no longer matches ",
        "ld_outlier_test()'s internal one; see this function's header.")

  ## Label each reconstructed group by its OWN span, the SAME construction
  ## the caller uses for test$regions -- both come from the same `pr$groups`
  ## here, so the strings agree whenever the two describe the same group.
  mp2 <- data.table::as.data.table(map)
  g_marker <- unlist(g$members, use.names = FALSE)
  g_idx <- match(g_marker, mp2$marker)
  g_id  <- rep.int(seq_len(nrow(g)), lengths(g$members))
  span <- data.table::data.table(grp_id = g_id, Chr = mp2$Chr[g_idx], Pos = mp2$Pos[g_idx])[
    , .(from = min(Pos), to = max(Pos)), by = .(grp_id, Chr)]
  span[, region_id := sprintf("%s:%.0f-%.0f", Chr, from, to)]

  merge(rep_marker[, .(unit_id, grp_id)], span[, .(grp_id, region_id)], by = "grp_id")[, .(unit_id, region_id)]
}

## @param dataset_dir Path to one dataset folder.
## @param engine Name for this p-value source (e.g. "emmax_unit", "emmax_simes",
##   or an external_pvalues/<engine>.rds basename). Determines output/stageB_<engine>/.
## @param p_obs,p_perm Optional. If NULL, resolved from
##   output/pvalues/<engine>/ (Stage A) or external_pvalues/<engine>.rds.
## @param statistic "unit" or "simes"; default guessed from `engine`
##   ("emmax_unit" -> "unit", anything else -> "simes", since external
##   p-values are always marker-aligned).
## @param annotation,chrom_lengths Optional bed-like overlap set for
##   ld_region_rotation(); if annotation is given and chrom_lengths is not,
##   chrom_lengths is derived from `map` (max Pos per Chr).
## @param cfg Resolved config; reads size_floor, alpha, assembly,
##   score_threshold, distance_threshold, gap, n_rotations, rotation_scheme, cores.
## @param force Rebuild even if the receipt says nothing changed.
## @return list(test, perm, rotation, snp_results, out_dir) -- perm/rotation NULL
##   when no p_perm/annotation was available.
run_stage_B <- function(dataset_dir, engine, p_obs = NULL, p_perm = NULL, p_display = NULL,
                        statistic = if (engine == "emmax_unit") "unit" else "simes",
                        cfg = resolve_config(dataset_dir),
                        annotation = NULL, chrom_lengths = NULL, force = FALSE) {
  s1 <- build_stage1(dataset_dir, cfg)
  map <- s1$map

  src_paths <- character()
  if (is.null(p_obs)) {
    stageA_dir <- file.path(dataset_dir, "output", "pvalues", engine)
    ext_f <- file.path(dataset_dir, "external_pvalues", paste0(engine, ".rds"))
    if (file.exists(file.path(stageA_dir, "p_obs.rds"))) {
      p_obs <- readRDS(file.path(stageA_dir, "p_obs.rds"))
      src_paths <- c(src_paths, file.path(stageA_dir, "p_obs.rds"))
      if (is.null(p_perm) && file.exists(file.path(stageA_dir, "p_perm.rds"))) {
        p_perm <- readRDS(file.path(stageA_dir, "p_perm.rds"))
        src_paths <- c(src_paths, file.path(stageA_dir, "p_perm.rds"))
      }
      if (is.null(p_display) && file.exists(file.path(stageA_dir, "p_obs_marker.rds"))) {
        p_display <- as.numeric(readRDS(file.path(stageA_dir, "p_obs_marker.rds"))[map$marker])
        src_paths <- c(src_paths, file.path(stageA_dir, "p_obs_marker.rds"))
      }
    } else if (file.exists(ext_f)) {
      p_named <- readRDS(ext_f)
      p_obs <- as.numeric(p_named[map$marker])
      src_paths <- c(src_paths, ext_f)
      statistic <- "simes"
    } else {
      stop("run_stage_B(): no p-values for engine \"", engine, "\" -- neither ",
           stageA_dir, "/p_obs.rds nor ", ext_f, " exists.")
    }
  }

  if (!is.null(annotation) && is.null(chrom_lengths))
    chrom_lengths <- map[, .(len = max(Pos)), by = Chr]

  out_dir <- file.path(dataset_dir, "output", paste0("stageB_", engine))
  ## Content digests, not just source FILE paths: when p_obs/p_perm/
  ## annotation/chrom_lengths are supplied directly (bypassing the
  ## src_paths auto-resolution above -- e.g. run_floor_sweep()'s reused
  ## "simes" p-values), stage_stale() had nothing to hash their actual
  ## VALUES against, so two calls with the same `engine` but genuinely
  ## different supplied p-values would false-positive as "up to date" and
  ## silently reuse the first call's result.
  params <- list(statistic = statistic, size_floor = cfg$size_floor, alpha = cfg$alpha,
                 assembly = cfg$assembly, score_threshold = cfg$score_threshold,
                 distance_threshold = cfg$distance_threshold, gap = cfg$gap,
                 n_rotations = cfg$n_rotations, rotation_scheme = cfg$rotation_scheme,
                 p_obs_digest = digest::digest(p_obs, algo = "sha256"),
                 p_perm_digest = if (!is.null(p_perm)) digest::digest(p_perm, algo = "sha256") else NA_character_,
                 annotation_digest = if (!is.null(annotation)) digest::digest(annotation, algo = "sha256") else NA_character_,
                 chrom_lengths_digest = if (!is.null(chrom_lengths)) digest::digest(chrom_lengths, algo = "sha256") else NA_character_)
  inputs <- c(src_paths, receipt_path(file.path(dataset_dir, "cache")))

  if (!force && !stage_stale(out_dir, inputs, params, label = paste0("stageB_", engine))) {
    return(list(test = readRDS(file.path(out_dir, "outlier_test.rds")),
               perm = if (file.exists(file.path(out_dir, "outlier_perm.rds")))
                        readRDS(file.path(out_dir, "outlier_perm.rds")) else NULL,
               rotation = if (file.exists(file.path(out_dir, "region_rotation.rds")))
                            readRDS(file.path(out_dir, "region_rotation.rds")) else NULL,
               snp_results = data.table::fread(file.path(out_dir, "snp_results.csv")),
               out_dir = out_dir))
  }

  say("    [stageB:%s] ld_outlier_test(statistic = \"%s\")\n", engine, statistic)
  test <- ld_outlier_test(s1$stage1, map, p_obs, statistic = statistic, size_floor = cfg$size_floor,
                          alpha = cfg$alpha, assembly = cfg$assembly, GTs = s1$genotypes,
                          LD_decay = s1$LD_decay, score_threshold = cfg$score_threshold,
                          distance_threshold = cfg$distance_threshold, gap = cfg$gap)
  print(test)

  perm <- NULL
  if (!is.null(p_perm)) {
    B <- if (is.matrix(p_perm)) ncol(p_perm) else length(p_perm)
    say("    [stageB:%s] ld_outlier_perm(): %d surrogates\n", engine, B)
    perm <- ld_outlier_perm(test, s1$stage1, map, p_perm, GTs = s1$genotypes, LD_decay = s1$LD_decay,
                            B = B, level = "units", cores = cfg$cores)
    print(perm)
  }

  rotation <- NULL
  if (!is.null(annotation) && nrow(test$regions)) {
    say("    [stageB:%s] ld_region_rotation(): %d draws\n", engine, cfg$n_rotations)
    rotation <- ld_region_rotation(data.table::copy(test$regions), annotation, chrom_lengths,
                                   scheme = cfg$rotation_scheme, n_rotations = cfg$n_rotations,
                                   seed = cfg$seed)
    print(rotation)
  }

  snp_results <- .build_snp_results(map, test, p_obs, statistic, stage1 = s1$stage1,
                                    size_floor = cfg$size_floor, GTs = s1$genotypes,
                                    LD_decay = s1$LD_decay, p_display = p_display)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(test, file.path(out_dir, "outlier_test.rds"))
  perm_f <- file.path(out_dir, "outlier_perm.rds")
  rot_f  <- file.path(out_dir, "region_rotation.rds")
  if (!is.null(perm)) saveRDS(perm, perm_f) else if (file.exists(perm_f)) file.remove(perm_f)
  if (!is.null(rotation)) saveRDS(rotation, rot_f) else if (file.exists(rot_f)) file.remove(rot_f)
  data.table::fwrite(snp_results, file.path(out_dir, "snp_results.csv"))
  data.table::fwrite(test$regions, file.path(out_dir, "region_table.csv"))
  ## Only the outputs THIS call actually wrote -- a leftover outlier_perm.rds/
  ## region_rotation.rds from an earlier call that supplied p_perm/annotation
  ## is removed above, not left for stage_stale()'s existence check to trip
  ## over, and not falsely promised here for a call that has neither.
  write_receipt(out_dir, inputs = inputs, params = params,
                outputs = file.path(out_dir, c("outlier_test.rds", "snp_results.csv", "region_table.csv",
                                               if (!is.null(perm)) "outlier_perm.rds",
                                               if (!is.null(rotation)) "region_rotation.rds")))

  list(test = test, perm = perm, rotation = rotation, snp_results = snp_results, out_dir = out_dir)
}
