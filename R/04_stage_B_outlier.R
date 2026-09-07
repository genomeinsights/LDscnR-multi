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
## not per marker -- turning that back into a per-marker table (snp_results.csv)
## uses a spatial join (marker position within a unit's/region's [from, to]
## span), the same pattern module_3sp/R_figures/figure_manhattan.R uses to
## colour points by region, rather than re-deriving cluster membership from
## stage1$clusters (which ld_outlier_test()'s internal unit numbering does not
## expose in the same row order it returns).
##
## For `statistic = "simes"` (marker-aligned p_obs, including every external
## engine), a marker's own p/q are its own -- q is BH over the full p_obs
## vector, independent of which clusters were tested, matching module_3sp's
## display convention. For `statistic = "unit"`, no per-marker p exists at all
## (p_obs is one value per cluster-summary-variable); a marker's p/q are its
## unit's aggregate values broadcast to every member. This is a display choice,
## documented rather than hidden: it is the evidence that drove that marker's
## cluster's significance, not independent per-SNP evidence.
## =============================================================================

## Turn an ld_outlier_test() result back into one row per marker.
.build_snp_results <- function(map, test, p_obs, statistic) {
  mp <- data.table::as.data.table(map)[, .(marker, Chr, Pos)]
  qx <- mp[, .(marker, Chr, from = Pos, to = Pos)]

  units <- data.table::copy(test$units)
  if (nrow(units)) {
    data.table::setkey(units, Chr, from, to)
    uj <- data.table::foverlaps(qx, units, by.x = c("Chr", "from", "to"),
                                type = "within", mult = "first", nomatch = NA)
  } else {
    uj <- data.table::data.table(unit_id = NA_integer_, p = NA_real_, q = NA_real_,
                                 significant = NA)[rep(1L, nrow(mp))]
  }

  regions <- data.table::copy(test$regions)
  if (nrow(regions)) {
    regions[, region_id := sprintf("%s:%.0f-%.0f", Chr, from, to)]
    data.table::setkey(regions, Chr, from, to)
    rj <- data.table::foverlaps(qx, regions[, .(Chr, from, to, region_id)],
                                by.x = c("Chr", "from", "to"), type = "within",
                                mult = "first", nomatch = NA)
    region_id <- rj$region_id
  } else region_id <- rep(NA_character_, nrow(mp))

  out <- data.table::data.table(marker = mp$marker, Chr = mp$Chr, Pos = mp$Pos,
                                statistic = statistic, unit_id = uj$unit_id, region_id = region_id)
  if (statistic == "simes") {
    if (length(p_obs) != nrow(map))
      stop("statistic = \"simes\": p_obs must have one value per marker, aligned to map (",
           length(p_obs), " vs ", nrow(map), ").")
    out[, p := as.numeric(p_obs)]
    out[, q := stats::p.adjust(p, method = "BH")]
    out[, significant := uj$significant]   # cluster-test significance, marker's own p/q
  } else {
    out[, p := uj$p][, q := uj$q][, significant := uj$significant]
  }
  out[, tested := !is.na(unit_id)]
  data.table::setcolorder(out, c("marker", "Chr", "Pos", "statistic", "p", "q",
                                 "unit_id", "region_id", "significant", "tested"))
  out[]
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
run_stage_B <- function(dataset_dir, engine, p_obs = NULL, p_perm = NULL,
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
  params <- list(statistic = statistic, size_floor = cfg$size_floor, alpha = cfg$alpha,
                 assembly = cfg$assembly, score_threshold = cfg$score_threshold,
                 distance_threshold = cfg$distance_threshold, gap = cfg$gap,
                 has_perm = !is.null(p_perm), has_rotation = !is.null(annotation),
                 n_rotations = cfg$n_rotations, rotation_scheme = cfg$rotation_scheme)
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

  snp_results <- .build_snp_results(map, test, p_obs, statistic)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(test, file.path(out_dir, "outlier_test.rds"))
  if (!is.null(perm)) saveRDS(perm, file.path(out_dir, "outlier_perm.rds"))
  if (!is.null(rotation)) saveRDS(rotation, file.path(out_dir, "region_rotation.rds"))
  data.table::fwrite(snp_results, file.path(out_dir, "snp_results.csv"))
  data.table::fwrite(test$regions, file.path(out_dir, "region_table.csv"))
  write_receipt(out_dir, inputs = inputs, params = params,
                outputs = file.path(out_dir, c("outlier_test.rds", "snp_results.csv", "region_table.csv")))

  list(test = test, perm = perm, rotation = rotation, snp_results = snp_results, out_dir = out_dir)
}
