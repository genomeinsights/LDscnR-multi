## =============================================================================
## LDscnR-multi/R/02_stage1_cluster.R
##
## GENOTYPE-ONLY LD CLUSTERING, SHARED BY STAGE A, STAGE B AND THE STRUCTURE-
## ALIGNMENT DIAGNOSTIC. No phenotype anywhere in this file --
## `ld_complexity_reduction()` clusters markers from genotypes alone, and the
## GRM built from its representatives (`SNPRelate::snpgdsGRM()`) is likewise
## phenotype-free: both `run_stage_A()` (unit matrix + EMMAX) and
## `check_structure_alignment()` (R/08_structure_alignment.R) need the same
## `stage1`/`GRM`, so both are built once here and cached, not recomputed by
## whichever caller runs first. (GRM used to be built inside `run_stage_A()`
## itself; hoisted here once the alignment diagnostic needed it too, since a
## user may want that diagnostic *before* ever running Stage A -- exactly how
## the manuscript's own pre-scan framing uses it.)
##
## Also computes `ld_w_095` (local-LD support at rho = 0.95) and attaches it
## to `map`. Until LDscnR commit 5f12cd2 this was REQUIRED --
## `ld_outlier_test(assembly = "stage2_discovered")` hardcoded a lookup of a
## column by that exact name; that commit replaced it with `n_loci` (always
## present, always positive -- a structural no-op for the same purpose), so
## it is no longer needed for testing at all. Kept/recomputed anyway as a
## genuinely useful diagnostic in its own right: a manhattan-style plot of
## `ld_w_095` (`ldm_manhattan(..., value = "ld_w_095")`) shows where local LD
## support is high (low recombination, structural variants, ...)
## independent of any association test, so a user can compare a significant
## region's location against local LD directly.
## =============================================================================

## Build (or load, if unchanged) a dataset's stage-1 LD clustering.
##
## @param dataset_dir Path to one dataset folder.
## @param cfg Resolved config (see resolve_config()); reads `cr_rho`,
##   `decay_args`, `grm_method`, `seed`, `cores`.
## @param force Rebuild even if the receipt says nothing changed.
## @return list(gds_path, LD_decay, stage1, map, genotypes, GRM) -- `gds_path`
##   is a closed file path (reopen with SNPRelate::snpgdsOpen() if needed).
build_stage1 <- function(dataset_dir, cfg = resolve_config(dataset_dir), force = FALSE) {
  cache_dir <- file.path(dataset_dir, "cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  gds_path   <- file.path(cache_dir, "genotypes.gds")
  stage1_rds <- file.path(cache_dir, "stage1.rds")

  d <- read_dataset(dataset_dir, need_phenotype = FALSE)
  inputs <- c(file.path(dataset_dir, "input", "genotypes.rds"),
              file.path(dataset_dir, "input", "map.rds"))
  params <- list(cr_rho = cfg$cr_rho, decay_args = cfg$decay_args,
                grm_method = cfg$grm_method, seed = cfg$seed)

  if (!force && !stage_stale(cache_dir, inputs, params, label = "stage1")) {
    x <- readRDS(stage1_rds)
    return(list(gds_path = gds_path, LD_decay = x$LD_decay, stage1 = x$stage1,
                map = x$map, genotypes = d$genotypes, GRM = x$GRM))
  }

  say("    [stage1] %d individuals x %s markers, rho = %.2f\n",
      nrow(d$genotypes), format(ncol(d$genotypes), big.mark = ","), cfg$cr_rho)

  if (file.exists(gds_path)) file.remove(gds_path)
  gds <- create_gds_from_geno(d$genotypes, d$map, gds_path)
  on.exit(SNPRelate::snpgdsClose(gds), add = TRUE)

  ## Deliberately NOT passing el_data_folder (nor keep_el = TRUE): at full-
  ## genome marker density, one chromosome's unfiltered edge list (every SNP
  ## against its next `slide` neighbours, no r2 floor) runs 4-5 GiB, and
  ## el_data_folder would persist every chromosome's simultaneously (compute_ld_w()
  ## below needs every chromosome's a_pred, so nothing is deleted until all are
  ## fitted) -- tens of GiB of temp disk for one dataset. Leaving `el` unset
  ## instead makes every downstream consumer rebuild it from `gds`, one
  ## chromosome at a time, discarded immediately after use:
  ##   - ld_complexity_reduction() below already has this fallback
  ##     (LDscnR's .chr_edge_list()), and for THAT caller specifically the
  ##     on-the-fly rebuild is floor-filtered at construction (el_floor =
  ##     the clustering r2 threshold) -- LDscnR's own comment on that path:
  ##     "the difference between ~15M rows and the few that clear the
  ##     threshold." Genuinely cheaper than reading a saved, unfiltered file.
  ##   - compute_ld_w() just below is given `gds` explicitly for the same
  ##     on-the-fly rebuild (it has no fallback of its own the way
  ##     ld_complexity_reduction() does).
  ## Clustering itself is unaffected by losing "long-range" edges beyond
  ## `slide`: ld_complexity_reduction() single-links markers by connected
  ## components (igraph::components()) on the thresholded graph, so a
  ## genuinely extended LD block (an inversion) still merges into one
  ## component through a chain of locally-adjacent above-threshold edges --
  ## no single pair at opposite ends of the block ever needs to be directly
  ## compared. See the LD-decay report (R/10_ld_decay_report.R) for tuning
  ## `slide` itself, the other lever on `el` size.
  decay_args <- utils::modifyList(cfg$decay_args, list(gds = gds, seed = cfg$seed))
  LD_decay <- do.call(compute_LD_decay, decay_args)

  say("    [stage1] ld_w (rho = 0.95) per marker\n")
  map <- data.table::copy(d$map)
  ## cores = 1, not cfg$cores: compute_ld_w() has no reopen-per-forked-worker
  ## logic (unlike ld_complexity_reduction()'s reopen_path below), so handing
  ## it this call's own already-open `gds` handle under parallel_apply's
  ## forking would share one live GDS connection across worker processes --
  ## exactly what ld_complexity_reduction() explicitly avoids doing. Forcing
  ## serial execution here keeps this correct rather than merely usually
  ## working; a smaller, calibrated `slide` keeps each chromosome's rebuild
  ## cheap enough that this doesn't dominate runtime.
  ldw <- compute_ld_w(LD_decay, rho = 0.95, cores = 1, gds = gds)   # single rho -> plain named vector
  map$ld_w_095 <- as.numeric(ldw[map$marker])

  stage1 <- ld_complexity_reduction(map = map, LD_decay = LD_decay, rho = cfg$cr_rho,
                                    cores = cfg$cores, gds = gds_path)

  grm_markers <- unique(stats::na.omit(stage1$pruned))
  say("    [stage1] GRM: %s markers (stage-1 representatives), method = %s\n",
      format(length(grm_markers), big.mark = ","), cfg$grm_method)
  ## autosome.only = FALSE: see R/03_stage_A_emmax.R's original note (moved
  ## here with the GRM build itself) -- SNPRelate::snpgdsGRM() otherwise
  ## silently excludes every marker whenever chromosome labels aren't in its
  ## human-centric numeric range.
  ##
  ## missing.rate = 1 (SNPRelate's own default is 0.01): the default excludes
  ## any marker with MORE than 1% of individuals missing a call -- with any
  ## genuinely missing genotypes (this pipeline's own real datasets have
  ## never been zero-missing, e.g. formica_hybrid's ~5.9%), and with n in the
  ## hundreds, the probability that a GIVEN marker has NO missing calls at
  ## all approaches zero, so the default filter silently discards nearly
  ## EVERY marker, not a handful -- confirmed directly: n=150, 6% missing,
  ## 2,000 candidate markers -> 1,998 excluded, GRM built from 2 markers,
  ## cor with the true (complete-data) GRM = NA. This is NOT a computation
  ## bug in snpgdsGRM() itself: with missing.rate = 1 (exclude nothing),
  ## snpgdsGRM()'s own GCTA formula on the SAME missing data matches an
  ## independently-implemented, pairwise-complete-observations R
  ## reimplementation (~/gitlab/pike_phenology/R/grm_functions.R's
  ## gcta_grm(), verified there to reproduce snpgdsGRM(GCTA) to machine
  ## precision on complete data) to machine precision too -- the bug is
  ## entirely the default marker-exclusion threshold, not the arithmetic.
  ## Diagnosed against real missing-genotype data 2026-09-30 (see also
  ## LDscnR-multi's own missing-data audit, same date).
  GRM_obj <- SNPRelate::snpgdsGRM(gds, snp.id = grm_markers, method = cfg$grm_method,
                                  autosome.only = FALSE, missing.rate = 1, verbose = FALSE)
  ## Self-verifying, not just trusting the parameter above forever: if a
  ## future SNPRelate version changes what missing.rate = 1 means, or this
  ## call is ever edited without noticing why, silently losing markers here
  ## again should be a loud, immediate error, not a repeat of the silent
  ## corruption this fix exists to close.
  if (length(GRM_obj$snp.id) != length(grm_markers))
    stop(sprintf(
      "snpgdsGRM() used %d of %d requested GRM markers -- missing.rate = 1 ",
      length(GRM_obj$snp.id), length(grm_markers)),
      "should make every requested marker contribute regardless of missing ",
      "genotypes; investigate before trusting this GRM.")
  GRM <- GRM_obj$grm

  saveRDS(list(stage1 = stage1, LD_decay = LD_decay, map = map, GRM = GRM), stage1_rds)
  write_receipt(cache_dir, inputs = inputs, params = params, outputs = stage1_rds)
  list(gds_path = gds_path, LD_decay = LD_decay, stage1 = stage1, map = map,
      genotypes = d$genotypes, GRM = GRM)
}
