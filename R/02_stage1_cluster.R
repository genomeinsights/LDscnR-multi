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
  el_dir     <- file.path(cache_dir, "edge_lists")

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

  dir.create(el_dir, recursive = TRUE, showWarnings = FALSE)
  ## Trailing separator, deliberately: compute_LD_decay(el_data_folder=)
  ## writes each chromosome's edge list via plain string concatenation
  ## (`paste0(el_data_folder, ch, ".el")`, R/compute_ld_structure.R), not
  ## file.path() -- passed a bare directory it writes a SIBLING file named
  ## by gluing the directory's own name onto the chromosome
  ## (".../cache/edge_listsChr1.el", leaving ".../cache/edge_lists/" empty),
  ## not a file inside it. A trailing "/" makes the same concatenation land
  ## correctly inside the directory instead.
  decay_args <- utils::modifyList(cfg$decay_args,
                                  list(gds = gds, el_data_folder = paste0(el_dir, "/"), seed = cfg$seed))
  LD_decay <- do.call(compute_LD_decay, decay_args)

  say("    [stage1] ld_w (rho = 0.95) per marker\n")
  map <- data.table::copy(d$map)
  ldw <- compute_ld_w(LD_decay, rho = 0.95, cores = cfg$cores)   # single rho -> plain named vector
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
  GRM <- SNPRelate::snpgdsGRM(gds, snp.id = grm_markers, method = cfg$grm_method,
                              autosome.only = FALSE, verbose = FALSE)$grm

  saveRDS(list(stage1 = stage1, LD_decay = LD_decay, map = map, GRM = GRM), stage1_rds)
  write_receipt(cache_dir, inputs = inputs, params = params, outputs = stage1_rds)
  list(gds_path = gds_path, LD_decay = LD_decay, stage1 = stage1, map = map,
      genotypes = d$genotypes, GRM = GRM)
}
