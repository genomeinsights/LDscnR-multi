## =============================================================================
## LDscnR-multi/R/02_stage1_cluster.R
##
## GENOTYPE-ONLY LD CLUSTERING, SHARED BY STAGE A AND STAGE B. No phenotype
## anywhere in this file -- `ld_complexity_reduction()` clusters markers from
## genotypes alone, and both `run_stage_A()` (GRM + unit matrix) and
## `run_stage_B()` (ld_outlier_test/perm) need the same `stage1` object, so it
## is built once here and cached, not recomputed by whichever stage runs first.
##
## Also computes `ld_w_095` and attaches it to `map` before clustering, because
## `ld_outlier_test(assembly = "stage2_discovered")` hardcodes a lookup of a
## column BY THAT NAME on `stage1$map_snp` (see R/ld_prune_and_eMLG.R's
## `ld_w_col = "ld_w_095"` inside R/ld_outlier_test.R) -- not optional, and not
## documented anywhere a first-time caller would find it before hitting the
## error. Getting this wrong is exactly the kind of thing this repo exists to
## get right once instead of per dataset.
## =============================================================================

## Build (or load, if unchanged) a dataset's stage-1 LD clustering.
##
## @param dataset_dir Path to one dataset folder.
## @param cfg Resolved config (see resolve_config()); only `cr_rho`,
##   `decay_args`, `seed`, `cores` are read.
## @param force Rebuild even if the receipt says nothing changed.
## @return list(gds_path, LD_decay, stage1, map, genotypes) -- `map` here
##   carries the `ld_w_095` column stage1 was built with; `gds_path` is a
##   closed file path (reopen with SNPRelate::snpgdsOpen() if needed).
build_stage1 <- function(dataset_dir, cfg = resolve_config(dataset_dir), force = FALSE) {
  cache_dir <- file.path(dataset_dir, "cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  gds_path   <- file.path(cache_dir, "genotypes.gds")
  stage1_rds <- file.path(cache_dir, "stage1.rds")
  el_dir     <- file.path(cache_dir, "edge_lists")

  d <- read_dataset(dataset_dir, need_phenotype = FALSE)
  inputs <- c(file.path(dataset_dir, "input", "genotypes.rds"),
              file.path(dataset_dir, "input", "map.rds"))
  params <- list(cr_rho = cfg$cr_rho, decay_args = cfg$decay_args, seed = cfg$seed)

  if (!force && !stage_stale(cache_dir, inputs, params, label = "stage1")) {
    x <- readRDS(stage1_rds)
    return(list(gds_path = gds_path, LD_decay = x$LD_decay, stage1 = x$stage1,
                map = x$map, genotypes = d$genotypes))
  }

  say("    [stage1] %d individuals x %s markers, rho = %.2f\n",
      nrow(d$genotypes), format(ncol(d$genotypes), big.mark = ","), cfg$cr_rho)

  if (file.exists(gds_path)) file.remove(gds_path)
  gds <- create_gds_from_geno(d$genotypes, d$map, gds_path)
  on.exit(SNPRelate::snpgdsClose(gds), add = TRUE)

  dir.create(el_dir, recursive = TRUE, showWarnings = FALSE)
  decay_args <- utils::modifyList(cfg$decay_args,
                                  list(gds = gds, el_data_folder = el_dir, seed = cfg$seed))
  LD_decay <- do.call(compute_LD_decay, decay_args)

  ## ld_w_095: required by ld_outlier_test(assembly = "stage2_discovered"); see
  ## the file header. compute_ld_w()'s return is a plain named vector when a
  ## single rho was requested, or a matrix (columns "rho_<value>") for more
  ## than one -- handle both rather than assume decay_args$ld_w_rho's length.
  map <- data.table::copy(d$map)
  ldw <- LD_decay$ld_ws
  map$ld_w_095 <- if (is.matrix(ldw)) ldw[map$marker, "rho_0.95"] else ldw[map$marker]

  stage1 <- ld_complexity_reduction(map = map, LD_decay = LD_decay, rho = cfg$cr_rho,
                                    cores = cfg$cores, gds = gds_path)

  saveRDS(list(stage1 = stage1, LD_decay = LD_decay, map = map), stage1_rds)
  write_receipt(cache_dir, inputs = inputs, params = params, outputs = stage1_rds)
  list(gds_path = gds_path, LD_decay = LD_decay, stage1 = stage1, map = map, genotypes = d$genotypes)
}
