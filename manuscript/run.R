## =============================================================================
## LDscnR-multi/manuscript/run.R
##
## Runs the LDscnR-multi pipeline (Stage 1, Stage A, Stage B and the floor
## profile) on the two full-genome panels built by manuscript/build_inputs.R,
## with the package defaults, which match the manuscript modules' settings.
## For 3sp, every engine's regions are also tested against the EcoPeaks
## (ld_region_rotation(), within-chromosome scheme, 10,000 relocations,
## seed 1), as in LDscnR-paper/module_3sp/R/03_EMMAX.R and 04_lfmm.R, and so
## are the floor-stability tiers, as in module_3sp/R/14_floor_tiers.R.
##
## Run from the repository root:  Rscript manuscript/run.R [3sp_full 9sp_full]
## Then compare with the manuscript:  Rscript manuscript/compare.R
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
source(file.path("R", "00_config.R"))
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)
invisible(check_ldscnr())

panels <- commandArgs(trailingOnly = TRUE)
if (!length(panels)) panels <- c("3sp_full", "9sp_full")

for (p in panels) {
  d <- file.path("manuscript", "data", p)
  stopifnot("run manuscript/build_inputs.R first" = dir.exists(d))
  say("\n########## %s ##########\n", p)
  t0 <- Sys.time()
  cfg <- resolve_config(d)
  say("    size_floor = %d (marker-density rule)\n", cfg$size_floor)
  s1 <- build_stage1(d, cfg)
  say("    [stage1] %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins")))
  run_stage_A(d, cfg)
  say("    [stage A] %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins")))
  ann_f <- file.path(d, "annotation.rds")
  ann <- if (file.exists(ann_f)) readRDS(ann_f) else list(annotation = NULL, chrom_lengths = NULL)
  engines <- c("emmax_unit", "emmax_simes", names(list_external_pvalues(d, s1$map)))
  for (eng in engines)
    run_stage_B(d, eng, cfg = cfg, annotation = ann$annotation, chrom_lengths = ann$chrom_lengths)
  fp <- run_floor_profile(d, cfg = cfg, annotation = ann$annotation, chrom_lengths = ann$chrom_lengths,
                          canonical_unsuffixed = TRUE)
  say("    [floor profile] %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins")))
  if (!is.null(ann$annotation)) {
    ## EcoPeak enrichment per floor-stability tier (module_3sp/R/14_floor_tiers.R)
    tier_rot <- function(reg, set) {
      if (!nrow(reg)) return(data.table(set = set, n = 0L))
      rr <- ld_region_rotation(reg[, .(Chr, from, to)], ann$annotation, ann$chrom_lengths,
                               scheme = cfg$rotation_scheme, n_rotations = cfg$n_rotations, seed = cfg$seed)
      data.table(set = set, n = nrow(reg), ecopeak = rr$observed, null_mean = rr$null_mean, fold = rr$fold, p = rr$p)
    }
    eco <- rbindlist(lapply(names(fp$regions), function(st) {
      ## the relocation draws depend on region order, so use the outlier
      ## test's own order (as module_3sp/R/14_floor_tiers.R does), not the
      ## profile table's order by tier
      tt <- readRDS(file.path(d, "output", paste0("stageB_emmax_", st), "outlier_test.rds"))$regions
      reg <- fp$regions[[st]][match(tt[, sprintf("%s:%d-%d", Chr, as.integer(from), as.integer(to))], region_id)]
      stopifnot(!anyNA(reg$region_id))
      rbindlist(list(tier_rot(reg, "all"), tier_rot(reg[n_ticks >= 2], ">=2/3"), tier_rot(reg[n_ticks == 3], "3/3"),
                     tier_rot(reg[n_ticks == 2], "2/3 only"), tier_rot(reg[n_ticks == 1], "1/3 only")),
                fill = TRUE)[, statistic := st]
    }), fill = TRUE)
    fwrite(eco, file.path(fp$out_dir, "ecopeak_by_tier.csv"))
  }
  say("    [%s] done in %.1f min\n", p, as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
