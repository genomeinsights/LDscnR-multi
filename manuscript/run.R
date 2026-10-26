## =============================================================================
## LDscnR-multi/manuscript/run.R
##
## Runs the LDscnR-multi pipeline (Stage 1, Stage A, Stage B) on the two
## full-genome panels built by manuscript/build_inputs.R, with the package
## defaults, which match the manuscript modules' settings. For 3sp, every
## engine's regions are also tested against the EcoPeaks (ld_region_rotation(),
## within-chromosome scheme, 10,000 relocations, seed 1), as in
## LDscnR-paper/module_3sp/R/03_EMMAX.R and 04_lfmm.R.
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
  say("    [%s] done in %.1f min\n", p, as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
