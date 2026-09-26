## examples/9sp_chr19_chr20/input/perm_fun.R
##
## Population-level permutation, held within structure_group.rds (lineage):
## same reasoning as examples/3sp_chr1_chr4/input/perm_fun.R -- ecotype is
## constant within every population in this panel too, so the shuffle
## reassigns WHICH population gets which value, restricted to populations
## sharing the same lineage, never mixing values within a population.
## Mirrors LDscnR-paper/module_9sp/R/03_EMMAX.R's within-lineage scheme
## (population-level, not this repo's covariate-adjustment step -- LDscnR-
## multi has no Covar mechanism; see R/03_stage_A_emmax.R).
population <- readRDS(file.path(.dataset_dir, "input", "population.rds"))
structure_group <- readRDS(file.path(.dataset_dir, "input", "structure_group.rds"))

perm_fun <- function(b, y) {
  pop <- population[names(y)]
  grp <- structure_group[names(y)]
  pops <- unique(pop)
  pop_group <- grp[match(pops, pop)]
  y_pop <- as.numeric(y[match(pops, pop)])
  y_pop_shuf <- stats::ave(y_pop, pop_group, FUN = sample)
  as.numeric(y_pop_shuf[match(pop, pops)])
}
