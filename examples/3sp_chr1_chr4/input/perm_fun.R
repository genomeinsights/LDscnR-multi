## examples/3sp_chr1_chr4/input/perm_fun.R
##
## Population-level permutation, held within structure_group.rds
## (pop_locality): ecotype is a population-level trait in this panel, not an
## individual one -- every individual in a population shares one value. The
## shuffle therefore reassigns WHICH population gets which value, restricted
## to populations sharing the same locality, and never mixes values within a
## population (an individual-level shuffle within locality does exactly
## that -- verified directly: after one such draw, ecotype varied within
## 14/35 populations, which never happens in the real data). Mirrors
## LDscnR-paper/module_3sp/R/03_EMMAX.R's own perm_regional().
population <- readRDS(file.path(.dataset_dir, "input", "population.rds"))
structure_group <- readRDS(file.path(.dataset_dir, "input", "structure_group.rds"))

perm_fun <- function(b, y) {
  pop <- population[names(y)]
  grp <- structure_group[names(y)]
  pops <- unique(pop)
  pop_group <- grp[match(pops, pop)]
  y_pop <- as.numeric(y[match(pops, pop)])          # one value per population (constant within)
  y_pop_shuf <- stats::ave(y_pop, pop_group, FUN = sample)   # shuffled within locality, at population resolution
  as.numeric(y_pop_shuf[match(pop, pops)])          # broadcast back to individuals
}
