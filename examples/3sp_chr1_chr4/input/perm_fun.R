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

## perm_within_groups() (R/01_dataset_io.R) visits populations and groups in
## order of first appearance, so the same seed gives the same surrogate as the
## manuscript modules' data.table shuffle.
perm_fun <- function(b, y) perm_within_groups(y, population, structure_group)
