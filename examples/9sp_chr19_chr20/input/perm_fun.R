## examples/9sp_chr19_chr20/input/perm_fun.R
##
## Population-level permutation, held within structure_group.rds (lineage):
## same reasoning as examples/3sp_chr1_chr4/input/perm_fun.R -- ecotype is
## constant within every population in this panel too, so the shuffle
## reassigns WHICH population gets which value, restricted to populations
## sharing the same lineage, never mixing values within a population.
## Mirrors LDscnR-paper/module_9sp/R/03_EMMAX.R's perm_lineage(). The
## surrogate is then residualised on input/covariate.rds (lineage) by Stage A,
## exactly as the observed phenotype is.
population <- readRDS(file.path(.dataset_dir, "input", "population.rds"))
structure_group <- readRDS(file.path(.dataset_dir, "input", "structure_group.rds"))

## perm_within_groups() (R/01_dataset_io.R) visits populations and groups in
## order of first appearance, so the same seed gives the same surrogate as the
## manuscript modules' data.table shuffle.
perm_fun <- function(b, y) perm_within_groups(y, population, structure_group)
