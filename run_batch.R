## =============================================================================
## LDscnR-multi/run_batch.R
##
## Rscript entry point: source R/*.R in order, then run_all() over every
## dataset folder passed on the command line (or every subfolder of a single
## root folder, if exactly one argument is given and it is not itself a
## dataset folder).
##
##   Rscript run_batch.R datasets/panelA datasets/panelB ...
##   Rscript run_batch.R datasets/                 # every subfolder of datasets/
##
## Interactive use doesn't need this file -- source R/*.R yourself and call
## run_dataset()/run_all() directly (see README.md).
## =============================================================================
HERE <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
if (!length(HERE) || !nzchar(HERE)) HERE <- "."
source(file.path(HERE, "R", "00_config.R"))   # defines LDSCNR_PIN, needed before loading LDscnR

## Use the installed GitHub package for public runs. Developers can opt into
## their local checkout with LDSCNR_DEV_LOAD_ALL=1; check_ldscnr() validates
## the package actually loaded in either mode before the first dataset runs.
if (identical(Sys.getenv("LDSCNR_DEV_LOAD_ALL"), "1")) {
  suppressMessages(devtools::load_all(LDSCNR_PIN$repo, quiet = TRUE))
} else {
  suppressPackageStartupMessages(library(LDscnR))
}

for (f in sort(list.files(file.path(HERE, "R"), pattern = "\\.R$", full.names = TRUE))) source(f)

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: Rscript run_batch.R <dataset_dir> [<dataset_dir> ...] | <root_folder>/")

is_dataset <- function(d) dir.exists(file.path(d, "input"))
dataset_dirs <- if (length(args) == 1L && dir.exists(args[1]) && !is_dataset(args[1])) {
  sub <- list.dirs(args[1], recursive = FALSE)
  Filter(is_dataset, sub)
} else args

if (!length(dataset_dirs)) stop("No dataset folders found (each needs an input/ subfolder).")
say("Running %d dataset(s):\n", length(dataset_dirs))
for (d in dataset_dirs) say("  - %s\n", d)

run_all(dataset_dirs, summary_path = file.path(HERE, "summary.csv"))
