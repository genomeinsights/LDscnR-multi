## =============================================================================
## LDscnR-multi/R/01_dataset_io.R
##
## THE DATASET-FOLDER CONTRACT, IN ONE PLACE. Every other stage reads a
## dataset through read_dataset() rather than touching input/*.rds directly, so
## the contract (what a folder must contain, what shape each file is) is
## checked once, loudly, before anything expensive runs -- same reasoning as
## module_3sp's 01_inputs.R: a missing/malformed input costs seconds here, and
## whatever stage first happens to touch it otherwise.
##
## Folder contract (see README.md for the full description):
##   <dataset_dir>/input/genotypes.rds   n x m numeric dosage matrix,
##                                       individuals (rows, rownames = IDs) x
##                                       markers (cols, colnames = marker IDs)
##   <dataset_dir>/input/map.rds          data.table/data.frame: marker, Chr, Pos
##                                       -- marker order must equal colnames(genotypes)
##   <dataset_dir>/input/phenotype.rds     optional (Stage A only): named numeric
##                                       vector (names = rownames(genotypes)) or
##                                       a 1-column data.frame with those rownames
##   <dataset_dir>/input/config.R           optional: DEFAULTS overrides
##   <dataset_dir>/input/perm_fun.R          REQUIRED for Stage A (run_stage_A()
##                                       errors without it). Defines
##                                       perm_fun(b, y): b is the permutation
##                                       index, y the observed phenotype
##                                       (individual-level, aligned to
##                                       genotypes' rownames); must return a
##                                       permuted phenotype of the same shape.
##                                       No generic default is offered on
##                                       purpose -- an earlier version fell
##                                       back to an unstratified sample(y), and
##                                       a later structure-aware default
##                                       (shuffle within structure_group.rds at
##                                       the INDIVIDUAL level) still broke the
##                                       sampling design (an external audit
##                                       found the phenotype, constant within
##                                       every population in both worked
##                                       examples, varied within 14/35 and
##                                       9/30 populations after one such
##                                       shuffle). What "no signal" means is
##                                       specific to each dataset's own
##                                       sampling design, so this repo asks
##                                       for it explicitly, per dataset,
##                                       rather than guess -- see examples/
##                                       */input/perm_fun.R for the population-
##                                       level, structure-preserving pattern
##                                       both worked examples actually use.
##                                       `.dataset_dir` (this dataset's own
##                                       path) is available to the sourced
##                                       script if it needs to read its own
##                                       population.rds/structure_group.rds.
##   <dataset_dir>/input/population.rds       optional (required only for
##                                       check_structure_alignment(), R/08_
##                                       structure_alignment.R): named vector,
##                                       individual -> population ID.
##   <dataset_dir>/input/structure_group.rds  optional: named vector, individual
##                                       -> coarser sampling group (locality/
##                                       lineage/etc); defaults to population.rds
##                                       when absent.
##   <dataset_dir>/external_pvalues/*.rds  optional: one named p-value vector
##                                       (names = map$marker) per file; the
##                                       basename (sans .rds) becomes the engine
## =============================================================================

## Read + validate the required inputs. Does not touch cache/ or output/.
##
## @param dataset_dir Path to one dataset folder.
## @param need_phenotype If TRUE, error when input/phenotype.rds is absent
##   (Stage A needs it; stage1 clustering and Stage B do not).
## @return list(genotypes, map, phenotype (or NULL), dataset_dir)
read_dataset <- function(dataset_dir, need_phenotype = FALSE) {
  in_dir <- file.path(dataset_dir, "input")
  gt_f  <- file.path(in_dir, "genotypes.rds")
  map_f <- file.path(in_dir, "map.rds")
  if (!file.exists(gt_f))  stop("missing required input: ", gt_f)
  if (!file.exists(map_f)) stop("missing required input: ", map_f)

  genotypes <- readRDS(gt_f)
  map <- data.table::as.data.table(readRDS(map_f))

  if (!is.matrix(genotypes) || !is.numeric(genotypes))
    stop(gt_f, ": must be a numeric matrix (individuals x markers).")
  if (is.null(rownames(genotypes)))
    stop(gt_f, ": rownames (individual IDs) are required.")
  if (is.null(colnames(genotypes)))
    stop(gt_f, ": colnames (marker IDs) are required.")
  req_map_cols <- c("marker", "Chr", "Pos")
  missing_cols <- setdiff(req_map_cols, names(map))
  if (length(missing_cols))
    stop(map_f, ": missing column(s): ", paste(missing_cols, collapse = ", "))
  if (!identical(colnames(genotypes), map$marker))
    stop(map_f, ": map$marker must equal colnames(genotypes), in the same order ",
         "(this is the alignment every downstream LDscnR call assumes).")

  phenotype <- NULL
  ph_f <- file.path(in_dir, "phenotype.rds")
  if (file.exists(ph_f)) {
    phenotype <- readRDS(ph_f)
    if (is.data.frame(phenotype)) {
      if (ncol(phenotype) != 1L) stop(ph_f, ": data.frame phenotype must have exactly one column.")
      nm <- rownames(phenotype)
      phenotype <- stats::setNames(phenotype[[1]], nm)
    }
    if (is.null(names(phenotype)))
      stop(ph_f, ": phenotype must be named (names = individual IDs matching genotypes rownames).")
    if (!identical(names(phenotype), rownames(genotypes)))
      phenotype <- phenotype[rownames(genotypes)]
    if (anyNA(names(phenotype)) || length(phenotype) != nrow(genotypes))
      stop(ph_f, ": phenotype names do not cover genotypes rownames.")
  } else if (need_phenotype) {
    stop("Stage A needs ", ph_f, ", which is absent.")
  }

  list(genotypes = genotypes, map = map, phenotype = phenotype, dataset_dir = dataset_dir)
}

## perm_fun(b, y): the caller's permutation scheme. REQUIRED -- no generic
## default. Two earlier defaults were tried and both broke the sampling
## design: an unstratified sample(y), and (once that was flagged) a
## structure_group.rds-based shuffle done at the INDIVIDUAL level -- an
## external audit found that in both worked examples the phenotype is
## constant within every population, and after one such shuffle it varied
## WITHIN 14/35 (3sp) and 9/30 (9sp) populations, which never happens in the
## real data. What "no signal" means -- which unit is permuted, what
## structure is held fixed -- is specific to each dataset's own sampling
## design (locality-stratified, lineage-stratified, population-level vs.
## individual-level, ...), so this repo asks for it explicitly, once per
## dataset, rather than guess a scheme that fits some designs and silently
## breaks others. See examples/*/input/perm_fun.R for the population-level,
## structure-preserving pattern both worked examples actually use.
resolve_perm_fun <- function(dataset_dir) {
  pf_file <- file.path(dataset_dir, "input", "perm_fun.R")
  if (!file.exists(pf_file))
    stop("run_stage_A() needs ", pf_file, ", which is absent -- there is no generic default ",
         "permutation scheme (see this function's own header for why). Define ",
         "`perm_fun(b, y)` there for this dataset's own sampling design.")
  env <- new.env(parent = globalenv())
  env$.dataset_dir <- dataset_dir   # available to the sourced script, e.g. to read its own population.rds
  sys.source(pf_file, envir = env)
  if (!is.function(env$perm_fun))
    stop(pf_file, " must define a function `perm_fun(b, y)`.")
  env$perm_fun
}

## Every external p-value engine dropped in external_pvalues/*.rds. Returns a
## named list of named numeric vectors (names = map$marker), keyed by the
## engine name (the file's basename without extension).
list_external_pvalues <- function(dataset_dir, map) {
  ext_dir <- file.path(dataset_dir, "external_pvalues")
  if (!dir.exists(ext_dir)) return(list())
  files <- list.files(ext_dir, pattern = "\\.rds$", full.names = TRUE)
  out <- stats::setNames(lapply(files, function(f) {
    p <- readRDS(f)
    if (is.null(names(p))) stop(f, ": external p-values must be a NAMED numeric vector (names = map$marker).")
    missing <- setdiff(map$marker, names(p))
    if (length(missing)) warning(sprintf(
      "%s: %d marker(s) in map have no p-value (will read as NA).", f, length(missing)))
    p
  }), tools::file_path_sans_ext(basename(files)))
  out
}
