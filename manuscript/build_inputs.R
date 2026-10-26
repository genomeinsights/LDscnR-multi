## =============================================================================
## LDscnR-multi/manuscript/build_inputs.R
##
## Builds the two full-genome dataset folders used to check that LDscnR-multi
## reproduces the manuscript's empirical analyses (manuscript/run.R,
## manuscript/compare.R). Inputs are taken from the manuscript modules'
## Stage-1 bundles in LDscnR-paper (module_3sp, module_9sp: out/02_bundle/
## bundle.rds), i.e. the same individuals, MAF > 0.1 marker set, map and
## phenotype the manuscript analysed. Nothing in LDscnR-paper is written.
##
##   manuscript/data/3sp_full/   EMMAX + LFMM (external p-values); EcoPeaks
##                               annotation for the relocation test
##   manuscript/data/9sp_full/   EMMAX with ecotype residualised on lineage
##                               (input/covariate.rds)
##
## manuscript/data/ is gitignored (hundreds of MB per panel).
## Run from the repository root:  Rscript manuscript/build_inputs.R
## =============================================================================
suppressMessages(library(data.table))
HERE  <- file.path(getwd(), "manuscript", "data")
PAPER <- path.expand(Sys.getenv("LDSCNR_PAPER", "~/gitlab/LDscnR-paper"))
stopifnot(dir.exists(PAPER))

write_panel <- function(out_dir, b, group_col, covariate_col = NULL) {
  GTs <- b$GTs
  ids <- paste0("ind_", seq_len(nrow(GTs)))
  rownames(GTs) <- ids
  stopifnot(identical(colnames(GTs), b$map$marker), nrow(b$pheno) == nrow(GTs))
  in_dir <- file.path(out_dir, "input")
  dir.create(in_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(GTs, file.path(in_dir, "genotypes.rds"))
  saveRDS(b$map[, .(marker, Chr, Pos)], file.path(in_dir, "map.rds"))
  saveRDS(setNames(as.numeric(b$eco), ids), file.path(in_dir, "phenotype.rds"))
  saveRDS(setNames(as.character(b$pheno$pop_ID), ids), file.path(in_dir, "population.rds"))
  saveRDS(setNames(as.character(b$pheno[[group_col]]), ids), file.path(in_dir, "structure_group.rds"))
  if (!is.null(covariate_col))
    saveRDS(setNames(as.character(b$pheno[[covariate_col]]), ids), file.path(in_dir, "covariate.rds"))
  writeLines(c(
    "## Population-level permutation within structure_group.rds, in the manuscript",
    "## modules' group order (see perm_within_groups(), R/01_dataset_io.R).",
    "population <- readRDS(file.path(.dataset_dir, \"input\", \"population.rds\"))",
    "structure_group <- readRDS(file.path(.dataset_dir, \"input\", \"structure_group.rds\"))",
    "perm_fun <- function(b, y) perm_within_groups(y, population, structure_group)"),
    file.path(in_dir, "perm_fun.R"))
  writeLines(sprintf("cores <- %dL", as.integer(Sys.getenv("LDSCNR_CORES", "4"))),
             file.path(in_dir, "config.R"))
  cat(sprintf("  wrote %s: %d individuals x %s markers\n", out_dir, nrow(GTs),
              format(ncol(GTs), big.mark = ",")))
}

## ---- three-spined ---------------------------------------------------------------
cat("[1] 3sp\n")
cfg3 <- new.env()
sys.source(file.path(PAPER, "module_3sp", "R", "00_config.R"), envir = cfg3)
b3 <- readRDS(file.path(PAPER, "module_3sp", "out", "02_bundle", "bundle.rds"))
d3 <- file.path(HERE, "3sp_full")
write_panel(d3, b3, group_col = "pop_locality")      # module_3sp perm_regional()

## LFMM: the manuscript's own precomputed F statistics (module_3sp/R/04_lfmm.R)
e <- new.env(); load(cfg3$PATHS$raw_3sp, envir = e)
keep <- e$map_3sp$maf > cfg3$MAF_KEEP
lfmm_F <- readRDS(cfg3$PATHS$lfmm_F)[keep]
stopifnot(identical(e$map_3sp$marker[keep], b3$map$marker))
lfmm_p <- setNames(stats::pf(lfmm_F, 1, nrow(b3$GTs) - 2, lower.tail = FALSE), b3$map$marker)
dir.create(file.path(d3, "external_pvalues"), showWarnings = FALSE)
saveRDS(lfmm_p, file.path(d3, "external_pvalues", "lfmm.rds"))

## EcoPeaks (Kingman et al. 2021, lifted to gasAcu1), read as in module_3sp/R/03_EMMAX.R
LIFT <- file.path(PAPER, "kingman2021", "data", "liftover")
ROMAN <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV",
           "XV","XVI","XVII","XVIII","XIX","XX","XXI")
.rd <- function(f) { x <- fread(file.path(LIFT, f), header = FALSE, col.names = c("chr","start","end","pv"))
  x[, chr_num := match(sub("^chr", "", chr), ROMAN)]
  x[!is.na(chr_num), .(chr = paste0("Chr", chr_num), start, end)] }
KG <- unique(rbindlist(lapply(cfg3$ECOPEAK_BEDS, .rd)))
LEN <- b3$map[, .(len = max(Pos)), by = Chr]; setnames(LEN, "Chr", "chr")
saveRDS(list(annotation = KG, chrom_lengths = LEN), file.path(d3, "annotation.rds"))
cat(sprintf("  EcoPeaks: %d intervals\n", nrow(KG)))
rm(b3, e); invisible(gc())

## ---- nine-spined ----------------------------------------------------------------
cat("[2] 9sp\n")
cfg9 <- new.env()
sys.source(file.path(PAPER, "module_9sp", "R", "00_config.R"), envir = cfg9)
b9 <- readRDS(file.path(PAPER, "module_9sp", "out", "02_bundle", "bundle.rds"))
write_panel(file.path(HERE, "9sp_full"), b9, group_col = cfg9$EMMAX_COVAR,
            covariate_col = cfg9$EMMAX_COVAR)              # module_9sp perm_lineage() + residualisation
cat("done.\n")
