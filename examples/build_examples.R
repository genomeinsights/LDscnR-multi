## =============================================================================
## LDscnR-multi/examples/build_examples.R
##
## BUILDS THE TWO WORKED EXAMPLES from already-computed LDscnR-paper bundles,
## rather than from raw data -- both bundles (117/149 individuals x ~800k/1.2M
## markers) already exist on this machine, and re-deriving them from scratch
## belongs to those modules, not to this repo. This script only SUBSETS: two
## chromosomes each, kept small enough to run the whole LDscnR-multi pipeline
## end to end in a few minutes.
##
##   examples/3sp_chr1_chr4/     three-spine stickleback, Chr1 + Chr4
##                                (module_3sp's own headline chromosomes:
##                                Chr4 carries Eda, Chr1 the marine/freshwater
##                                inversion -- see module_3sp/R_figures/
##                                figure_manhattan.R). Also carries an
##                                external_pvalues/lfmm.rds, aligned the same
##                                way module_3sp/R/04_lfmm.R aligns it, to
##                                demonstrate Stage B on a non-EMMAX engine.
##   examples/9sp_chr19_chr20/   nine-spine stickleback, Chr19 + Chr20 -- a
##                                second species, different sample size and
##                                population structure (see module_9sp/R/
##                                00_config.R), no external p-value engine.
##
## NOT RUN as part of the pipeline itself -- run once to (re)generate the
## example inputs, checked in under examples/ (small: two chromosomes, not
## the genome).
## =============================================================================
suppressMessages(library(data.table))

HERE <- path.expand("~/gitlab/LDscnR-multi/examples")

## `pheno`/`structure_col`: optional. When given, also writes
## input/population.rds (from pheno$pop_ID) and input/structure_group.rds
## (from pheno[[structure_col]]) for the structure-alignment diagnostic --
## see R/08_structure_alignment.R. `pheno` must be in the same row order as
## GTs/eco (true of both bundles here; never re-sorted).
write_example <- function(out_dir, GTs, map, eco, chrs, pheno = NULL, structure_col = NULL) {
  keep <- map$Chr %in% chrs
  map_sub <- map[keep, .(marker, Chr, Pos)]
  GT_sub <- GTs[, keep, drop = FALSE]
  rownames(GT_sub) <- paste0("ind_", seq_len(nrow(GT_sub)))
  stopifnot(identical(colnames(GT_sub), map_sub$marker))

  in_dir <- file.path(out_dir, "input")
  dir.create(in_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(GT_sub, file.path(in_dir, "genotypes.rds"))
  saveRDS(map_sub, file.path(in_dir, "map.rds"))
  saveRDS(stats::setNames(as.numeric(eco), rownames(GT_sub)), file.path(in_dir, "phenotype.rds"))
  cat(sprintf("  wrote %s: %d individuals x %s markers (%s)\n", out_dir,
              nrow(GT_sub), format(ncol(GT_sub), big.mark = ","), paste(chrs, collapse = "+")))

  if (!is.null(pheno)) {
    stopifnot(nrow(pheno) == nrow(GT_sub))
    saveRDS(stats::setNames(as.character(pheno$pop_ID), rownames(GT_sub)),
           file.path(in_dir, "population.rds"))
    saveRDS(stats::setNames(as.character(pheno[[structure_col]]), rownames(GT_sub)),
           file.path(in_dir, "structure_group.rds"))
    cat(sprintf("  wrote population.rds (%d pops) + structure_group.rds (from %s, %d groups)\n",
                length(unique(pheno$pop_ID)), structure_col, length(unique(pheno[[structure_col]]))))
  }
  map_sub
}

## ---- 3sp: Chr1 + Chr4 ----------------------------------------------------------
cat("[1] 3sp (Chr1 + Chr4)\n")
b3 <- readRDS(path.expand("~/gitlab/LDscnR-paper/module_3sp/out/02_bundle/bundle.rds"))
map3_sub <- write_example(file.path(HERE, "3sp_chr1_chr4"), b3$GTs, b3$map, b3$eco, c("Chr1", "Chr4"),
                          pheno = b3$pheno, structure_col = "pop_locality")

## external p-value engine: LFMM, aligned to map3_sub exactly as
## module_3sp/R/04_lfmm.R aligns it to the full bundle's map (see that file's
## comment for why this two-step mask -- maf filter, then this subset -- is
## the only way to recover correct alignment; lfmm_F.rds is over the FULL
## PRE-MAF map, in map_3sp's own row order).
cfg3 <- new.env()
sys.source(path.expand("~/gitlab/LDscnR-paper/module_3sp/R/00_config.R"), envir = cfg3)
e <- new.env(); load(cfg3$PATHS$raw_3sp, envir = e)
keep_maf <- e$map_3sp$maf > cfg3$MAF_KEEP
lfmm_F_full <- readRDS(cfg3$PATHS$lfmm_F)
stopifnot(length(lfmm_F_full) == length(keep_maf))
lfmm_F <- lfmm_F_full[keep_maf]
stopifnot(identical(e$map_3sp$marker[keep_maf], b3$map$marker))
names(lfmm_F) <- b3$map$marker
lfmm_p_sub <- stats::pf(lfmm_F[map3_sub$marker], 1, nrow(b3$GTs) - 2, lower.tail = FALSE)
ext_dir <- file.path(HERE, "3sp_chr1_chr4", "external_pvalues")
dir.create(ext_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(lfmm_p_sub, file.path(ext_dir, "lfmm.rds"))
cat(sprintf("  wrote %s: %s LFMM p-values\n", file.path(ext_dir, "lfmm.rds"),
            format(length(lfmm_p_sub), big.mark = ",")))
rm(b3, e); gc()

## ---- 9sp: Chr19 + Chr20 ---------------------------------------------------------
cat("\n[2] 9sp (Chr19 + Chr20)\n")
b9 <- readRDS(path.expand("~/gitlab/LDscnR-paper/module_9sp/out/02_bundle/bundle.rds"))
write_example(file.path(HERE, "9sp_chr19_chr20"), b9$GTs, b9$map, b9$eco, c("Chr19", "Chr20"),
             pheno = b9$pheno, structure_col = "lineage")
rm(b9); gc()

cat("\ndone.\n")
