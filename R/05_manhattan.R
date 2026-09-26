## =============================================================================
## LDscnR-multi/R/05_manhattan.R
##
## FLEXIBLE, MULTI-TRACK MANHATTAN PLOTTING. A thin composition layer over
## LDscnR::ld_manhattan(), not a reimplementation of its per-chromosome
## faceting/colouring -- the same pattern already proven in
## module_3sp/R_figures/figure_manhattan.R (stack one ld_manhattan()-shaped
## panel per engine with patchwork). "Flexible" here means: any number of
## engines (any stageB_<engine> the dataset has), any of ld_manhattan()'s
## colouring modes, and no per-dataset custom plotting code required -- not a
## new plotting grammar.
##
## Two things this file adds on top of ld_manhattan() itself:
##
## - It plots `p_display`/`q_display`, not the tested `p`/`q` (see
##   R/04_stage_B_outlier.R's header) -- for the "unit" statistic this is a
##   genuine per-marker EMMAX scan, not the tested cluster's value broadcast
##   flat across every member, so a "unit" engine's panel shows real per-SNP
##   variation.
## - `colour_by = "region"` colours by a genome-wide, CROSS-ENGINE-CONSISTENT
##   set of "loci" (.master_loci(), below) rather than letting each panel's
##   own ld_manhattan() call assign colours independently -- so the same
##   physical region reads as the same colour in every stacked panel, not a
##   coincidence of that panel's own region ordering.
##
## Grey-vs-colour draw order is NOT handled here: ld_manhattan()'s own
## group/regions mode already draws the "ns" (grey) layer before the coloured
## one (R/ld_manhattan.R in the LDscnR package), so every mode this file uses
## inherits that ordering for free -- significant/coloured points are always
## the top layer and are never covered by the (now much denser, since
## p_display fills in a real value for markers the "unit" test itself never
## individually assessed) grey background.
## =============================================================================

## Build a genome-wide, cross-engine set of "loci": the union of every
## requested engine's significant regions, merged where they overlap
## physically. Different engines' own stage-2 assembly rarely produces
## identical bounds for "the same" region, but the same underlying signal
## typically does overlap -- merging on overlap, not on an exact match, is
## what lets a shared colour survive that. One persistent colour per locus
## (assigned by the caller) is what makes a region read as the same colour in
## every panel.
.master_loci <- function(dataset_dir, engines) {
  regs <- data.table::rbindlist(lapply(engines, function(eng) {
    f <- file.path(dataset_dir, "output", paste0("stageB_", eng), "region_table.csv")
    if (!file.exists(f)) return(NULL)
    r <- data.table::fread(f)
    if (!nrow(r)) return(NULL)
    r[, .(Chr, from, to)]
  }), use.names = TRUE, fill = TRUE)
  if (is.null(regs) || !nrow(regs))
    return(data.table::data.table(Chr = character(), from = numeric(), to = numeric(), locus_id = character()))
  data.table::setorder(regs, Chr, from)
  ## Chromosome-safe interval merge (by = Chr): a global cummax(to) that never
  ## resets at a chromosome boundary silently over-merges every chromosome
  ## after the first sizeable one -- see LDscnR's own R/ld_outlier_internal.R
  ## header for the exact failure this `by = Chr` avoids.
  regs[, grp := cumsum(c(TRUE, from[-1] - cummax(to)[-.N] > 0)), by = Chr]
  loci <- regs[, .(from = min(from), to = max(to)), by = .(Chr, grp)][, grp := NULL]
  loci[, locus_id := sprintf("%s:%.0f-%.0f", Chr, from, to)]
  loci[]
}

## @param dataset_dir Path to one dataset folder (must have output/stageB_*/snp_results.csv).
## @param engines Character vector of engine names to plot, one panel each, in
##   order. Default: every output/stageB_<engine> present.
## @param colour_by "region" (every locus its own persistent colour, shared
##   across panels -- see .master_loci() above), "significant" (two-colour
##   group: significant vs not), or "none".
## @param value "q" or "p" -- plots -log10(q_display)/-log10(p_display).
## @param alpha Draws the dashed significance reference at -log10(alpha)
##   (ignored if `hline` is set explicitly).
## @param hline Override the reference line's y value, or NULL for none.
## @param highlight Character vector of markers to mark with crosses (passed
##   to every panel's `qtn =`).
## @param ncol Panels per row when stacking (patchwork::wrap_plots(ncol=)).
## @param point_size Passed through to every panel.
## @return A single ggplot (one engine) or a patchwork object (multiple engines).
ldm_manhattan <- function(dataset_dir, engines = NULL,
                          colour_by = c("region", "significant", "none"),
                          value = c("q", "p"), alpha = 0.05, hline = -log10(alpha),
                          highlight = NULL, ncol = 1, point_size = 1.2) {
  colour_by <- match.arg(colour_by)
  value <- match.arg(value)
  value_col <- paste0(value, "_display")

  if (is.null(engines)) {
    out_dir <- file.path(dataset_dir, "output")
    dd <- list.dirs(out_dir, recursive = FALSE, full.names = FALSE)
    engines <- sub("^stageB_", "", grep("^stageB_", dd, value = TRUE))
  }
  if (!length(engines)) stop("ldm_manhattan(): no stageB_<engine> output under ", dataset_dir, "/output.")

  loci <- pal <- NULL
  if (colour_by == "region") {
    loci <- .master_loci(dataset_dir, engines)
    if (nrow(loci)) {
      pal <- stats::setNames(rep(LDscnR::default_cluster_colours(), length.out = nrow(loci)), loci$locus_id)
      data.table::setkey(loci, Chr, from, to)
    }
  }

  panels <- lapply(engines, function(eng) {
    f <- file.path(dataset_dir, "output", paste0("stageB_", eng), "snp_results.csv")
    if (!file.exists(f)) stop("ldm_manhattan(): missing ", f,
                              " -- run_stage_B(dataset_dir, \"", eng, "\") first.")
    d <- data.table::fread(f)
    map_e <- d[, .(marker, Chr, Pos)]
    yv <- -log10(pmax(d[[value_col]], .Machine$double.xmin))
    names(yv) <- d$marker

    args <- list(map = map_e, value = yv, value_label = sprintf("-log10(%s)", value),
                title = eng, hline = hline, qtn = highlight, point_size = point_size)
    if (colour_by == "region") {
      if (!is.null(loci) && nrow(loci)) {
        qx <- d[, .(marker, Chr, from = Pos, to = Pos)]
        lj <- data.table::foverlaps(qx, loci, by.x = c("Chr", "from", "to"),
                                    type = "within", mult = "first", nomatch = NA)
        grp <- stats::setNames(lj$locus_id, d$marker)
        grp <- grp[!is.na(grp)]
        if (length(grp)) { args$group <- grp; args$group_colours <- pal[unique(grp)] }
      }
    } else if (colour_by == "significant") {
      tested <- d[tested == TRUE]
      if (nrow(tested)) args$group <- stats::setNames(
        ifelse(tested$significant, "significant", "not significant"), tested$marker)
    }
    do.call(LDscnR::ld_manhattan, args)
  })
  names(panels) <- engines

  if (length(panels) == 1L) return(panels[[1L]])
  ## guides = "collect": each panel's own ld_manhattan() call builds its own
  ## legend; since colour_by = "region" now hands every panel the SAME
  ## group_colours palette, those legends share entries rather than merely
  ## resembling each other -- collect them into one so the same locus_id
  ## doesn't print once per panel that happens to contain it.
  patchwork::wrap_plots(panels, ncol = ncol, guides = "collect")
}
