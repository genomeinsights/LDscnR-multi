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
## =============================================================================

## @param dataset_dir Path to one dataset folder (must have output/stageB_*/snp_results.csv).
## @param engines Character vector of engine names to plot, one panel each, in
##   order. Default: every output/stageB_<engine> present.
## @param colour_by "region" (each significant region its own colour, via
##   ld_manhattan(regions=)), "significant" (two-colour group: significant vs
##   not, via ld_manhattan(group=)), or "none".
## @param value "q" or "p" -- which column of snp_results.csv to plot as
##   -log10(value) on the y-axis.
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

  if (is.null(engines)) {
    out_dir <- file.path(dataset_dir, "output")
    dd <- list.dirs(out_dir, recursive = FALSE, full.names = FALSE)
    engines <- sub("^stageB_", "", grep("^stageB_", dd, value = TRUE))
  }
  if (!length(engines)) stop("ldm_manhattan(): no stageB_<engine> output under ", dataset_dir, "/output.")

  panels <- lapply(engines, function(eng) {
    f <- file.path(dataset_dir, "output", paste0("stageB_", eng), "snp_results.csv")
    if (!file.exists(f)) stop("ldm_manhattan(): missing ", f,
                              " -- run_stage_B(dataset_dir, \"", eng, "\") first.")
    d <- data.table::fread(f)
    map_e <- d[, .(marker, Chr, Pos)]
    yv <- -log10(pmax(d[[value]], .Machine$double.xmin))
    names(yv) <- d$marker

    args <- list(map = map_e, value = yv, value_label = sprintf("-log10(%s)", value),
                title = eng, hline = hline, qtn = highlight, point_size = point_size)
    if (colour_by == "region") {
      sig <- d[!is.na(region_id)]
      if (nrow(sig)) args$regions <- unname(split(sig$marker, sig$region_id))
    } else if (colour_by == "significant") {
      tested <- d[tested == TRUE]
      if (nrow(tested)) args$group <- stats::setNames(
        ifelse(tested$significant, "significant", "not significant"), tested$marker)
    }
    do.call(LDscnR::ld_manhattan, args)
  })
  names(panels) <- engines

  if (length(panels) == 1L) return(panels[[1L]])
  patchwork::wrap_plots(panels, ncol = ncol)
}
