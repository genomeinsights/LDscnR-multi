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
## - It plots `p_display`/`q_display`, not the tested `unit_p`/`unit_q` (see
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
## Colouring is driven by snp_results.csv's OWN `region_id` (correct,
## membership-based -- see R/04_stage_B_outlier.R's header) mapped to the
## shared locus that contains it, via a REGION-level lookup (this engine's own
## region_table.csv against `loci`, a handful of rows) -- never a fresh
## marker-position join against `loci` here. An earlier version did exactly
## that fresh join, which (an external audit caught) coloured every marker
## physically inside a region's bounding span regardless of whether it was
## ever actually a tested member of that region -- the same coordinate-vs-
## membership mistake `.build_snp_results()`'s header describes, reproduced a
## second time in the plotting code.
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
## @param value "q" or "p" (plots -log10(q_display)/-log10(p_display)), or
##   "ld_w_095" -- local-LD support at rho = 0.95 (R/02_stage1_cluster.R),
##   plotted RAW (no -log10 transform; already a bounded [0,1] LD statistic),
##   with `hline`/hline-based defaults skipped since -log10(alpha) isn't a
##   meaningful reference on this scale. Otherwise identical: same panels,
##   same faceting, same colour_by. Useful to compare a significant region's
##   location against local LD directly, independent of any association test.
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
                          value = c("q", "p", "ld_w_095"), alpha = 0.05,
                          hline = if (value == "ld_w_095") NULL else -log10(alpha),
                          highlight = NULL, ncol = 1, point_size = 1.2) {
  colour_by <- match.arg(colour_by)
  value <- match.arg(value)
  value_col <- if (value == "ld_w_095") "ld_w_095" else paste0(value, "_display")

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
    yv <- if (value == "ld_w_095") as.numeric(d[[value_col]])
          else -log10(pmax(d[[value_col]], .Machine$double.xmin))
    names(yv) <- d$marker
    value_label <- if (value == "ld_w_095") "ld_w (rho = 0.95)" else sprintf("-log10(%s)", value)

    args <- list(map = map_e, value = yv, value_label = value_label,
                title = eng, hline = hline, qtn = highlight, point_size = point_size)
    if (colour_by == "region") {
      if (!is.null(loci) && nrow(loci)) {
        rt_f <- file.path(dataset_dir, "output", paste0("stageB_", eng), "region_table.csv")
        rt <- if (file.exists(rt_f)) data.table::fread(rt_f) else data.table::data.table()
        if (nrow(rt)) {
          ## same region_id string .build_snp_results() used, so this engine's
          ## own d$region_id values match these rows exactly.
          rt[, region_id := sprintf("%s:%.0f-%.0f", Chr, from, to)]
          data.table::setkey(rt, Chr, from, to)
          rl <- data.table::foverlaps(rt[, .(Chr, from, to, region_id)], loci,
                                      by.x = c("Chr", "from", "to"), type = "within",
                                      mult = "first", nomatch = NA)
          region_to_locus <- stats::setNames(rl$locus_id, rl$region_id)
          sig <- d[!is.na(region_id) & region_id %chin% names(region_to_locus)]
          if (nrow(sig)) {
            grp <- stats::setNames(region_to_locus[sig$region_id], sig$marker)
            grp <- grp[!is.na(grp)]
            if (length(grp)) { args$group <- grp; args$group_colours <- pal[unique(grp)] }
          }
        }
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
