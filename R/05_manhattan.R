## =============================================================================
## LDscnR-multi/R/05_manhattan.R
##
## GENOME-WIDE MANHATTAN PLOTTING, styled directly after the ACTUAL plotting
## code (not just the rendered image) of
## ~/gitlab/formica_hybrid/module_manuscript_rho05/module_BayPass/R/
## formica_region_concordance.R's make_panel()/build_block(), which produces
## figures/formica_local_score_vs_ld_regions.png: a single continuous genome
## axis (chromosomes concatenated, not faceted into separate per-chromosome
## panels), grey background points (colour = "grey73", size = 0.45, alpha =
## 0.45), each significant locus as its own coloured vertical streak (size =
## point_size [default 1.6], alpha = 0.90), and the most significant loci
## labelled directly on the plot -- a WHITE label box with bold, per-locus
## COLOURED text and border (fill = "white", colour = locus colour,
## label.size = 0.2), not a solid-colour box, with a leader line to its
## peak via ggrepel::geom_label_repel() (that package IS the
## collision-avoiding label placement the reference figure uses; reused
## here rather than reimplemented, with the reference's own nudge_y/
## box.padding/point.padding/max.iter/seed tuning). Colours themselves come
## from `.assign_locus_colours()`, porting the reference's own algorithm:
## LDscnR::default_cluster_colours(), luminance-filtered to drop near-white
## entries, walked in genomic order picking the available colour with
## greatest minimum CIE Lab distance from the last few already assigned, so
## physically adjacent loci never get confusably similar hues. Stacked
## panels (one per engine) share one x-axis, shown only on the bottom
## panel, and one persistent colour+label per physical locus across every
## panel (`.master_loci()`), so the same peak reads as the same colour and
## the same short ID ("R1", "R2", ...) in every engine's track. This is a
## from-scratch ggplot2 build, NOT a wrapper over LDscnR::ld_manhattan()
## any more -- that function's own per-chromosome facet_wrap() is a
## genuinely different layout from the reference's continuous axis, not a
## parameterisation of it.
##
## Labelling is deliberately selective, matching the reference: every
## significant marker is coloured, but only the `max_labels` most significant
## LOCI per panel get a text label (chosen by that panel's own best q/p
## value among the locus's members) -- labelling every one of dozens of
## significant regions would be unreadable, and the reference itself labels
## only a curated subset of its own coloured streaks.
## =============================================================================

## Build a genome-wide, cross-engine set of "loci": the union of every
## requested engine's significant regions, merged where they overlap
## physically, each assigned a persistent colour AND a short label ("R1",
## "R2", ... in genomic order). Different engines' own stage-2 assembly
## rarely produces identical bounds for "the same" region, but the same
## underlying signal typically does overlap -- merging on overlap, not on an
## exact match, is what lets a shared colour/label survive that.
.master_loci <- function(dataset_dir, engines) {
  regs <- data.table::rbindlist(lapply(engines, function(eng) {
    f <- file.path(dataset_dir, "output", paste0("stageB_", eng), "region_table.csv")
    if (!file.exists(f)) return(NULL)
    r <- data.table::fread(f)
    if (!nrow(r)) return(NULL)
    r[, .(Chr, from, to)]
  }), use.names = TRUE, fill = TRUE)
  if (is.null(regs) || !nrow(regs))
    return(data.table::data.table(Chr = character(), from = numeric(), to = numeric(),
                                  locus_id = character(), label = character()))
  data.table::setorder(regs, Chr, from)
  ## Chromosome-safe interval merge (by = Chr): a global cummax(to) that never
  ## resets at a chromosome boundary silently over-merges every chromosome
  ## after the first sizeable one -- see LDscnR's own R/ld_outlier_internal.R
  ## header for the exact failure this `by = Chr` avoids.
  regs[, grp := cumsum(c(TRUE, from[-1] - cummax(to)[-.N] > 0)), by = Chr]
  loci <- regs[, .(from = min(from), to = max(to)), by = .(Chr, grp)][, grp := NULL]
  loci[, locus_id := sprintf("%s:%.0f-%.0f", Chr, from, to)]
  ## Label order follows genomic order (natural chromosome sort, then
  ## position), so "R1" is always the first physical locus, "R2" the next,
  ## regardless of which engine happened to discover it first. A temporary
  ## column, not setorder(loci, ..ord, from) -- setorder() takes real column
  ## names, and `..ord` (data.table's "from the calling environment" prefix)
  ## is not valid there.
  loci[, .sort_key := .chr_sort_order(Chr)]
  data.table::setorder(loci, .sort_key, from)
  loci[, .sort_key := NULL]
  loci[, label := paste0("R", seq_len(.N))]
  loci[]
}

## Locus colour assignment, ported verbatim (algorithm, not code) from the
## reference figure's own R/formica_region_concordance.R
## (~/gitlab/formica_hybrid/module_manuscript_rho05/module_BayPass):
## LDscnR::default_cluster_colours(), filtered to drop near-white entries
## (invisible against a white panel background) by standard luminance, then
## walked in genomic order (the order `loci` is already sorted in) picking
## at each step the available colour with the greatest minimum CIE Lab
## distance from the last `lookback` colours already assigned -- so two
## PHYSICALLY NEARBY loci (the ones most likely to sit close enough on the
## page to be confused) are never handed visually similar hues by chance,
## even though the underlying palette repeats once exhausted.
.assign_locus_colours <- function(n, lookback = 5L) {
  pal <- unique(LDscnR::default_cluster_colours())
  rgb <- t(grDevices::col2rgb(pal)) / 255
  luminance <- 0.299 * rgb[, 1] + 0.587 * rgb[, 2] + 0.114 * rgb[, 3]
  pal <- pal[luminance <= 0.75]
  lab <- grDevices::convertColor(t(grDevices::col2rgb(pal)) / 255, from = "sRGB", to = "Lab")
  dist_mat <- as.matrix(stats::dist(lab))

  chosen <- integer(n)
  available <- seq_along(pal)
  for (i in seq_len(n)) {
    if (!length(available)) available <- seq_along(pal)
    if (i == 1L) {
      pick <- available[1L]
    } else {
      previous <- utils::tail(chosen[seq_len(i - 1L)], lookback)
      separation <- vapply(available, function(j) min(dist_mat[j, previous]), numeric(1))
      pick <- available[which.max(separation)]
    }
    chosen[i] <- pick
    available <- setdiff(available, pick)
  }
  pal[chosen]
}

## Natural chromosome sort key ("Chr2" before "Chr10"): the numeric part of
## the chromosome name if every chromosome has one and they're distinct,
## else a plain alphabetical fallback.
.chr_sort_order <- function(chr) {
  n <- suppressWarnings(as.numeric(gsub("[^0-9]", "", chr)))
  if (anyNA(n) || anyDuplicated(n)) match(chr, sort(unique(chr))) else n
}

## Cumulative genome coordinate: chromosomes concatenated in natural order,
## separated by `chr_gap` (as a fraction of the largest chromosome's own
## length), so a genome-wide x-axis can be built without per-chromosome
## facets. Returns a data.table (Chr, len, offset, mid) keyed on Chr; `mid`
## is where that chromosome's axis tick label goes. NOTE: setkey() below
## re-sorts the table by Chr's own (lexicographic) order for O(1) keyed
## lookup elsewhere (`gc[Chr, offset]`) -- callers that need NATURAL
## genomic order back (e.g. the alternating chromosome-shading bands in
## ldm_manhattan()) must re-derive it via .chr_sort_order(), not assume
## this return value is still in the order it was built in.
.genome_coords <- function(map, chr_gap = 0.03) {
  chr_len <- data.table::as.data.table(map)[, .(len = max(Pos)), by = Chr]
  chr_len <- chr_len[order(.chr_sort_order(Chr))]
  gap <- max(chr_len$len) * chr_gap
  chr_len[, offset := cumsum(data.table::shift(len, fill = 0)) + (seq_len(.N) - 1) * gap]
  chr_len[, mid := offset + len / 2]
  data.table::setkey(chr_len, Chr)
  chr_len[]
}

## @param dataset_dir Path to one dataset folder (must have output/stageB_*/snp_results.csv).
## @param engines Character vector of engine names to plot, one panel each, in
##   order. Default: every output/stageB_<engine> present.
## @param colour_by "region" (every locus its own persistent colour AND label,
##   shared across panels -- see .master_loci() above), "significant"
##   (two-colour group: significant vs not, no labels), or "none".
## @param value "q" or "p" (plots -log10(q_display)/-log10(p_display)), or
##   "ld_w_095" -- local-LD support at rho = 0.95 (R/02_stage1_cluster.R),
##   plotted RAW (no -log10 transform), with no significance line by default.
##   Otherwise identical: same panels, same genome axis, same colouring.
## @param alpha Draws the dashed significance reference at -log10(alpha)
##   (ignored if `hline` is set explicitly).
## @param hline Override the reference line's y value, or NULL for none.
## @param highlight Character vector of markers to mark with black crosses.
## @param max_labels Label at most this many of the most significant loci PER
##   PANEL (by that panel's own best q/p among each locus's members) -- every
##   significant marker is still coloured regardless; only text labels are
##   capped, matching the reference figure's own curated (not exhaustive)
##   labelling. Ignored when colour_by != "region".
## @param ncol Panels per column when stacking.
## @param point_size Passed through to every panel; default 1.6 matches the
##   reference figure's own coloured-point size.
## @param chr_gap Space between adjacent chromosomes on the genome axis, as
##   a fraction of the largest chromosome's own length.
## @param chr_shade Alternating light-grey background bands, one per
##   chromosome, so a chromosome's start/end is visible even where it has
##   no significant (coloured) markers of its own to break up the grey
##   point cloud. Set FALSE to turn off.
## @return A single ggplot (one engine) or a patchwork object (multiple engines).
ldm_manhattan <- function(dataset_dir, engines = NULL,
                          colour_by = c("region", "significant", "none"),
                          value = c("q", "p", "ld_w_095"), alpha = 0.05,
                          hline = if (value == "ld_w_095") NULL else -log10(alpha),
                          highlight = NULL, max_labels = 20L, ncol = 1, point_size = 1.6,
                          chr_gap = 0.03, chr_shade = TRUE) {
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
      ## loci is already sorted in genomic order (.master_loci()), which is
      ## exactly the order .assign_locus_colours()'s local-contrast walk
      ## needs -- physically adjacent loci get maximally separated hues.
      pal <- stats::setNames(.assign_locus_colours(nrow(loci)), loci$locus_id)
      data.table::setkey(loci, Chr, from, to)
    }
  }

  ## One genome coordinate system for every panel (same underlying map/
  ## markers across all engines in this pipeline), so peaks line up
  ## vertically when panels are stacked.
  gc <- .genome_coords(data.table::fread(file.path(dataset_dir, "output", paste0("stageB_", engines[1]),
                                                    "snp_results.csv"), select = c("Chr", "Pos")),
                       chr_gap = chr_gap)

  ## Alternating chromosome-shading bands, EVERY OTHER chromosome in NATURAL
  ## genomic order -- re-derived here, not read off `gc`'s own row order,
  ## because .genome_coords()'s setkey() leaves `gc` sorted lexicographically
  ## by Chr (needed for its OWN O(1) keyed lookups elsewhere), which is not
  ## genomic order (".chr_sort_order()"'s own header explains why "Chr2"
  ## would otherwise sort after "Chr10").
  chr_bands <- NULL
  if (isTRUE(chr_shade)) {
    chr_bands <- data.table::copy(gc)
    chr_bands[, .sort_key := .chr_sort_order(Chr)]
    data.table::setorder(chr_bands, .sort_key)
    chr_bands <- chr_bands[seq(2L, .N, by = 2L)]   ## every other chromosome, starting with the 2nd
    chr_bands[, `:=`(xmin = offset, xmax = offset + len)]
  }

  panels <- lapply(engines, function(eng) {
    f <- file.path(dataset_dir, "output", paste0("stageB_", eng), "snp_results.csv")
    if (!file.exists(f)) stop("ldm_manhattan(): missing ", f,
                              " -- run_stage_B(dataset_dir, \"", eng, "\") first.")
    ## na.strings = c("NA", ""): fwrite() writes NA_character_ as an empty
    ## CSV field; fread()'s own default na.strings = "NA" doesn't read that
    ## back as NA -- region_id/marker_p etc. would silently become "" instead
    ## of NA for every untested marker. See R/04_stage_B_outlier.R's own
    ## matching comment (its cache-hit path has the identical read).
    d <- data.table::fread(f, na.strings = c("NA", ""))
    d[, gpos := Pos + gc[Chr, offset]]
    d[, yv := if (value == "ld_w_095") as.numeric(get(value_col))
             else -log10(pmax(get(value_col), .Machine$double.xmin))]
    value_label <- if (value == "ld_w_095") "ld_w (rho = 0.95)" else sprintf("-log10(%s)", value)

    d[, locus_id := NA_character_]
    if (colour_by == "region" && !is.null(loci) && nrow(loci)) {
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
        d[!is.na(region_id) & region_id %chin% names(region_to_locus),
         locus_id := region_to_locus[region_id]]
      }
    } else if (colour_by == "significant") {
      d[tested == TRUE, locus_id := ifelse(significant, "significant", "not significant")]
    }

    grey <- d[is.na(locus_id) & is.finite(yv)]
    coloured <- d[!is.na(locus_id) & is.finite(yv)]

    ## Headroom for labels, and the axis range itself: from EVERY finite
    ## point in the panel (grey included), matching the reference figure's
    ## own `hi <- max(pdat$stat); lo <- min(pdat$stat)` -- not just the
    ## coloured subset, so a label nudged above the highest coloured peak
    ## still has defined headroom even when the true max is a grey point.
    yrange <- d$yv[is.finite(d$yv)]
    hi <- if (length(yrange)) max(yrange) else 1
    lo <- if (length(yrange)) min(yrange) else 0
    span <- hi - lo
    if (!is.finite(span) || span == 0) span <- 1

    p <- ggplot2::ggplot()
    if (!is.null(chr_bands) && nrow(chr_bands))
      ## ymin/ymax = -Inf/Inf, not `lo`/`hi`: the panel's own visible range
      ## is set below by coord_cartesian(), and an infinite rect is cropped
      ## to whatever that ends up being -- so the shading always reaches
      ## the true top/bottom of the panel regardless of this panel's own
      ## data range, without needing to duplicate that ylim math here.
      p <- p + ggplot2::geom_rect(data = chr_bands,
                                  ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                                  inherit.aes = FALSE, fill = "grey95")
    p <- p +
      ggplot2::geom_point(data = grey, ggplot2::aes(x = gpos, y = yv),
                         colour = "grey73", size = 0.45, alpha = 0.45) +
      ggplot2::scale_x_continuous(breaks = gc$mid, labels = sub("^Chr", "", gc$Chr),
                                  expand = ggplot2::expansion(mult = c(0.005, 0.005))) +
      ggplot2::coord_cartesian(ylim = c(lo - 0.04 * span, hi + 0.42 * span), clip = "off") +
      ggplot2::labs(x = "Chromosome", y = value_label, title = eng) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())

    if (nrow(coloured)) {
      if (colour_by == "region") {
        p <- p + ggplot2::geom_point(data = coloured, ggplot2::aes(x = gpos, y = yv, colour = locus_id),
                                     size = point_size, alpha = 0.90, show.legend = FALSE) +
          ggplot2::scale_colour_manual(values = pal)

        top_loci <- coloured[, .(peak_y = max(yv), gpos = gpos[which.max(yv)],
                                 region_id = region_id[which.max(yv)]), by = locus_id]
        data.table::setorder(top_loci, -peak_y)
        top_loci <- utils::head(top_loci, max_labels)
        ## match(), not loci[locus_id, label]: `loci` is keyed on (Chr, from,
        ## to) for the foverlaps() join above, not on locus_id, and indexing
        ## a keyed data.table by a non-key vector silently returns NA rather
        ## than erroring -- caught directly, not assumed.
        top_loci[, label := loci$label[match(locus_id, loci$locus_id)]]
        ## Floor-stability tier, appended to the label text ("R5 (3/3)") when
        ## run_floor_profile() (R/11_floor_profile.R) has been run for this
        ## engine's own statistic -- keyed on THIS engine's own region_id
        ## (the peak marker's, picked above), not the cross-engine `locus_id`
        ## .master_loci() may have widened by merging several engines'
        ## regions together. Silently absent (no label suffix) when the
        ## floor profile hasn't been run, or this engine isn't "emmax_unit"/
        ## "emmax_simes" (external engines have no floor-profile concept).
        stat_for_eng <- if (eng == "emmax_unit") "unit" else if (eng == "emmax_simes") "simes" else NA
        stab_f <- if (!is.na(stat_for_eng))
          file.path(dataset_dir, "output", "floor_profile", sprintf("canonical_regions_%s.csv", stat_for_eng))
        if (!is.na(stat_for_eng) && !is.null(stab_f) && file.exists(stab_f)) {
          stab <- data.table::fread(stab_f, select = c("region_id", "stability_tier"))
          top_loci[stab, on = "region_id", tier := i.stability_tier]
          top_loci[!is.na(tier), label := paste0(label, " (", tier, ")")]
        }
        if (nrow(top_loci))
          ## White label fill + coloured (per-locus) bold text/border, not a
          ## solid-colour box -- matches the reference figure's own
          ## geom_label_repel() call exactly (fill = "white", colour =
          ## region_id), including its nudge/collision-avoidance tuning.
          p <- p + ggrepel::geom_label_repel(
            data = top_loci, ggplot2::aes(x = gpos, y = peak_y, label = label, colour = locus_id),
            fill = "white", alpha = 0.96, fontface = "bold", size = 3, label.size = 0.2,
            box.padding = 0.25, point.padding = 0.1, nudge_y = 0.19 * span,
            min.segment.length = 0, max.overlaps = Inf, max.time = 3, max.iter = 20000,
            seed = 100, show.legend = FALSE)
      } else {
        pal2 <- c(significant = "firebrick", `not significant` = "steelblue3")
        p <- p + ggplot2::geom_point(data = coloured, ggplot2::aes(x = gpos, y = yv, colour = locus_id),
                                     size = point_size, alpha = 0.90) +
          ggplot2::scale_colour_manual(values = pal2, name = NULL)
      }
    }

    if (!is.null(hline)) p <- p + ggplot2::geom_hline(yintercept = hline, linetype = "dashed", linewidth = 0.3)
    if (!is.null(highlight)) {
      hi <- d[marker %chin% highlight & is.finite(yv)]
      if (nrow(hi)) p <- p + ggplot2::geom_point(data = hi, ggplot2::aes(x = gpos, y = yv),
                                                 shape = 3, size = 3, stroke = 1, colour = "black")
    }
    p
  })
  names(panels) <- engines

  ## Chromosome axis text/title only on the bottom panel when stacking --
  ## every panel shares the same genome coordinate system (built once,
  ## above), so repeating it on each panel would be redundant.
  if (length(panels) > 1L) {
    for (i in seq_len(length(panels) - 1L)) {
      panels[[i]] <- panels[[i]] + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                                   axis.title.x = ggplot2::element_blank(),
                                                   axis.ticks.x = ggplot2::element_blank())
    }
  }

  if (length(panels) == 1L) return(panels[[1L]])
  patchwork::wrap_plots(panels, ncol = ncol)
}
