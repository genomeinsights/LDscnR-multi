## =============================================================================
## LDscnR-multi/R/07_run_all.R
##
## THE BATCH ENTRY POINT: loop run_dataset() over as many folders as exist,
## isolate one dataset's failure from the rest (a bad panel should not abort a
## run of hundreds), and roll up TWO summary tables -- the actual deliverable
## for screening many datasets:
##   summary.csv             one row per (dataset, engine): what was tested,
##                           how much was significant, how long it took, plus
##                           the two headline structure-alignment numbers
##                           (repeated per engine row, so a screen of many
##                           datasets can sort/filter on either concern at
##                           once).
##   alignment_summary.csv   one row per dataset (only those with
##                           input/population.rds): every observed alignment
##                           measure plus each null scheme's comparison stats
##                           -- the detailed companion to summary.csv's two
##                           headline columns. See R/08_structure_alignment.R.
## =============================================================================

## Wide, one-row-per-dataset detail: every observed measure from
## check_structure_alignment(), plus each (scheme, metric) pair's null
## comparison, flattened out of null_summary's long format. NULL (zero-row)
## when the dataset has no alignment result.
.alignment_detail_row <- function(dataset_dir, alignment) {
  if (is.null(alignment)) return(NULL)
  obs <- alignment$observed
  row <- data.table::data.table(dataset = basename(dataset_dir))
  row <- cbind(row, obs[, .(n_individuals, n_populations, n_structure_groups,
                            grm_axes_r2_population, grm_axes_adjusted_r2_population,
                            structure_r2_population, structure_r2_individual_weighted)])
  ns <- alignment$null_summary
  for (i in seq_len(nrow(ns))) {
    prefix <- paste0(gsub("[^A-Za-z0-9]+", "_", ns$scheme[i]), "__", ns$metric[i])
    row[[paste0(prefix, "__null_mean")]] <- ns$null_mean[i]
    row[[paste0(prefix, "__observed_within_null95")]] <- ns$observed_within_null95[i]
    row[[paste0(prefix, "__proportion_null_ge_observed")]] <- ns$proportion_null_ge_observed[i]
  }
  row
}

## The two headline alignment numbers (observed value + one-sided null
## comparison, "within-group permutation" scheme -- the manuscript's primary
## null), repeated onto every summary.csv row for this dataset. All-NA when
## there is no alignment result, so summary.csv's schema doesn't depend on
## which datasets happened to carry input/population.rds.
.alignment_headline <- function(alignment) {
  if (is.null(alignment)) return(data.table::data.table(
    align_grm_axes_r2 = NA_real_, align_grm_axes_p = NA_real_,
    align_structure_r2 = NA_real_, align_structure_p = NA_real_))
  ns <- alignment$null_summary[scheme == "within-group permutation"]
  g <- ns[metric == "grm_axes_r2_population"]
  s <- ns[metric == "structure_r2_population"]
  data.table::data.table(
    align_grm_axes_r2 = g$observed[1], align_grm_axes_p = g$proportion_null_ge_observed[1],
    align_structure_r2 = s$observed[1], align_structure_p = s$proportion_null_ge_observed[1])
}

.summarise_one <- function(dataset_dir, res, runtime_s, status, err_msg) {
  align_hl <- .alignment_headline(if (identical(status, "ok")) res$alignment else NULL)
  if (!identical(status, "ok")) {
    return(cbind(data.table::data.table(
      dataset = basename(dataset_dir), engine = NA_character_, n_markers = NA_integer_,
      n_units_tested = NA_integer_, n_significant = NA_integer_, n_regions = NA_integer_,
      runtime_s = runtime_s, status = status, error = err_msg), align_hl))
  }
  n_markers <- nrow(res$stage1$map)
  if (!length(res$results)) {
    return(cbind(data.table::data.table(
      dataset = basename(dataset_dir), engine = NA_character_, n_markers = n_markers,
      n_units_tested = NA_integer_, n_significant = NA_integer_, n_regions = NA_integer_,
      runtime_s = runtime_s, status = "no_engine", error = NA_character_), align_hl))
  }
  data.table::rbindlist(lapply(names(res$results), function(eng) {
    r <- res$results[[eng]]
    cbind(data.table::data.table(
      dataset = basename(dataset_dir), engine = eng, n_markers = n_markers,
      n_units_tested = nrow(r$test$units), n_significant = sum(r$test$units$significant),
      n_regions = nrow(r$test$regions), runtime_s = runtime_s, status = "ok", error = NA_character_),
      align_hl)
  }))
}

## @param dataset_dirs Character vector of dataset folder paths.
## @param summary_path Where to write the per-(dataset,engine) roll-up CSV
##   (NULL to skip writing; always returned as a data.table either way).
## @param alignment_summary_path Where to write the per-dataset alignment
##   detail CSV (NULL to skip writing). Default derives from
##   `summary_path`'s directory when given, else "alignment_summary.csv".
## @param cores Datasets to process in parallel (parallel::mclapply, Unix only).
## @param stop_on_error If TRUE, a dataset's error propagates and stops the
##   batch; default FALSE isolates it as one "error" row instead.
## @param ... Passed through to run_dataset() for every dataset (engines, cfg,
##   annotation, chrom_lengths, plot, force).
## @return A list(summary, alignment_summary), each a data.table (also
##   written to CSV unless the corresponding path is NULL).
run_all <- function(dataset_dirs, summary_path = "summary.csv",
                    alignment_summary_path = if (!is.null(summary_path))
                      file.path(dirname(summary_path), "alignment_summary.csv") else NULL,
                    cores = 1L, stop_on_error = FALSE, ...) {
  results <- vector("list", length(dataset_dirs))
  run_one <- function(dd) {
    say("\n>>> %s\n", dd)
    t0 <- Sys.time()
    res <- tryCatch(run_dataset(dd, ...), error = function(e) e)
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (inherits(res, "error")) {
      if (stop_on_error) stop(res)
      say("    [ERROR] %s\n", conditionMessage(res))
      return(list(summary = .summarise_one(dd, NULL, dt, "error", conditionMessage(res)),
                 alignment = NULL))
    }
    list(summary = .summarise_one(dd, res, dt, "ok", NA_character_),
        alignment = .alignment_detail_row(dd, res$alignment))
  }
  rows <- if (cores > 1L) parallel::mclapply(dataset_dirs, run_one, mc.cores = cores)
          else lapply(dataset_dirs, run_one)

  summary <- data.table::rbindlist(lapply(rows, `[[`, "summary"), fill = TRUE)
  alignment_summary <- data.table::rbindlist(lapply(rows, `[[`, "alignment"), fill = TRUE)

  if (!is.null(summary_path)) {
    data.table::fwrite(summary, summary_path)
    say("\n=== wrote %s (%d row(s)) ===\n", summary_path, nrow(summary))
  }
  if (!is.null(alignment_summary_path) && nrow(alignment_summary)) {
    data.table::fwrite(alignment_summary, alignment_summary_path)
    say("=== wrote %s (%d row(s)) ===\n", alignment_summary_path, nrow(alignment_summary))
  }
  list(summary = summary[], alignment_summary = alignment_summary[])
}
