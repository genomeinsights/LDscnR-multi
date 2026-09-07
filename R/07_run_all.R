## =============================================================================
## LDscnR-multi/R/07_run_all.R
##
## THE BATCH ENTRY POINT: loop run_dataset() over as many folders as exist,
## isolate one dataset's failure from the rest (a bad panel should not abort a
## run of hundreds), and roll up a summary table -- the actual deliverable for
## screening many datasets: one row per (dataset, engine) with what was tested,
## how much was significant, and how long it took.
## =============================================================================

.summarise_one <- function(dataset_dir, res, runtime_s, status, err_msg) {
  if (!identical(status, "ok")) {
    return(data.table::data.table(
      dataset = basename(dataset_dir), engine = NA_character_, n_markers = NA_integer_,
      n_units_tested = NA_integer_, n_significant = NA_integer_, n_regions = NA_integer_,
      runtime_s = runtime_s, status = status, error = err_msg))
  }
  n_markers <- nrow(res$stage1$map)
  if (!length(res$results)) {
    return(data.table::data.table(
      dataset = basename(dataset_dir), engine = NA_character_, n_markers = n_markers,
      n_units_tested = NA_integer_, n_significant = NA_integer_, n_regions = NA_integer_,
      runtime_s = runtime_s, status = "no_engine", error = NA_character_))
  }
  data.table::rbindlist(lapply(names(res$results), function(eng) {
    r <- res$results[[eng]]
    data.table::data.table(
      dataset = basename(dataset_dir), engine = eng, n_markers = n_markers,
      n_units_tested = nrow(r$test$units), n_significant = sum(r$test$units$significant),
      n_regions = nrow(r$test$regions), runtime_s = runtime_s, status = "ok", error = NA_character_)
  }))
}

## @param dataset_dirs Character vector of dataset folder paths.
## @param summary_path Where to write the roll-up CSV (NULL to skip writing;
##   the summary is always returned as a data.table either way).
## @param cores Datasets to process in parallel (parallel::mclapply, Unix only).
## @param stop_on_error If TRUE, a dataset's error propagates and stops the
##   batch; default FALSE isolates it as one "error" row instead.
## @param ... Passed through to run_dataset() for every dataset (engines, cfg,
##   annotation, chrom_lengths, plot, force).
## @return A data.table, one row per (dataset, engine): n_markers,
##   n_units_tested, n_significant, n_regions, runtime_s, status, error.
run_all <- function(dataset_dirs, summary_path = "summary.csv", cores = 1L,
                    stop_on_error = FALSE, ...) {
  run_one <- function(dd) {
    say("\n>>> %s\n", dd)
    t0 <- Sys.time()
    res <- tryCatch(run_dataset(dd, ...), error = function(e) e)
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (inherits(res, "error")) {
      if (stop_on_error) stop(res)
      say("    [ERROR] %s\n", conditionMessage(res))
      return(.summarise_one(dd, NULL, dt, "error", conditionMessage(res)))
    }
    .summarise_one(dd, res, dt, "ok", NA_character_)
  }
  rows <- if (cores > 1L) parallel::mclapply(dataset_dirs, run_one, mc.cores = cores)
          else lapply(dataset_dirs, run_one)
  summary <- data.table::rbindlist(rows, fill = TRUE)
  if (!is.null(summary_path)) {
    data.table::fwrite(summary, summary_path)
    say("\n=== wrote %s (%d row(s)) ===\n", summary_path, nrow(summary))
  }
  summary[]
}
