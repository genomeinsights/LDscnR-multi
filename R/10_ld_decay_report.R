## =============================================================================
## LDscnR-multi/R/10_ld_decay_report.R
##
## REPORTING LAYER over the LD_decay object build_stage1() already computes
## (compute_LD_decay(), R/02_stage1_cluster.R) and caches in cache/stage1.rds
## -- no recomputation happens here, this just writes out what's already
## there in a form a user can inspect without loading the .rds. As with
## Stage B's report.txt (R/04_stage_B_outlier.R), this reuses LDscnR's own
## reporting methods rather than reimplementing them:
##   - decay_summary.csv : LD_decay$decay_sum, one row per chromosome (fitted
##     decay parameters a/c/b, size-predicted a_pred/c_pred, recommended
##     slide per rho target -- see print.ld_decay()/compute_LD_decay() in
##     LDscnR for what each column means).
##   - decay_print.txt   : print(LD_decay)'s own console summary, captured.
##   - decay_curves.pdf  : LDscnR's own plot.ld_decay() -- one genome-wide
##     "summary" page (decay rate vs chromosome size, rho covered by the
##     current slide, informative windows per chromosome, recommended slide
##     sizes), then one "chr" page per chromosome (window-wise a/c, LD
##     contrast across windows, and the fitted decay curves themselves).
## =============================================================================

## @param dataset_dir Path to one dataset folder.
## @param cfg Resolved config (build_stage1()'s own default).
## @param force Passed through to build_stage1() -- rebuild stage1 first if
##   its receipt is stale, then report from the (possibly rebuilt) result.
## @return LD_decay, invisibly.
report_ld_decay <- function(dataset_dir, cfg = resolve_config(dataset_dir), force = FALSE) {
  s1 <- build_stage1(dataset_dir, cfg, force = force)
  LD_decay <- s1$LD_decay

  out_dir <- file.path(dataset_dir, "output", "ld_decay")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  data.table::fwrite(LD_decay$decay_sum, file.path(out_dir, "decay_summary.csv"))
  writeLines(utils::capture.output(print(LD_decay)), file.path(out_dir, "decay_print.txt"))

  grDevices::pdf(file.path(out_dir, "decay_curves.pdf"), width = 9, height = 7)
  ## on.exit, not a bare dev.off() below: a plot.ld_decay() error partway
  ## through the per-chromosome loop would otherwise leave this pdf() device
  ## open, silently catching whatever base-graphics call runs next (in this
  ## call or a later one in the same session) instead of its intended target.
  on.exit(grDevices::dev.off(), add = TRUE)
  plot(LD_decay, type = "summary")
  for (ch in names(LD_decay$by_chr)) plot(LD_decay, type = "chr", chr = ch)

  say("    [ld_decay] %s\n", out_dir)
  invisible(LD_decay)
}
