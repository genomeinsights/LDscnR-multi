## =============================================================================
## LDscnR-multi/R/00_config.R
##
## DEFAULTS FOR EVERY DATASET, AND THE RECEIPT MACHINERY THAT TRACKS WHETHER A
## STAGE NEEDS RERUNNING. Not a repo-wide constant file the way module_3sp's
## 00_config.R is -- there is no single dataset here -- but the same discipline
## applies one level up: every tunable is named once, in DEFAULTS, and a
## dataset overrides it explicitly (input/config.R) rather than a stage script
## growing its own copy of a threshold.
##
## Modelled directly on LDscnR-paper/module_3sp/R/00_config.R's receipt
## machinery and check_ldscnr() pin, generalised from one hardcoded module path
## to one `stage_dir` per (dataset, stage).
## =============================================================================
suppressMessages({library(data.table); library(digest)})

## ---- 1. DEFAULTS -------------------------------------------------------------
## Every tunable, named once. A dataset's input/config.R may reassign any of
## these names; run_dataset()/run_stage_A()/run_stage_B() argument defaults
## resolve to DEFAULTS, so a call-site argument can override a dataset config,
## which overrides this.
DEFAULTS <- list(
  ## stage1 clustering (genotype-only)
  cr_rho              = 0.5,      # ld_complexity_reduction rho
  decay_args          = list(min_maf_decay = 0.1, q = 0.95, n_sub_bg = 5000,
                              n_win_decay = 20, overlap = 0.5, max_SNPs_decay = Inf,
                              prob_robust = 0.95, max_pairs = 5000, ld_method = "corr",
                              n_strata = 20, slide = 1000, cores = 1),
  ## Stage A: EMMAX
  statistics          = c("unit", "simes"),   # which arm(s) to compute; either or both
  unit_repr           = "consensus_dosage",   # ld_unit_matrix() representation for the "unit" arm
  grm_method           = "GCTA",               # SNPRelate::snpgdsGRM method
  b_unit               = 1000L,                 # permutations for the "unit" arm (cheap)
  b_simes              = 200L,                  # permutations for the "simes" arm (rescans every marker)
  ## Stage B: outlier test / permutation null / region rotation
  size_floor           = 8L,
  alpha                = 0.05,
  assembly              = "stage2_discovered",  # or "physical"
  score_threshold       = 0.80,
  distance_threshold    = 1e5,
  gap                   = 3e5,
  n_rotations           = 10000L,
  rotation_scheme        = "within",
  ## misc
  seed                  = 1L,
  cores                 = 1L
)

## ---- 2. LDscnR VERSION PIN ---------------------------------------------------
## Same reasoning as module_3sp/R/00_config.R: `packageVersion("LDscnR")` never
## changes across commits of a 0.0.0.9000 package, so the only check that
## actually catches a stale install is a content hash over R/*.R at the commit
## this repo was validated against. A batch of hundreds of datasets is exactly
## where a silent stale install would be most expensive to discover late.
LDSCNR_PIN <- list(
  repo    = path.expand("~/gitlab/LDscnR"),
  ## outlier-scan (the branch this was first pinned to) was fast-forwarded
  ## into main and development has continued there since -- outlier-scan
  ## itself is now a stale ancestor, not a separate line.
  branch  = "main",
  sha     = "2993ea80e6c8",
  src_sha = "dc89dbdc504fb660"
)

check_ldscnr <- function(stop_on_fail = !nzchar(Sys.getenv("LDSCNR_LAX"))) {
  v <- as.character(utils::packageVersion("LDscnR"))
  g <- function(...) tryCatch(system2("git", c("-C", LDSCNR_PIN$repo, ...),
                                      stdout = TRUE, stderr = FALSE), error = function(e) character())
  head_sha <- substr(paste(g("rev-parse", "HEAD"), collapse = ""), 1, 12)
  dirty <- length(g("status", "--porcelain", "--untracked-files=no")) > 0
  src <- sort(g("ls-files", "R/"))
  cur <- if (!length(src)) NA_character_ else {
    fp <- file.path(LDSCNR_PIN$repo, src)
    substr(digest::digest(paste(vapply(fp[file.exists(fp)],
             function(f) paste(readLines(f, warn = FALSE), collapse = "\n"), ""), collapse = "\n"),
             algo = "sha256", serialize = FALSE), 1, 16)
  }
  ok <- !dirty && identical(cur, LDSCNR_PIN$src_sha)
  cat(sprintf("  LDscnR %s | commit %s (%s)%s\n", v, head_sha, LDSCNR_PIN$branch,
              if (dirty) " [TRACKED CHANGES]" else ""))
  cat(sprintf("  source hash %s %s pin %s\n", cur, if (identical(cur, LDSCNR_PIN$src_sha)) "==" else "!=",
              LDSCNR_PIN$src_sha))
  if (!ok) {
    m <- paste("LDscnR source does not match the pin recorded in R/00_config.R.",
               "Reinstall from", LDSCNR_PIN$repo, "and update LDSCNR_PIN$src_sha (see check_ldscnr()).",
               "Set LDSCNR_LAX=1 to proceed anyway.")
    if (stop_on_fail) stop(m) else warning(m)
  }
  invisible(list(version = v, sha = head_sha, dirty = dirty, src_sha = cur, ok = ok))
}

## ---- 3. RECEIPT MACHINERY -----------------------------------------------------
## Per (dataset, stage), not per module: `stage_dir` is a path the caller
## supplies (e.g. <dataset>/output/pvalues/emmax_unit), so the exact same
## staleness contract module_3sp uses works unmodified across every dataset a
## batch touches.
receipt_path <- function(stage_dir) file.path(stage_dir, "_receipt.rds")

sha <- function(f) if (file.exists(f)) digest::digest(f, algo = "sha256", file = TRUE) else NA_character_

## Every input path is normalised before it's used as a hash-table key (by
## write_receipt() when storing, by stage_stale() when comparing) -- a raw
## path string is not, e.g. "examples/x/input/map.rds" and
## "/Users/.../examples/x/input/map.rds" name the same file but are different
## strings, so a receipt written from one and checked from the other looked
## "changed" and triggered a needless rebuild (found running the exact same
## dataset via a relative-path caller after an earlier absolute-path run).
## normalizePath(mustWork = FALSE) so a genuinely missing input still shows up
## as "changed" (sha() returns NA for it) rather than erroring here.
.norm_path <- function(p) if (!length(p)) character() else normalizePath(p, mustWork = FALSE)

write_receipt <- function(stage_dir, inputs = character(), params = list(), outputs = character()) {
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  inputs <- .norm_path(inputs)
  saveRDS(list(when = Sys.time(),
               ldscnr = tryCatch(check_ldscnr(stop_on_fail = FALSE), error = function(e) NA),
               inputs = data.table(path = inputs, sha256 = vapply(inputs, sha, "")),
               params = params, outputs = outputs), receipt_path(stage_dir))
  invisible(TRUE)
}

## TRUE when the stage must (re)run: no receipt, a changed parameter, or a
## changed input. Reports why, so a rerun a user didn't expect is explainable.
stage_stale <- function(stage_dir, inputs = character(), params = list(), label = basename(stage_dir)) {
  rp <- receipt_path(stage_dir)
  if (!file.exists(rp)) { message("    [", label, "] no receipt -- will run"); return(TRUE) }
  r <- readRDS(rp)
  if (!identical(params, r$params)) { message("    [", label, "] parameters changed -- will run"); return(TRUE) }
  inputs <- .norm_path(inputs)
  now <- vapply(inputs, sha, "")
  old <- stats::setNames(r$inputs$sha256, r$inputs$path)
  ch <- names(now)[is.na(old[names(now)]) | old[names(now)] != now]
  if (length(ch)) {
    message("    [", label, "] inputs changed: ", paste(basename(ch), collapse = ", "), " -- will run")
    return(TRUE)
  }
  message("    [", label, "] up to date (", format(r$when, "%Y-%m-%d %H:%M"), ")")
  FALSE
}

say <- function(...) { cat(sprintf(...)); flush(stdout()) }

## ---- 4. CONFIG RESOLUTION -----------------------------------------------------
## Merge order: DEFAULTS < dataset input/config.R < explicit call-site args
## (handled by each stage function's own formals, which default to
## `cfg[[name]]` after this merge -- see R/01_dataset_io.R).
resolve_config <- function(dataset_dir) {
  cfg <- DEFAULTS
  cfg_file <- file.path(dataset_dir, "input", "config.R")
  if (file.exists(cfg_file)) {
    env <- new.env(parent = emptyenv())
    sys.source(cfg_file, envir = env)
    override <- as.list(env)
    unknown <- setdiff(names(override), names(DEFAULTS))
    if (length(unknown)) warning(sprintf(
      "%s assigns %d name(s) not in DEFAULTS (typo?): %s",
      cfg_file, length(unknown), paste(unknown, collapse = ", ")))
    cfg[names(override)] <- override
  }
  cfg
}
