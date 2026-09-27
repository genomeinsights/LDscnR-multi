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
  ## slide = 450, not compute_LD_decay()'s own default of 1000: `slide` sizes
  ## the raw pairwise edge list (every SNP against its next `slide` SNPs, no
  ## r2 floor) that both ld_complexity_reduction() and compute_ld_w() build
  ## per chromosome -- at full-genome marker density, slide=1000 runs 4-5 GiB
  ## PER CHROMOSOME, which is why el_data_folder is no longer used at all
  ## (R/02_stage1_cluster.R) and every consumer rebuilds it on the fly
  ## instead. 450 comes from a calibration run on the full 3sp genome (20
  ## chromosomes, LDscnR-paper/module_3sp): a cheap max_SNPs_decay = 5000 scan
  ## put the tool's own "suggested slide for rho = 0.99" at up to 397 SNPs
  ## (Chr4, the chromosome carrying the ecotype-associated Eda region --
  ## genuinely needs the most). A follow-up max_SNPs_decay sweep (5k/10k/
  ## 20k/40k) at that slide showed the fitted decay rate `a` still rising at
  ## every step, not yet converged -- but monotonically, EVERY chromosome,
  ## always in the direction of a subsampled (sparser) fit reading a SLOWER
  ## decay than the true one. A slower apparent decay implies a LARGER
  ## required window, so the 5000-subsample estimate this 450 is based on is
  ## systematically conservative (oversized), never undersized, relative to
  ## the true full-density value -- confirmed directly on Chr17, whose true
  ## marker count (34,714) is itself below the 40k sweep step, so its value
  ## there already IS the unbiased full-density fit and still followed the
  ## same direction. Safe to use without paying for a full, unsampled
  ## calibration pass.
  decay_args          = list(min_maf_decay = 0.1, q = 0.95, n_sub_bg = 5000,
                              n_win_decay = 20, overlap = 0.5, max_SNPs_decay = Inf,
                              prob_robust = 0.95, max_pairs = 5000, ld_method = "corr",
                              n_strata = 20, slide = 450, cores = 1),
  ## Stage A: EMMAX
  statistics          = c("unit", "simes"),   # which arm(s) to compute; either or both
  unit_repr           = "consensus_dosage",   # ld_unit_matrix() representation for the "unit" arm
  grm_method           = "GCTA",               # SNPRelate::snpgdsGRM method
  b_unit               = 1000L,                 # permutations for the "unit" arm (cheap)
  b_simes              = 200L,                  # permutations for the "simes" arm (rescans every marker)
  ## Stage B: outlier test / permutation null / region rotation
  ## size_floor: NULL means "derive from this dataset's own marker count" --
  ## see resolve_config()'s DERIVED SIZE_FLOOR section below. Set to a number
  ## in DEFAULTS or a dataset's input/config.R to pin it instead.
  size_floor           = NULL,
  size_floor_per_markers = 1e5,   # 1 unit of floor per this many assayed markers
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
## A content hash over R/*.R at the commit this repo was validated against,
## not the package Version string alone: LDscnR is now versioned (0.9.0, was
## 0.0.0.9000), but a version bump is not guaranteed for every commit, so the
## hash remains the check that actually catches a stale install. A batch of
## hundreds of datasets is exactly where that would be most expensive to
## discover late. LDscnR-multi's own repo shares this working tree with other,
## sometimes-concurrent sessions (module_manuscript_rho05, vignette
## reorganisation) -- confirm any diff since the last pin is documentation/
## non-functional (or deliberately intended) before moving this pin, not just
## that check_ldscnr() currently passes.
LDSCNR_PIN <- list(
  repo    = path.expand("~/gitlab/LDscnR"),
  ## outlier-scan (the branch this was first pinned to) was fast-forwarded
  ## into main and development has continued there since -- outlier-scan
  ## itself is now a stale ancestor, not a separate line.
  branch  = "main",
  sha     = "2fc1d44a2552",
  src_sha = "3060402d53eecdb2"
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

## TRUE when the stage must (re)run: no receipt, a changed parameter, a
## changed input, a missing recorded output, or LDscnR's own source having
## moved since the receipt was written. Reports why, so a rerun a user didn't
## expect is explainable.
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
  ## A receipt claiming "up to date" is worthless if what it promises exists
  ## no longer does -- e.g. an output manually deleted, or a stage that wrote
  ## fewer optional outputs on its last run than an earlier receipt recorded.
  missing_out <- r$outputs[nzchar(r$outputs) & !file.exists(r$outputs)]
  if (length(missing_out)) {
    message("    [", label, "] recorded output(s) missing: ",
           paste(basename(missing_out), collapse = ", "), " -- will run")
    return(TRUE)
  }
  ## The receipt recorded which LDscnR source built this result; if that
  ## source has since moved (a real code change, not just a version bump the
  ## package may never make -- see LDSCNR_PIN's own header), the result is
  ## not provably reproducible from the CURRENT package regardless of how
  ## unchanged this stage's own inputs/params are. Recorded, not merely
  ## checked at call time, so this survives being read back by a later run.
  old_src <- tryCatch(r$ldscnr$src_sha, error = function(e) NA_character_)
  if (!is.null(old_src) && !is.na(old_src)) {
    cur <- tryCatch(check_ldscnr(stop_on_fail = FALSE), error = function(e) NULL)
    if (!is.null(cur) && !is.na(cur$src_sha) && !identical(cur$src_sha, old_src)) {
      message("    [", label, "] LDscnR source changed since this receipt (",
             old_src, " -> ", cur$src_sha, ") -- will run")
      return(TRUE)
    }
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
    ## parent = baseenv(), NOT emptyenv(): a config.R sourced into an
    ## environment with no path at all to base R cannot resolve `<-` itself
    ## (verified directly: even `x <- 5` throws "could not find function
    ## '<-'") -- every dataset override was silently failing to execute.
    ## baseenv() gives base functions (assignment, arithmetic, c(), list())
    ## without inheriting from globalenv() or attached packages, keeping the
    ## override sandboxed but actually able to run.
    env <- new.env(parent = baseenv())
    sys.source(cfg_file, envir = env)
    override <- as.list(env)
    unknown <- setdiff(names(override), names(DEFAULTS))
    if (length(unknown)) warning(sprintf(
      "%s assigns %d name(s) not in DEFAULTS (typo?): %s",
      cfg_file, length(unknown), paste(unknown, collapse = ", ")))
    cfg[names(override)] <- override
  }

  ## ---- DERIVED SIZE_FLOOR ------------------------------------------------
  ## 1 unit of floor per size_floor_per_markers assayed markers (default
  ## 1e5): module_3sp's own SIZE_FLOOR = 8, hand-derived as "2x its median
  ## stage-1 cluster size", was fitted to that panel's ~790,578 markers --
  ## 790578 / 1e5 ~= 7.9 ~= 8, so this recovers that value on a full-size
  ## panel from a rule rather than a per-panel judgement call. On a small
  ## worked example (e.g. this repo's own two-chromosome subsets, ~120k
  ## markers) it derives a SMALLER floor than a full dataset would get --
  ## expected, not a bug: real use of this pipeline assumes full datasets,
  ## and a subset's floor being smaller than the full panel's is just what
  ## "fewer markers" means under this rule. Only runs when `size_floor`
  ## wasn't set explicitly (DEFAULTS or a dataset's own config.R); reads
  ## map.rds's row count rather than loading the full genotype matrix.
  if (is.null(cfg$size_floor)) {
    map_f <- file.path(dataset_dir, "input", "map.rds")
    n_markers <- if (file.exists(map_f)) nrow(readRDS(map_f)) else NA_integer_
    cfg$size_floor <- if (is.na(n_markers)) 8L
                      else round(n_markers / cfg$size_floor_per_markers)
  }
  ## Hard floor of 2, unconditionally -- a singleton "cluster" (size_floor =
  ## 1) is the single most effect-size-inflated, least reliable unit this
  ## pipeline can test (one marker's own sampling noise, with nothing to
  ## average it against), and removing singletons has the single biggest
  ## effect on result quality of anything size_floor controls. Applied AFTER
  ## the derivation above, and to an explicit override from DEFAULTS/a
  ## dataset's own config.R too, so "regardless of how many SNPs are
  ## analysed" is a genuine invariant, not just this rule's own default.
  cfg$size_floor <- max(2L, as.integer(cfg$size_floor))
  cfg
}
