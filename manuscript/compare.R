## =============================================================================
## LDscnR-multi/manuscript/compare.R
##
## Compares LDscnR-multi's results on the full-genome panels (manuscript/run.R)
## with the manuscript modules' stored results in LDscnR-paper:
##   module_3sp/out/02_bundle/bundle.rds   Stage 1, LD decay, GRM
##   module_3sp/out/03_EMMAX/scan.rds      consensus + Simes: test, null, EcoPeak relocation
##   module_3sp/out/04_lfmm/lfmm.rds       LFMM Simes: test, EcoPeak relocation
##   module_9sp/out/02_bundle/bundle.rds, module_9sp/out/03_EMMAX/scan.rds
## Optionally also with a reference re-run of the module scans under the
## installed LDscnR (--reference-dir=DIR holding ref_3sp.rds / ref_9sp.rds),
## which separates LDscnR-version effects from pipeline differences.
##
## Writes manuscript/comparison.{csv,md} (or --out=STEM).
## Run from the repository root:  Rscript manuscript/compare.R [--reference-dir=DIR]
## =============================================================================
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
REF_DIR <- sub("^--reference-dir=", "", grep("^--reference-dir=", args, value = TRUE))
PAPER <- path.expand(Sys.getenv("LDSCNR_PAPER", "~/gitlab/LDscnR-paper"))
## --stored-3sp=DIR / --stored-9sp=DIR: read the module's scan results (03_EMMAX/,
## 04_lfmm/) from DIR instead of <module>/out, e.g. a backup of earlier outputs.
## The Stage-1 bundle is always read from <module>/out/02_bundle.
stored_dir <- function(panel, module) {
  a <- sub(sprintf("^--stored-%s=", panel), "", grep(sprintf("^--stored-%s=", panel), args, value = TRUE))
  if (length(a)) normalizePath(a, mustWork = TRUE) else file.path(PAPER, module, "out")
}
OUT_STEM <- sub("^--out=", "", grep("^--out=", args, value = TRUE))
if (!length(OUT_STEM)) OUT_STEM <- file.path("manuscript", "comparison")

members_key <- function(cl) {
  cl <- as.data.table(cl)
  vapply(cl$members, function(m) paste(sort(m), collapse = ","), "")
}
sig_key  <- function(t) t$units[significant == TRUE, paste(Chr, from, to)]
reg_key  <- function(t) t$regions[, paste(Chr, from, to)]
maxdiff  <- function(a, b) if (length(a) != length(b)) NA_real_ else max(abs(a - b), na.rm = TRUE)

rows <- list()
add <- function(panel, arm, metric, module, multi, reference = NA) {
  fmt <- function(x) if (is.null(x) || !length(x)) NA_character_ else
    if (is.numeric(x)) format(signif(x, 6), scientific = FALSE, trim = TRUE) else as.character(x)
  rows[[length(rows) + 1L]] <<- data.table(panel, arm, metric, module = fmt(module), multi = fmt(multi),
                                           reference = fmt(reference))
}
arm_rows <- function(panel, arm, mod, mul, ref = NULL) {
  tm <- mod$test; tl <- mul$test; tr <- ref$test
  add(panel, arm, "units tested", nrow(tm$units), nrow(tl$units), if (!is.null(tr)) nrow(tr$units))
  add(panel, arm, "significant units", sum(tm$units$significant), sum(tl$units$significant),
      if (!is.null(tr)) sum(tr$units$significant))
  add(panel, arm, "reported regions", nrow(tm$regions), nrow(tl$regions), if (!is.null(tr)) nrow(tr$regions))
  add(panel, arm, "same significant units", "-", identical(sig_key(tm), sig_key(tl)),
      if (!is.null(tr)) identical(sig_key(tm), sig_key(tr)))
  add(panel, arm, "same region spans", "-", identical(reg_key(tm), reg_key(tl)),
      if (!is.null(tr)) identical(reg_key(tm), reg_key(tr)))
  add(panel, arm, "max |unit p difference|", "-", maxdiff(tm$units$p, tl$units$p),
      if (!is.null(tr)) maxdiff(tm$units$p, tr$units$p))
  if (!is.null(mod$null)) {
    add(panel, arm, "permutation: observed", mod$null$observed, mul$perm$observed, ref$null$observed)
    add(panel, arm, "permutation: null mean", mean(mod$null$surrogates), mean(mul$perm$surrogates),
        if (!is.null(ref$null)) mean(ref$null$surrogates))
    add(panel, arm, "permutation: p", mod$null$p, mul$perm$p, ref$null$p)
    add(panel, arm, "permutation: same surrogate counts", "-",
        identical(as.numeric(mod$null$surrogates), as.numeric(mul$perm$surrogates)),
        if (!is.null(ref$null)) identical(as.numeric(mod$null$surrogates), as.numeric(ref$null$surrogates)))
  }
  if (!is.null(mod$rotation)) {
    add(panel, arm, "EcoPeaks: regions overlapping", mod$rotation$observed, mul$rotation$observed, ref$rotation$observed)
    add(panel, arm, "EcoPeaks: relocation null mean", mod$rotation$null_mean, mul$rotation$null_mean, ref$rotation$null_mean)
    add(panel, arm, "EcoPeaks: fold", mod$rotation$fold, mul$rotation$fold, ref$rotation$fold)
    add(panel, arm, "EcoPeaks: p", mod$rotation$p, mul$rotation$p, ref$rotation$p)
  }
}
stage_b <- function(d, eng) {
  o <- file.path("manuscript", "data", d, "output", paste0("stageB_", eng))
  rd <- function(f) { f <- file.path(o, f); if (file.exists(f)) readRDS(f) }
  list(test = rd("outlier_test.rds"), perm = rd("outlier_perm.rds"), rotation = rd("region_rotation.rds"))
}

for (pp in list(list(panel = "3sp", module = "module_3sp", lfmm = TRUE),
                list(panel = "9sp", module = "module_9sp", lfmm = FALSE))) {
  d <- paste0(pp$panel, "_full")
  cat("[", pp$panel, "]\n")
  ref <- if (length(REF_DIR)) { f <- file.path(REF_DIR, sprintf("ref_%s.rds", pp$panel)); if (file.exists(f)) readRDS(f) }

  ## Stage 1, LD decay and GRM against the module bundle
  b  <- readRDS(file.path(PAPER, pp$module, "out", "02_bundle", "bundle.rds"))
  s1 <- readRDS(file.path("manuscript", "data", d, "cache", "stage1.rds"))
  ds_m <- as.data.table(b$LD_decay$decay_sum); ds_l <- as.data.table(s1$LD_decay$decay_sum)
  num <- intersect(names(ds_m)[vapply(ds_m, is.numeric, TRUE)], names(ds_l))
  add(pp$panel, "Stage 1", "LD decay: max |difference| (decay_sum)", "-",
      max(vapply(num, function(v) maxdiff(ds_m[[v]], ds_l[[v]]), 0)))
  add(pp$panel, "Stage 1", "clusters", nrow(b$stage1$clusters), nrow(s1$stage1$clusters))
  add(pp$panel, "Stage 1", "same cluster membership", "-",
      identical(sort(members_key(b$stage1$clusters)), sort(members_key(s1$stage1$clusters))))
  add(pp$panel, "Stage 1", "GRM markers", length(b$grm_markers), length(unique(na.omit(s1$stage1$pruned))))
  add(pp$panel, "Stage 1", "GRM: max |difference|", "-", maxdiff(unname(b$GRM), unname(s1$GRM)))
  rm(b); invisible(gc())

  sd <- stored_dir(pp$panel, pp$module)
  cat("  module scan results from", sd, "\n")
  scan <- readRDS(file.path(sd, "03_EMMAX", "scan.rds"))
  arm_rows(pp$panel, "EMMAX consensus", scan$consensus, stage_b(d, "emmax_unit"), ref$consensus)
  arm_rows(pp$panel, "EMMAX Simes", scan$simes, stage_b(d, "emmax_simes"), ref$simes)
  if (pp$lfmm) {
    lf <- readRDS(file.path(sd, "04_lfmm", "lfmm.rds"))
    arm_rows(pp$panel, "LFMM Simes", lf, stage_b(d, "lfmm"), ref$lfmm)
  }
}
out <- rbindlist(rows)
out[, match := fifelse(module == "-", multi %in% c("TRUE"), module == multi)]
out[grepl("max \\|", metric), match := !is.na(as.numeric(multi)) & as.numeric(multi) == 0]
fwrite(out, paste0(OUT_STEM, ".csv"))

md <- c("| panel | arm | metric | module (stored) | LDscnR-multi | module re-run, installed LDscnR | match |",
        "|---|---|---|---|---|---|---|",
        out[, sprintf("| %s | %s | %s | %s | %s | %s | %s |", panel, arm, metric, module, multi,
                      fifelse(is.na(reference), "", reference), fifelse(match, "yes", "**no**"))])
writeLines(md, paste0(OUT_STEM, ".md"))
print(out, nrows = 200)
cat(sprintf("\n%d of %d rows match the stored module results\n", sum(out$match), nrow(out)))
