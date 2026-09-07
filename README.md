# LDscnR-multi

Batch driver for [LDscnR](https://github.com/genomeinsights/LDscnR)'s
association + LD-aware outlier-region pipeline: run it across as many
self-contained dataset folders as exist on disk, write the results back into
each folder, and roll up a cross-dataset summary. Built for a comparative
project that will screen hundreds of GT panels; this repo is the looping
code, not a package -- it depends on `LDscnR` being present locally (see
"LDscnR version" below) and sources plain `R/*.R` files.

Two stages, per dataset:

- **Stage A** (`run_stage_A()`): genotypes + map + phenotype -> EMMAX
  observed and permuted-null p-values (two arms: `emmax_unit`, one test per
  LD-cluster; `emmax_simes`, one test per marker aggregated per cluster).
- **Stage B** (`run_stage_B()`): *any* p-value vector -> LDscnR's
  `ld_outlier_test()`/`ld_outlier_perm()` -> per-marker results and a region
  table. Engine-agnostic: Stage A's own arms and any external p-values you
  drop in (LFMM, or anything else that produces one p-value per marker) go
  through the identical call.

Plus a flexible, multi-track Manhattan plot (`ldm_manhattan()`) and a batch
loop (`run_all()`) that isolates one dataset's failure from the rest and
writes a `summary.csv`.

## Quick start

```r
devtools::load_all("~/gitlab/LDscnR")   # see "LDscnR version" below
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

run_dataset("examples/3sp_chr1_chr4")
run_all(c("examples/3sp_chr1_chr4", "examples/9sp_chr19_chr20"))
```

or from the shell:

```bash
Rscript run_batch.R examples/3sp_chr1_chr4 examples/9sp_chr19_chr20
Rscript run_batch.R examples/                 # every subfolder with an input/
```

## Dataset folder contract

```
<dataset_dir>/
  input/
    genotypes.rds     required. n x m numeric dosage matrix (0/1/2, NA ok),
                       individuals in rows (rownames = individual IDs),
                       markers in columns (colnames = marker IDs).
    map.rds             required. data.table/data.frame with marker, Chr, Pos.
                       map$marker must equal colnames(genotypes), same order.
    phenotype.rds        optional -- only needed to run Stage A. A named
                       numeric vector (names = genotypes' rownames), or a
                       1-column data.frame with those rownames.
    config.R              optional. Plain R assigning any DEFAULTS name from
                       R/00_config.R (e.g. `cr_rho <- 0.35`) to override it
                       for this dataset only.
    perm_fun.R              optional. Defines `perm_fun(b, y)` returning a
                       permuted phenotype for permutation index b. Default
                       (file absent): a plain label permutation, sample(y).
  external_pvalues/       optional. One .rds per extra p-value engine, each a
                       NAMED numeric vector (names = map$marker). The
                       filename (sans .rds) becomes the engine name -- e.g.
                       external_pvalues/lfmm.rds -> engine "lfmm".
  cache/                   gitignored, regenerable: the GDS file, LD-decay fit,
                       stage-1 clustering. Rebuilt automatically when input/
                       changes (content-hash receipts, not timestamps).
  output/                   gitignored, regenerable:
    pvalues/emmax_unit/, pvalues/emmax_simes/   Stage A's p_obs.rds/p_perm.rds
    stageB_<engine>/     outlier_test.rds, outlier_perm.rds (if a permuted
                       null was available), region_rotation.rds (if an
                       annotation was given), snp_results.csv, region_table.csv
    figures/manhattan.{png,pdf}
```

`snp_results.csv` columns: `marker, Chr, Pos, statistic, p, q, unit_id,
region_id, significant, tested`. For `statistic = "simes"` (marker-aligned
p-values, including every external engine), `p`/`q` are the marker's own
values. For `statistic = "unit"`, no per-marker p exists -- `p`/`q` are the
marker's cluster's aggregate value, broadcast to every member (documented,
not hidden: it's the evidence that drove that cluster's significance, not
independent per-SNP evidence).

Config precedence: `R/00_config.R`'s `DEFAULTS` < dataset `input/config.R` <
an explicit `cfg =` argument at the call site.

## Every stage is receipt-gated

Every stage (`build_stage1()`, `run_stage_A()`, each `run_stage_B()` engine)
hashes its own inputs and parameters into `_receipt.rds` next to its output,
and skips recomputing when nothing has changed -- content hashes, not
timestamps, so a `git checkout` or a copied file can't produce a false
"unchanged". Calling `run_dataset()`/`run_all()` again on unchanged datasets
is a no-op walk over receipts. Pass `force = TRUE` to rebuild regardless.

## LDscnR version

`R/00_config.R`'s `check_ldscnr()` pins `LDscnR` to a source-content hash
(`LDSCNR_PIN$src_sha`), not a version number -- the package's `Version` field
never changes across commits, so a version check alone cannot catch a stale
install, and a stale `LDscnR` install partway through a batch of hundreds of
datasets is exactly the kind of thing worth catching loudly and immediately.
`run_batch.R` and the examples above use `devtools::load_all()`, not
`library()`, matching every `LDscnR-paper` module -- `LDscnR` is still under
active development on the `outlier-scan` branch. If you update `LDscnR` and
this check fails, either update `LDSCNR_PIN` in `R/00_config.R` (after
confirming the change is intentional) or set `LDSCNR_LAX=1` to proceed
anyway.

## Examples

`examples/3sp_chr1_chr4` and `examples/9sp_chr19_chr20` (see
`examples/build_examples.R`) are real two-chromosome subsets of the
`LDscnR-paper` `module_3sp`/`module_9sp` datasets -- small enough to run the
whole pipeline in a few minutes, real enough to reproduce known biology (the
3sp Chr1 region lands almost exactly on the marine/freshwater inversion
`module_3sp/R_figures/figure_manhattan.R` annotates by hand at
Chr1:21.49-21.94 Mb). `3sp_chr1_chr4` also carries an
`external_pvalues/lfmm.rds`, demonstrating Stage B on a non-EMMAX engine.
9sp has no external engine and a plain-character (not factor) `Chr` column,
which is what surfaced the `SNPRelate::snpgdsGRM()` `autosome.only` fix in
`R/03_stage_A_emmax.R` -- kept as the second example specifically because it
does not look like the first one.

## Scope

LDscnR-only: Stage 2 is engine-agnostic on *p-values*, not on detection
method -- there is no plugin interface here for non-LDscnR outlier callers
(OutFLANK, BayeScan, pcadapt, ...) that a wider comparative project might
also use; those are a separate concern from this repo's job of batch-running
LDscnR itself.
