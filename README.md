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

Plus a flexible, multi-track Manhattan plot (`ldm_manhattan()`), a pre-scan
diagnostic for phenotype/structure confounding (`check_structure_alignment()`),
a `size_floor` sensitivity sweep (`run_floor_sweep()`), and a batch loop
(`run_all()`) that isolates one dataset's failure from the rest and writes
`summary.csv` + `alignment_summary.csv`.

`ldm_manhattan()` stacks one panel per engine (`ggplot2`/`patchwork`, via
`LDscnR::ld_manhattan()`). Non-significant markers are always the grey
background layer, drawn before -- never over -- the coloured ones. With the
default `colour_by = "region"`, colour is assigned to a genome-wide set of
loci built once across every panel being plotted (the union of all engines'
significant regions, merged where they physically overlap -- different
engines' own region assembly rarely agrees on exact bounds), so the same
underlying region reads as the same colour in every stacked panel, with one
shared, deduplicated legend.

## Quick start

```r
devtools::load_all("~/gitlab/LDscnR")   # see "LDscnR version" below
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

run_dataset("examples/3sp_chr1_chr4")
result <- run_all(c("examples/3sp_chr1_chr4", "examples/9sp_chr19_chr20"))
result$summary            # one row per (dataset, engine)
result$alignment_summary  # one row per dataset with input/population.rds
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
    perm_fun.R              REQUIRED for Stage A -- run_stage_A() errors
                       without it, with no generic fallback. Defines
                       `perm_fun(b, y)` returning a permuted phenotype for
                       permutation index b. There is deliberately no default:
                       an unstratified sample(y), and later a
                       structure_group.rds-based shuffle done at the
                       INDIVIDUAL level, were both tried and both broke the
                       sampling design -- the second was caught by an
                       external audit reproducing, on both worked examples,
                       a shuffle that made a phenotype constant-within-
                       population vary WITHIN populations. What "no signal"
                       means is specific to each dataset's own sampling
                       design, so this repo asks for it explicitly, per
                       dataset (`.dataset_dir` is available to the sourced
                       script). See examples/*/input/perm_fun.R for the
                       population-level, structure-preserving pattern both
                       worked examples actually use.
    population.rds          optional (required only for
                       check_structure_alignment()). Named vector, individual
                       -> population ID.
    structure_group.rds       optional. Named vector, individual -> coarser
                       sampling group (locality/lineage/etc). Defaults to
                       population.rds when absent.
  external_pvalues/       optional. One .rds per extra p-value engine, each a
                       NAMED numeric vector (names = map$marker). The
                       filename (sans .rds) becomes the engine name -- e.g.
                       external_pvalues/lfmm.rds -> engine "lfmm".
  cache/                   gitignored, regenerable: the GDS file, LD-decay fit,
                       stage-1 clustering, GRM. Rebuilt automatically when
                       input/ or size_floor/decay/GRM settings change
                       (content-hash receipts, not timestamps).
  output/                   gitignored, regenerable:
    alignment/                observed.csv, null_draws.csv, null_summary.csv
                       -- see "Structure-alignment diagnostic" below.
    pvalues/emmax_unit/     p_obs.rds/p_perm.rds (the unit-level test) plus
                       p_obs_marker.rds -- a genuine per-marker EMMAX scan
                       (same GRM, same phenotype), computed purely so the
                       "unit" engine's plot can show real per-SNP variation
                       (see snp_results.csv below).
    pvalues/emmax_simes/    p_obs.rds/p_perm.rds
    stageB_<engine>/     outlier_test.rds, outlier_perm.rds (if a permuted
                       null was available), region_rotation.rds (if an
                       annotation was given), snp_results.csv, region_table.csv
    figures/manhattan.{png,pdf}                    -- q-value, per engine
    figures/manhattan_ld_w.{png,pdf}                 -- same panels/faceting/
                       colouring, ld_w_095 (local-LD support, no test
                       involved) on the y-axis instead
    floor_sweep_summary.csv, figures/manhattan_floor_sweep.{png,pdf}
                       -- only after run_floor_sweep().
```

`snp_results.csv` columns: `marker, Chr, Pos, statistic, marker_p, marker_q,
unit_p, unit_q, p_display, q_display, unit_id, region_id, significant,
tested, ld_w_095`. Assignment to `unit_id`/`region_id` is by EXACT
MEMBERSHIP at both levels, never by whether a marker's or unit's coordinate
falls inside a bounding span -- Stage-1 clusters AND the Stage-2 groups
`assembly = "stage2_discovered"` merges them into are both membership-based,
not necessarily contiguous, so a wide span can physically contain markers/
units that actually belong to a different, interleaved cluster or group.
Two related bugs from this, both caught by an external audit and fixed (see
`R/04_stage_B_outlier.R`'s header): joining markers to units by coordinate
mislabelled thousands of markers; separately, even after that fix, joining
a unit's span to a region's span stayed ambiguous whenever assembled
regions overlap physically (quantified: several significant units in each
worked example fit inside more than one region's bounds). `region_id` for
`assembly = "stage2_discovered"` is now recovered by rerunning
`ld_outlier_test()`'s own internal `ld_prune_and_eMLG()` call
(`.stage2_unit_to_region()`) to get the real group membership, not a span
join; `assembly = "physical"` keeps the span join, which is exact there
(its regions are gap-merged directly from the units' own spans, so they
cannot overlap). Every column is also assigned either directly onto a
freshly-`map`-ordered table (before anything can reorder it) or via a
name-keyed update-join, never via a row-order-dependent `merge()` +
manual resort -- an external audit reproduced a marker/p-value swap on a
two-marker case under the earlier approach, which never actually required
`map` to be position-sorted despite silently assuming it was.

`marker_p`/`marker_q` are the marker's own raw value (only meaningful when
p_obs is marker-aligned -- `statistic = "simes"` or any external engine; NA
for `statistic = "unit"`, where no per-marker p exists at all). `unit_p`/
`unit_q` are the actual tested unit's value -- what `significant` is
decided from, for BOTH statistics, broadcast to every member. `p_display`/
`q_display` are what `ldm_manhattan()` plots: identical to `marker_p`/
`marker_q` when those exist; for `statistic = "unit"`, a genuine per-marker
EMMAX scan run purely for display resolution (`p_obs_marker.rds` above),
distinct from the flat broadcast `unit_p`. `ld_w_095` is local-LD support at
rho = 0.95 (computed once per dataset in `build_stage1()`, carried onto
every engine's snp_results.csv for convenience) -- not needed for testing
any more (see `R/02_stage1_cluster.R`), kept as a diagnostic in its own
right: `ldm_manhattan(dataset_dir, value = "ld_w_095")` plots it directly,
same panels/faceting/colouring as the q-value plot, letting a significant
region be compared against local LD independent of any association test.

Config precedence: `R/00_config.R`'s `DEFAULTS` < dataset `input/config.R` <
an explicit `cfg =` argument at the call site.

`size_floor` defaults to a DERIVED value -- 1 unit of floor per 1e5 assayed
markers (`size_floor_per_markers` in `DEFAULTS`), recovering module_3sp's own
hand-derived `SIZE_FLOOR = 8` on its ~790k-marker full panel from a rule
rather than a per-panel judgement call. On this repo's own small worked
examples (~120k markers) it derives a smaller floor than a full dataset
would get -- expected, not a bug: real use of this pipeline assumes full
datasets. Set `size_floor` explicitly in `DEFAULTS` or a dataset's own
`config.R` to pin it instead of deriving it.

## Structure-alignment diagnostic

`check_structure_alignment(dataset_dir)` -- generalised from
`LDscnR-paper/module_9sp/R/02b_env_structure_alignment.R`, the script that
actually produces the manuscript's empirical alignment numbers (confirmed
via `LDscnR_manuscript/generate_empirical_values.R`) -- asks, before any
association test runs, how much of the tested phenotype is explained by (a)
`structure_group.rds` (the sampling groups a permutation null would shuffle
within) and (b) the leading axes of the population-level GRM. A WARNING
diagnostic, not a gate: per the manuscript's own framing, it "reveals when
phenotype and ancestry provide nearly interchangeable explanations," it does
not itself invalidate a result. Runs automatically inside `run_dataset()`
whenever `input/population.rds` exists. Verified to reproduce the
manuscript's published structure-based numbers exactly on both worked
examples (3sp: 5.4%/11.1%; 9sp: 54.8%/72.2%) -- the GRM-axis numbers differ
from the manuscript's own (computed genome-wide) since these examples' GRM
is built from only two chromosomes.

## size_floor sensitivity sweep

`run_floor_sweep(dataset_dir, factors = c(0.5, 1, 2))` reruns Stage A/B at
each of `cfg$size_floor`'s half/default/double (deduplicated by resolved
integer value -- on a small panel where the floor already derives to 1, the
"half" point collapses onto the default rather than wasting a repeat). Only
the "unit" arm is genuinely expensive to sweep (`ld_unit_matrix()`'s columns
depend on `size_floor`, so its permutations must be redone per floor); the
"simes" arm's Stage-A scan doesn't depend on `size_floor` at all, so it's
computed once and reused, receipt-deduplicated, across every floor. Not
wired into `run_dataset()`'s default flow -- call it separately when you
want the extra robustness check.

## Every stage is receipt-gated

Every stage (`build_stage1()`, `run_stage_A()`, each `run_stage_B()` engine,
`check_structure_alignment()`) hashes its own inputs and parameters into
`_receipt.rds` next to its output, and skips recomputing when nothing has
changed -- content hashes, not timestamps, so a `git checkout` or a copied
file can't produce a false "unchanged". Calling `run_dataset()`/`run_all()`
again on unchanged datasets is a no-op walk over receipts. Pass `force =
TRUE` to rebuild regardless. A receipt is also treated as stale (not just
"unchanged params/inputs") when: a recorded output file no longer exists on
disk (an output manually deleted, or a stage that wrote fewer optional
outputs than an earlier receipt promised); or LDscnR's own source has moved
since the receipt was written (recorded in every receipt, not just checked
at call time -- see "LDscnR version" below). `run_stage_A()`'s own receipt
additionally depends on `build_stage1()`'s receipt, and `run_stage_B()`'s
params include a content digest of any directly-supplied `p_obs`/`p_perm`/
`annotation`/`chrom_lengths` (not just their source file, when they have
one) -- both closing gaps an external audit found: a changed LD-decay
setting could otherwise rebuild Stage 1 while Stage A kept stale p-values,
and two calls passing different p-values under the same engine name could
otherwise false-positive as "up to date".

## LDscnR version

`R/00_config.R`'s `check_ldscnr()` pins `LDscnR` to a source-content hash
(`LDSCNR_PIN$src_sha`), not a version number -- the package's `Version` field
is not guaranteed to change on every commit, so a version check alone cannot
catch a stale install, and a stale `LDscnR` install partway through a batch
of hundreds of datasets is exactly the kind of thing worth catching loudly
and immediately. `run_dataset()` calls it with its own default
(`stop_on_fail = TRUE` unless `LDSCNR_LAX` is set) -- a batch stops on the
FIRST call against a mismatched or dirty install, rather than warning and
continuing through the rest against results that may not be reproducible.
`run_batch.R` and the examples above use `devtools::load_all()`, not
`library()`, matching every `LDscnR-paper` module -- `LDscnR` is still under
active development (currently on `main`; the `outlier-scan` branch this was
first pinned to was fast-forwarded into it and is now a stale ancestor, not
a separate line -- `LDSCNR_PIN$branch` tracks whichever is current), and its
working tree is sometimes concurrently edited by other sessions -- confirm
any diff since the last pin is documentation/non-functional (or deliberately
intended) before moving the pin, not just that `check_ldscnr()` currently
passes. USE `devtools::load_all(LDSCNR_PIN$repo)`, NOT `library(LDscnR)`:
the hash check itself always scans the repo's actual source tree regardless
of how the package was attached, but `library(LDscnR)` attaches whatever is
currently INSTALLED -- possibly a stale build with a different `Version` --
so the check can report a clean pin match while the code actually executing
in that session doesn't reflect it. If you update `LDscnR` and this check fails, either update
`LDSCNR_PIN` in `R/00_config.R` (after confirming the change is intentional)
or set `LDSCNR_LAX=1` to proceed anyway.

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

## Known limitations (not yet fixed)

An external audit of this repo (2026-09-26, uncommitted work at the time)
found several correctness/robustness issues; the higher-severity ones
(marker-to-unit/region assignment by coordinate instead of membership, the
unstratified default permutation null, `config.R` overrides not executing at
all, several receipt/provenance gaps) are fixed above. Two it raised are
deliberately deferred, not silently dropped:

- **Memory/plotting at genome scale.** `ldm_manhattan()` plots every marker
  by default, and a full-size dataset's permutation matrices (`p_perm`) can
  be large (the audit estimated ~1.8 GiB for 200 permutations x 1.2M
  markers). Fine on this repo's small worked examples; worth a point-
  thinning/rasterisation pass and a memory-conscious `p_perm` representation
  before running this at genome scale.
- **One dataset folder per (genotype panel, phenotype) pair.** Screening
  many phenotypes against the SAME genotype panel currently means one
  folder per phenotype, each independently rebuilding Stage 1 (genotype-
  only, and so identical across them) and duplicating `genotypes.rds`/
  `map.rds` on disk. Fine when "dataset" means a genuinely different
  genomic panel (this repo's own two examples); a real cost when it means
  many phenotypes on one panel. Splitting genotype-panel-level shared state
  from phenotype-level per-test state is a real architecture change, not a
  bolt-on fix -- worth its own design pass before hundreds of phenotypes on
  one panel is the actual use case.
