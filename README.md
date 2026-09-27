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
a full-permutation `size_floor` sensitivity sweep (`run_floor_sweep()`), a
cheaper observed-only floor-profile alternative for batch use
(`run_floor_profile()`), and a batch loop (`run_all()`) that isolates one
dataset's failure from the rest and writes `summary.csv` +
`alignment_summary.csv`.

`ldm_manhattan()` is a from-scratch `ggplot2`/`patchwork` build (no longer a
wrapper over `LDscnR::ld_manhattan()`'s per-chromosome `facet_wrap()`),
styled directly after the ACTUAL plotting code (not just the rendered
image) of
`~/gitlab/formica_hybrid/module_manuscript_rho05/module_BayPass/R/formica_region_concordance.R`'s
`make_panel()`, which produces
`figures/formica_local_score_vs_ld_regions.png`: one continuous genome-wide
axis (chromosomes concatenated, separated by `chr_gap` [default 0.03, as a
fraction of the largest chromosome's own length], not faceted into
separate panels), with alternating light-grey bands (`chr_shade = TRUE` by
default, one band every other chromosome) marking where each chromosome
starts and ends -- LDscnR-multi's own addition, not in the reference
figure, added because a chromosome with no significant markers of its own
is otherwise invisible in the grey point cloud. Stacked one panel per
engine sharing that axis (shown only on the bottom panel). Non-significant
markers are always the grey background layer (`colour = "grey73"`,
`size = 0.45`, `alpha = 0.45`), drawn before -- never over -- the coloured
ones (`size = point_size` [default 1.6], `alpha = 0.90`). With the default
`colour_by = "region"`, colour is assigned to a genome-wide set of loci
built once across every panel being
plotted (the union of all engines' significant regions, merged where they
physically overlap -- different engines' own region assembly rarely agrees
on exact bounds), so the same underlying region reads as the same colour
in every stacked panel. Colours themselves come from
`.assign_locus_colours()`, porting the reference script's own algorithm:
`LDscnR::default_cluster_colours()`, filtered by standard luminance to
drop near-white entries (invisible against a white background), then
walked in genomic order picking, at each step, the available colour with
the greatest minimum CIE Lab distance from the last few already assigned
-- so two physically nearby loci (the ones most likely to be confused) are
never handed similar hues by chance. Each locus also gets a short
persistent ID ("R1", "R2", ... in genomic order); rather than a legend,
the `max_labels` most significant loci per panel are labelled directly on
the plot via `ggrepel::geom_label_repel()`, using the reference's own
label styling exactly: a WHITE box (`fill = "white"`, `alpha = 0.96`) with
bold, per-locus COLOURED text and border (`colour = locus colour`,
`label.size = 0.2`) -- not a solid-colour box -- and the reference's own
nudge/collision-avoidance tuning (`nudge_y` proportional to the panel's
own y-range, `box.padding = 0.25`, `point.padding = 0.1`,
`min.segment.length = 0`, `max.overlaps = Inf`). The y-axis itself
reserves headroom for labels above the highest point the same way
(`coord_cartesian(ylim = c(lo - 0.04*span, hi + 0.42*span), clip = "off")`).
Labelling is deliberately selective (matching the reference figure): every
significant marker is coloured, but labelling dozens of regions at once
would be unreadable, so each panel labels only its own top `max_labels`
loci by that panel's best q/p value.

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
                       (content-hash receipts, not timestamps). No edge-list
                       cache here (no cache/edge_lists/) -- at full-genome
                       marker density a chromosome's raw pairwise edge list
                       runs several GiB, so it is never persisted; every
                       consumer (ld_complexity_reduction(), compute_ld_w())
                       rebuilds it from the GDS on the fly, one chromosome at
                       a time, discarded right after use. See the `slide`
                       comment in R/00_config.R's DEFAULTS for how its default
                       was calibrated down from compute_LD_decay()'s own 1000
                       to keep each on-the-fly rebuild small.
  output/                   gitignored, regenerable:
    ld_decay/                 decay_summary.csv, decay_print.txt, decay_curves.pdf
                       -- see "LD-decay report" below. Genotype-only, always
                       written (like stage1 itself).
    alignment/                observed.csv, null_draws.csv, null_summary.csv
                       -- see "Structure-alignment diagnostic" below.
    pvalues/emmax_unit/     p_obs.rds/p_perm.rds (the unit-level test) plus
                       p_obs_marker.rds -- a genuine per-marker EMMAX scan
                       (same GRM, same phenotype), computed purely so the
                       "unit" engine's plot can show real per-SNP variation
                       (see snp_results.csv below).
    pvalues/emmax_simes/    p_obs.rds (every marker) + p_perm_compact.rds
                       (permutations, restricted to markers whose cluster
                       clears size_floor -- see "The simes arm's compact
                       permutation scan" below)
    stageB_<engine>/     outlier_test.rds, outlier_perm.rds (if a permuted
                       null was available), region_rotation.rds (if an
                       annotation was given), snp_results.csv, region_table.csv,
                       report.txt -- the console text LDscnR's own
                       print.ld_outlier_test()/print.ld_outlier_perm()/
                       print.ld_region_rotation() methods already produce,
                       captured verbatim (no separate summary logic of this
                       repo's own).
    figures/manhattan.{png,pdf}                    -- q-value, per engine
    figures/manhattan_ld_w.{png,pdf}                 -- same panels/layout/
                       colouring, ld_w_095 (local-LD support, no test
                       involved) on the y-axis instead
    floor_sweep_summary.csv, figures/manhattan_floor_sweep.{png,pdf}
                       -- only after run_floor_sweep().
    floor_profile/            grid.csv, summary.csv, canonical_regions_unit.csv,
                       canonical_regions_simes.csv -- only after
                       run_floor_profile(); see "Floor profile" below.
```

`snp_results.csv` -- one row per marker in `map`:

| column | meaning |
| --- | --- |
| `marker` | Marker ID, from `map$marker`. |
| `Chr`, `Pos` | Chromosome and position, from `map`. |
| `statistic` | `"unit"` or `"simes"` (or whatever `statistic` an external engine was run with) -- which arm produced this row. |
| `marker_p`, `marker_q` | The marker's own raw p-value / BH-adjusted q-value. Only meaningful when `p_obs` is marker-aligned (`statistic = "simes"`, or any external engine); `NA` for `statistic = "unit"`, where no per-marker p-value exists at all. |
| `unit_p`, `unit_q` | The marker's tested unit's p-value / q-value, broadcast to every member marker. This is what `significant` is decided from, for BOTH statistics. |
| `p_display`, `q_display` | What `ldm_manhattan()` actually plots. Identical to `marker_p`/`marker_q` when those exist; for `statistic = "unit"`, a genuine per-marker EMMAX scan run purely for display resolution (`p_obs_marker.rds`, see `pvalues/emmax_unit/` above) -- distinct from the flat broadcast `unit_p`, so the "unit" engine's plot still shows real per-SNP variation instead of one flat value per cluster. |
| `unit_id` | Which Stage-1 cluster (that cleared `size_floor`) this marker belongs to, or `NA` if the marker's cluster didn't clear the floor (never tested). |
| `region_id` | Which Stage-2 (or physical-merge) region this marker's unit was assembled into, or `NA` if its unit wasn't significant / assembled into any region. String form `"<Chr>:<from>-<to>"`. |
| `significant` | Whether this marker's unit cleared `alpha` after BH correction. `NA` when `tested` is `FALSE`. |
| `tested` | Whether this marker's unit cleared `size_floor` and was actually tested (`!is.na(unit_id)`). |
| `ld_w_095` | Local-LD support at rho = 0.95 (see the LD-decay report section below) -- carried from `map`, identical across every engine's `snp_results.csv` for this dataset (it isn't a property of any one test). |

Assignment to `unit_id`/`region_id` is by EXACT
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

`region_table.csv` -- one row per assembled region for this engine (the
`test$regions` LDscnR's `ld_outlier_test()` returns, written as-is):

| column | meaning |
| --- | --- |
| `Chr`, `from`, `to` | The region's genomic span (bp). |
| `n_units` | How many of this engine's significant units were assembled into this region. |
| `n_markers` | How many markers belong to this region's own assembled units (or group, for `assembly = "stage2_discovered"`). |
| `occupancy` | `n_markers` divided by how many markers on `map` physically fall within `[from, to]`. 1.0 means every marker in the region's own physical span actually belongs to it; lower values mean other, unrelated markers are interleaved within the same bounding box -- a lower-occupancy region is a looser, less contiguous signal. |

Config precedence: `R/00_config.R`'s `DEFAULTS` < dataset `input/config.R` <
an explicit `cfg =` argument at the call site.

`size_floor` defaults to a DERIVED value -- 1 unit of floor per 1e5 assayed
markers (`size_floor_per_markers` in `DEFAULTS`), recovering module_3sp's own
hand-derived `SIZE_FLOOR = 8` on its ~790k-marker full panel from a rule
rather than a per-panel judgement call. On this repo's own small worked
examples (~120k markers) it derives a smaller floor than a full dataset
would get -- expected, not a bug: real use of this pipeline assumes full
datasets. Never below **2** regardless of that derivation, though, or of an
explicit override in `DEFAULTS`/a dataset's own `config.R`: a singleton
cluster (`size_floor = 1`) is a single marker's own sampling noise with
nothing to average it against, the single most effect-size-inflated unit
this pipeline could test, and is excluded unconditionally. Set `size_floor`
explicitly to pin it instead of deriving it (still clamped to `>= 2`).

## LD-decay report

`report_ld_decay(dataset_dir)` writes out the `LD_decay` object
`build_stage1()` already computes and caches (`compute_LD_decay()`) -- purely
a reporting layer, no recomputation. Genotype-only, so it runs automatically
inside `run_dataset()` unconditionally, the same as stage1 itself:

- `decay_summary.csv` -- `LD_decay$decay_sum`, one row per chromosome:

  | column | meaning |
  | --- | --- |
  | `Chr` | Chromosome. |
  | `chr_size` | Chromosome span (bp), `max(Pos)` on this chromosome. |
  | `n_snp_chr` | Number of markers on this chromosome. |
  | `bp_per_snp` | `chr_size / (n_snp_chr - 1)` -- average marker spacing, used to convert every SNP-count column below to/from bp. |
  | `a`, `c`, `b` | This chromosome's OWN fitted decay-curve parameters for `r2 ~ b + (c - b) / (1 + a * d)`: `a` is the decay rate (larger = faster decay), `c` the short-range asymptote, `b` the genome-wide background LD floor (same value every row). |
  | `a_pred`, `c_pred` | The genome-wide, chromosome-size-predicted values (a robust regression of `log(a)` on `log(chr_size)` across every chromosome) -- what `ld_w_095` and `decay_curves.pdf`'s "size-predicted" curve actually use, not this chromosome's own noisier `a`/`c`. |
  | `n_w_used` | How many of this chromosome's sliding windows cleared the "structured" contrast-over-background gate and were used to fit `a`/`c`. |
  | `rho_slide_raw`, `rho_slide_pred` | The LD retention (`rho`) the CURRENT `slide` setting actually covers, using this chromosome's own `a` / the size-predicted `a_pred`. |
  | `slide_bp_rho_<X>`, `slide_snp_rho_<X>` | For each target in `rho_targets` (default 0.90/0.95/0.99): the physical window (bp) / marker-count window (SNPs) needed to retain that much LD, using `a_pred`. This is the number the `slide` default in `R/00_config.R` was calibrated against (see its comment there). |

- `decay_print.txt` -- `print(LD_decay)`'s own console output, captured
  verbatim (LDscnR already has `print.ld_decay()`; this repo doesn't
  reimplement a summary, just persists what that method already reports).
- `decay_curves.pdf` -- LDscnR's own `plot.ld_decay()` (not a custom plot):
  one genome-wide `type = "summary"` page (decay rate vs. chromosome size,
  rho covered by the current slide, informative windows per chromosome,
  recommended slide sizes), then one `type = "chr"` page per chromosome
  (window-wise `a`, window-wise `c`, LD contrast across windows, and the
  fitted decay curves themselves against the pooled per-window points). A
  chromosome whose own decay disagrees with what its size alone would
  predict -- an inversion, a region of suppressed recombination -- is
  visible directly on its own page.

The same "own print.* method, captured, not reimplemented" pattern applies
to Stage B: `run_stage_B()` already called `print(test)`/`print(perm)`/
`print(rotation)` (LDscnR's `print.ld_outlier_test()`/`print.ld_outlier_perm()`/
`print.ld_region_rotation()`) for console progress; that same text is now
also captured into `output/stageB_<engine>/report.txt`.

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

Three files, `output/alignment/`:

`observed.csv` -- one row, the observed alignment measures:

| column | meaning |
| --- | --- |
| `n_individuals`, `n_populations`, `n_structure_groups` | Sample sizes: individuals, distinct `population.rds` values, distinct `structure_group.rds` values. |
| `grm_axes_r2_population` | R² of population-mean phenotype ~ the top `n_axes` eigenvectors of the population-level GRM. |
| `grm_axes_adjusted_r2_population` | Same, adjusted R² (penalised for `n_axes`). |
| `structure_r2_population` | R² of population-mean phenotype ~ `structure_group`, population-weighted (one point per population). |
| `structure_r2_individual_weighted` | Same, individual-weighted (one point per individual -- populations with more individuals count more). |

`null_draws.csv` -- one row per (scheme, draw): `scheme` (`"within-group
permutation"` or `"GRM-matched continuous"`), `draw` (1..`n_draws`), and that
draw's own `grm_axes_r2_population`/`structure_r2_population`/
`structure_r2_individual_weighted` (same definitions as `observed.csv`,
recomputed on the null-shuffled/null-surrogate phenotype).

`null_summary.csv` -- one row per (scheme, metric), `observed.csv`'s three
measures compared against their null distribution:

| column | meaning |
| --- | --- |
| `scheme`, `metric` | Which null, and which of the three measures. |
| `observed` | The observed value (from `observed.csv`). |
| `null_mean`, `null_median`, `null_sd` | That null distribution's own summary. |
| `null_q025`, `null_q975` | The null's 2.5th/97.5th percentile. |
| `observed_within_null95` | Whether `observed` falls inside `[null_q025, null_q975]` -- if `TRUE`, the observed alignment is not distinguishable from this null. |
| `proportion_null_ge_observed` | One-sided p-value-like quantity: fraction of null draws at least as large as `observed`. Low = observed alignment is unusually strong relative to this null. |

## The "simes" arm's compact permutation scan

`run_stage_A()`'s "simes" arm tests every marker directly (unlike "unit",
which tests one summary variable per cluster) -- but a marker only ever
contributes to a tested unit's Simes p-value once its own cluster clears
`size_floor` (LDscnR's `.ld_outlier_units()`/`.ld_outlier_tested_units()`).
On a small panel that is nearly every marker; at full-genome scale it can be
a tiny fraction -- on the full 3sp panel (`LDscnR-paper/module_3sp`), 16,890
of 790,578 markers (2.1%) sit in a cluster clearing `size_floor = 8`.
Rescanning the other 97.9% on every one of `b_simes` permutations was pure
waste, and the dominant memory/runtime cost at full-genome scale: found by
timing a full-genome run, where it drove free memory below 200 MB via
`emmax_setup()`'s per-permutation working object (multiplied further by
`cores > 1`'s forked workers, each holding their own copy).

The fix, verified to give bit-for-bit identical results to the old
unrestricted scan (`stage1$clusters`'s membership provably never lets an
ineligible marker's value be read): the **observed** scan
(`p_obs.rds`) still covers every marker, unchanged -- it is also
`marker_p`/`p_display`'s source, needed for full genome-wide manhattan
resolution regardless of which markers are ever tested. Only the
**permutation** scan is restricted, to markers whose cluster clears
`size_floor`, and saved compact (`p_perm_compact.rds`, rows = only those
markers, `rownames()` giving each one). `run_stage_B()` expands it back to
the full-length, `NA`-padded form `ld_outlier_perm()` requires only
transiently, per surrogate (`.simes_perm_accessor()`) -- a dense
markers-by-`B` matrix is never held in memory or written to disk. `NA`
outside the eligible rows is safe and never actually read: exactly the
markers a tested unit's members can be drawn from.

This makes the "simes" arm's permutation scan floor-DEPENDENT where it
previously wasn't (`simes_floor`, a parameter to `run_stage_A()` separate
from `cfg$size_floor`, defaulting to it) -- it is part of that arm's own
receipt params, so a later call at a lower floor (needing a strictly larger
eligible set) correctly invalidates an earlier, narrower scan.
`run_floor_sweep()` accounts for this: every floor in a sweep requests
`simes_floor = min(floors)`, since a lower floor only ever makes MORE
markers eligible, so one compact scan at the sweep's own minimum floor is a
superset covering every higher floor tested afterward too -- preserving the
"simes is cheap to reuse across a floor sweep" design this arm has always
had.

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

`floor_sweep_summary.csv` -- one row per (floor, statistic):

| column | meaning |
| --- | --- |
| `floor` | The resolved integer `size_floor` this row used. |
| `statistic` | `"unit"` or `"simes"`. |
| `n_units_tested` | How many units cleared this floor and were tested. |
| `n_significant` | How many of those were significant at `alpha`. |
| `n_regions` | How many regions the significant units assembled into. |

## Floor profile: a lightweight, observed-only sensitivity check

`run_floor_sweep()` above reruns full permutations at every floor it tests --
fine for a one-off robustness check, not for hundreds of datasets.
`run_floor_profile(dataset_dir, targets = c(0.995, 0.997, 0.999))` is a
separate, cheaper alternative: it runs the pipeline's own, unmodified Stage
A/B functions exactly once **with** permutations (the *canonical* floor), and
twice more **observed-only** (no permutations, no `ld_region_rotation()`) at
two side floors -- all three chosen from the dataset's actual Stage-1
cluster-size distribution, not from marker count. This does not change what
`size_floor` means in `resolve_config()`'s own derivation, or in
`run_floor_sweep()` -- both are untouched.

**`run_dataset()` calls this automatically (`floor_profile = TRUE` by
default)** and adopts its canonical floor as `cfg$size_floor` for that call's
own Stage A/B and manhattan figures, in place of `resolve_config()`'s own
marker-count-based derivation -- the two side floors add only a fast,
permutation-free scan each (see "How the three floors are chosen" below), so
this is cheap even at full-genome scale (this session's own full-genome
timing: canonical Stage A/B already ran in ~10s per call; the two extra
observed-only passes are a small fraction of that, since they skip the
permutation loop entirely). Pass `floor_profile = FALSE` to
`run_dataset()` to keep `resolve_config()`'s own floor unchanged, exactly as
before this feature existed. `run_floor_sweep()` and a direct
`run_floor_profile()` call are both unaffected by this default either way.

**How the three floors are chosen.** For a candidate integer floor `f`,
`reduction(f) = 1 - (clusters with n_loci >= f) / (assayed markers)` -- the
fraction of per-marker tests that clustering-plus-floor avoids running
outright. For each target reduction, the floor whose achieved reduction is
closest is selected (ties -- only possible when no cluster has the exact size
that would separate two candidate floors -- broken toward the **larger**
floor). Verified against the real cluster-size data: 99.7% resolves to floor
8 on `examples/3sp_chr1_chr4`, floor 12 on `examples/9sp_chr19_chr20`. If two
targets resolve to the same integer floor, Stage A/B for it only ever runs
once -- both targets' rows in `grid.csv` still show it, honestly, rather than
hiding the coincidence.

**Only the middle (99.7%, "canonical") floor gets the full treatment**:
permutations, `ld_outlier_perm()`, and (if `annotation` is given)
`ld_region_rotation()` -- the structure-aware null this pipeline's primary
region reporting relies on. The two side floors get a fast, floor-specific
**observed** EMMAX-consensus scan and Stage-2 assembly (`run_stage_A(permute
= FALSE)`), and `ld_outlier_test()`'s own BH correction and assembly, but
never a permutation-based test (`run_stage_B(observed_only = TRUE)`) -- that
mode does not merely skip using a cached permutation file if one is missing,
it never attempts to auto-load one from disk at all, specifically so a
canonical floor's cached `p_perm`/`p_perm_compact.rds` can never be picked up
by an observed-only call by accident. The "simes" arm's genome-wide observed
marker scan is floor-independent and is reused directly from the canonical
floor's own run rather than recomputed.

**Regions are matched across floors by Stage-1 cluster identity, not by
physical span overlap.** A unit's own `unit_id` is a position index into that
floor's own size_floor-filtered cluster subset, and is *not* comparable
across floors -- two different clusters can both be "unit_id 5" at two
different floors. Each cluster's `core_snp` (its representative marker,
assigned once during clustering, independent of any floor) is stable, and is
what "is this canonical region's cluster also significant at the other floor"
is actually asked about. Physical span overlap is used only as a secondary
sanity check on a recovered hit.

Four output files under `output/floor_profile/`:

`grid.csv` -- one row per target: `target, floor, achieved_reduction,
n_units_tested, marker_coverage, is_canonical`.

`summary.csv` -- one compact row per method (`unit`/`simes`): each target's
`floor_<target*1000>`/`reduction_<target*1000>` (e.g. `floor_997`,
`reduction_997`), plus the canonical floor's `n_units_tested_canonical`,
`n_significant_canonical`, `n_regions_canonical`.

`canonical_regions_unit.csv`/`canonical_regions_simes.csv` -- one row per
**canonical-floor** region (never a region found only at a side floor -- this
is about how robust the primary, reported findings are, not a census of
every floor's own discoveries):

| column | meaning |
| --- | --- |
| `region_id`, `Chr`, `from`, `to`, `span` | The canonical region's own identity and physical extent. |
| `tick_<floor>` | TRUE if this region's cluster(s) are *also* a significant unit at that side floor. A cluster whose `n_loci` doesn't clear a side floor is simply not tested there, which correctly reads as not-recovered -- this is not a bug, it's what floor selection means. |
| `overlap_check_<floor>` | Secondary sanity check, only meaningful when the tick above is TRUE: does that side floor's own region (containing the recovered cluster) physically overlap this canonical region's span? Almost always TRUE; a FALSE here is a genuine anomaly worth a human look, not expected behaviour. |
| `n_ticks`, `stability_tier` | How many of the 3 floors (canonical + 2 side) found this region -- `"3/3"`, `"2/3"`, or `"1/3"` (never `"0/3"`, since every row starts from a canonical, already-found region). |

**This count is a descriptive robustness label, not a significance measure**
-- never a C-score, p-value, FDR estimate, or independent replication. Only
the canonical floor's regions get the structure-aware permutation null; a
"3/3" region has not been separately tested against a null at the other two
floors, it has simply also cleared BH correction there. Claiming a 3/3 region
is *statistically* more trustworthy than a 1/3 one would require repeating
the full permutation-based ranking at all three floors -- exactly the cost
this function exists to avoid for a batch of hundreds of datasets.

`run_floor_profile()` itself is not wired into `run_all()`; a user wanting a
multi-dataset profile calls it per dataset the same way `run_all()` calls
`run_dataset()`, and `rbindlist()`s the returned `summary` tables. (Its
*canonical floor*, though, is what `run_dataset()` uses by default -- see above.)

**Manhattan labels show the stability tier directly** (`ldm_manhattan()`,
R/05_manhattan.R): whenever `output/floor_profile/canonical_regions_<statistic>.csv`
exists for an engine's own statistic ("emmax_unit" -> "unit", "emmax_simes"
-> "simes"), that region's label gets the tier appended -- `"R5 (3/3)"`
instead of plain `"R5"`. Matched on that engine's own `region_id` (the peak
labelled marker's), not the cross-engine `locus_id` `.master_loci()` may have
widened by merging several engines' regions together. External engines (e.g.
an `lfmm.rds` in `external_pvalues/`) have no floor-profile concept and keep
plain, unsuffixed labels regardless. Nothing needs to be passed to
`ldm_manhattan()` for this -- it looks for the file itself, so it appears
automatically once `run_floor_profile()` (directly, or via `run_dataset()`'s
own default) has run for that dataset.

**A known, deliberately unresolved discrepancy**: applying this same
99.7%-target rule to the *full* 3sp genome panel (not the small worked
example) resolves to floor **7**, not the floor of **8** the manuscript's
existing full-panel analysis actually used (8 corresponds to a ~99.83%
reduction there, closer to a 99.8% target). `run_floor_profile()` reports
this kind of thing; it never silently revises an existing analysis to match.

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

## Batch summary (`run_all()`)

`run_all(dataset_dirs)` loops `run_dataset()` over every folder, isolates one
dataset's failure from the rest (`status = "error"`, not an abort -- see
`stop_on_error`), and writes two roll-up CSVs.

`summary.csv` -- one row per (dataset, engine):

| column | meaning |
| --- | --- |
| `dataset` | `basename(dataset_dir)`. |
| `engine` | Which p-value engine this row is for (`emmax_unit`, `emmax_simes`, or an external engine's name). `NA` if the dataset had no engine to run, or errored before reaching one. |
| `n_markers` | Total markers in this dataset's `map`. |
| `n_units_tested`, `n_significant`, `n_regions` | Same meaning as `floor_sweep_summary.csv`'s columns, for this engine's actual (non-swept) run. |
| `runtime_s` | Wall-clock seconds for this dataset's whole `run_dataset()` call (shared across every engine row for that dataset -- it isn't per-engine timing). |
| `status` | `"ok"`, `"no_engine"` (ran, but had no p-value source), or `"error"`. |
| `error` | The error message, when `status = "error"`; `NA` otherwise. |
| `align_grm_axes_r2`, `align_grm_axes_p`, `align_structure_r2`, `align_structure_p` | The two headline structure-alignment numbers (see below), repeated onto every engine row for this dataset so a screen of many datasets can sort/filter on association results and confounding risk at once. All `NA` when the dataset has no `input/population.rds`. |

The four `align_*` columns are the **"within-group permutation" scheme**
only (the manuscript's primary null) for two of `null_summary.csv`'s
metrics: `align_grm_axes_r2`/`align_structure_r2` are that metric's
`observed` value, and `align_grm_axes_p`/`align_structure_p` are its
`proportion_null_ge_observed` -- low values flag a phenotype whose signal is
hard to distinguish from ancestry/sampling structure alone.

`alignment_summary.csv` -- one row per dataset with `input/population.rds`
(the detailed companion to `summary.csv`'s four headline columns): `dataset`,
then every `observed.csv` column (see "Structure-alignment diagnostic"
above), then three columns per (scheme, metric) pair from `null_summary.csv`,
named `<scheme>__<metric>__null_mean`, `<scheme>__<metric>__observed_within_null95`,
`<scheme>__<metric>__proportion_null_ge_observed` (`scheme` with any
non-alphanumeric characters collapsed to `_`, e.g.
`within_group_permutation__structure_r2_population__null_mean`).

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

- **Plotting at genome scale.** `ldm_manhattan()` plots every marker by
  default; fine on this repo's small worked examples, worth a point-
  thinning/rasterisation pass before running this at genome scale. (The
  memory-conscious `p_perm` representation this bullet used to ask for is
  done -- see "The simes arm's compact permutation scan" above -- but the
  plotting side of the concern is still open.)
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
