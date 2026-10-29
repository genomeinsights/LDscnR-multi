# LDscnR-multi

LDscnR-multi is a batch workflow for
[LDscnR](https://github.com/genomeinsights/LDscnR), the R package for genome
scans that account for linkage disequilibrium (LD). It runs the LDscnR
association and outlier-region analysis on any number of self-contained
dataset folders, writes the results back into each folder and combines them
into a summary across data sets. It works for one data set or for hundreds.
This repository is a set of R scripts, not a package. It needs LDscnR
installed (see [LDscnR version](#ldscnr-version)) and sources the plain
`R/*.R` files.

LDscnR-multi version 0.9.1 matches LDscnR v0.9.1, the version used in the
manuscript (see [How to cite](#how-to-cite)).

On one core of an Apple M3 Pro laptop (36 GB memory, R 4.6.1), a complete run
of `run_dataset()` took about 19 minutes for the full three-spined stickleback
panel (790,578 markers) and about 36 minutes for the full nine-spined
stickleback panel (1,195,557 markers), including the Manhattan plots. The
manual breaks these times down by step. These are measured examples, not
promised run times. The association scans are fast because `emmax_fast()`
reuses the fitted relationship model, whereas the Stage-1 LD calculations and
the floor profile take most of the time.

## What the workflow does

The manuscript describes the method in two stages. Stage 1 groups markers
into Stage-1 units from genotype LD alone, without the phenotype. Stage 2
assembles the significant Stage-1 units into outlier regions. The workflow
runs these steps for each data set:

1. **Stage 1** (`build_stage1()`) builds the Stage-1 units, the LD-decay fit
   and the genomic relationship matrix (GRM).
2. **Association scan** (`run_stage_A()`) uses EMMAX (efficient mixed-model
   association expedited; Kang et al. 2010) to compute observed p-values and
   p-values for permuted phenotypes. It has two tests. The *consensus* test
   fits one association per Stage-1 unit. The *Simes* test fits one
   association per marker and combines the marker p-values within each
   Stage-1 unit with Simes' method.
3. **Unit testing and Stage-2 assembly** (`run_stage_B()`) decides which
   Stage-1 units are significant and assembles them into outlier regions. It
   accepts the association scan's results or marker-wise p-values from
   another method, such as LFMM2 (latent factor mixed model; Caye et al.
   2019).

The function names `run_stage_A()` and `run_stage_B()` are kept for
compatibility. They are not the manuscript's Stage 1 and Stage 2. The table
below lists the other names that differ between the code and the manuscript.

| In the code and output files | In the manuscript and this README |
| --- | --- |
| `run_stage_A()` | association scan |
| `run_stage_B()`, `stageB_<engine>/` | unit testing and Stage-2 assembly |
| statistic `unit`, engine `emmax_unit` | consensus test |
| statistic `simes`, engine `emmax_simes` | Simes test |
| `ld_region_rotation()`, `n_rotations`, `region_rotation.rds` | relocation test, number of relocations |
| `canonical_regions_*.csv`, "canonical floor" | outlier regions at the reported floor |
| cluster, `unit_id`, `core_snp` | Stage-1 unit, its index and its representative marker |
| `size_floor` | minimum Stage-1 unit size (the floor) |

The workflow also provides:

- a multi-track Manhattan plot (`ldm_manhattan()`);
- a check of how closely the phenotype follows population structure, run
  before the association scan (`check_structure_alignment()`);
- floor-stability labels for every outlier region (`run_floor_profile()`);
- a sensitivity sweep over the floor with full permutations
  (`run_floor_sweep()`);
- a batch loop (`run_all()`) that keeps one data set's failure from stopping
  the others and writes `summary.csv` and `alignment_summary.csv`.

The Manhattan plots show every method on one genome-wide axis. Grey points
are background markers, and coloured points are markers in outlier regions.
Nearby regions get contrasting colours, and the strongest regions are
labelled directly. The manual describes the plotting options.

## Installation

LDscnR-multi needs R (tested with 4.6.1) and the packages below. The minimum
versions are the versions it was tested with. Earlier versions may work but
have not been tested.

| Package | Minimum version | Source |
| --- | --- | --- |
| LDscnR | 0.9.1 (GitHub tag `v0.9.1`) | GitHub |
| data.table | 1.18.2 | CRAN |
| digest | 0.6.39 | CRAN |
| ggplot2 | 4.0.3 | CRAN |
| ggrepel | 0.9.8 | CRAN |
| patchwork | 1.3.2 | CRAN |
| igraph | 2.3.1 | CRAN |
| remotes | 2.5.0 | CRAN |
| SNPRelate | 1.46.0 | Bioconductor |
| gdsfmt | 1.48.0 | Bioconductor |

The base packages `parallel`, `stats`, `tools`, `utils` and `grDevices` come
with R. `devtools` is needed only by developers who load an uninstalled LDscnR
checkout.

```r
install.packages(c("remotes", "BiocManager", "data.table", "digest",
                   "ggplot2", "ggrepel", "patchwork", "igraph"))
BiocManager::install(c("SNPRelate", "gdsfmt"))
remotes::install_github("genomeinsights/LDscnR@v0.9.1")
```

Then download or clone this repository. LDscnR-multi itself is not installed.

## Quick start

The same `run_dataset()` call works for one data set, whereas `run_all()` or
`run_batch.R` loops over many folders. The full guide is the manual,
[`doc/manual.pdf`](doc/manual.pdf).

```r
library(LDscnR)
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

run_dataset("examples/3sp_chr1_chr4")
result <- run_all(c("examples/3sp_chr1_chr4", "examples/9sp_chr19_chr20"))
result$summary            # one row per data set and engine
result$alignment_summary  # one row per data set with input/population.rds
```

From the shell:

```bash
Rscript run_batch.R examples/3sp_chr1_chr4 examples/9sp_chr19_chr20
Rscript run_batch.R examples/                 # every subfolder with an input/
```

`run_batch.R` writes `summary.csv` and `alignment_summary.csv` to the
repository root, whereas an interactive `run_all()` writes them to the
working directory.

### Very dense panels

For very dense panels, consider thinning to about one million well-spread
markers before you create the input files. Use a recorded random seed and
do not use the phenotype. Sample across chromosomes and short physical bins,
and keep several markers in each bin. Do not keep only one marker per Stage-1
unit, because that defeats unit-level testing. Apply the same marker subset
to the map, the genotype matrix and any external p-value vectors. One million
markers is a practical run-time target, **not** a guarantee of detection,
because thinning can change the Stage-1 units and remove a real signal. If
possible, compare one dense and one thinned analysis.

## Examples

`examples/3sp_chr1_chr4` and `examples/9sp_chr19_chr20` are two-chromosome
subsets of the three-spined and nine-spined stickleback panels analysed in
the manuscript (data from Fang et al. 2021). They are small enough to run the
whole workflow in a few minutes. The three-spined example reproduces known
biology, because its strongest Chr1 outlier region (Chr1:21.49-21.94 Mb) lies
on the known marine-freshwater inversion. It also includes
`external_pvalues/lfmm.rds`, which shows unit testing on p-values from a
method other than EMMAX. The nine-spined example has no external p-values,
but it has `input/covariate.rds` (lineage), as in the manuscript's
nine-spined analysis.

`examples/build_examples.R` regenerates the example inputs from the
manuscript modules' bundles. Only maintainers need it (see its header).

## Dataset folder contract

```
<dataset_dir>/
  input/
    genotypes.rds        required. Numeric dosage matrix (0, 1, 2; NA allowed),
                         individuals in rows (row names = individual IDs) and
                         markers in columns (column names = marker IDs).
    map.rds              required. data.table or data.frame with columns marker,
                         Chr and Pos. map$marker must equal colnames(genotypes),
                         in the same order.
    phenotype.rds        needed for the association scan. Named numeric vector
                         (names = row names of genotypes), or a one-column
                         data.frame with those row names.
    perm_fun.R           required for the association scan. Defines
                         perm_fun(b, y), which returns one permuted phenotype
                         for permutation b. There is no default (see below).
    population.rds       needed for the structure-alignment check. Named
                         vector from individual to population ID.
    structure_group.rds  optional. Named vector from individual to a coarser
                         sampling group (locality or lineage, for example).
                         Defaults to population.rds when absent.
    covariate.rds        optional. Named vector from individual to a structure
                         covariate, such as lineage. When present, the
                         phenotype is replaced by resid(lm(y ~ covariate))
                         before EMMAX, and every permuted phenotype is
                         residualised in the same way.
    config.R             optional. Plain R that assigns any DEFAULTS name from
                         R/00_config.R (for example cr_rho <- 0.35) to
                         override it for this data set only.
  external_pvalues/      optional. One .rds per extra p-value engine, each a
                         named numeric vector (names = map$marker). The file
                         name without .rds is the engine name, so
                         external_pvalues/lfmm.rds is engine "lfmm".
  cache/                 written by the workflow and regenerable (see below).
  output/                written by the workflow and regenerable (see below).
```

**The permutation function.** There is no default permutation, because what
"no signal" means depends on each data set's sampling design. A shuffle of
individuals can, for example, create variation within populations that never
occurs in the real data. Each data set therefore defines its own
`perm_fun.R`. The sourced script can use `.dataset_dir`. Both examples use
`perm_within_groups()` (in `R/01_dataset_io.R`), which reassigns whole
populations' values among populations in the same structure group. This is
the structure-preserving permutation used in the manuscript.

**`cache/`** holds the genotype-only results: the SNPRelate genotype file in
GDS (Genomic Data Structure) format, the LD-decay fit, the Stage-1 units and
the GRM. They are rebuilt automatically when the inputs or the relevant
settings change. Pairwise LD edge lists are not stored, because at
full-genome density one chromosome's edge list can take several gigabytes.
Every step that needs them rebuilds them from the GDS file, one chromosome at
a time. The `slide` setting in `R/00_config.R` sets their size.

**`output/`** holds:

| Folder or file | Contents |
| --- | --- |
| `ld_decay/` | `decay_summary.csv`, `decay_print.txt`, `decay_curves.pdf` (see [LD-decay report](#ld-decay-report)) |
| `alignment/` | `observed.csv`, `null_draws.csv`, `null_summary.csv` (see [Structure-alignment check](#structure-alignment-check)) |
| `pvalues/emmax_unit/` | `p_obs.rds` and `p_perm.rds` for the consensus test, plus `p_obs_marker.rds`, a marker-wise EMMAX scan used only for plotting |
| `pvalues/emmax_simes/` | `p_obs.rds` for every marker, and `p_perm_compact.rds`, the permutations for markers in tested Stage-1 units only |
| `stageB_<engine>/` | `outlier_test.rds`, `outlier_perm.rds` (with a permutation null), `region_rotation.rds` (with an annotation), `snp_results.csv`, `region_table.csv` and `report.txt` |
| `floor_profile/` | `grid.csv`, `summary.csv`, `canonical_regions_unit.csv`, `canonical_regions_simes.csv` (see [Floor profile](#floor-profile)) |
| `figures/` | `manhattan.{png,pdf}` (q-values) and `manhattan_ld_w.{png,pdf}` (local LD support) |
| `floor_sweep_summary.csv`, `figures/manhattan_floor_sweep.{png,pdf}` | only after `run_floor_sweep()` |

`report.txt` is the console summary that LDscnR's own print methods produce,
saved as it is.

### `snp_results.csv`

One row per marker in `map`:

| Column | Meaning |
| --- | --- |
| `marker` | Marker ID, from `map$marker`. |
| `Chr`, `Pos` | Chromosome and position, from `map`. |
| `statistic` | `"unit"` (consensus) or `"simes"`, the test that produced this row. |
| `marker_p`, `marker_q` | The marker's own p-value and Benjamini-Hochberg (BH) q-value. Defined only for marker-wise p-values (the Simes test or an external engine). `NA` for the consensus test, which has no marker-wise p-value. |
| `unit_p`, `unit_q` | The p-value and q-value of the marker's Stage-1 unit, copied to every member marker. Significance is decided from these for both tests. |
| `p_display`, `q_display` | The values that `ldm_manhattan()` plots. They equal `marker_p` and `marker_q` when those exist. For the consensus test they come from the separate marker-wise EMMAX scan (`p_obs_marker.rds`), so the plot shows variation between markers instead of one flat value per unit. |
| `unit_id` | The marker's tested Stage-1 unit, or `NA` if its unit is smaller than the floor and was not tested. |
| `region_id` | The marker's outlier region, as `"<Chr>:<from>-<to>"`, or `NA` if its unit was not significant. |
| `significant` | Whether the marker's unit is significant at `alpha` after BH correction. `NA` when `tested` is `FALSE`. |
| `tested` | Whether the marker's unit reached the floor and was tested (`!is.na(unit_id)`). |
| `ld_w_095` | Local LD support at rho = 0.95 (see [LD-decay report](#ld-decay-report)). It depends only on the genotypes, so it is the same in every engine's file. |

Markers are assigned to `unit_id` and `region_id` by Stage-1 and Stage-2
membership, not by coordinates. Stage-1 units and outlier regions are not
always physically contiguous, so a region's span can contain markers that
belong to another unit or region.

### `region_table.csv`

One row per outlier region for this engine, as returned by LDscnR's
`ld_outlier_test()`:

| Column | Meaning |
| --- | --- |
| `Chr`, `from`, `to` | The region's span in base pairs. |
| `n_units` | The number of significant Stage-1 units assembled into the region. |
| `n_markers` | The number of markers in those units. |
| `occupancy` | `n_markers` divided by the number of markers in `map` between `from` and `to`. A value of 1 means every marker in the span belongs to the region. Lower values mean that other markers are interleaved, so the signal is less contiguous. |

## Configuration and the floor

Settings are resolved in this order, with later ones taking precedence:
`DEFAULTS` in `R/00_config.R`, then the data set's `input/config.R`, then a
`cfg =` argument in the function call. The manual lists every setting.

The minimum Stage-1 unit size (the floor, `size_floor`) is derived from
marker density unless it is set explicitly. The floor is the number of
markers per Mb divided by 250 (`size_floor_markers_per_mb`), rounded. Density
is the number of markers divided by the summed span of the chromosomes in
`map.rds`. Because it is a density rather than a count, a subset of
chromosomes gets the same floor as the full panel. The floor is 8 for the
three-spined panel (about 2,080 markers per Mb) and for its two-chromosome
example, and 12 for the nine-spined panel (about 2,960 markers per Mb). The
floor is never below two, even when it is set explicitly. A single-marker
unit has nothing to average its sampling noise against, so its effect size
is the most inflated.

## LD-decay report

`report_ld_decay(dataset_dir)` writes out the LD-decay fit that
`build_stage1()` computes with LDscnR's `compute_LD_decay()`. It recomputes
nothing and depends only on the genotypes, so `run_dataset()` always runs it.

`decay_summary.csv` has one row per chromosome:

| Column | Meaning |
| --- | --- |
| `Chr` | Chromosome. |
| `chr_size` | Chromosome span in base pairs (`max(Pos)`). |
| `n_snp_chr` | Number of markers on the chromosome. |
| `bp_per_snp` | Average marker spacing, `chr_size / (n_snp_chr - 1)`, used to convert marker counts to base pairs. |
| `a`, `c`, `b` | The chromosome's fitted parameters of the decay curve r² = b + (c - b) / (1 + a·d), where d is distance. `a` is the decay rate (larger means faster decay), `c` the short-range level, and `b` the genome-wide background LD (the same in every row). |
| `a_pred`, `c_pred` | The values predicted from chromosome size by a robust regression of log(a) on log(chr_size) across chromosomes. `ld_w_095` and the size-predicted curve in `decay_curves.pdf` use these. With fewer than three distinct chromosome sizes, as in the two-chromosome examples, the regression cannot be fitted and the median of `a` is used. |
| `n_w_used` | The number of windows on the chromosome that passed the contrast-over-background check and were used in the fit. |
| `rho_slide_raw`, `rho_slide_pred` | The share of LD (rho) that the current `slide` setting retains, using `a` and `a_pred`, respectively. |
| `slide_bp_rho_<X>`, `slide_snp_rho_<X>` | For each target in `rho_targets` (default 0.99), the window in base pairs and in markers needed to retain that share of LD, using `a_pred`. |

`decay_print.txt` is the console output of LDscnR's `print.ld_decay()`.
`decay_curves.pdf` is LDscnR's `plot.ld_decay()`. Its first page summarises
the whole genome (decay rate against chromosome size, LD retained by the
current `slide`, informative windows per chromosome and recommended `slide`
sizes). Each chromosome then gets a page with its window-wise fits. A
chromosome whose decay differs from what its size predicts, because of an
inversion or suppressed recombination, stands out on its page.

## Structure-alignment check

`check_structure_alignment(dataset_dir)` asks, before any association test,
how much of the phenotype is explained by (a) `structure_group.rds`, the
sampling groups within which the permutation null shuffles, and (b) the
leading axes of the GRM computed between populations. It follows the
manuscript's alignment analysis. It is a warning, not a gate. As the
manuscript puts it, the check shows when phenotype and ancestry offer nearly
interchangeable explanations, but it does not invalidate a result.
`run_dataset()` runs it when both `input/population.rds` and
`input/phenotype.rds` exist. It uses the raw phenotype, not the phenotype
adjusted for `covariate.rds`. On both two-chromosome examples it reproduces
the manuscript's values for structure groups (three-spined: 5.4% weighted by
population and 11.1% weighted by individual; nine-spined: 54.8% and 72.2%).
The values for the GRM axes differ from the manuscript's, because the
examples' GRM comes from two chromosomes only.

It writes three files to `output/alignment/`.

`observed.csv` has one row with the observed measures:

| Column | Meaning |
| --- | --- |
| `n_individuals`, `n_populations`, `n_structure_groups` | The numbers of individuals, populations and structure groups. |
| `grm_axes_r2_population` | R² of the population-mean phenotype regressed on the top `n_axes` eigenvectors (default five) of the between-population GRM. |
| `grm_axes_adjusted_r2_population` | The same, as adjusted R². |
| `structure_r2_population` | R² of the population-mean phenotype regressed on `structure_group`, with one point per population. |
| `structure_r2_individual_weighted` | The same with one point per individual, so populations with more individuals count more. |

`null_draws.csv` has one row per null scheme and draw (`n_draws`, default
1,000). It holds the same three measures computed on a permuted phenotype
(`"within-group permutation"`) or on a continuous surrogate matched to the
GRM (`"GRM-matched continuous"`).

`null_summary.csv` compares each observed measure with its null:

| Column | Meaning |
| --- | --- |
| `scheme`, `metric` | The null scheme and the measure. |
| `observed` | The observed value. |
| `null_mean`, `null_median`, `null_sd` | Summary of the null distribution. |
| `null_q025`, `null_q975` | The 2.5th and 97.5th percentiles of the null. |
| `observed_within_null95` | `TRUE` if the observed value lies between those percentiles, so it cannot be distinguished from the null. |
| `proportion_null_ge_observed` | The fraction of null draws at least as large as the observed value. Low values mean the observed alignment is unusually strong. |

## Floor profile

`run_floor_profile(dataset_dir, factors = c(0.5, 2))` labels each outlier
region found at the reported floor by whether it is also recovered at half
and at twice that floor. The labels are `"3/3"` (recovered at both),
`"2/3"` (at one) and `"1/3"` (at neither). For the three-spined panel, with a
floor of 8, the comparison floors are 4 and 16. The low comparison floor may
be 1 (single markers included), but it is used only for comparison.

`run_dataset()` runs the floor profile by default. It never changes the
reported floor, so `floor_profile = FALSE` only skips the labels.

**Recovery rule.** A region is recovered at a comparison floor when at least
one of its significant Stage-1 units is also significant at that floor and
the region overlaps an outlier region assembled at that floor. This is the
manuscript's rule. Units are matched across floors by their representative
marker (`core_snp`), because `unit_id` is only an index into one floor's own
list of tested units.

**Why half and twice the floor.** In the manuscript's simulations (600 runs
with known causal variants and a floor of two), the precision of outlier
regions was about 0.2 for regions labelled 1/3, about 0.5 for 2/3 and about
0.8 for 3/3. Dropping 1/3 regions cost little recall. Recovery at the lower
floor, where more units are tested and the BH threshold is harder to reach,
was the more informative check. In the three-spined panel, the EcoPeak
enrichment (the overlap with the marine-freshwater differentiated regions of
Kingman et al. 2021) was carried by regions recovered at a second floor.

**Cost.** The reported floor runs the association scan and unit testing
exactly as `run_dataset()` does, with permutations and, if an annotation is
given, the relocation test. The comparison floors are tested on the observed
data only (`run_stage_B(observed_only = TRUE)`). The consensus test is fitted
once more, at the low floor. The units of the high floor are a subset of those
of the low floor, and a unit's EMMAX p-value does not depend on which other
units are tested. The high floor therefore reuses those p-values, and only
the BH correction and Stage-2 assembly are repeated. The Simes test reuses
the marker-wise scan at every floor.

It writes four files to `output/floor_profile/`.

`grid.csv` has one row per floor: `role` (low, reported or high), `floor`,
`achieved_reduction` (the fraction of marker-wise tests avoided),
`n_units_tested` and `marker_coverage`.

`summary.csv` has one row per test (`unit` or `simes`): the three floors and
their reductions, the number of units tested, significant units and outlier
regions at the reported floor, and the number of regions in each tier
(`n_tier_1`, `n_tier_2`, `n_tier_3`).

`canonical_regions_unit.csv` and `canonical_regions_simes.csv` have one row
per outlier region at the reported floor. Regions found only at a comparison
floor are not listed, because the profile describes how robust the reported
regions are.

| Column | Meaning |
| --- | --- |
| `region_id`, `Chr`, `from`, `to`, `span` | The region and its span. |
| `tick_<floor>` | `TRUE` if the region is recovered at that comparison floor. A unit smaller than a comparison floor is not tested there, so it cannot be recovered there. |
| `n_ticks`, `stability_tier` | The number of floors, out of three, at which the region is found, as a number and as `"3/3"`, `"2/3"` or `"1/3"`. |

**The tier is a robustness label, not a measure of significance.** It is not
a p-value, a false discovery rate (FDR) or an independent replication,
because every floor analyses the same data and only the reported floor gets
the structure-preserving permutation null. Its meaning comes from the
simulations above.

`run_all()` gets the labels through `run_dataset()`. To profile many data
sets without the rest of the workflow, call `run_floor_profile()` for each
and combine the returned `summary` tables with `rbindlist()`.

`ldm_manhattan()` adds the tier to the labels of the consensus and Simes
panels whenever `output/floor_profile/canonical_regions_<statistic>.csv`
exists, for example `"R5 (3/3)"` instead of `"R5"`. External engines, such as
LFMM, have no floor profile and keep plain labels.

## Simes permutations at full-genome scale

The Simes test computes a p-value for every marker, but a marker only
contributes to a tested unit when its Stage-1 unit reaches the floor. In the
full three-spined panel, 16,890 of 790,578 markers (2.1%) are in units of at
least eight markers. The observed scan (`p_obs.rds`) still covers every
marker, because the Manhattan plot needs them. The permutation scan, however,
is restricted to the markers in tested units and saved in compact form
(`p_perm_compact.rds`, with one row per eligible marker). `run_stage_B()`
expands it one permutation at a time, so a dense matrix of all markers by all
permutations is never held in memory. The results are identical to an
unrestricted scan.

The permutation scan therefore depends on the floor. `run_stage_A()` takes
`simes_floor` (by default `cfg$size_floor`), and a call with a lower floor
rebuilds the scan. `run_floor_sweep()` and `run_floor_profile()` request the
lowest floor they need once and reuse that scan for the higher floors.

## Floor sensitivity sweep

`run_floor_sweep(dataset_dir, factors = c(0.5, 1, 2))` reruns the association
scan and unit testing, with full permutations, at half, one and two times the
floor. Each floor is at least two, and factors that give the same floor are
run once. The consensus test is expensive to sweep, because its tested units,
and therefore its permutations, change with the floor. The Simes test reuses
one marker-wise scan. `run_dataset()` does not run the sweep, so call it
separately when you want this more thorough check.

`floor_sweep_summary.csv` has one row per floor and test, with the columns
`floor`, `statistic`, `n_units_tested`, `n_significant` and `n_regions`.

## Receipts: every step is skipped when nothing changed

Every step (`build_stage1()`, `run_stage_A()`, each `run_stage_B()` engine
and `check_structure_alignment()`) saves a hash of its inputs and settings in
`_receipt.rds` next to its output. It skips the work when nothing has
changed. The floor profile and the floor sweep are made of these steps, so
they use the same receipts. The hashes are of file contents, not timestamps, so
a `git checkout` or a copied file cannot make a step look unchanged. Running
`run_dataset()` or `run_all()` again on unchanged data sets only reads the
receipts. Pass `force = TRUE` to rebuild anyway.

A receipt is also out of date when:

- one of its recorded output files is missing;
- the LDscnR revision has changed since it was written (see
  [LDscnR version](#ldscnr-version));
- for the association scan, the Stage-1 receipt has changed;
- for unit testing, p-values or annotations passed directly to
  `run_stage_B()` differ in content from those recorded.

## Batch summary (`run_all()`)

`run_all(dataset_dirs)` runs `run_dataset()` on every folder. A failing data
set is recorded with `status = "error"` and the others still run (see
`stop_on_error`). It writes two summary files.

`summary.csv` has one row per data set and engine:

| Column | Meaning |
| --- | --- |
| `dataset` | `basename(dataset_dir)`. |
| `engine` | The p-value engine (`emmax_unit`, `emmax_simes` or an external engine). `NA` if the data set had no engine or failed before reaching one. |
| `n_markers` | The number of markers in the data set's `map`. |
| `n_units_tested`, `n_significant`, `n_regions` | Units tested, significant units and outlier regions for this engine at the reported floor. |
| `runtime_s` | Wall-clock seconds for the whole `run_dataset()` call. It is the same on every row of a data set. |
| `status` | `"ok"`, `"no_engine"` (the data set ran but had no p-values to test) or `"error"`. |
| `error` | The error message when `status = "error"`, and `NA` otherwise. |
| `align_grm_axes_r2`, `align_grm_axes_p`, `align_structure_r2`, `align_structure_p` | The two main structure-alignment measures and their tail proportions, repeated on every row of the data set, so that many data sets can be sorted on association results and on the risk of confounding at once. `NA` without `input/population.rds`. |

The four `align_*` columns use the within-group permutation null, the
manuscript's primary null. `align_grm_axes_r2` and `align_structure_r2` are
the observed values, and `align_grm_axes_p` and `align_structure_p` are the
corresponding `proportion_null_ge_observed`. Low values flag a phenotype
whose signal is hard to separate from ancestry or sampling structure.

`alignment_summary.csv` has one row per data set with `input/population.rds`.
It holds `dataset`, every column of `observed.csv` and three columns for each
scheme and measure in `null_summary.csv`, named
`<scheme>__<metric>__null_mean`, `<scheme>__<metric>__observed_within_null95`
and `<scheme>__<metric>__proportion_null_ge_observed`. Characters other than
letters and digits in the scheme name become `_`, for example
`within_group_permutation__structure_r2_population__null_mean`.

## LDscnR version

`R/00_config.R` records the validated LDscnR GitHub commit (tag `v0.9.1`).
`check_ldscnr()` checks the installed package's GitHub revision before a data
set runs, and an older or non-GitHub installation stops the batch. If that
happens, reinstall the revision from [Installation](#installation). A version
number alone is not enough, because package code can change without a
version bump.

Developers who need an uninstalled local checkout can set
`LDSCNR_DEV_LOAD_ALL=1`, and `LDSCNR_REPO` to the checkout's path, before
`Rscript run_batch.R`. This uses `devtools::load_all()` and checks the local
source hash instead. Ordinary users need neither an LDscnR checkout nor
`devtools`.

## Reproducing the manuscript

`manuscript/` runs this workflow on the two full stickleback panels of the
manuscript and compares every result with the manuscript's analysis modules
(`module_3sp` and `module_9sp` in the LDscnR-paper repository). The defaults
match the modules' settings:

| Setting | Value |
| --- | --- |
| Stage-1 LD threshold | rho = 0.50 |
| LD-decay window (`slide`) | 500 markers |
| GRM | GCTA (Yang et al. 2011), from the representative marker of every Stage-1 unit |
| Floor | marker-density rule (8 for three-spined, 12 for nine-spined) |
| Permutations | 1,000 (consensus) and 200 (Simes) |
| EcoPeak relocation test | 10,000 relocations within chromosomes |
| Random seed | 1 |

The nine-spined panel adds `input/covariate.rds` (lineage). Both panels
permute populations within groups with `perm_within_groups()`, which draws
random numbers in the same order as the modules.

The scripts need the modules' Stage-1 bundles and data from a local
LDscnR-paper checkout (environment variable `LDSCNR_PAPER`). Without them,
`manuscript/` documents the reproduction but cannot rerun it.

```bash
Rscript manuscript/build_inputs.R    # data set folders from the module bundles
Rscript manuscript/run.R             # all steps, including the floor profile
Rscript manuscript/compare.R --out=manuscript/comparison_vs_v091
```

`manuscript/run.R` calls the workflow's steps directly, without the Manhattan
plots and the structure-alignment check. For a new data set, use
`run_dataset()` or `run_batch.R` instead.

| Manuscript result | Produced by | Compared in |
| --- | --- | --- |
| Stage-1 units, GRM, LD decay | `build_stage1()` | `comparison_vs_v091.md`, "Stage 1" rows |
| Units tested, significant units and outlier regions (consensus, Simes, LFMM) | `run_stage_A()`, `run_stage_B()` | rows per test |
| Permutation tests (null mean and p) | `run_stage_B()` | "permutation" rows |
| EcoPeak relocation tests | `run_stage_B()` with the EcoPeak annotation | "EcoPeaks" rows |
| Floor-stability tiers and EcoPeaks per tier | `run_floor_profile()`, `manuscript/run.R` | "floor tier" rows |

Result (LDscnR v0.9.1): all 108 checks in
`manuscript/comparison_vs_v091.md` match the modules' results. They cover the
Stage-1 units, the GRM, the LD-decay fit, every unit p-value, the significant
units, the outlier-region spans, the permutation counts, the EcoPeak
relocation tests, the floor-stability tiers and the EcoPeak tests per tier.
`manuscript/comparison_vs_stored.md` compares with the modules' earlier
outputs from LDscnR commits before v0.9.1. It matches except for the EcoPeak
relocation null mean, fold and p of the three three-spined tests (9 of 108
checks), which changed because v0.9.1 keeps relocated regions inside their
chromosome.

The LFMM2 p-values are not computed here. `manuscript/build_inputs.R` converts
the manuscript's stored LFMM2 F statistics (`LEA::lfmm2.test()`, K = 7) to
p-values with `pf(F, 1, n - 2)`, exactly as the modules do. These F
statistics are not divided by the genomic inflation factor.

## Scope

LDscnR-multi accepts p-values from any method, but it runs only LDscnR's own
outlier analysis. It has no interface for other outlier methods (OutFLANK,
BayeScan or pcadapt, for example).

## Known limitations

- **Plotting at genome scale.** `ldm_manhattan()` plots every marker. This
  works for the full panels used here, but a much larger panel would benefit
  from thinning or rasterising the points.
- **One data set folder per genotype panel and phenotype.** Screening many
  phenotypes on the same genotype panel needs one folder per phenotype. Each
  folder then rebuilds the same Stage 1 and stores its own copy of
  `genotypes.rds` and `map.rds`. Sharing Stage 1 between phenotypes would need
  a change of design.

## How to cite

If you use LDscnR-multi, please cite the manuscript and the LDscnR package:

- Kemppainen P, Guillaume F (2026). LDscnR: LD-aware complexity reduction
  for genome scans. *bioRxiv*.
  <https://doi.org/10.64898/2026.10.27.762116>. Submitted to *Molecular
  Ecology Resources*.
- Kemppainen P. LDscnR, version 0.9.1. R package.
  <https://github.com/genomeinsights/LDscnR>

`CITATION.cff` holds the same information in machine-readable form.

The example data and the manuscript reproduction use:

- Fang B, Kemppainen P, Momigliano P, Merilä J (2021). Population structure
  limits parallel evolution in sticklebacks. *Molecular Biology and
  Evolution* 38: 4205-4221. <https://doi.org/10.1093/molbev/msab144>.
  Genotype data: <https://doi.org/10.5281/zenodo.4722879> (CC BY 4.0). The
  example genotypes are two-chromosome subsets of these data.
- Kingman GAR et al. (2021). Predicting future from past: the genomic basis
  of recurrent and rapid stickleback evolution. *Science Advances* 7:
  eabg5285. <https://doi.org/10.1126/sciadv.abg5285>. EcoPeaks: FigShare
  project 162634 (<https://figshare.com/projects/kingsleyAssemblyHub/162634>).
  Used only by the manuscript reproduction and not distributed here.

Methods referred to above:

- Kang HM et al. (2010). Variance component model to account for sample
  structure in genome-wide association studies. *Nature Genetics* 42:
  348-354.
- Caye K, Jumentier B, Lepeule J, François O (2019). LFMM 2: fast and accurate
  inference of gene-environment associations in genome-wide studies.
  *Molecular Biology and Evolution* 36: 852-860.
- Yang J, Lee SH, Goddard ME, Visscher PM (2011). GCTA: a tool for
  genome-wide complex trait analysis. *American Journal of Human Genetics*
  88: 76-82.

## Licence

MIT (see [`LICENSE`](LICENSE)), the same licence as LDscnR.
