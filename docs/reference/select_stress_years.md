# Select k representative stress-test years with good sample coverage

Picks a small subset of `k` years from the candidate set that
collectively represent the full distribution of stress severity. Three
selection strategies are available:

- `"stratified"` (default):

  Divides the `composite_score` range into `k` equal-count bins and
  returns the year whose score is closest to each bin's median.
  Guarantees one representative per severity level — from near-miss to
  worst-case.

- `"diverse"`:

  Greedy maximin selection in normalised metric space. Starts from the
  most extreme year, then iteratively adds the year that is furthest
  from all already-selected years. Maximises spread in the
  multi-dimensional metric space.

- `"top"`:

  Returns the `k` years with the highest `composite_score`. Useful when
  the goal is a set of purely high-severity scenarios.

**Portfolio mode vs per-location mode:** When `metrics` has no
`location` column (the typical output of
[`aggregate_stress_metrics()`](https://tanerumit.github.io/ipdcstorm/reference/aggregate_stress_metrics.md)
with multiple locations), all rows are treated as a single pool and `k`
years are selected for the full portfolio. When a `location` column is
present (single-location input or a pre-aggregated table), selection is
performed independently within each location and the result is sorted by
`location` then `selection_rank`.

If `composite_score` is absent from `metrics`,
[`aggregate_stress_metrics()`](https://tanerumit.github.io/ipdcstorm/reference/aggregate_stress_metrics.md)
is called internally using the supplied `weights` and `metrics_used`.

## Usage

``` r
select_stress_years(
  metrics,
  k,
  method = c("stratified", "diverse", "top"),
  metrics_used = NULL,
  weights = NULL
)
```

## Arguments

- metrics:

  Tibble from
  [`compute_stress_year_metrics()`](https://tanerumit.github.io/ipdcstorm/reference/compute_stress_year_metrics.md)
  or
  [`aggregate_stress_metrics()`](https://tanerumit.github.io/ipdcstorm/reference/aggregate_stress_metrics.md).

- k:

  Positive integer: number of years to select (portfolio-wide, or per
  location when a `location` column is present). If the candidate pool
  has \\\le k\\ rows, all are returned.

- method:

  Character scalar selection strategy; one of `"stratified"` (default),
  `"diverse"`, or `"top"`.

- metrics_used:

  Forwarded to
  [`aggregate_stress_metrics()`](https://tanerumit.github.io/ipdcstorm/reference/aggregate_stress_metrics.md)
  if `composite_score` is not already present.

- weights:

  Forwarded to
  [`aggregate_stress_metrics()`](https://tanerumit.github.io/ipdcstorm/reference/aggregate_stress_metrics.md)
  if `composite_score` is not already present.

## Value

Subset of `metrics` rows for the selected years with an additional
integer column `selection_rank` (1 = first/most important pick, as
determined by the method). In portfolio mode rows are sorted by
`selection_rank`; in per-location mode by `location` then
`selection_rank`.

## See also

[`compute_stress_year_metrics`](https://tanerumit.github.io/ipdcstorm/reference/compute_stress_year_metrics.md),
[`aggregate_stress_metrics`](https://tanerumit.github.io/ipdcstorm/reference/aggregate_stress_metrics.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Portfolio mode — k years for the full location set (typical use)
scored   <- aggregate_stress_metrics(metrics)       # no location column
selected <- select_stress_years(scored, k = 10)

# Diversity-maximising selection with custom scoring
selected <- select_stress_years(
  scored,
  k            = 8,
  method       = "diverse",
  metrics_used = c("peak_wind_kt", "cum_damage", "n_events"),
  weights      = c(2, 1, 0.5)
)
} # }
```
