# Compute a weighted composite stress score from per-year metrics

Aggregates metrics across locations (portfolio mode), normalises the
chosen metric columns to \\\[0, 1\]\\ using min-max scaling, then
combines them into a single `composite_score` via a weighted mean.

**Portfolio mode (default, multiple locations):** When the input
contains more than one location, metric values are first aggregated
across locations for each `sim_year` using `location_agg` (`"mean"`,
`"max"`, or `"sum"`). The result has one row per `sim_year` with no
`location` column. This is the recommended path when a single
representative year-set is needed for the full location portfolio.

**Single-location mode:** When only one location is present the input
rows are scored directly and the `location` column is retained.

Normalisation is relative to the rows supplied, not the full ensemble.
Pass the full-ensemble metrics table if you need ensemble-relative
scores.

## Usage

``` r
aggregate_stress_metrics(
  metrics,
  metrics_used = NULL,
  weights = NULL,
  location_agg = c("mean", "max", "sum")
)
```

## Arguments

- metrics:

  Tibble returned by
  [`compute_stress_year_metrics()`](https://tanerumit.github.io/ipdcstorm/reference/compute_stress_year_metrics.md).

- metrics_used:

  Character vector of metric column names to include in the composite.
  Defaults to all numeric columns except `sim_year` and any previously
  computed score/rank columns.

- weights:

  Optional numeric vector specifying relative importance.

  - **Named vector** — names must match `metrics_used`; missing names
    receive weight 0.

  - **Unnamed vector** — must have the same length as `metrics_used`;
    assigned in order.

  - `NULL` (default) — uniform weights across all selected metrics.

  Weights are automatically rescaled to sum to 1.

- location_agg:

  Character scalar controlling how metrics are aggregated across
  locations before scoring. One of `"mean"` (default), `"max"`, or
  `"sum"`. Ignored when only one location is present.

## Value

Tibble with two additional columns:

- `composite_score`:

  Weighted mean of min-max-normalised metric values; 0 = least extreme,
  1 = most extreme within the provided set.

- `composite_rank`:

  Integer rank (1 = highest score).

When multiple locations are present (portfolio mode) the returned tibble
has columns `sim_year`, the metric columns, `composite_score`, and
`composite_rank` — no `location` column.

## See also

[`compute_stress_year_metrics`](https://tanerumit.github.io/ipdcstorm/reference/compute_stress_year_metrics.md),
[`select_stress_years`](https://tanerumit.github.io/ipdcstorm/reference/select_stress_years.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Portfolio mode — single year-set across all locations (default)
scored <- aggregate_stress_metrics(metrics)

# Up-weight intensity and damage, use max across locations
scored <- aggregate_stress_metrics(
  metrics,
  metrics_used = c("peak_wind_kt", "cum_damage", "n_events"),
  weights      = c(peak_wind_kt = 2, cum_damage = 2, n_events = 1),
  location_agg = "max"
)
} # }
```
