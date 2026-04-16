# Plot hindcast validation figures

Creates the retained hindcast return-level comparison figure.

## Usage

``` r
plot_hindcast_validation(
  val,
  cfg = NULL,
  out_dir = NULL,
  base_size = NULL,
  save = TRUE
)
```

## Arguments

- val:

  Output from
  [`run_validation_suite()`](https://tanerumit.github.io/ipdcstorm/reference/run_validation_suite.md).

- cfg:

  Optional `validation_cfg` object. When provided, `out_dir` and
  `base_size` are read from the config (explicit arguments override).

- out_dir:

  Directory to save plots in. Overrides `cfg$out_dir`.

- base_size:

  Base font size for plots. Overrides `cfg$advanced$base_size`.

## Value

Validation plot object with optional saved path metadata, or `NULL`.
