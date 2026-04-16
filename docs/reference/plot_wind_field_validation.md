# Plot wind-field validation figures

Plot wind-field validation figures

## Usage

``` r
plot_wind_field_validation(
  val,
  out = NULL,
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

- out:

  Model output list from
  [`run_hazard_model()`](https://tanerumit.github.io/ipdcstorm/reference/run_hazard_model.md).

- cfg:

  Optional `validation_cfg` object.

- out_dir:

  Directory to save plots in. Overrides `cfg$out_dir`.

- base_size:

  Base font size for plots. Overrides `cfg$advanced$base_size`.

## Value

Validation plot object with optional saved path metadata, or `NULL`.
