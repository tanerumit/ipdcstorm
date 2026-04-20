# Save standard hazard visualization plots as PNG files

Builds a standard set of five hazard plots, applies consistent titles,
subtitles, and base theme sizing, and saves each plot to `output_dir` as
a PNG file.

## Usage

``` r
save_hazard_viz_plots(
  daily,
  output_dir,
  location_name,
  width = 9,
  height = 6,
  dpi = 300,
  base_size = 11,
  thr_tc = 34,
  thr_hur = 64
)
```

## Arguments

- daily:

  A data frame/tibble or named list of daily tables. Inputs must contain
  the columns required by
  [`plot_monthly_events()`](https://tanerumit.github.io/ipdcstorm/reference/plot_monthly_events.md),
  [`plot_annual_counts()`](https://tanerumit.github.io/ipdcstorm/reference/plot_annual_counts.md),
  [`plot_wind_timeseries()`](https://tanerumit.github.io/ipdcstorm/reference/plot_wind_timeseries.md),
  [`plot_wind_distribution()`](https://tanerumit.github.io/ipdcstorm/reference/plot_wind_distribution.md),
  and
  [`plot_intensity_duration()`](https://tanerumit.github.io/ipdcstorm/reference/plot_intensity_duration.md).
  Named lists are flattened using the list names as `location`.

- output_dir:

  Character scalar path to the folder where PNG files are saved. The
  directory is created if it does not already exist.

- location_name:

  Character scalar used in standardized plot titles.

- width:

  Numeric scalar plot width in inches passed to
  [`ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

- height:

  Numeric scalar plot height in inches passed to
  [`ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

- dpi:

  Numeric scalar resolution passed to
  [`ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

- base_size:

  Numeric scalar base font size applied via `plot_theme()`.

- thr_tc:

  Numeric scalar; tropical-cyclone threshold in knots.

- thr_hur:

  Numeric scalar; hurricane threshold in knots.

## Value

A named list with components `plots` and `paths`.

## Examples

``` r
if (FALSE) { # \dontrun{
d <- data.frame(
  date = as.Date("2001-01-01") + 0:364,
  wind_kt = pmax(1, 20 + rnorm(365, sd = 8)),
  event_id = sample(c(NA, 1:8), 365, replace = TRUE),
  location = "A",
  sim_year = 2001
)
save_hazard_viz_plots(d, tempdir(), "A")
} # }
```
