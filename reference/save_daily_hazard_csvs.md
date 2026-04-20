# Save daily hazard output to per-location CSV files

Writes a named list of daily hazard tables to disk, creating one CSV per
location. Filenames include both the scenario label and the location
name.

## Usage

``` r
save_daily_hazard_csvs(daily, scenario, out_dir = file.path("output", "raw"))
```

## Arguments

- daily:

  Named list of per-location tibbles, typically returned by
  [`generate_daily_hazard_impact_spatial()`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_hazard_impact_spatial.md).

- scenario:

  Character scalar scenario label to embed in each filename.

- out_dir:

  Character scalar output directory. Created if it does not already
  exist.

## Value

Named character vector of file paths, one per location.

## Examples

``` r
daily <- list(
  Saba = tibble::tibble(
    sim_year = 1L,
    doy = 1L,
    wind_kt = 0,
    surge_m = NA_real_,
    event_id = NA_character_,
    pressure_hpa = NA_real_,
    pressure_deficit_hpa = NA_real_,
    rmw_km = NA_real_,
    damage_intensity = 0,
    damage_rate = 0,
    cum_damage = 0
  )
)
save_daily_hazard_csvs(daily, scenario = "baseline", out_dir = tempdir())
#>                                                                          Saba 
#> "C:\\Users\\taner\\AppData\\Local\\Temp\\RtmpYZVK5j/daily_baseline__Saba.csv" 
```
