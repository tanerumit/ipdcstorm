# Compute site winds for storm track points.

Computes distance, bearing, storm-motion diagnostics, resolved wind
radii, and site wind estimates for each track point relative to a fixed
target location.

## Usage

``` r
compute_site_winds_full(df, target_lat, target_lon)
```

## Arguments

- df:

  Track-point data frame containing the fields required for wind-field
  reconstruction.

- target_lat:

  Numeric scalar target latitude in decimal degrees.

- target_lon:

  Numeric scalar target longitude in decimal degrees.

## Value

Data frame with the original track-point fields plus derived wind-field
and site-wind columns.

## Note

Input data must already contain the expected IBTrACS-derived fields used
in the internal wind solver.

## See also

[`compute_storm_heading`](https://tanerumit.github.io/ipdcstorm/reference/compute_storm_heading.md),
[`estimate_R34_climo`](https://tanerumit.github.io/ipdcstorm/reference/estimate_R34_climo.md),
[`estimate_RMW_knaff`](https://tanerumit.github.io/ipdcstorm/reference/estimate_RMW_knaff.md)

## Examples

``` r
df <- data.frame(
  SID = c("AL012000", "AL012000"),
  iso_time = as.POSIXct(c("2000-08-01 00:00:00", "2000-08-01 06:00:00"), tz = "UTC"),
  lat = c(18, 18.2),
  lon = c(-63, -63.2),
  dist_km = c(NA_real_, NA_real_),
  wind_kt = c(60, 65),
  rmw_km = c(25, 25),
  r34_ne_nm = c(60, 60), r34_se_nm = c(55, 55), r34_sw_nm = c(50, 50), r34_nw_nm = c(55, 55),
  r50_ne_nm = c(30, 30), r50_se_nm = c(25, 25), r50_sw_nm = c(20, 20), r50_nw_nm = c(25, 25),
  r64_ne_nm = c(15, 15), r64_se_nm = c(10, 10), r64_sw_nm = c(10, 10), r64_nw_nm = c(10, 10),
  storm_speed_kt = c(12, 12)
)
compute_site_winds_full(df, target_lat = 18.05, target_lon = -63.05)
#> # A tibble: 2 × 44
#>   SID      iso_time              lat   lon dist_km wind_kt rmw_km r34_ne_nm
#>   <chr>    <dttm>              <dbl> <dbl>   <dbl>   <dbl>  <dbl>     <dbl>
#> 1 AL012000 2000-08-01 00:00:00  18   -63      7.68      60     25        60
#> 2 AL012000 2000-08-01 06:00:00  18.2 -63.2   23.0       65     25        60
#> # ℹ 36 more variables: r34_se_nm <dbl>, r34_sw_nm <dbl>, r34_nw_nm <dbl>,
#> #   r50_ne_nm <dbl>, r50_se_nm <dbl>, r50_sw_nm <dbl>, r50_nw_nm <dbl>,
#> #   r64_ne_nm <dbl>, r64_se_nm <dbl>, r64_sw_nm <dbl>, r64_nw_nm <dbl>,
#> #   storm_speed_kt <dbl>, bearing_to_target <dbl>, quadrant <chr>, nq34 <dbl>,
#> #   R34_nm_dir <dbl>, R34_nm_mean <dbl>, R34_nm <dbl>, R50_nm <dbl>,
#> #   R64_nm <dbl>, R34_km <dbl>, R50_km <dbl>, R64_km <dbl>, R34_mean_km <dbl>,
#> #   R50_mean_km <dbl>, R64_mean_km <dbl>, Vmax_kt <dbl>, heading_deg <dbl>, …
```
