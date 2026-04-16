# Compute SST anomalies relative to a climatological baseline

Computes `delta_SST_t = SST_t - SST_clim`, where `SST_clim` is the mean
SST over the specified baseline period (default: 1991-2020, the current
WMO standard climatological normal).

## Usage

``` r
compute_sst_anomaly(sst_df, baseline_years = 1991L:2020L)
```

## Arguments

- sst_df:

  Tibble with columns: year, sst_mdr_aso.

- baseline_years:

  Integer vector of years defining the climatological baseline (default:
  1991:2020).

## Value

The input tibble with added columns:

- sst_clim:

  Climatological mean SST (degC) over the baseline period.

- sst_anomaly:

  SST anomaly (degC) relative to baseline.

## Examples

``` r
sst <- get_mdr_sst_builtin()
sst_anom <- compute_sst_anomaly(sst, baseline_years = 1991:2020)
head(sst_anom)
#> # A tibble: 6 × 4
#>    year sst_mdr_aso sst_clim sst_anomaly
#>   <int>       <dbl>    <dbl>       <dbl>
#> 1  1970        26.9     27.1      -0.231
#> 2  1971        26.6     27.1      -0.451
#> 3  1972        26.7     27.1      -0.371
#> 4  1973        26.8     27.1      -0.291
#> 5  1974        26.6     27.1      -0.471
#> 6  1975        26.7     27.1      -0.411
```
