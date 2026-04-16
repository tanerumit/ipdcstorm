# Look up a time-slice scenario SST shift

Returns a scalar `delta_sst` for a future climate target by
interpolating the scenario targets in `sst_scenario_info("all")` at
`target_year`.

## Usage

``` r
get_scenario_delta(scenario, target_year, baseline_years = 1991L:2020L)
```

## Arguments

- scenario:

  Character scalar naming a scenario in `sst_scenario_info("all")`.

- target_year:

  Numeric scalar target year used to derive `delta_sst`.

- baseline_years:

  Integer vector of climatological reference years.

## Value

Numeric scalar `delta_sst` in degC.

## Examples

``` r
get_scenario_delta("ssp585", target_year = 2050)
#> [1] 1
```
