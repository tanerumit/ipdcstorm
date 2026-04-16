# Retrieve available SST scenario definitions

Returns a combined scenario table for all available climate-information
sources. If KNMI'23 support is available (knmi_scenario_info()), those
scenarios are included as well.

## Usage

``` r
sst_scenario_info(source = c("all", "ipcc_ar6", "knmi23"))
```

## Arguments

- source:

  Character; one of "all", "ipcc_ar6", "knmi23".

## Value

Tibble with at least: scenario, source, delta_sst_2050, delta_sst_2100.
