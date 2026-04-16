# IPCC AR6 / CMIP6 SSP MDR SST anomaly targets

Returns a tibble describing default MDR SST anomaly targets (degC
relative to the 1991-2020 baseline) used by the built-in SSP scenario
generator.

These are deliberately simple, stakeholder-facing targets used for the
piecewise-linear SST anomaly generator (to 2050 and 2100, then held
constant). They are not intended to reproduce a particular CMIP6 model
member.

## Usage

``` r
ipcc_ar6_sst_scenario_info()
```

## Value

Tibble with columns: scenario, source, delta_sst_2050, delta_sst_2100,
description.
