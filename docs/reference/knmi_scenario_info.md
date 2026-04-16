# KNMI'23 scenario reference table

Returns a tibble describing the four KNMI'23 climate scenarios for the
Dutch Caribbean, including their SSP equivalents, MDR SST anomaly
targets, and recommended Level 3 precipitation scaling modifiers.

The temperature axis (H/L) determines the emissions pathway and SST
warming. The moisture axis (d/n) affects precipitation response only
(Level 3).

MDR SST anomalies are relative to the 1991-2020 baseline, derived from:

- KNMI'23 global warming levels (IPCC AR6 constrained estimates)

- Tropical Atlantic SST scaling factor beta ~ 0.71 K/K

- Cross-validated against Hibbert et al. (2025) Caribbean CMIP6 SST
  study

## Usage

``` r
knmi_scenario_info()
```

## Value

Tibble with columns: scenario, ssp, variant, delta_sst_2050,
delta_sst_2100, precip_scale, air_temp_2050, air_temp_2100, description.
