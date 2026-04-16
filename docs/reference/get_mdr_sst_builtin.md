# Built-in MDR SST annual means from NOAA ERSST v5

Returns a tibble of annual-mean SST (degC) averaged over the Main
Development Region (MDR: 10-20N, 80-20W) for the Atlantic hurricane
season (Aug-Oct, the peak months most predictive of activity).

Values are derived from NOAA ERSST v5 monthly data, spatially averaged
over the MDR box, then temporally averaged over ASO (Aug-Sep-Oct) each
year.

These are provided as a built-in fallback so the model can run without
requiring users to download and process NetCDF files.

Source: NOAA/NCEI ERSST v5, accessed 2024. Reference: Huang et al.
(2017), J. Climate, 30, 8179-8205.

## Usage

``` r
get_mdr_sst_builtin()
```

## Value

Tibble with columns:

- `year`: integer calendar year.

- `sst_mdr_aso`: MDR ASO seasonal mean SST (degC).
