# Read MDR SST from ERSST v5 NetCDF (optional)

Reads monthly ERSST v5 NetCDF data, subsets to the MDR box (10-20N,
80-20W), averages spatially, then computes ASO (Aug-Sep-Oct) seasonal
means per year.

Requires the `ncdf4` package. If not available, falls back to the
built-in reference data with a warning.

## Usage

``` r
read_mdr_sst_ersst(
  nc_path,
  mdr_lat = c(10, 20),
  mdr_lon = c(-80, -20),
  aso_months = 8L:10L
)
```

## Arguments

- nc_path:

  Character; path to ERSST v5 NetCDF file (e.g., "sst.mnmean.nc").

- mdr_lat:

  Range of latitudes for MDR (default: c(10, 20)).

- mdr_lon:

  Range of longitudes for MDR (default: c(-80, -20)). Note: ERSST uses
  0-360 longitude convention; this function converts internally.

- aso_months:

  Integer vector of months for seasonal average (default: 8:10).

## Value

Tibble with columns: year, sst_mdr_aso.
