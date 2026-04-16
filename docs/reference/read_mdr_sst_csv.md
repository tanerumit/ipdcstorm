# Read MDR SST from a CSV file

Reads a user-supplied CSV containing annual MDR SST values. The CSV must
have at minimum a `year` column and an SST column (name specified by
`sst_col`). This allows users to supply their own ERSST v5 extractions
or alternative SST products.

## Usage

``` r
read_mdr_sst_csv(csv_path, sst_col = "sst_mdr_aso", year_col = "year")
```

## Arguments

- csv_path:

  Character; path to the CSV file.

- sst_col:

  Character; name of the SST column (default: "sst_mdr_aso").

- year_col:

  Character; name of the year column (default: "year").

## Value

Tibble with columns: year, sst_mdr_aso.
