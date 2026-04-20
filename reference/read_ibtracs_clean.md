# Read and clean IBTrACS CSV (USA fields + pressure/structure when present)

Reads an IBTrACS CSV and returns a cleaned track-point tibble using USA
best-track variables. Adds central pressure `pres_hpa` (USA_PRES with
WMO_PRES fallback) and optional structure fields
(USA_POCI/USA_ROCI/USA_RMW) when present.

## Usage

``` r
read_ibtracs_clean(
  ibtracs_csv,
  basin = "NA",
  season = NULL,
  keep_all = FALSE,
  verbose = TRUE
)
```

## Arguments

- ibtracs_csv:

  Path to an IBTrACS CSV file.

- basin:

  Character vector or NULL; filter BASIN if provided (e.g. "NA").

- season:

  Integer vector or NULL; filter SEASON if provided.

- keep_all:

  Logical; keep original columns if TRUE.

- verbose:

  Logical; print read diagnostics if TRUE.

## Value

Tibble of cleaned track points.
