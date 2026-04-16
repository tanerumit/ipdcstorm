# Parse IBTrACS wind radii fields (nautical miles) with caps and sentinels

Parse IBTrACS wind radii fields (nautical miles) with caps and sentinels

## Usage

``` r
.to_radii_nm(x, cap_nm)
```

## Arguments

- x:

  A vector (character/numeric) of radii values.

- cap_nm:

  Numeric scalar; values greater than this are treated as NA.

## Value

Numeric vector of radii (nm) with invalid/sentinel values set to NA.
