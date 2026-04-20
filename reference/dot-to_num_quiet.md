# Quiet numeric parsing for messy IBTrACS string fields

Converts character inputs to numeric robustly by stripping common text
decorations (e.g., "deg", "degrees_north") and treating known
placeholders as missing.

## Usage

``` r
.to_num_quiet(x)
```

## Arguments

- x:

  A vector (character/numeric) to parse as numeric.

## Value

Numeric vector with non-parsable values set to NA.
