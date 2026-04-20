# Add hazard intensity and damage forcing from daily wind

Maps daily sustained wind speed to a bounded hazard intensity index and
a bounded daily damage rate for downstream impact calculations.

## Usage

``` r
add_damage_forcing(daily, V0 = 34, V1 = 120, p = 3, dmax = 0.02)
```

## Arguments

- daily:

  Tibble/data.frame with at least a `wind_kt` column.

- V0:

  Numeric scalar threshold wind in kt below which intensity is zero.

- V1:

  Numeric scalar wind in kt at which intensity saturates at one.

- p:

  Numeric scalar nonlinearity exponent.

- dmax:

  Numeric scalar maximum daily damage fraction.

## Value

Tibble with added `damage_intensity` and `damage_rate` columns.

## See also

[`damage_rate_from_wind`](https://tanerumit.github.io/ipdcstorm/reference/damage_rate_from_wind.md)

## Examples

``` r
daily <- tibble::tibble(wind_kt = c(20, 40, 80))
add_damage_forcing(daily)
#> Error in add_damage_forcing(daily): could not find function "add_damage_forcing"
```
