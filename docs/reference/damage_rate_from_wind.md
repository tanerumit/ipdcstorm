# Bounded power-law damage rate from wind speed

Converts sustained wind speed to a bounded daily damage fraction using a
thresholded power-law response calibrated at `V_ref`.

## Usage

``` r
damage_rate_from_wind(
  wind_kt,
  thr = 34,
  V_ref = 80,
  d_ref = 0.03,
  p = 3,
  d_max = 0.1
)
```

## Arguments

- wind_kt:

  Numeric vector of sustained wind speeds in kt.

- thr:

  Numeric scalar threshold wind in kt below which damage is zero.

- V_ref:

  Numeric scalar reference wind in kt where damage equals `d_ref`.

- d_ref:

  Numeric scalar damage fraction at `V_ref`.

- p:

  Numeric scalar exponent controlling nonlinearity.

- d_max:

  Numeric scalar upper cap on daily damage fraction.

## Value

Numeric vector of daily damage rates.

## See also

[`add_damage_forcing`](https://tanerumit.github.io/ipdcstorm/reference/add_damage_forcing.md)

## Examples

``` r
damage_rate_from_wind(c(20, 40, 80))
#> Error in damage_rate_from_wind(c(20, 40, 80)): could not find function "damage_rate_from_wind"
```
