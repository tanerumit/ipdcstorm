# Estimate climatological R34.

Provides an empirical estimate of the 34-kt wind radius for storms
without observed radii data using a North Atlantic climatological
relationship.

## Usage

``` r
estimate_R34_climo(Vmax_kt, lat = 18)
```

## Arguments

- Vmax_kt:

  Numeric vector of maximum sustained wind speeds in knots.

- lat:

  Numeric vector of storm latitudes in decimal degrees north.

## Value

Numeric vector of estimated 34-kt wind radii in kilometres, with `NA`
where `Vmax_kt < 34` or inputs are non-finite.

## Note

This is used as a fallback when directional 34-kt radii are unavailable.

## See also

[`estimate_RMW_knaff`](https://tanerumit.github.io/ipdcstorm/reference/estimate_RMW_knaff.md),
[`compute_site_winds_full`](https://tanerumit.github.io/ipdcstorm/reference/compute_site_winds_full.md)

## Examples

``` r
estimate_R34_climo(Vmax_kt = c(30, 50, 100), lat = c(18, 18, 20))
#> [1]       NA 157.2867 217.2915
```
