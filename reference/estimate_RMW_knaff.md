# Estimate climatological radius of maximum wind.

Estimates radius of maximum wind from storm intensity and latitude using
a Knaff and Zehr climatological relationship.

## Usage

``` r
estimate_RMW_knaff(Vmax_kt, lat = 18)
```

## Arguments

- Vmax_kt:

  Numeric vector of maximum sustained wind speeds in knots.

- lat:

  Numeric vector of storm latitudes in decimal degrees north.

## Value

Numeric vector of estimated radius of maximum wind values in kilometres.

## See also

[`estimate_R34_climo`](https://tanerumit.github.io/ipdcstorm/reference/estimate_R34_climo.md),
[`compute_site_winds_full`](https://tanerumit.github.io/ipdcstorm/reference/compute_site_winds_full.md)

## Examples

``` r
estimate_RMW_knaff(Vmax_kt = c(40, 80, 120), lat = c(15, 18, 22))
#> [1] 97.27667 96.43383 75.50000
```
