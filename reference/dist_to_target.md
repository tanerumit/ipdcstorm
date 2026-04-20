# Compute great-circle distance from storm center to target.

Computes Haversine distance in kilometres from each storm track point to
a fixed target coordinate.

## Usage

``` r
dist_to_target(lat, lon, t_lat, t_lon)
```

## Arguments

- lat:

  Numeric vector of storm latitudes in decimal degrees.

- lon:

  Numeric vector of storm longitudes in decimal degrees.

- t_lat:

  Numeric scalar target latitude in decimal degrees.

- t_lon:

  Numeric scalar target longitude in decimal degrees.

## Value

Numeric vector of distances in kilometres with the same length as `lat`
and `lon`.

## See also

[`calculate_bearing`](https://tanerumit.github.io/ipdcstorm/reference/calculate_bearing.md)

## Examples

``` r
dist_to_target(lat = c(18, 18.1), lon = c(-63, -63.1), t_lat = 18.05, t_lon = -63.05)
#> [1] 7.680746 7.679709
```
