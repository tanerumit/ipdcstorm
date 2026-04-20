# Compute bearing from storm center to target.

Computes great-circle initial bearing in degrees from each storm track
point to a fixed target coordinate.

## Usage

``` r
calculate_bearing(lat, lon, t_lat, t_lon)
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

Numeric vector of bearings in degrees with the same length as `lat` and
`lon`.

## See also

[`dist_to_target`](https://tanerumit.github.io/ipdcstorm/reference/dist_to_target.md)

## Examples

``` r
calculate_bearing(lat = c(18, 18.1), lon = c(-63, -63.1), t_lat = 18.05, t_lon = -63.05)
#> [1] -43.72507 136.26768
```
