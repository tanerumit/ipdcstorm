# Resolve track-point radius of maximum wind.

Resolves the radius of maximum wind at each track point by prioritising
valid observations, then radii-based estimates, then climatological
estimates.

## Usage

``` r
.resolve_trackpoint_rmw_km(
  rmw_obs_km,
  R64_mean_km,
  R50_mean_km,
  R34_mean_km,
  Vmax_kt,
  lat
)
```

## Arguments

- rmw_obs_km:

  Numeric vector of observed radius of maximum wind values in
  kilometres.

- R64_mean_km:

  Numeric vector of mean 64-kt radii in kilometres.

- R50_mean_km:

  Numeric vector of mean 50-kt radii in kilometres.

- R34_mean_km:

  Numeric vector of mean 34-kt radii in kilometres.

- Vmax_kt:

  Numeric vector of maximum sustained wind speeds in knots.

- lat:

  Numeric vector of storm latitudes in decimal degrees north.

## Value

Numeric vector of resolved radius of maximum wind values in kilometres.
