# Estimate site wind with Holland profile.

Computes site wind speed from storm intensity, storm size, and site
distance using a Holland-type radial wind profile with internal fallback
logic for missing radii.

## Usage

``` r
.estimate_site_wind_holland(
  Vmax_kt,
  r_km,
  R34_km,
  R50_km = NA,
  R64_km = NA,
  RMW_km,
  Pn = 1013,
  Pc = NA,
  lat = 18
)
```

## Arguments

- Vmax_kt:

  Numeric vector of maximum sustained wind speeds in knots.

- r_km:

  Numeric vector of distances from storm centre to site in kilometres.

- R34_km:

  Numeric vector of 34-kt radii in kilometres.

- R50_km:

  Numeric vector of 50-kt radii in kilometres.

- R64_km:

  Numeric vector of 64-kt radii in kilometres.

- RMW_km:

  Numeric vector of radius of maximum wind values in kilometres.

- Pn:

  Numeric vector of ambient pressures in hPa.

- Pc:

  Numeric vector of central pressures in hPa.

- lat:

  Numeric vector of storm latitudes in decimal degrees north.

## Value

Numeric vector of estimated sustained site winds in knots.
