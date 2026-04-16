# Add forward-motion asymmetry.

Applies a translation-speed asymmetry adjustment to a symmetric site
wind estimate using the angle between storm heading and site bearing.

## Usage

``` r
.add_forward_motion_asymmetry(
  V_site_base_kt,
  storm_speed_kt,
  r_km,
  bearing_to_target,
  storm_heading,
  RMW_km = 40
)
```

## Arguments

- V_site_base_kt:

  Numeric vector of symmetric site wind estimates in knots.

- storm_speed_kt:

  Numeric vector of storm translation speeds in knots.

- r_km:

  Numeric vector of distances from storm centre to site in kilometres.

- bearing_to_target:

  Numeric vector of bearings from storm centre to site in degrees.

- storm_heading:

  Numeric vector of storm motion headings in degrees.

- RMW_km:

  Numeric vector of radius of maximum wind values in kilometres.

## Value

Numeric vector of asymmetry-adjusted site winds in knots.
