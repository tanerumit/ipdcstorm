# Select directional wind radius.

Returns the quadrant-specific wind radius matching the bearing-derived
quadrant.

## Usage

``` r
.get_directional_radius(quadrant, r_ne, r_se, r_sw, r_nw)
```

## Arguments

- quadrant:

  Character vector of quadrant labels such as `"NE"`, `"SE"`, `"SW"`, or
  `"NW"`.

- r_ne:

  Numeric vector or scalar radius for the northeast quadrant.

- r_se:

  Numeric vector or scalar radius for the southeast quadrant.

- r_sw:

  Numeric vector or scalar radius for the southwest quadrant.

- r_nw:

  Numeric vector or scalar radius for the northwest quadrant.

## Value

Numeric vector of selected radii.
