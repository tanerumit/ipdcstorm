# Get KNMI'23 Level 3 cc_params adjusted for scenario variant

Returns the default cc_params list with the precipitation scaling factor
adjusted for the dry (d) or wet (n) variant of the KNMI'23 scenarios.

## Usage

``` r
knmi_cc_params(scenario, base_params = NULL)
```

## Arguments

- scenario:

  Character; one of "knmi_Ld", "knmi_Ln", "knmi_Hd", "knmi_Hn".

- base_params:

  Optional named list; base cc_params to modify. If NULL, uses
  [`default_cc_params()`](https://tanerumit.github.io/ipdcstorm/reference/default_cc_params.md).

## Value

Named list of cc_params with adjusted precip_scale.
