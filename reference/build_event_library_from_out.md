# Convenience wrapper: build an event library from run_hazard_model() output

Extracts one location's track-point and event tables from
[`run_hazard_model()`](https://tanerumit.github.io/ipdcstorm/reference/run_hazard_model.md)
output and forwards them to
[`build_event_library()`](https://tanerumit.github.io/ipdcstorm/reference/build_event_library.md).

## Usage

``` r
build_event_library_from_out(out, location, ..., seed = NULL)
```

## Arguments

- out:

  List returned by
  [`run_hazard_model()`](https://tanerumit.github.io/ipdcstorm/reference/run_hazard_model.md).

- location:

  Character scalar target location name.

- ...:

  Additional arguments passed to
  [`build_event_library()`](https://tanerumit.github.io/ipdcstorm/reference/build_event_library.md).

- seed:

  Optional integer seed for deterministic library construction.

## Value

List in the format returned by
[`build_event_library()`](https://tanerumit.github.io/ipdcstorm/reference/build_event_library.md).

## See also

[`build_event_library`](https://tanerumit.github.io/ipdcstorm/reference/build_event_library.md)

## Examples

``` r
if (FALSE) { # \dontrun{
lib <- build_event_library_from_out(out, location = "Saba")
lib$events
} # }
```
