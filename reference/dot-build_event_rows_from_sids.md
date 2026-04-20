# Build event rows from shared SIDs using a location-specific library

For each SID in `sids` that exists in `lib$events`, constructs one event
row in the schema produced by
[`sample_events_for_year_extended()`](https://tanerumit.github.io/ipdcstorm/reference/sample_events_for_year_extended.md).
SIDs absent from the location library are silently skipped, so the
number of rows returned may be less than `length(sids)`.

## Usage

``` r
.build_event_rows_from_sids(lib, sids, sev, year, counter_start)
```

## Arguments

- lib:

  Event library from
  [`build_event_library_from_out()`](https://tanerumit.github.io/ipdcstorm/reference/build_event_library_from_out.md).

- sids:

  Character vector of SIDs to look up.

- sev:

  Character scalar severity; one of `"TS"` or `"HUR"`.

- year:

  Integer scalar calendar year for date construction and `event_id`
  encoding.

- counter_start:

  Integer scalar starting value for the within-year event counter (used
  to produce unique `event_id` values across TS and HUR draws for the
  same year and location).

## Value

Named list with elements:

- `rows`:

  Tibble of event rows (zero or more).

- `n_filled`:

  Integer count of rows produced.

- `counter`:

  Updated event counter.
