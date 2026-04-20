# Observation quality bias thresholds for wind field validation

Returns acceptable absolute bias percentages per observation quality
tier. Quality tiers reflect measurement uncertainty:

- **A** (Direct): Station measurement, instrument survived. Threshold:
  15%.

- **B** (Converted): 10-min reading with averaging-period conversion.
  Threshold: 25%.

- **C** (Estimated): NHC best-track or indirect estimate. Threshold:
  35%.

## Usage

``` r
.wf_quality_thresholds()
```

## Value

Named numeric vector: quality tier -\> acceptable \|bias\| (%).
